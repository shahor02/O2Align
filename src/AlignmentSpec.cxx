// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <cmath>
#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <memory>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include <Eigen/Dense>
#include <GblTrajectory.h>
#include <GblData.h>
#include <GblPoint.h>
#include <GblMeasurement.h>
#include <MilleBinary.h>
#include <nlohmann/json.hpp>
#include <boost/interprocess/sync/named_semaphore.hpp>
#include <TFile.h>
#include <TMethodCall.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include "Headers/DataHeader.h"
#include "Framework/CCDBParamSpec.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "Framework/TimingInfo.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DetectorsVertexing/PVertexer.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "Steer/MCKinematicsReader.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "CommonUtils/StringUtils.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "MathUtils/Utils.h"
#include "ITStracking/IOUtils.h"
#include "ITSMFTTracking/MathUtils.h"
#include "ITS3Reconstruction/IOUtils.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "TPCCalibration/VDriftHelper.h"
#include "TPCFastTransformPOD.h"
#include "GPUO2ExternalUser.h"
#include "GPUParam.h"
#include "O2Align/TrackFit.h"
#include "O2Align/AlignmentSpec.h"
#include "O2Align/Params.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/Volume.h"
#include "O2Align/MisalignmentUtils.h"
#include "O2Align/SensorITS.h"
#include "O2Align/DetectorITS.h"
#include "O2Align/DetectorPVT.h"
#include "O2Align/DetectorTPC.h"
#include "O2Align/DetectorTRD.h"
#include "O2Align/DetectorTOF.h"

namespace o2::alignrs
{
using namespace o2::framework;
using DetID = o2::detectors::DetID;
using DataRequest = o2::globaltracking::DataRequest;
using PVertex = o2::dataformats::PrimaryVertex;
using V2TRef = o2::dataformats::VtxTrackRef;
using VTIndex = o2::dataformats::VtxTrackIndex;
using GTrackID = o2::dataformats::GlobalTrackID;
using GlobalIDSet = std::array<GTrackID, GTrackID::NSources>;
using TrackD = o2::track::TrackParCovD;
using ClusterD = o2::BaseCluster<double>;
using ClusterF = o2::BaseCluster<float>;

namespace
{
DerivativeContext makeDerivativeContext(const FrameInfoExt& frame, const TrackD& trk)
{
  const auto slopes = TrackSlopes::computeTrackSlopes(trk.getSnp(), trk.getTgl());
  const bool isITS3 = o2::its3::constants::detID::isDetITS3(frame.cluster.getSensorID());
  return {.detID = static_cast<int>(frame.label.det()),
          .volID = static_cast<int>(frame.label.id()),
          .sensorID = isITS3 ? o2::its3::constants::detID::getSensorID(frame.cluster.getSensorID()) : -1,
          .layerID = isITS3 ? o2::its3::constants::detID::getDetID2Layer(frame.cluster.getSensorID()) : -1,
          .measX = frame.x,
          .measAlpha = frame.alpha,
          .measZ = frame.cluster.getZ(),
          .trkY = trk.getY(),
          .trkZ = trk.getZ(),
          .snp = trk.getSnp(),
          .tgl = trk.getTgl(),
          .dydx = slopes.dydx,
          .dzdx = slopes.dzdx};
}

Matrix26 getRigidBodyBaseDerivatives(const DerivativeContext& ctx)
{
  static const RigidBodyDOFSet sRigidBodyBasis;
  Eigen::MatrixXd dyn(2, sRigidBodyBasis.nDOFs());
  sRigidBodyBasis.fillDerivatives(ctx, dyn);
  return dyn;
}

/// d(prediction)/d(rigid-body DOFs of the measurement leaf, in its local frame), thesis A.57-A.58:
/// the base derivative (A.57) is evaluated in the tracking frame wrt rotations about the point
/// (x,0,0) of the measurement plane, then transformed to the local frame of the leaf with the
/// jacobian J_L2T for that pivot.
Matrix26 getLeafRigidBodyDerivatives(const Volume& leaf, const FrameInfoExt& frame, const DerivativeContext& ctx)
{
  // the tracking frame of an ITS3 tile depends on the point, the other sensors have a fixed one
  TGeoHMatrix t2lPoint;
  const TGeoHMatrix* t2l = &leaf.getT2L();
  if (!leaf.hasFixedTrackingFrame()) {
    t2lPoint = leaf.computeT2L(frame.alpha);
    t2l = &t2lPoint;
  }
  const double posTrk[3] = {frame.x, 0., 0.};
  double posLoc[3];
  t2l->LocalToMaster(posTrk, posLoc);
  Matrix66 jacL2T;
  Volume::computeJacobianL2T(*t2l, posLoc, jacL2T);
  return getRigidBodyBaseDerivatives(ctx) * jacL2T;
}

/// precision matrix of a 2D measurement in the tracking frame, accounting for the Y-Z correlation
/// (e.g. the pad tilt of the TRD). Returns false if the covariance is not positive definite.
bool getMeasurementPrecision(const o2::BaseCluster<float>& cluster, Eigen::Matrix2d& prec)
{
  Eigen::Matrix2d cov;
  cov << cluster.getSigmaY2(), cluster.getSigmaYZ(), cluster.getSigmaYZ(), cluster.getSigmaZ2();
  if (cov(0, 0) <= 0. || cov.determinant() <= 0.) {
    return false;
  }
  prec = cov.inverse();
  return true;
}

/// precision of the multiple-scattering kinks of the GBL local slopes, which are the ALICE (snp, tgl):
/// a scattering angle theta0 in each of the two orthogonal directions transverse to the track gives
/// var(snp) = theta0^2 (1-snp^2)(1+tgl^2) and var(tgl) = theta0^2 (1+tgl^2)^2, as in the material
/// correction of the KF (TrackParCov::correctForMaterial)
Eigen::Vector2d getScatteringPrecision(const TrackD& trk, double theta0)
{
  const double theta2 = theta0 * theta0, snp = trk.getSnp(), tgl2 = trk.getTgl() * trk.getTgl();
  return {1. / (theta2 * (1. - snp * snp) * (1. + tgl2)), 1. / (theta2 * (1. + tgl2) * (1. + tgl2))};
}

/// reorder the ALICE (Y,Z,Snp,Tgl,Q/Pt) transport jacobian to the GBL (Q/Pt,Snp,Tgl,Y,Z) convention
gbl::Matrix5d toGBLOrder(const gbl::Matrix5d& jacALICE)
{
  constexpr int perm[5] = {4, 2, 3, 0, 1};
  gbl::Matrix5d jacGBL;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      jacGBL(i, j) = jacALICE(perm[i], perm[j]);
    }
  }
  return jacGBL;
}
} // namespace

class AlignmentSpec final : public Task
{
 public:

  struct ProcStat {
    enum {
      kInput,
      kAccepted,
      kNStatCl
    };
    enum {
      kVertices,
      kTracks,
      kTracksWithVertex,
      kCosmic,
      kMaxStat
    };
    std::array<std::array<size_t, kMaxStat>, kNStatCl> data{};
    void print() const;
  };

  /// statistics of the GBL trajectories construction and fit
  struct GBLStat {
    int failedProp{0};       // tracks lost on the propagation between the points
    int construct{0};        // trajectories which could not be constructed
    int fit{0};              // successfully fitted trajectories
    int fitFail{0};          // failed fits
    int chi2Rej{0};          // fits rejected by the chi2/ndf cut
    double chi2Sum{0};       // sum of the chi2 of the accepted fits
    double lostWeightSum{0}; // sum of the lost weights of the accepted fits
    int ndfSum{0};           // sum of the ndf of the accepted fits
    void print() const
    {
      LOGP(info, "\tGBL SUMMARY: fitted {}, construction failed {}, fit failed {}, chi2Ndf rejected {}, propagation failed {}",
           fit, construct, fitFail, chi2Rej, failedProp);
      LOGP(info, "\t\tGBL Chi2/Ndf = {}, LostWeight = {}", ndfSum ? chi2Sum / ndfSum : -1., lostWeightSum);
    }
  };


  ~AlignmentSpec() final = default;
  AlignmentSpec(const AlignmentSpec&) = delete;
  AlignmentSpec(AlignmentSpec&&) = delete;
  AlignmentSpec& operator=(const AlignmentSpec&) = delete;
  AlignmentSpec& operator=(AlignmentSpec&&) = delete;
  AlignmentSpec(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr, GTrackID::mask_t src, DetID::mask_t dmask, bool useMC, bool withITS3, o2::alignrs::OutputEnum out)
    : mDataRequest(dr), mGGCCDBRequest(gr), mTracksSrcMask(src), mUseMC(useMC), mIsITS3(withITS3), mOutOpt(out)
  {
    if (dmask[DetID::ITS]) {
      mITS = std::make_unique<DetectorITS>(withITS3);
      mDetectors.push_back(mITS.get());
    }
    if (dmask[DetID::TPC]) {
      mTPC = std::make_unique<DetectorTPC>();
      mDetectors.push_back(mTPC.get());
    }
    if (dmask[DetID::TRD]) {
      mTRD = std::make_unique<DetectorTRD>();
      mDetectors.push_back(mTRD.get());
    }
    if (dmask[DetID::TOF]) {
      mTOF = std::make_unique<DetectorTOF>();
      mDetectors.push_back(mTOF.get());
    }
    mPVT = std::make_unique<DetectorPVT>(); // virtual, always created: it has no data of its own
    mDetectors.push_back(mPVT.get());
  }

  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;
  void process();

 private:
  void updateTimeDependentParams(ProcessingContext& pc);
  void buildHierarchy();

  // calculate the transport jacobian for points FROM and TO numerically via ridder's method
  // this assumes the track is already at point FROM and will be extrapolated to TO's x (xTo)
  // method does not modify the original track
  bool getTransportJacobian(const TrackD& track, double xTo, double alphaTo, gbl::Matrix5d& jac, gbl::Matrix5d& err);

  // refit the primary vertex vtxOrig with those tracks of resTracks which are its contributors and
  // were successfully refitted in the current alignment, the result is put to vtxRefit.
  // Returns false if the vertex has too few refitted contributors or the refit failed.
  bool refitPV(const PVertex& vtxOrig, const std::vector<Track>& resTracks, PVertex& vtxRefit);

  // fill the GBL points of the track frames from the ipStart slot outward
  // mvPriorCovScale > 0: impose on the vertex point, besides the refitted vertex, the mean vertex prior
  // with the luminous region covariance scaled by this factor (per-track PV constraint mode)
  bool fillGBLPoints(Track& resTrack, int ipStart, bool skipFirstMeas, std::vector<gbl::GblPoint>& points, double mvPriorCovScale = 0.);

  // fit the constructed trajectory and store it to gblTraj if the Mille data is requested
  bool fitGBLTrajectory(gbl::GblTrajectory& traj, float kfChi2Ndf, std::vector<gbl::GblTrajectory>& gblTraj, FitInfo& fitOut);

  // build and fit the GBL trajectory of a single track, accounting its frames from the ipStart slot outward
  bool buildGBLTrack(Track& resTrack, int ipStart, std::vector<gbl::GblTrajectory>& gblTraj, double mvPriorCovScale = 0.);

  // d(local offsets at the vertex point) / d(vertex position) of the track, which is at the same
  // time the derivative of the local position of the mean vertex wrt its global position
  static Eigen::MatrixXd computeVertexTransformation(const Track& resTrack);
  static Eigen::MatrixXd makeZeroFieldInnerTransformation(const Eigen::MatrixXd& vtxTrans, int iTrk, int nTrk);

  // impose the prior of the mean interaction point on the vertex point of one track of a collision
  void addMeanVertexPrior(const Track& resTrack, const Eigen::MatrixXd& trans, gbl::GblPoint& vtxPoint, double covScale = 1.);

  // build and fit a single composed GBL trajectory for all tracks of one collision, with their
  // common vertex position as parameters shared by all of them
  bool buildGBLVertex(const std::vector<Track*>& contributors, std::vector<gbl::GblTrajectory>& gblTraj);

  /// counters of the loop over the vertices of a TF
  struct VertexLoopStat {
    int nVtx{0}, nVtxAcc{0}, nTrc{0}, nTrcAcc{0};
  };
  // steps of process(), see there
  void collectVertexTracks(const V2TRef& trackRef, bool useVertexConstraint, std::unordered_map<GTrackID, bool>& ambigTable, std::vector<Track>& resTracks);
  void refitTracks(std::vector<Track>& resTracks, bool useVertexConstraint);
  bool constrainWithVertex(const PVertex& vtx, int ivref, std::vector<Track>& resTracks);
  void buildVertexTrajectories(std::vector<Track>& resTracks, bool useCommonVertex, std::vector<gbl::GblTrajectory>& gblTraj, VertexLoopStat& stat);
  void writeMilleRecords(std::vector<gbl::GblTrajectory>& gblTraj);

  // build track to vertex association
  void buildT2V();

  // apply some misalignment on inner ITS3 layers
  // it can happen that a measurement is pushed outside of
  // ITS3 acceptance so false is to discard track
  bool applyMisalignment(Eigen::Vector2d& res, const FrameInfoExt& frame, const TrackD& wTrk, size_t iTrk);

  // combine the fitted rigid-body corrections with the initial alignment and write them as AlignParam vector
  void writeAlignParams(const std::map<uint32_t, double>& labelToValue) const;

  /// global (alignment + calibration) derivatives of one measured point, with their labels
  struct PointGlobals {
    std::vector<int> labels;
    Eigen::MatrixXd der;
  };
  PointGlobals buildPointGlobals(const FrameInfoExt& frame, const TrackD& wTrk) const;

  // steps of updateTimeDependentParams
  void initOnFirstTF();
  void initVertexer();
  void initMisalignment();
  void updateCalibrationSlots();
  void updateTPCCalibration(ProcessingContext& pc);
  void loadConfigMacro();
  void executeConfigMacro();

  long mTimeStamp{0}; // current TF time stamp in ms
  o2::framework::TimingInfo mTimeInfo;
  o2::alignrs::OutputEnum mOutOpt;
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDBGOut;
  std::unique_ptr<gbl::MilleBinary> mMille; // Millepede binary, open for the whole run
  std::vector<dataformats::VertexBase> mPVMC;
  std::vector<int> mT2PVMC;
  bool mIsITS3{true};
  o2::globaltracking::RecoContainer* mRecoData = nullptr;
  std::unique_ptr<steer::MCKinematicsReader> mcReader;
  o2::vertexing::PVertexer mVertexer; // used to refit the PV with the tracks refitted in the current alignment

  std::unique_ptr<DetectorPVT> mPVT; // virtual detector of the mean vertex, owner of its calibration state
  std::unique_ptr<DetectorITS> mITS;
  std::unique_ptr<DetectorTPC> mTPC;
  std::unique_ptr<DetectorTRD> mTRD;
  std::unique_ptr<DetectorTOF> mTOF;
  std::vector<Detector*> mDetectors; // all created detectors, in the order of the detector index

  std::unique_ptr<o2::gpu::GPUParam> mTPCParam; // TPC cluster error parametrization, field-dependent
  o2::tpc::VDriftHelper mTPCVDriftHelper{};     // drift calibration accounted by the correction maps
  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::unique_ptr<Volume> mHierarchy; // single tree-hierarchy of all detectors, rooted in a virtual volume
  Volume::SensorMapping mChip2Hiearchy; // global label mapping to leaves in the tree
  ProcStat mStat{};     // processing statistics
  GBLStat mGBLStat{};   // GBL construction and fit statistics
  bool mFieldOFF{false}; // set per TF
  bool mUseMC{false};
  bool mUsePVConstraint{false}; // use PV as additional constraint in a given track refit  //RSTODO: this should be a per-track decision, not global, should not be datamember
  GTrackID::mask_t mTracksSrcMask;
  std::vector<int> mTrackSources;
  int mNThreads{1};
  std::string mConfMacro{};                  // optional user macro configuring the detectors
  std::unique_ptr<TMethodCall> mUsrConfMethod; // its entry point, loaded by init
  const Params* mParams{nullptr};
  MisalignmentModel mMisalignment;
};

void AlignmentSpec::init(InitContext& ic)
{
  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);
  mNThreads = ic.options().get<int>("nthreads");
  if (mOutOpt) {
    LOG(info) << mOutOpt.pstring();
    mDBGOut = std::make_unique<o2::utils::TreeStreamRedirector>("its3_debug_alg.root", "recreate");
  }
  if (mUseMC) {
    mcReader = std::make_unique<steer::MCKinematicsReader>("collisioncontext.root");
  }
  for (int src = GTrackID::NSources; src--;) {
    if (mTracksSrcMask[src]) {
      mTrackSources.push_back(src);
    }
  }
  mConfMacro = ic.options().get<std::string>("config-macro");
  if (!mConfMacro.empty() && mConfMacro != "none") {
    loadConfigMacro();
  }
}

// Compile and load the user configuration macro, whose entry point is the function named after the
// macro file. It is executed later, once the hierarchy is built, see executeConfigMacro.
void AlignmentSpec::loadConfigMacro()
{
  if (!std::filesystem::exists(mConfMacro)) {
    LOGP(fatal, "Requested user macro {} does not exist", mConfMacro);
  }
  const std::string tmpMacro = mConfMacro + "+";
  TString cmd = gSystem->GetMakeSharedLib();
  cmd += " -O0 -g -ggdb";
  { // protect the compilation by a semaphore, to avoid clashes between the pipelined devices
    const auto semName = "align_macro_" + std::to_string(std::hash<std::string>{}(mConfMacro)).substr(0, 16);
    std::unique_ptr<boost::interprocess::named_semaphore> sem;
    try {
      sem = std::make_unique<boost::interprocess::named_semaphore>(boost::interprocess::open_or_create_t{}, semName.c_str(), 1);
    } catch (const std::exception& e) {
      LOGP(error, "Exception {} during the compilation semaphore setup of {}", e.what(), tmpMacro);
    }
    if (sem) {
      sem->wait(); // wait until we can enter (no one else there)
    }
    gSystem->SetMakeSharedLib(cmd.Data());
    const auto res = gROOT->LoadMacro(tmpMacro.c_str());
    if (sem) {
      sem->post();
      if (sem->try_wait()) { // nobody else is waiting: remove the semaphore resource
        sem->post();
        sem.reset();
        boost::interprocess::named_semaphore::remove(semName.c_str());
      }
    }
    if (res) {
      LOGP(fatal, "Failed to load the user macro {}", tmpMacro);
    }
  }
  const auto funcName = std::filesystem::path(mConfMacro).stem().native();
  mUsrConfMethod = std::make_unique<TMethodCall>();
  mUsrConfMethod->InitWithPrototype(funcName.c_str(), "std::vector<o2::alignrs::Detector*>*, int");
  if (!mUsrConfMethod->IsValid()) {
    LOGP(fatal, "Did not find {}(std::vector<o2::alignrs::Detector*>*, int) in the user macro {}", funcName, mConfMacro);
  }
}

// Hand the created detectors to the user macro, if any. Called once the hierarchy is built and
// configured, so that the macro has the last word on the setup.
void AlignmentSpec::executeConfigMacro()
{
  if (!mUsrConfMethod) {
    return;
  }
  int dummyPar = 0, ret = -1;
  auto* detectors = &mDetectors;
  const void* args[2] = {&detectors, &dummyPar};
  mUsrConfMethod->Execute(nullptr, args, 2, &ret);
  if (ret != 0) {
    LOGP(fatal, "Execution of the user config macro {} failed with {}", mConfMacro, ret);
  }
}

void AlignmentSpec::run(ProcessingContext& pc)
{
  if (mOutOpt[o2::alignrs::OutputOpt::MilleRes]) {
    updateTimeDependentParams(pc);
    const auto fitted = Volume::readMillepedeResults(mParams->milleResFile);
    Volume::writeMillepedeResults(mHierarchy.get(), fitted, mParams->milleResOutJson, mParams->misAlgJson);
    if (!mParams->algParamsOutFile.empty()) {
      writeAlignParams(fitted);
    }
  } else {
    o2::globaltracking::RecoContainer recoData;
    mRecoData = &recoData;
    mRecoData->collectData(pc, *mDataRequest);
    updateTimeDependentParams(pc);
    process();
  }
  mRecoData = nullptr;
}

void AlignmentSpec::process() // collisions
{
  for (auto* det : mDetectors) {
    det->prepareData(mRecoData);
  }
  if (mParams->usePVConstraintMinTracks > 0) {
    buildT2V(); // RSTODO for data
  }
  LOGP(info, "Starting fits with {} threads", mNThreads);
  std::vector<gbl::GblTrajectory> gblTraj;
  std::vector<Track> resTracks;
  const auto primVertices = mRecoData->getPrimaryVertices();
  const auto primVer2TRefs = mRecoData->getPrimaryVertexMatchedTrackRefs();
  std::unordered_map<GTrackID, bool> ambigTable;
  const int nvRefs = primVer2TRefs.size();
  VertexLoopStat stat;
  for (int ivref = 0; ivref < nvRefs; ivref++) {
    // the last reference holds the tracks not attached to any vertex
    const PVertex* vtx = (ivref < nvRefs - 1) ? &primVertices[ivref] : nullptr;
    bool useVertexConstraint = vtx && mParams->usePVConstraintMinTracks > 0 && vtx->getNContributors() >= mParams->usePVConstraintMinTracks;
    if (useVertexConstraint) {
      mStat.data[ProcStat::kInput][ProcStat::kVertices]++;
    }
    if (mParams->verbose > 1) {
      LOGP(info, "processing vtref {} of {} with {} tracks, {}", ivref, nvRefs, primVer2TRefs[ivref].getEntries(), vtx ? vtx->asString() : std::string{});
    }
    stat.nVtx++;
    collectVertexTracks(primVer2TRefs[ivref], useVertexConstraint, ambigTable, resTracks);
    refitTracks(resTracks, useVertexConstraint);
    if (useVertexConstraint && !constrainWithVertex(*vtx, ivref, resTracks)) {
      useVertexConstraint = false; // the tracks of this vertex are fitted w/o the vertex point
    }
    if (useVertexConstraint) {
      stat.nVtxAcc++;
    }
    buildVertexTrajectories(resTracks, useVertexConstraint && mParams->useMultyTrackPVConstraint, gblTraj, stat);
  }
  LOGP(info, "Fitted {} of {} tracks of {} of {} vertices", stat.nTrcAcc, stat.nTrc, stat.nVtxAcc, stat.nVtx);
  mGBLStat.print();
  writeMilleRecords(gblTraj);
}

// Book the tracks of one vertex reference for the refit, with an extra frame prebooked for the
// vertex point if the vertex constraint is requested.
void AlignmentSpec::collectVertexTracks(const V2TRef& trackRef, bool useVertexConstraint, std::unordered_map<GTrackID, bool>& ambigTable, std::vector<Track>& resTracks)
{
  const auto primVerGIs = mRecoData->getPrimaryVertexMatchedTracks();
  resTracks.clear();
  for (int src : mTrackSources) {
    const int start = trackRef.getFirstEntryOfSource(src), end = start + trackRef.getEntriesOfSource(src);
    for (int ti = start; ti < end; ti++) {
      const auto trackIndex = primVerGIs[ti];
      if (trackIndex.isAmbiguous()) {
        auto& ambSeen = ambigTable[trackIndex]; // an ambiguous track is accounted only with the 1st vertex it is attached to
        if (ambSeen) {
          continue;
        }
        ambSeen = true;
      }
      const auto& trPar = mRecoData->getTrackParam(trackIndex);
      if (!mFieldOFF && trPar.getPt() < mParams->minPt) {
        continue;
      }
      mStat.data[ProcStat::kInput][ProcStat::kTracks]++;
      auto& tr = resTracks.emplace_back();
      tr.gid = trackIndex;
      tr.track = convertTrack<double>(trPar);
      resetTrackCovariance(tr.track);
      tr.kfFit.chi2 = 0.f;       // the detectors accumulate the chi2 of their own points into it
      if (useVertexConstraint) { // reserve a frame for the eventual vertex point
        tr.info.emplace_back();
      }
    }
  }
}

// Let every contributing detector add its points to the tracks and refit them inward, to the
// innermost point (or the slot of the vertex point). A failed track has its gid cleared.
void AlignmentSpec::refitTracks(std::vector<Track>& resTracks, bool useVertexConstraint)
{
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(mNThreads)
#endif
  for (int itr = 0; itr < (int)resTracks.size(); itr++) {
    auto& track = resTracks[itr];
    const auto contributorsGID = mRecoData->getSingleDetectorRefs(track.gid);
    bool ok = true;
    for (auto* det : mDetectors) { // in the order of the radial position
      const int detID = det->getO2DetID();
      if (detID >= 0 && track.gid.includesDet(detID) && !det->prepareTrack(mRecoData, contributorsGID, track)) {
        ok = false; // the whole track is abandoned if any of its detectors fails
        break;
      }
    }
    if (!ok || !track.fitTrack((int)track.info.size() - 1, useVertexConstraint ? 1 : 0, false, true)) {
      track.gid.clear(); // mark as failed
    }
  }
}

// Refit the vertex with the refitted tracks and add the refitted vertex as the point prebooked in
// the info[0] slot of its contributors. Returns false if the vertex refit failed.
bool AlignmentSpec::constrainWithVertex(const PVertex& vtx, int ivref, std::vector<Track>& resTracks)
{
  PVertex vtxRefit{};
  if (!refitPV(vtx, resTracks, vtxRefit)) {
    return false;
  }
  mStat.data[ProcStat::kAccepted][ProcStat::kVertices]++;
  if (mParams->verbose > 1) {
    LOGP(info, "refitted vtref {}: {} (original: {})", ivref, vtxRefit.asString(), vtx.asString());
  }
  int nTrcWithPV = 0;
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(mNThreads) reduction(+ : nTrcWithPV)
#endif
  for (int itr = 0; itr < (int)resTracks.size(); itr++) {
    auto& track = resTracks[itr];
    if (!track.gid.isPVContributor()) {
      continue;
    }
    if (!track.updateWithVertex(vtxRefit)) {
      LOGP(debug, "Failed to update track {} with {}", track.gid.asString(), vtxRefit.asString());
      continue;
    }
    nTrcWithPV++;
  }
  mStat.data[ProcStat::kAccepted][ProcStat::kTracksWithVertex] += nTrcWithPV;
  return true;
}

// Create the GBL input of the tracks of one vertex: 1st, if the common vertex constraint is
// requested, a single composed trajectory of the tracks whose prebooked info[0] slot holds the
// vertex point. Then the remaining tracks, each separately: those carrying the vertex point are
// fitted from it on, treating it as an ordinary measured point of this track only, the others from
// their 1st measured point.
void AlignmentSpec::buildVertexTrajectories(std::vector<Track>& resTracks, bool useCommonVertex, std::vector<gbl::GblTrajectory>& gblTraj, VertexLoopStat& stat)
{
  if (useCommonVertex) {
    std::vector<Track*> contributors;
    for (auto& resTrack : resTracks) {
      if (resTrack.gid.isIndexSet() && !resTrack.info.empty() && resTrack.info.front().isVertex()) {
        contributors.push_back(&resTrack);
      }
    }
    stat.nTrc += (int)contributors.size();
    if (buildGBLVertex(contributors, gblTraj)) { // the collision is a single Millepede local fit object
      stat.nTrcAcc += (int)contributors.size();
    }
  }
  // per-track PV constraint: every track carrying the vertex point gets the mean vertex prior,
  // optionally deweighted by their number for the prior to count once per collision
  double mvPriorCovScale = 0.;
  if (!useCommonVertex) {
    const auto nWithVertex = std::count_if(resTracks.begin(), resTracks.end(), [](const Track& t) { return t.gid.isIndexSet() && !t.info.empty() && t.info.front().isVertex(); });
    if (nWithVertex > 0) {
      mvPriorCovScale = mParams->scaleMVPriorWithNTracks ? static_cast<double>(nWithVertex) : 1.;
    }
  }
  for (auto& resTrack : resTracks) {
    if (!resTrack.gid.isIndexSet() || resTrack.info.empty()) {
      continue; // failed track
    }
    if (useCommonVertex && resTrack.info.front().isVertex()) {
      continue; // already accounted with the common vertex constraint
    }
    stat.nTrc++;
    const bool hasVertex = resTrack.info.front().isVertex();
    if (buildGBLTrack(resTrack, resTrack.info.front().isValid() ? 0 : 1, gblTraj, hasVertex ? mvPriorCovScale : 0.)) {
      stat.nTrcAcc++;
    }
  }
}

// Append the trajectories of this TF to the Millepede binary, which is kept open for the whole run
void AlignmentSpec::writeMilleRecords(std::vector<gbl::GblTrajectory>& gblTraj)
{
  if (!mOutOpt[o2::alignrs::OutputOpt::MilleData]) {
    return;
  }
  if (!mMille) {
    mMille = std::make_unique<gbl::MilleBinary>(mParams->milleBinFile, true);
  }
  for (auto& traj : gblTraj) {
    traj.milleOut(*mMille);
  }
  LOGP(info, "Wrote {} trajectories to {}", gblTraj.size(), mParams->milleBinFile);
}

void AlignmentSpec::updateTimeDependentParams(ProcessingContext& pc)
{
  mTimeInfo = pc.services().get<o2::framework::TimingInfo>();
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  mFieldOFF = std::abs(o2::base::PropagatorD::Instance()->getNominalBz()) < 0.1;
  mTimeStamp = (o2::base::GRPGeomHelper::instance().getOrbitResetTimeMUS() +  static_cast<long>(mTimeInfo.firstTForbit * o2::constants::lhc::LHCOrbitMUS)) * 1e-3;
  if (static bool initOnce{false}; !initOnce) {
    initOnce = true;
    initOnFirstTF();
  }
  updateCalibrationSlots();
  if (!mOutOpt[o2::alignrs::OutputOpt::MilleRes] && Params::Instance().usePVConstraintMinTracks > 0) {
    pc.inputs().get<o2::dataformats::MeanVertexObject*>("meanvtx"); // triggers finaliseCCDB
  }
  if (mTRD && !mOutOpt[o2::alignrs::OutputOpt::MilleRes]) {
    pc.inputs().get<o2::trd::CalVdriftExB*>("calvdexb"); // triggers finaliseCCDB
  }
  if (mTPC && !mOutOpt[o2::alignrs::OutputOpt::MilleRes]) {;
    updateTPCCalibration(pc); // must precede initOnFirstTF: the drift calibration DOFs need the maps
  }
}

void AlignmentSpec::updateTPCCalibration(ProcessingContext& pc)
{
  mTPCVDriftHelper.extractCCDBInputs(pc);
  const auto* raw = pc.inputs().get<const char*>("corrMap");
  mTPC->setCorrMaps(&o2::gpu::TPCFastTransformPOD::get(raw));
  if (mTPCVDriftHelper.isUpdated()) {
    const auto& vd = mTPCVDriftHelper.getVDriftObject();
    LOGP(info, "TPC VDrift factor {} wrt reference {} and time offset correction {} wrt {} from source {}",
         vd.corrFact, vd.refVDrift, vd.timeOffsetCorr, vd.refTimeOffset, mTPCVDriftHelper.getSourceName());
    mTPC->setVDrift(vd);
    mTPCVDriftHelper.acknowledgeUpdate();
  }
}

// One-time initialisation, deferred to the 1st TF since it needs the geometry
void AlignmentSpec::initOnFirstTF()
{
  mParams = &Params::Instance();
  mParams->printKeyValues(true, true);
  if (mITS) {
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G, o2::math_utils::TransformType::T2G));
  }
  buildHierarchy();
  if (mTPC) {
    mTPCParam = std::make_unique<o2::gpu::GPUParam>();
    mTPCParam->SetDefaults(o2::base::PropagatorD::Instance()->getNominalBz(), false);
    mTPC->setTPCParam(mTPCParam.get());
    mTPC->initCalib();
  }
  if (mTRD) {
    mTRD->initCalib();
  }
  initVertexer();
  initMisalignment();
}

void AlignmentSpec::initVertexer()
{
  if (mParams->usePVConstraintMinTracks <= 0) {
    return;
  }
  o2::conf::ConfigurableParam::updateFromString("pvertexer.useTimeInChi2=false;"); // the PV refit does not use the track time
  mVertexer.init();
  mVertexer.setMeanVertex(&mPVT->getMeanVertexPrior());
  mVertexer.initMeanVertexConstraint();
}

void AlignmentSpec::initMisalignment()
{
  if (!(mParams->doMisalignmentLeg || mParams->doMisalignmentRB || mParams->doMisalignmentInex)) {
    return;
  }
  mMisalignment = loadMisalignmentModel(mParams->misAlgJson); // empty model for an empty path
}

// Assign every detector with a time-sliced calibration to the slot covering this TF. The detectors
// re-point their own slot-dependent labels and priors; only the mean vertex has an effect outside
// of its own state, the prior of the PV refit having to be re-imposed on the vertexer.
void AlignmentSpec::updateCalibrationSlots()
{
  for (auto* det : mDetectors) {
    if (!det->setTimeStamp(mTimeStamp)) { // nothing changed for this detector
      continue;
    }
    if (det == mPVT.get()) {
      // copies the object and re-inits the XY constraint, no init() needed
      mVertexer.setMeanVertex(&mPVT->getMeanVertexPrior());
      mVertexer.initMeanVertexConstraint();
    }
  }
}

// Derivatives of the prediction of one measured point (i.e. minus those of the residual, as Millepede
// expects) wrt all the global parameters it depends on:
// the rigid-body DOFs of the measurement leaf and of all its ancestors up to (excluding) the common
// root, followed by the calibration DOFs of the direct parent of the leaf.
// The base derivative is computed in the tracking frame while the alignment is done in the local
// one, dr/da_(LOC) = dr/da_(TRK) * da_(TRK)/da_(LOC), and is then transported level by level with
// the local-to-parent jacobian of each child.
AlignmentSpec::PointGlobals AlignmentSpec::buildPointGlobals(const FrameInfoExt& frame, const TrackD& wTrk) const
{
  const auto volIt = mChip2Hiearchy.find(frame.label); // as assigned by the detector owning the frame
  if (volIt == mChip2Hiearchy.end()) {
    LOGP(fatal, "Cannot find global label: {}", frame.label.asString());
  }
  const auto* tileVol = volIt->second;
  const auto derCtx = makeDerivativeContext(frame, wTrk);

  // count rigid body columns: only volumes with real DOFs (not DOFPseudo).
  // The chain walks up to the common root of all detectors, which owns no DOFs: the top
  // volume of the detector is therefore included, its DOFs being free or fixed as configured.
  int nColRB{0};
  for (const auto* v = tileVol; v && !v->isRoot(); v = v->getParent()) {
    if (v->getRigidBody()) {
      nColRB += v->getRigidBody()->nDOFs();
    }
  }

  // count calibration columns
  const auto* sensorVol = tileVol->getParent();
  const auto* calibSet = sensorVol ? sensorVol->getCalib() : nullptr;
  const int nCalib = calibSet ? calibSet->nDOFs() : 0;

  PointGlobals globals{{}, Eigen::MatrixXd::Zero(2, nColRB + nCalib)};
  globals.labels.reserve(nColRB + nCalib);
  Eigen::Index curCol{0};

  // 1) tile: TRK -> LOC via precomputed T2L and J_L2T
  Matrix26 der = getLeafRigidBodyDerivatives(*tileVol, frame, derCtx);
  if (tileVol->getRigidBody()) {
    const int nd = tileVol->getRigidBody()->nDOFs();
    for (int iDOF = 0; iDOF < nd; ++iDOF) {
      globals.labels.push_back(tileVol->getLabel().rawGBL(iDOF));
    }
    globals.der.middleCols(curCol, nd) = der;
    curCol += nd;
  }

  // 2) chain through parents: child's J_L2P
  for (const auto* child = tileVol; child->getParent() && !child->getParent()->isRoot(); child = child->getParent()) {
    der *= child->getJL2P();
    const auto* parent = child->getParent();
    if (parent->getRigidBody()) {
      const int nd = parent->getRigidBody()->nDOFs();
      for (int iDOF = 0; iDOF < nd; ++iDOF) {
        globals.labels.push_back(parent->getLabel().rawGBL(iDOF));
      }
      globals.der.middleCols(curCol, nd) = der;
      curCol += nd;
    }
  }

  // 3) calibration derivatives (apply directly on the whole sensor, not on individual tiles).
  // The label is that of the calibration time slot of the processed TF, the DOFs of every slot
  // being independent parameters of the fit.
  if (calibSet) {
    const int nd = calibSet->nDOFs();
    Eigen::MatrixXd calDer(2, nd);
    calibSet->fillDerivatives(derCtx, calDer);
    const auto& calibLbl = sensorVol->getActiveCalibLabel();
    for (int iDOF = 0; iDOF < nd; ++iDOF) {
      globals.labels.push_back(calibLbl.rawGBL(iDOF));
    }
    globals.der.middleCols(curCol, nd) = calDer;
    curCol += nd;
  }
  if (frame.zFromTrack) { // the Z "measurement" follows the track, whatever the alignment
    globals.der.row(1).setZero();
  }
  return globals;
}

void AlignmentSpec::buildHierarchy()
{
  // all detectors share a single hierarchy: the top volume of every detector is a branch of a
  // fictitious root. The root owns no rigid-body DOFs, hence the branches, which may belong to
  // detectors as unrelated as the mean vertex and the TPC, are not mutually constrained.
  mHierarchy = Volume::makeRoot();
  for (auto* det : mDetectors) {
    if (det->getDetIdx() == Detector::DetPVT && mParams->usePVConstraintMinTracks <= 0) {
      continue; // w/o the PV constraint the mean vertex is neither used nor aligned
    }
    // the calibration slots must be known before the hierarchy is built: a detector calibrated per
    // time slot prepares the labels of all its slots on the volume owning the DOFs
    det->loadTimeSlots();
    det->attachTo(mHierarchy.get(), mChip2Hiearchy);
  }
  const bool withPVT = mPVT->getTopVolume() != nullptr;
  if (withPVT) {
    LOGP(info, "Mean vertex prior: {}", mPVT->getMeanVertexCCDB().asString());
  }

  if (!mParams->dofConfigJson.empty()) {
    // a rule may e.g. fix the mean vertex position (calib type meanvertex), as well as fix or free the
    // top volume of any detector: no DOF is fixed or freed by the hierarchy construction itself
    Volume::applyDOFConfig(mHierarchy.get(), mParams->dofConfigJson);
  }
  executeConfigMacro(); // the user macro overrides whatever the JSON configuration has set

  mHierarchy->finalise();
  if (withPVT) {
    mPVT->updatePositionLabels();
    if (mPVT->getPositionLabels().empty()) {
      LOGP(info, "Mean vertex position is fixed, it is imposed as a prior w/o being aligned");
    }
  }
  if (mOutOpt[o2::alignrs::OutputOpt::MilleSteer]) {
    std::ofstream tree(mParams->milleTreeFile);
    mHierarchy->writeTree(tree);
    std::ofstream cons(mParams->milleConFile);
    mHierarchy->writeRigidBodyConstraints(cons);
    std::ofstream par(mParams->milleParamFile);
    mHierarchy->writeParameters(par);
  }
}

bool AlignmentSpec::getTransportJacobian(const TrackD& track, double xTo, double alphaTo, gbl::Matrix5d& jac, gbl::Matrix5d& err)
{
  auto prop = o2::base::PropagatorD::Instance();
  const auto minStep = std::sqrt(std::numeric_limits<double>::epsilon());
  const gbl::Vector5d x0(track.getParams());

  auto propagate = [&](gbl::Vector5d& p) -> bool {
    TrackD tmp(track);
    for (int i{0}; i < track::kNParams; ++i) {
      tmp.setParam(p[i], i);
    }
    // the reference is propagated together with the track, hence it must restart from the
    // unperturbed state at every call
    o2::track::TrackParD ref(track), *refLin = mParams->useStableRef ? &ref : nullptr;
    if (!prop->propagateToAlphaX(tmp, refLin, alphaTo, xTo, false, mParams->maxSnp, mParams->maxStep, 1, mParams->corrType)) {
      return false;
    }
    p = gbl::Vector5d(tmp.getParams());
    return true;
  };

  for (int iPar{0}; iPar < track::kNParams; ++iPar) {
    // step size
    double h = std::min(mParams->ridderMaxIniStep[iPar], std::max(minStep, std::abs(track.getParam(iPar)) * mParams->ridderRelIniStep[iPar]) * std::pow(mParams->ridderShrinkFac, mParams->ridderMaxExtrap / 2));
    ;
    // romberg tableu
    Eigen::MatrixXd cur(track::kNParams, mParams->ridderMaxExtrap);
    Eigen::MatrixXd pre(track::kNParams, mParams->ridderMaxExtrap);
    double normErr = std::numeric_limits<double>::max();
    gbl::Vector5d bestDeriv = gbl::Vector5d::Constant(std::numeric_limits<double>::max());
    for (int iExt{0}; iExt < mParams->ridderMaxExtrap; ++iExt) {
      gbl::Vector5d xPlus = x0, xMinus = x0;
      xPlus(iPar) += h;
      xMinus(iPar) -= h;
      if (!propagate(xPlus) || !propagate(xMinus)) {
        return false;
      }
      cur.col(0) = (xPlus - xMinus) / (2.0 * h);
      if (!iExt) {
        bestDeriv = cur.col(0);
      }
      // shrink step in next iteration
      h /= mParams->ridderShrinkFac;
      // richardson extrapolation
      double fac = mParams->ridderShrinkFac * mParams->ridderShrinkFac;
      for (int k{1}; k <= iExt; ++k) {
        cur.col(k) = (fac * cur.col(k - 1) - pre.col(k - 1)) / (fac - 1.0);
        fac *= mParams->ridderShrinkFac * mParams->ridderShrinkFac;
        double e = std::max((cur.col(k) - cur.col(k - 1)).norm(), (cur.col(k) - pre.col(k - 1)).norm());
        if (e <= normErr) {
          normErr = e;
          bestDeriv = cur.col(k);
          if (normErr < mParams->ridderEps) {
            break;
          }
        }
      }
      if (normErr < mParams->ridderEps) {
        break;
      }
      // check stability
      if (iExt > 0) {
        double tableauErr = (cur.col(iExt) - pre.col(iExt - 1)).norm();
        if (tableauErr >= 2.0 * normErr) {
          break;
        }
      }
      std::swap(cur, pre);
    }
    if (bestDeriv.isApproxToConstant(std::numeric_limits<double>::max())) {
      return false;
    }
    jac.col(iPar) = bestDeriv;
    err.col(iPar) = gbl::Vector5d::Constant(normErr);
  }

  if (jac.isIdentity(1e-8)) {
    LOGP(error, "Near jacobian idendity for taking track from {} to {}", track.getX(), xTo);
    return false;
  }

  return true;
}

bool AlignmentSpec::refitPV(const PVertex& vtxOrig, const std::vector<Track>& resTracks, PVertex& vtxRefit)
{
  // Refit the vertex with the tracks refitted in the current alignment. Only the successfully refitted
  // contributors of this vertex participate: a track whose refit failed has its gid cleared, which also
  // resets the PVContributor flag. Hence the refitted vertex may differ from the original one not only
  // because of the alignment but also because of the reduced number of contributors.
  std::vector<o2::track::TrackParCov> tracks;
  tracks.reserve(resTracks.size());
  for (const auto& resTrack : resTracks) {
    if (resTrack.gid.isIndexSet() && resTrack.gid.isPVContributor()) {
      tracks.push_back(convertTrack<float>(resTrack.track)); // the track is at its innermost update point
    }
  }
  if (static_cast<int>(tracks.size()) < mParams->usePVConstraintMinTracks) {
    LOGP(debug, "Abandon the refit of {}: only {} of {} contributors were refitted", vtxOrig.asString(), tracks.size(), vtxOrig.getNContributors());
    return false;
  }
  if (!mVertexer.prepareVertexRefit(tracks, vtxOrig)) {
    LOGP(debug, "Failed to prepare the refit of {} with {} tracks", vtxOrig.asString(), tracks.size());
    return false;
  }
  // the pool contains only the tracks we want to use, hence the empty useTrack mask;
  // refitVertexFull (rather than refitVertex) is used since the alignment may shift the vertex
  vtxRefit = mVertexer.refitVertexFull({}, vtxOrig);
  if (vtxRefit.getChi2() < 0.f) {
    LOGP(debug, "Failed to refit {} with {} tracks", vtxOrig.asString(), tracks.size());
    return false;
  }
  return true;
}

// Fill the GBL points of the track frames from the ipStart slot outward (resTrack.track must be the
// state at the ipStart frame). The measurement residuals are stored in resTrack.points.
// skipFirstMeas: add no measurement on the 1st accounted point, used for the common vertex point of
// a composed trajectory, whose position enters via the inner transformation instead.
bool AlignmentSpec::fillGBLPoints(Track& resTrack, int ipStart, bool skipFirstMeas, std::vector<gbl::GblPoint>& points, double mvPriorCovScale)
{
  auto prop = o2::base::PropagatorD::Instance();
  const int np = (int)resTrack.info.size();
  auto wTrk = resTrack.track; // working copy, the seed must be preserved
  o2::track::TrackParD trkRef, *refLin = nullptr;
  if (mParams->useStableRef) {
    refLin = &(trkRef = wTrk);
  }
  // the MC index of the ITS part, needed only by the MC-based misalignment simulation
  const auto itsGID = mRecoData->getITSContributorGID(resTrack.gid);
  const size_t iTrk = itsGID.isIndexSet() ? itsGID.getIndex() : 0;
  points.clear();
  points.reserve(np - ipStart);
  resTrack.points.clear();
  resTrack.points.reserve(np - ipStart);
  track::TrackLTIntegral lt;
  lt.setTimeNotNeeded();
  for (int ip = ipStart; ip < np; ++ip) {
    const auto& frame = resTrack.info[ip];
    if (!frame.isValid()) {
      continue;
    }
    gbl::Matrix5d err = gbl::Matrix5d::Identity(), jacALICE = gbl::Matrix5d::Identity(), jacGBL = gbl::Matrix5d::Identity();
    float msErr = 0.f;
    if (!points.empty()) { // not the 1st accounted point: step to it
      // numerically calculates the transport jacobian from prev. point to this point
      // then we actually do the step to the point and accumulate the material
      if (!getTransportJacobian(wTrk, frame.x, frame.alpha, jacALICE, err) ||
          !prop->propagateToAlphaX(wTrk, refLin, frame.alpha, frame.x, false, mParams->maxSnp, mParams->maxStep, 1, mParams->corrType, &lt)) {
        ++mGBLStat.failedProp;
        return false;
      }
      msErr = its::math_utils::MSangle(wTrk.getPID().getMass(), wTrk.getP(), lt.getX2X0());
      jacGBL = toGBLOrder(jacALICE);
    }

    // wTrk is now in the measurment frame
    gbl::GblPoint point(jacGBL);
    // measurement
    const auto& cluster = frame.cluster;
    Eigen::Vector2d res = Eigen::Vector2d::Zero();
    Eigen::Matrix2d prec = Eigen::Matrix2d::Zero();
    const bool addMeas = !(skipFirstMeas && points.empty());
    if (addMeas) {
      res << cluster.getY() - wTrk.getY(), cluster.getZ() - wTrk.getZ();

      // here we can apply some misalignment on the measurment
      if (!applyMisalignment(res, frame, wTrk, iTrk)) {
        return false;
      }
      if (!getMeasurementPrecision(cluster, prec)) {
        LOGP(warn, "Abandon track {}: non-positive covariance of the point {}", resTrack.gid.asString(), frame.asString());
        return false;
      }
      // measured directly in the tracking frame (identity projection); GBL diagonalizes the precision
      // matrix and transforms the residuals and the global derivatives accordingly
      point.addMeasurement(res, prec);
      auto& meas = resTrack.points.emplace_back();
      meas.dy = res[0];
      meas.dz = res[1];
      meas.sig2y = cluster.getSigmaY2();
      meas.sig2z = cluster.getSigmaZ2();
      meas.z = wTrk.getZ();
      meas.phi = wTrk.getPhi();
      o2::math_utils::bringTo02Pid(meas.phi);
    }
    if (msErr > mParams->minMS && ip < np - 1) {
      // the material crossed since the previous point is lumped into a thin scatterer at this one
      point.addScatterer(Eigen::Vector2d::Zero(), getScatteringPrecision(wTrk, msErr));
      lt.clearFast(); // clear if accounted
    }

    if (!frame.isVertex()) { // the vertex point has no alignable volume behind it
      const auto globals = buildPointGlobals(frame, wTrk);
      point.addGlobals(globals.labels, globals.der);
    } else if (addMeas && mvPriorCovScale > 0.) {
      // Per-track PV constraint: the track passes through the true vertex V, measured (i) by the
      // refitted vertex with its fit covariance, added above w/o global derivatives: it does not move
      // with the mean vertex, and (ii) by the mean vertex mu with the luminous region covariance,
      // V ~ N(mu, Sigma_lumi), which carries the derivatives wrt mu. Integrating V out, the
      // information on mu comes from (V_refit - mu) ~ N(0, Sigma_fit + Sigma_lumi).
      addMeanVertexPrior(resTrack, computeVertexTransformation(resTrack), point, mvPriorCovScale);
    }

    if (mOutOpt[o2::alignrs::OutputOpt::VerboseGBL]) {
      static Eigen::IOFormat fmt(4, 0, ", ", "\n", "[", "]");
      LOGP(info, "WORKING-POINT {}", ip);
      LOGP(info, "Track: {}", wTrk.asString());
      LOGP(info, "FrameInfo: {}", frame.asString());
      std::cout << "jacALICE:\n"
                << jacALICE.format(fmt) << '\n';
      std::cout << "jacGBL:\n"
                << jacGBL.format(fmt) << '\n';
      LOGP(info, "residual: dy={} dz={}", res[0], res[1]);
      LOGP(info, "precision: precYY={} precYZ={} precZZ={}", prec(0, 0), prec(0, 1), prec(1, 1));
      point.printPoint(5);
    }
    points.push_back(point);
    if (refLin) { // displace the reference to the measurement, as done in the KF refit
      refLin->setY(cluster.getY());
      refLin->setZ(cluster.getZ());
    }
  }
  return true;
}

// Fit the constructed trajectory, account the statistics and store it for the Mille output.
// kfChi2Ndf (if >0) is the chi2/ndf of the KF refit of the same data, used for diagnostics only.
// On success the fit result is put to fitOut.
bool AlignmentSpec::fitGBLTrajectory(gbl::GblTrajectory& traj, float kfChi2Ndf, std::vector<gbl::GblTrajectory>& gblTraj, FitInfo& fitOut)
{
  if (!traj.isValid()) {
    ++mGBLStat.construct;
    return false;
  }
  double chi2 = NAN, lostWeight = NAN;
  int ndf = 0;
  if (auto ierr = traj.fit(chi2, ndf, lostWeight); ierr) {
    ++mGBLStat.fitFail;
    return false;
  }
  if (mOutOpt[o2::alignrs::OutputOpt::VerboseGBL]) {
    LOGP(info, "GBL FIT chi2 {} ndf {}", chi2, ndf);
    traj.printTrajectory(5);
  }
  if (!ndf || chi2 / ndf > mParams->maxChi2Ndf) {
    if (mGBLStat.chi2Rej++ < 10) {
      LOGP(error, "GBL fit exceeded red chi2 {} (ndf {})", ndf ? chi2 / ndf : -1., ndf);
      if (kfChi2Ndf > 0 && std::abs(kfChi2Ndf - 1) < 0.02) {
        LOGP(error, "\tGBL is far away from good KF fit!!!!");
      }
    }
    return false;
  }
  ++mGBLStat.fit;
  mGBLStat.chi2Sum += chi2;
  mGBLStat.lostWeightSum += lostWeight;
  mGBLStat.ndfSum += ndf;
  if (mOutOpt[o2::alignrs::OutputOpt::MilleData]) {
    gblTraj.push_back(traj);
  }
  fitOut = {.chi2Ndf = (float)(chi2 / ndf), .chi2 = (float)chi2, .ndf = ndf};
  return true;
}

// Build and fit the GBL trajectory of a single refitted track, accounting its frames from the
// ipStart slot outward. The GBL fit result is stored in resTrack.gblFit.
bool AlignmentSpec::buildGBLTrack(Track& resTrack, int ipStart, std::vector<gbl::GblTrajectory>& gblTraj, double mvPriorCovScale)
{
  std::vector<gbl::GblPoint> points;
  if (!fillGBLPoints(resTrack, ipStart, false, points, mvPriorCovScale)) {
    return false;
  }
  gbl::GblTrajectory traj(points, !mFieldOFF); // no curvature w/o field
  return fitGBLTrajectory(traj, resTrack.kfFit.chi2Ndf, gblTraj, resTrack.gblFit);
}

// d(offsets at the vertex point) / d(vertex position): the vertex is displaced in the global frame
// at fixed track direction and momentum, hence the offsets of the track in the tracking frame of
// the vertex point change by the displacement rotated to this frame, corrected for the shift of the
// reference X along the track. The track state must be the one at the vertex point.
Eigen::MatrixXd AlignmentSpec::computeVertexTransformation(const Track& resTrack)
{
  double ca{0}, sa{0};
  o2::math_utils::sincosd(resTrack.info.front().alpha, sa, ca);
  const auto slopes = TrackSlopes::computeTrackSlopes(resTrack.track.getSnp(), resTrack.track.getTgl());
  Eigen::MatrixXd trans(2, 3);
  trans(0, 0) = -sa - slopes.dydx * ca; // dY / dVx
  trans(0, 1) = ca - slopes.dydx * sa;  // dY / dVy
  trans(0, 2) = 0.;                     // dY / dVz
  trans(1, 0) = -slopes.dzdx * ca;      // dZ / dVx
  trans(1, 1) = -slopes.dzdx * sa;      // dZ / dVy
  trans(1, 2) = 1.;                     // dZ / dVz
  return trans;
}

// Inner transformation of the iTrk-th of nTrk tracks of a composed trajectory w/o magnetic field. The
// geometric constraint (2 rows) would add a curvature per track, which is not defined at B=0, hence
// the kinematic constraint (5 rows, the full local parameters q/pt, snp, tgl, Y, Z at the vertex point)
// is used instead, with the external parameters being the vertex position followed by the (snp, tgl)
// of every track. The q/pt row is null, i.e. the curvature is fixed. vtxTrans is the
// d(Y,Z)/d(vertex) transformation from computeVertexTransformation.
Eigen::MatrixXd AlignmentSpec::makeZeroFieldInnerTransformation(const Eigen::MatrixXd& vtxTrans, int iTrk, int nTrk)
{
  Eigen::MatrixXd trans = Eigen::MatrixXd::Zero(5, 3 + 2 * nTrk);
  trans(1, 3 + 2 * iTrk) = 1.;     // snp of this track
  trans(2, 3 + 2 * iTrk + 1) = 1.; // tgl of this track
  trans.block(3, 0, 2, 3) = vtxTrans;
  return trans;
}

// Impose the prior on the vertex of the collision: it must agree with the mean interaction point
// within the sigmas of the luminous region. Being a property of the collision, the prior is
// accounted only once, as a measurement on the vertex point of a single track of the composed
// trajectory (adding it per track would multiply its weight by their number). In the per-track PV
// constraint mode it is added to every track, with the covariance scaled by covScale to compensate.
// The measurement carries the derivatives wrt the position of the mean vertex, which is a global
// alignment parameter (the dummy volume of the virtual PVT detector), unless it is kept fixed.
// trans is the d(local offsets)/d(vertex position) transformation of the hosting track, which is at
// the same time the derivative of the local position of the mean vertex wrt its global position.
void AlignmentSpec::addMeanVertexPrior(const Track& resTrack, const Eigen::MatrixXd& trans, gbl::GblPoint& vtxPoint, double covScale)
{
  const auto& frame = resTrack.info.front();
  double ca{0}, sa{0};
  o2::math_utils::sincosd(frame.alpha, sa, ca);
  const auto slopes = TrackSlopes::computeTrackSlopes(resTrack.track.getSnp(), resTrack.track.getTgl());
  const auto& mvPrior = mPVT->getMeanVertexPrior(); // frozen for the whole calibration slot
  // the mean vertex in the tracking frame of the vertex point, brought to the plane of this point
  const double muX = mvPrior.getX() * ca + mvPrior.getY() * sa; // along the local X
  const double dX = muX - frame.x;
  const double muY = -mvPrior.getX() * sa + mvPrior.getY() * ca - slopes.dydx * dX;
  const double muZ = mvPrior.getZ() - slopes.dzdx * dX;
  Eigen::Vector2d res;
  res << muY - resTrack.track.getY(), muZ - resTrack.track.getZ();
  // covariance of the luminous region rotated to this frame: the transformation of the vertex
  // position to the local offsets is the same as for the common parameters of the trajectory
  Eigen::Matrix3d covGlo;
  covGlo << mvPrior.getSigmaX2(), mvPrior.getSigmaXY(), mvPrior.getSigmaXZ(),
    mvPrior.getSigmaXY(), mvPrior.getSigmaY2(), mvPrior.getSigmaYZ(),
    mvPrior.getSigmaXZ(), mvPrior.getSigmaYZ(), mvPrior.getSigmaZ2();
  const Eigen::Matrix2d cov = covScale * trans * covGlo * trans.transpose();
  if (cov.determinant() < 1e-16) {
    LOGP(warn, "Skipping the mean vertex prior: singular projected covariance of {}", mvPrior.asString());
    return;
  }
  vtxPoint.addMeasurement(res, Eigen::Matrix2d(cov.inverse()));
  if (!mPVT->getPositionLabels().empty()) { // the mean vertex position is aligned as well
    // Millepede expects the derivatives of the prediction, i.e. minus those of the residual. Here it
    // is the measurement, the mean vertex, which moves with the parameter: d(res)/d(MV) = +trans.
    vtxPoint.addGlobals(mPVT->getPositionLabels(), Eigen::MatrixXd(-trans));
  }
}

// Build and fit a single composed GBL trajectory for all tracks of one collision, with the 3
// coordinates of their common vertex as parameters shared by all of them: the vertex constraint is
// then exact by parameterization, while every track keeps its own curvature and scattering
// parameters. Follows the GBL composed trajectory with geometric constraint (GBL exampleComposedGeo).
// W/o magnetic field the kinematic constraint w/o curvature is used, see makeZeroFieldInnerTransformation.
// Every contributor must carry the vertex point in its info[0] slot, with its track state there,
// as left by Track::updateWithVertex.
bool AlignmentSpec::buildGBLVertex(const std::vector<Track*>& contributors, std::vector<gbl::GblTrajectory>& gblTraj)
{
  std::vector<std::pair<std::vector<gbl::GblPoint>, Eigen::MatrixXd>> pointsAndTrans;
  std::vector<Track*> used;
  pointsAndTrans.reserve(contributors.size());
  used.reserve(contributors.size());
  for (auto* trc : contributors) {
    std::vector<gbl::GblPoint> points;
    // the vertex point is the 1st point of every sub-trajectory and carries no measurement of its
    // own: the common vertex position enters via the inner transformation below
    if (!fillGBLPoints(*trc, 0, true, points) || points.size() < 2) {
      continue;
    }
    Eigen::MatrixXd innerTrans = computeVertexTransformation(*trc);
    if (used.empty()) { // impose the prior on the vertex of this collision
      addMeanVertexPrior(*trc, innerTrans, points.front());
    }
    pointsAndTrans.emplace_back(std::move(points), std::move(innerTrans));
    used.push_back(trc);
  }
  if (used.size() < 2) { // a single track does not define the common vertex
    return false;
  }
  if (mFieldOFF) { // the number of the external parameters depends on the number of tracks
    for (size_t i = 0; i < pointsAndTrans.size(); i++) {
      pointsAndTrans[i].second = makeZeroFieldInnerTransformation(pointsAndTrans[i].second, (int)i, (int)pointsAndTrans.size());
    }
  }
  gbl::GblTrajectory traj(pointsAndTrans);
  FitInfo fit{};
  bool res = fitGBLTrajectory(traj, -1.f, gblTraj, fit);
  if (mVerbose > 1) {
    LOGP(info, "GBL vertex fit of {} tracks (out of {}): chi2Ndf={} chi2={} ndf={} -> {}", 
      used.size(), contributors.size(), fit.chi2Ndf, fit.chi2, fit.ndf, res ? "success" : "failure");
  }
  if (!res) {
    return false;
  }
  for (auto* trc : used) { // the composed fit is common for all the tracks of the collision
    trc->gblFit = fit;
  }
  return true;
}

void AlignmentSpec::buildT2V()
{
  const auto& itsTracks = mRecoData->getITSTracks();
  mT2PVMC.clear();
  mT2PVMC.resize(itsTracks.size(), -1);
  if (mUseMC) {
    mPVMC.reserve(mcReader->getNEvents(0));
    for (int iEve{0}; iEve < mcReader->getNEvents(0); ++iEve) {
      const auto& eve = mcReader->getMCEventHeader(0, iEve);
      dataformats::VertexBase vtx;
      constexpr float err{22e-4f};
      vtx.setX((float)eve.GetX());
      vtx.setY((float)eve.GetY());
      vtx.setZ((float)eve.GetZ());
      vtx.setSigmaX(err);
      vtx.setSigmaY(err);
      vtx.setSigmaZ(err);
      mPVMC.push_back(vtx);
    }
    const auto& mcLbls = mRecoData->getITSTracksMCLabels();
    for (size_t iTrk{0}; iTrk < mcLbls.size(); ++iTrk) {
      const auto& lbl = mcLbls[iTrk];
      if (!lbl.isValid() || !lbl.isCorrect()) {
        continue;
      }
      const auto& mcTrk = mcReader->getTrack(lbl);
      if (mcTrk->isPrimary()) {
        mT2PVMC[iTrk] = lbl.getEventID();
      }
    }
  } else {
    LOGP(warn, "Data PV to track TODO");
  }
}

bool AlignmentSpec::applyMisalignment(Eigen::Vector2d& res, const FrameInfoExt& frame, const TrackD& wTrk, size_t iTrk)
{
  if (!o2::its3::constants::detID::isDetITS3(frame.cluster.getSensorID())) {
    return true;
  }

  const int sensorID = o2::its3::constants::detID::getSensorID(frame.cluster.getSensorID());
  const int layerID = o2::its3::constants::detID::getDetID2Layer(frame.cluster.getSensorID());
  const MisalignmentFrame misFrame{
    .sensorID = sensorID,
    .layerID = layerID,
    .x = frame.x,
    .alpha = frame.alpha,
    .z = frame.cluster.getZ()};

  // --- Legendre deformation (non-rigid-body) ---
  if (mParams->doMisalignmentLeg && mIsITS3 && mUseMC) {
    const auto prop = o2::base::PropagatorD::Instance();

    const auto lbl = mRecoData->getITSTracksMCLabels()[iTrk];
    if (lbl.isFake()) {
      return false;
    }
    const auto mcTrk = mcReader->getTrack(lbl);
    if (!mcTrk) {
      return false;
    }
    std::array<double, 3> xyz{mcTrk->GetStartVertexCoordinatesX(), mcTrk->GetStartVertexCoordinatesY(), mcTrk->GetStartVertexCoordinatesZ()};
    std::array<double, 3> pxyz{mcTrk->GetStartVertexMomentumX(), mcTrk->GetStartVertexMomentumY(), mcTrk->GetStartVertexMomentumZ()};
    TParticlePDG* pPDG = TDatabasePDG::Instance()->GetParticle(mcTrk->GetPdgCode());
    if (!pPDG) {
      return false;
    }
    o2::track::TrackParD mcPar(xyz, pxyz, TMath::Nint(pPDG->Charge() / 3), false);

    auto mcAtCl = mcPar;
    if (!mcAtCl.rotate(frame.alpha) || !prop->PropagateToXBxByBz(mcAtCl, frame.x)) {
      return false;
    }

    const auto shift = evaluateLegendreShift(mMisalignment[sensorID], misFrame, TrackSlopes::computeTrackSlopes(mcAtCl.getSnp(), mcAtCl.getTgl()));
    if (!shift.accepted) {
      return false;
    }

    res[0] += shift.dy;
    res[1] += shift.dz;
  }

  // --- Rigid body misalignment ---
  // Must use the same derivative chain as GBL:
  //   dres/da_parent = dres/da_TRK * J_L2T_tile * J_L2P_tile
  // The tile is a pseudo-volume; Millepede fits at the halfBarrel (parent) level.
  if (mParams->doMisalignmentRB) {
    // use the label assigned by the detector owning the frame: it carries the right detector index
    const auto volIt = mChip2Hiearchy.find(frame.label);
    if (volIt == mChip2Hiearchy.end()) {
      return true; // sensor not in hierarchy, skip
    }
    const auto* tileVol = volIt->second;

    // TRK -> tile LOC -> halfBarrel LOC (same chain as the GBL hierarchy walk)
    const Matrix26 der = getLeafRigidBodyDerivatives(*tileVol, frame, makeDerivativeContext(frame, wTrk)) * tileVol->getJL2P();
    // apply: delta_res = der * delta_a_halfBarrel
    const Eigen::Vector2d shift = der * Eigen::Map<const Eigen::Matrix<double, 6, 1>>(mMisalignment[sensorID].rigidBody.data());
    res[0] += shift[0]; // dy
    res[1] += shift[1]; // dz
  }

  // --- In-extensional deformation ---
  // displacement field u(phi,z) = (u_phi, u_z, u_r)
  //   dy = -u_phi + y' * u_r,  dz = -u_z + z' * u_r
  if (mParams->doMisalignmentInex) {
    const auto shift = evaluateInextensionalShift(mMisalignment[sensorID], misFrame, TrackSlopes::computeTrackSlopes(wTrk.getSnp(), wTrk.getTgl()));
    res[0] += shift.dy;
    res[1] += shift.dz;
  }

  if (mOutOpt[o2::alignrs::OutputOpt::MisRes] && mDBGOut) {    
    (*mDBGOut) << "mis"
               << "dy=" << res[0]
               << "dz=" << res[1]
               << "sens=" << sensorID
               << "lay=" << layerID
               << "z=" << frame.cluster.getZ()
               << "phi=" << frame.alpha
               << "\n";
  }

  return true;
}

// For every attached detector with a geometry counterpart: read its initial alignment (the CCDB
// object the fit started from, given as DET:file in algParamsInitial), combine it with the fitted
// rigid-body corrections of its branch (Detector::MP2AlignParams) and write the result, in the CCDB
// format, to <DET>_<algParamsOutFile>. A detector without initial file is assumed to have been ideal.
void AlignmentSpec::writeAlignParams(const std::map<uint32_t, double>& labelToValue) const
{
  using AlgParVec = std::vector<o2::detectors::AlignParam>;
  std::map<std::string, std::string> initialFiles; // detector name -> file
  for (const auto& item : o2::utils::Str::tokenize(mParams->algParamsInitial, ',')) {
    const auto pos = item.find(':');
    if (pos == std::string::npos) {
      LOGP(fatal, "algParamsInitial entry '{}' is not of the form DET:file", item);
    }
    initialFiles[item.substr(0, pos)] = item.substr(pos + 1);
  }

  for (const auto* det : mDetectors) {
    if (det->getTopVolume() == nullptr) {
      continue; // branch discarded
    }
    AlgParVec initial;
    if (auto it = initialFiles.find(det->getDetName()); it != initialFiles.end()) {
      TFile fin(it->second.c_str());
      std::unique_ptr<AlgParVec> pars{fin.IsZombie() ? nullptr : fin.Get<AlgParVec>("ccdb_object")};
      if (pars == nullptr) {
        LOGP(fatal, "Failed to read the initial alignment of {} from {}", det->getDetName(), it->second);
      }
      initial = std::move(*pars);
      LOGP(info, "Read {} initial AlignParam objects of {} from {}", initial.size(), det->getDetName(), it->second);
    } else {
      LOGP(warn, "No initial alignment provided for {}: the geometry the fit was done on is assumed to be the ideal one", det->getDetName());
    }
    auto result = det->MP2AlignParams(labelToValue, initial, mParams->writeLocalAlignParams);
    if (result.empty()) {
      continue; // nothing to align, e.g. a detector made only of virtual volumes and never aligned
    }
    const auto outName = std::format("{}_{}", det->getDetName(), mParams->algParamsOutFile);
    TFile fout(outName.c_str(), "RECREATE");
    fout.WriteObjectAny(&result, "std::vector<o2::detectors::AlignParam>", "ccdb_object");
    LOGP(info, "Wrote {} AlignParam objects of {} to {}", result.size(), det->getDetName(), outName);
  }
}

void AlignmentSpec::endOfStream(EndOfStreamContext& /*ec*/)
{
  if (mDBGOut) {
    mDBGOut->Close();
    mDBGOut.reset();
  }
  mMille.reset(); // flushes and closes the binary
}

void AlignmentSpec::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    LOG(info) << "its cluster dictionary updated";
    if (mITS) { // after the 1st TF the detector holds its own pointers
      mITS->setTopologyDictionaries((const o2::itsmft::TopologyDictionary*)obj, nullptr);
    }
    return;
  }
  if (matcher == ConcreteDataMatcher("IT3", "CLUSDICT", 0)) {
    LOG(info) << "it3 cluster dictionary updated";
    if (mITS) {
      mITS->setTopologyDictionaries(nullptr, (const o2::its3::TopologyDictionary*)obj);
    }
    return;
  }
  if (matcher == ConcreteDataMatcher("TRD", "CALVDRIFTEXB", 0)) {
    LOG(info) << "TRD CalVdriftExB updated";
    mTRD->setCalVdriftExB((const o2::trd::CalVdriftExB*)obj);
    return;
  }
  if (matcher == ConcreteDataMatcher("GLO", "MEANVERTEX", 0)) {
    mPVT->setMeanVertexCCDB(*(const o2::dataformats::MeanVertexObject*)obj);
    return;
  }
  if (mTPCVDriftHelper.accountCCDBInputs(matcher, obj)) {
    return;
  }
}

DataProcessorSpec getAlignmentSpec(GTrackID::mask_t srcTracks, GTrackID::mask_t srcClusters, bool useMC, bool withITS3, bool requestCTPLumi, o2::alignrs::OutputEnum out)
{
  auto dataRequest = std::make_shared<DataRequest>();
  std::shared_ptr<o2::base::GRPGeomRequest> ggRequest{nullptr};
  auto detMask = GTrackID::getSourcesDetectorsMask(srcClusters);

  if (!out[o2::alignrs::OutputOpt::MilleRes]) {
    dataRequest->requestTracks(srcTracks, useMC);
    if (withITS3) {
      dataRequest->requestIT3Clusters(useMC);
    } else {
      dataRequest->requestClusters(srcClusters, useMC);
    }
    if (requestCTPLumi) {
      dataRequest->inputs.emplace_back("lumiCTP", o2::header::gDataOriginCTP, "LUMICTP", 0, Lifetime::Timeframe);
    }

    dataRequest->requestPrimaryVertices(useMC);
    // the mean vertex is the prior of the primary vertex of every collision and the starting point
    // of the alignment of the virtual PVT detector
    dataRequest->inputs.emplace_back("meanvtx", "GLO", "MEANVERTEX", 0, Lifetime::Condition, ccdbParamSpec("GLO/Calib/MeanVertex", {}, 1));
    if (detMask[DetID::TRD]) { // the tracklet transformation needs the drift velocity and ExB
      dataRequest->inputs.emplace_back("calvdexb", "TRD", "CALVDRIFTEXB", 0, Lifetime::Condition, ccdbParamSpec("TRD/Calib/CalVdriftExB"));
    }
    if (detMask[DetID::TPC]) {
      // the cluster transformation: the maps are delivered per TF, the drift calibration they
      // account for is read from the CCDB just to be able to interpret the fitted drift correction
      o2::tpc::VDriftHelper::requestCCDBInputs(dataRequest->inputs);
      dataRequest->inputs.emplace_back("corrMap", o2::header::gDataOriginTPC, "TPCCORRMAP", 0, Lifetime::Timeframe);
    }
    ggRequest = std::make_shared<o2::base::GRPGeomRequest>(true,                              // orbitResetTime
                                                           false,                             // GRPECS=true
                                                           true,                              // GRPLHCIF
                                                           true,                              // GRPMagField
                                                           true,                              // askMatLUT
                                                           o2::base::GRPGeomRequest::Aligned, // geometry
                                                           dataRequest->inputs,               // inputs
                                                           true,                              // askOnce
                                                           true);                             // propagatorD
  } else {
    dataRequest->inputs.emplace_back("dummy", "GLO", "DUMMY_OUT", 0);
    ggRequest = std::make_shared<o2::base::GRPGeomRequest>(true,                              // orbitResetTime
                                                           false,                             // GRPECS=true
                                                           false,                             // GRPLHCIF
                                                           false,                             // GRPMagField
                                                           false,                             // askMatLUT
                                                           o2::base::GRPGeomRequest::Aligned, // geometry
                                                           dataRequest->inputs);
  }

  Options opts{
    {"nthreads", VariantType::Int, 1, {"number of threads"}},
    {"config-macro", VariantType::String, "none", {"configuration macro with signature (std::vector<o2::alignrs::Detector*>*, int) to execute from init"}},
  };

  return DataProcessorSpec{
    .name = "barrel-alignment",
    .inputs = dataRequest->inputs,
    .outputs = {},
    .algorithm = AlgorithmSpec{adaptFromTask<AlignmentSpec>(dataRequest, ggRequest, srcTracks, detMask, useMC, withITS3, out)},
    .options = opts};
}
} // namespace o2::alignrs
