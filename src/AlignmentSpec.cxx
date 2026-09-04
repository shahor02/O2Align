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
#include <fstream>
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

#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "Steer/MCKinematicsReader.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "ITStracking/MathUtils.h"
#include "ITStracking/IOUtils.h"
#include "ITS3Reconstruction/IOUtils.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "O2Align/TrackFit.h"
#include "O2Align/AlignmentSpec.h"
#include "O2Align/Params.h"
#include "O2Align/AlignmentTypes.h"
#include "O2Align/Volume.h"
#include "O2Align/MisalignmentUtils.h"
#include "O2Align/DetectorITS.h"
#include "O2Align/SensorITS.h"

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
  return {.sensorID = isITS3 ? o2::its3::constants::detID::getSensorID(frame.cluster.getSensorID()) : -1,
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

Matrix36 getRigidBodyBaseDerivatives(const DerivativeContext& ctx)
{
  static const RigidBodyDOFSet sRigidBodyBasis;
  Eigen::MatrixXd dyn(3, sRigidBodyBasis.nDOFs());
  sRigidBodyBasis.fillDerivatives(ctx, dyn);
  return dyn;
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


  ~AlignmentSpec() final = default;
  AlignmentSpec(const AlignmentSpec&) = delete;
  AlignmentSpec(AlignmentSpec&&) = delete;
  AlignmentSpec& operator=(const AlignmentSpec&) = delete;
  AlignmentSpec& operator=(AlignmentSpec&&) = delete;
  AlignmentSpec(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr, GTrackID::mask_t src, bool useMC, bool withITS, o2::alignrs::OutputEnum out)
    : mDataRequest(dr), mGGCCDBRequest(gr), mTracksSrcMask(src), mUseMC(useMC), mIsITS3(!withITS), mOutOpt(out), mITS(std::make_unique<DetectorITS>(!withITS))
  {
  }

  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;
  void process();

 private:
  void updateTimeDependentParams(ProcessingContext& pc);
  void buildHierarchy();

  bool processITSPart(Track& resTrack, const GlobalIDSet& contributorsGID);
  bool processTPCPart(Track& resTrack, const GlobalIDSet& contributorsGID);
  bool processTRDPart(Track& resTrack, const GlobalIDSet& contributorsGID);
  bool processTOFPart(Track& resTrack, const GlobalIDSet& contributorsGID);

  // calculate the transport jacobian for points FROM and TO numerically via ridder's method
  // this assumes the track is already at point FROM and will be extrapolated to TO's x (xTo)
  // method does not modify the original track
  bool getTransportJacobian(const TrackD& track, double xTo, double alphaTo, gbl::Matrix5d& jac, gbl::Matrix5d& err);

  // refit ITS track with inward/outward fit (opt. impose pv as additional constraint)
  // after this we have the refitted track at the innermost update point
  bool prepareITSTrack(int iTrk, const o2::its::TrackITS& itsTrack, Track& resTrack);

  // prepare ITS measuremnt points
  // build track to vertex association
  void buildT2V();

  // apply some misalignment on inner ITS3 layers
  // it can happen that a measurement is pushed outside of
  // ITS3 acceptance so false is to discard track
  bool applyMisalignment(Eigen::Vector2d& res, const FrameInfoExt& frame, const TrackD& wTrk, size_t iTrk);

  o2::alignrs::OutputEnum mOutOpt;
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDBGOut;
  std::vector<dataformats::VertexBase> mPVMC;
  std::vector<int> mT2PVMC;
  bool mIsITS3{true};
  const o2::itsmft::TopologyDictionary* mITSDict{nullptr};
  const o2::its3::TopologyDictionary* mIT3Dict{nullptr};
  o2::globaltracking::RecoContainer* mRecoData = nullptr;
  std::unique_ptr<steer::MCKinematicsReader> mcReader;

  std::unique_ptr<DetectorITS> mITS;

  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::unique_ptr<Volume> mHierarchy;   // tree-hiearchy
  Volume::SensorMapping mChip2Hiearchy; // global label mapping to leaves in the tree
  ProcStat mStat{}; // processing statistics
  bool mUseMC{false};
  bool mUsePVConstraint{false}; // use PV as additional constraint in a given track refit  //RSTODO: this should be a per-track decision, not global, should not be datamember
  GTrackID::mask_t mTracksSrcMask;
  std::vector<int> mTrackSources;
  int mNThreads{1};
  const Params* mParams{nullptr};
  MisalignmentModel mMisalignment;
  std::array<Eigen::Matrix<double, 6, 1>, 6> mRigidBodyParams; // (dx,dy,dz,rx,ry,rz) in LOC per sensorID
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
}

void AlignmentSpec::run(ProcessingContext& pc)
{
  if (mOutOpt[o2::alignrs::OutputOpt::MilleRes]) {
    updateTimeDependentParams(pc);
    Volume::writeMillepedeResults(mHierarchy.get(), mParams->milleResFile, mParams->milleResOutJson, mParams->misAlgJson); // RSTODO loop over possible top volumes (detectors)
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
  auto prop = o2::base::PropagatorD::Instance();
  const auto bz = prop->getNominalBz();
  const bool fieldON = std::abs(bz) > 0.1;
  std::span<const o2::MCCompLabel> mcLbls;
  if (mUseMC) {
    mcLbls = mRecoData->getITSTracksMCLabels();
  }
  // prepare detector data
  mITS->prepareData(mRecoData);

  if (mParams->usePVConstraintMinTrackes > 0) {
    buildT2V(); // RSTODO for data
  }

  if (mNThreads > 1 && !(mParams->misAlgJson.empty())) {
    LOGP(warn, "Applying misalignment works only single-threaded, forcing to 1");
    mNThreads = 1;
  }
  LOGP(info, "Starting fits with {} threads", mNThreads);
  std::vector<gbl::GblTrajectory> gblTraj;
  std::vector<Track> resTracks;

  const auto primVertices = mRecoData->getPrimaryVertices();
  const auto primVer2TRefs = mRecoData->getPrimaryVertexMatchedTrackRefs();
  const auto primVerGIs = mRecoData->getPrimaryVertexMatchedTracks();
  std::unordered_map<GTrackID, bool> ambigTable;
  // process vertices with contributor tracks
  const int nvRefs = primVer2TRefs.size();

  auto timeStart = std::chrono::high_resolution_clock::now();
  int cFailedRefit{0}, cFailedProp{0}, cSelected{0}, cGBLFit{0}, cGBLFitFail{0}, cGBLChi2Rej{0}, cGBLConstruct{0};
  double chi2Sum{0}, lostWeightSum{0};
  int ndfSum{0};

  int nVtx = 0, nVtxAcc = 0, nTrc = 0, nTrcAcc = 0;
  for (int ivref = 0; ivref < nvRefs; ivref++) {
    const o2::dataformats::PrimaryVertex* vtx = (ivref < nvRefs - 1) ? &primVertices[ivref] : nullptr;
    const auto& trackRef = primVer2TRefs[ivref];
    bool useVertexConstraint = false;
    if (vtx && mParams->usePVConstraintMinTrackes > 0 && vtx->getNContributors() >= mParams->usePVConstraintMinTrackes) {
      useVertexConstraint = true;
      mStat.data[ProcStat::kInput][ProcStat::kVertices]++;
    }
    if (mParams->verbose > 1) {
      LOGP(info, "processing vtref {} of {} with {} tracks, {}", ivref, nvRefs, trackRef.getEntries(), vtx ? vtx->asString() : std::string{});
    }
    nVtx++;
    for (int src : mTrackSources) {      
      int start = trackRef.getFirstEntryOfSource(src), end = start + trackRef.getEntriesOfSource(src);
      for (int ti = start; ti < end; ti++) {
        const auto trackIndex = primVerGIs[ti];
        if (trackIndex.isAmbiguous()) {
          auto& ambSeen = ambigTable[trackIndex]; // ambiguous track does not use PV constraint, no need to account it again
          if (ambSeen) { // processed
            continue;
          }
          ambSeen = true;
        }
        const auto& trPar = mRecoData->getTrackParam(trackIndex);
        if (fieldON && trPar.getPt() < mParams->minPt) {
          continue;
        }
        mStat.data[ProcStat::kInput][ProcStat::kTracks]++;
        // account preliminary
        resTracks.emplace_back().gid = trackIndex;
      }
    }
    // fit the tracks, prepare for GLB
#ifdef WITH_OPENMP
#pragma omp parallel num_threads(mNThreads)
#endif
    for (size_t itr = 0; itr < (int)resTracks.size(); itr++) {
      auto &track = resTracks[itr];
      auto contributorsGID = mRecoData->getSingleDetectorRefs(track.gid);
      if ( track.gid.includesDet(DetID::ITS) && !processITSPart(track, contributorsGID) ) {
        track.gid.clear(); // mark as failed
        continue;
      }
      if ( track.gid.includesDet(DetID::TPC) && !processTPCPart(track, contributorsGID) ) { // do we want to abandont the track if TPC fails? or just continue with ITS?
        track.gid.clear(); // mark as failed
        continue;
      }
      if ( track.gid.includesDet(DetID::TRD) && !processTRDPart(track, contributorsGID) ) { // do we want to abandont the track if TRD fails? or just continue with ITS?
        track.gid.clear(); // mark as failed
        continue;
      }
      if ( track.gid.includesDet(DetID::TOF) && !processTOFPart(track, contributorsGID) ) { // do we want to abandont the track if TOF fails? or just continue with ITS?
        track.gid.clear(); // mark as failed
        continue;
      }
    }
  }

/*
  // Data

#ifdef WITH_OPENMP
#pragma omp parallel num_threads(mNThreads) \
  reduction(+ : cFailedRefit)               \
  reduction(+ : cFailedProp)                \
  reduction(+ : cSelected)                  \
  reduction(+ : cGBLFit)                    \
  reduction(+ : cGBLFitFail)                \
  reduction(+ : cGBLChi2Rej)                \
  reduction(+ : cGBLConstruct)              \
  reduction(+ : chi2Sum)                    \
  reduction(+ : lostWeightSum)              \
  reduction(+ : ndfSum)
#endif
  {
#ifdef WITH_OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    auto& gblTrajSlot = gblTrajSlots[tid];
    auto& resTrackSlot = resTrackSlots[tid];

#ifdef WITH_OPENMP
#pragma omp for schedule(dynamic)
#endif
    for (size_t iTrk = 0; iTrk < (int)itsTracks.size(); ++iTrk) {
      const auto& trk = itsTracks[iTrk];
      if (trk.getNClusters() < mParams->minITSCls ||
          (trk.getChi2() / ((float)trk.getNClusters() * 2 - 5)) >= mParams->maxITSChi2Ndf ||
          trk.getPt() < mParams->minPt ||
          (mUseMC && (!mcLbls[iTrk].isValid() || !mcLbls[iTrk].isCorrect()))) {
        continue;
      }
      ++cSelected;
      Track& resTrack = resTrackSlot.emplace_back();
      if (!prepareITSTrack((int)iTrk, trk, resTrack)) {
        ++cFailedRefit;
        resTrackSlot.pop_back();
        continue;
      }

      o2::track::TrackParD* refLin = nullptr;
      if (mParams->useStableRef) {
        refLin = &resTrack.track;
      }

      // outward stepping from track IU
      auto wTrk = resTrack.track;
      const bool hasPV = resTrack.info[0].lr == -1;
      std::vector<gbl::GblPoint> points;
      bool failed = false;
      const int np = (int)resTrack.points.size();
      track::TrackLTIntegral lt;
      lt.setTimeNotNeeded();
      constexpr int perm[5] = {4, 2, 3, 0, 1}; // ALICE->GBL: Q/Pt,Snp,Tgl,Y,Z
      for (int ip{0}; ip < np; ++ip) {
        const auto& frame = resTrack.info[ip];
        gbl::Matrix5d err = gbl::Matrix5d::Identity(), jacALICE = gbl::Matrix5d::Identity(), jacGBL;
        float msErr = 0.f;
        if (ip) {
          // numerically calculates the transport jacobian from prev. point to this point
          // then we actually do the step to the point and accumulate the material
          if (!getTransportJacobian(wTrk, frame.x, frame.alpha, jacALICE, err) ||
              !prop->propagateToAlphaX(wTrk, refLin, frame.alpha, frame.x, false, mParams->maxSnp, mParams->maxStep, 1, mParams->corrType, &lt)) {
            ++cFailedProp;
            failed = true;
            break;
          }
          msErr = its::math_utils::MSangle(trk.getPID().getMass(), trk.getP(), lt.getX2X0());
          // after computing jac, reorder to GBL convention
          for (int i = 0; i < 5; i++) {
            for (int j = 0; j < 5; j++) {
              jacGBL(i, j) = jacALICE(perm[i], perm[j]);
            }
          }
        }

        // wTrk is now in the measurment frame
        gbl::GblPoint point(jacGBL);
        // measurement
        Eigen::Vector2d res, prec;
        res << frame.positionTrackingFrame[0] - wTrk.getY(), frame.positionTrackingFrame[1] - wTrk.getZ();

        // here we can apply some misalignment on the measurment
        if (!applyMisalignment(res, frame, wTrk, iTrk)) {
          failed = true;
          break;
        }

        prec << 1. / resTrack.points[ip].sig2y, 1. / resTrack.points[ip].sig2z;
        // the projection matrix is in the tracking frame the idendity so no need to diagonalize it
        point.addMeasurement(res, prec);
        if (msErr > mParams->minMS && ip < np - 1) {
          Eigen::Vector2d scat(0., 0.), scatPrec = Eigen::Vector2d::Constant(1. / (msErr * msErr));
          point.addScatterer(scat, scatPrec);
          lt.clearFast(); // clear if accounted
        }

        if (frame.lr >= 0) {
          Label lbl(0, frame.sens, true);
          if (mChip2Hiearchy.find(lbl) == mChip2Hiearchy.end()) {
            LOGP(fatal, "Cannot find global label: {}", lbl.asString());
          }

          // derivatives for all sensitive volumes and their parents
          // this is the derivative in TRK but we want to align in LOC
          // so dr/da_(LOC) = dr/da_(TRK) * da_(TRK)/da_(LOC)
          const auto* tileVol = mChip2Hiearchy.at(lbl);
          const auto derCtx = makeDerivativeContext(frame, wTrk);
          Matrix36 der = getRigidBodyBaseDerivatives(derCtx);

          // count rigid body columns: only volumes with real DOFs (not DOFPseudo)
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
          const int nCol = nColRB + nCalib;

          std::vector<int> gLabels;
          gLabels.reserve(nCol);
          Eigen::MatrixXd gDer(3, nCol);
          gDer.setZero();
          Eigen::Index curCol{0};

          // 1) tile: TRK -> LOC via precomputed T2L and J_L2T
          const double posTrk[3] = {frame.x, 0., 0.};
          double posLoc[3];
          tileVol->getT2L().LocalToMaster(posTrk, posLoc);
          Matrix66 jacL2T;
          tileVol->computeJacobianL2T(posLoc, jacL2T);
          der *= jacL2T;
          if (tileVol->getRigidBody()) {
            const int nd = tileVol->getRigidBody()->nDOFs();
            for (int iDOF = 0; iDOF < nd; ++iDOF) {
              gLabels.push_back(tileVol->getLabel().rawGBL(iDOF));
            }
            gDer.middleCols(curCol, nd) = der;
            curCol += nd;
          }

          // 2) chain through parents: child's J_L2P
          for (const auto* child = tileVol; child->getParent() && !child->getParent()->isRoot(); child = child->getParent()) {
            der *= child->getJL2P();
            const auto* parent = child->getParent();
            if (parent->getRigidBody()) {
              const int nd = parent->getRigidBody()->nDOFs();
              for (int iDOF = 0; iDOF < nd; ++iDOF) {
                gLabels.push_back(parent->getLabel().rawGBL(iDOF));
              }
              gDer.middleCols(curCol, nd) = der;
              curCol += nd;
            }
          }

          // 3) calibration derivatives (apply directly on the whole sensor, not on individual tiles)
          if (calibSet) {
            const int nd = calibSet->nDOFs();
            Eigen::MatrixXd calDer(3, nd);
            calibSet->fillDerivatives(derCtx, calDer);
            for (int iDOF = 0; iDOF < nd; ++iDOF) {
              gLabels.push_back(sensorVol->getLabel().asCalib().rawGBL(iDOF));
            }
            gDer.middleCols(curCol, nd) = calDer;
            curCol += nd;
          }
          point.addGlobals(gLabels, gDer);
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
          LOGP(info, "Point {}: GBL res=({}, {}), KF stored res=({}, {})",
               ip, res[0], res[1], resTrack.points[ip].dy, resTrack.points[ip].dz);
          LOGP(info, "residual: dy={} dz={}", res[0], res[1]);
          LOGP(info, "precision: precY={} precZ={}", prec[0], prec[1]);
          point.printPoint(5);
        }
        points.push_back(point);
      }
      if (!failed) {
        gbl::GblTrajectory traj(points, std::abs(bz) > 0.01);
        if (traj.isValid()) {
          double chi2 = NAN, lostWeight = NAN;
          int ndf = 0;
          if (auto ierr = traj.fit(chi2, ndf, lostWeight); !ierr) {
            if (mOutOpt[o2::alignrs::OutputOpt::VerboseGBL]) {
              LOGP(info, "GBL FIT chi2 {} ndf {}", chi2, ndf);
              traj.printTrajectory(5);
            }
            if (chi2 / ndf > mParams->maxChi2Ndf && cGBLChi2Rej++ < 10) {
              LOGP(error, "GBL fit exceeded red chi2 {}", chi2 / ndf);
              if (std::abs(resTrack.kfFit.chi2Ndf - 1) < 0.02) {
                LOGP(error, "\tGBL is far away from good KF fit!!!!");
                continue;
              }
            } else {
              ++cGBLFit;
              chi2Sum += chi2;
              lostWeightSum += lostWeight;
              ndfSum += ndf;
              if (mOutOpt[o2::alignrs::OutputOpt::MilleData]) {
                gblTrajSlot.push_back(traj);
              }
              FitInfo fit;
              fit.ndf = ndf;
              fit.chi2 = (float)chi2;
              fit.chi2Ndf = (float)chi2 / (float)ndf;
              resTrack.gblFit = fit;
            }
          } else {
            ++cGBLFitFail;
          }
        } else {
          ++cGBLConstruct;
        }
      }
    }
  }
  auto timeEnd = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart);
  LOGP(info, "Fitted {} tracks out of {} (selected {}) in {} sec", cGBLFit, itsTracks.size(), cSelected, duration.count() / 1e3);
  LOGP(info, "\tRefit failed for {} tracks; Failed prop for {} tracks", cFailedRefit, cFailedProp);
  LOGP(info, "\tGBL SUMMARY:");
  LOGP(info, "\t\tGBL construction failed {}", cGBLConstruct);
  LOGP(info, "\t\tGBL fit failed {}", cGBLFitFail);
  LOGP(info, "\t\tGBL chi2Ndf rejected {}", cGBLChi2Rej);
  if (!ndfSum) {
    LOGP(info, "\t\tGBL Chi2/Ndf = NDF IS 0");
  } else {
    LOGP(info, "\t\tGBL Chi2/Ndf = {}", chi2Sum / ndfSum);
  }
  LOGP(info, "\t\tGBL LostWeight = {}", lostWeightSum);
  LOGP(info, "Streaming results to output");
  if (mOutOpt[o2::alignrs::OutputOpt::MilleData]) {
    gbl::MilleBinary mille(mParams->milleBinFile, true);
    for (auto& slot : gblTrajSlots) {
      for (auto& traj : slot) {
        traj.milleOut(mille);
      }
    }
  }
  if (mOutOpt[o2::alignrs::OutputOpt::Debug]) {
    for (auto& slot : resTrackSlots) {
      for (auto& res : slot) {
        (*mDBGOut) << "res"
                   << "trk=" << res
                   << "\n";
      }
    }
  }
  */
}

void AlignmentSpec::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);  
  if (static bool initOnce{false}; !initOnce) {
    initOnce = true;
    mITS->setTopologyDictionaries(mITSDict, mIT3Dict);
    auto geom = o2::its::GeometryTGeo::Instance();
    o2::its::GeometryTGeo::Instance()->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G, o2::math_utils::TransformType::T2G));
    mParams = &Params::Instance();
    mParams->printKeyValues(true, true);
    buildHierarchy();

    if (mParams->doMisalignmentLeg || mParams->doMisalignmentRB || mParams->doMisalignmentInex) {
      mMisalignment = {};
      for (auto& rb : mRigidBodyParams) {
        rb.setZero();
      }
      if (!mParams->misAlgJson.empty()) {
        mMisalignment = loadMisalignmentModel(mParams->misAlgJson);
        if (mParams->doMisalignmentRB) {
          using json = nlohmann::json;
          std::ifstream f(mParams->misAlgJson);
          auto data = json::parse(f);
          for (const auto& item : data) {
            int id = item["id"].get<int>();
            if (!item.contains("rigidBody")) {
              continue;
            }
            auto rb = item["rigidBody"].get<std::vector<double>>();
            for (int k = 0; k < 6 && k < static_cast<int>(rb.size()); ++k) {
              mRigidBodyParams[id](k) = rb[k];
            }
          }
        }
      }
    }
  }
}

void AlignmentSpec::buildHierarchy()
{
  mITS->buildHierarchy(mChip2Hiearchy);

  if (!mParams->dofConfigJson.empty()) {
    Volume::applyDOFConfig(mHierarchy.get(), mParams->dofConfigJson); // RSTODO loop over detectors
  }

  mHierarchy->finalise();
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
  const auto bz = prop->getNominalBz();
  const auto minStep = std::sqrt(std::numeric_limits<double>::epsilon());
  const gbl::Vector5d x0(track.getParams());
  auto trackC = track;
  o2::track::TrackParD* refLin{nullptr};
  if (mParams->useStableRef) {
    refLin = &trackC;
  }

  auto propagate = [&](gbl::Vector5d& p) -> bool {
    TrackD tmp(track);
    for (int i{0}; i < track::kNParams; ++i) {
      tmp.setParam(p[i], i);
    }
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

bool AlignmentSpec::processITSPart(Track& resTrack, const GlobalIDSet& contributorsGID)
{
  return true;
}

bool AlignmentSpec::processTPCPart(Track& resTrack, const GlobalIDSet& contributorsGID)
{
  return true;
}


bool AlignmentSpec::processTRDPart(Track& resTrack, const GlobalIDSet& contributorsGID)
{
  return true;
}

bool AlignmentSpec::processTOFPart(Track& resTrack, const GlobalIDSet& contributorsGID)
{
  return true;
}

bool AlignmentSpec::prepareITSTrack(int iTrk, const o2::its::TrackITS& itsTrack, Track& resTrack)
{
  const auto itsClRefs = mRecoData->getITSTracksClusterRefs();
  auto trFit = convertTrack<double>(itsTrack.getParamOut()); // take outer track fit as start of refit
  auto prop = o2::base::PropagatorD::Instance();
  auto geom = o2::its::GeometryTGeo::Instance();
  const auto bz = prop->getNominalBz();
  std::array<const FrameInfoExt*, 8> frameArr{};
  o2::track::TrackParD trkOut, *refLin = nullptr;
  if (mParams->useStableRef) {
    refLin = &(trkOut = trFit);
  }

  auto accountCluster = [&](int i, TrackD& tr, float& chi2, Measurement& meas, o2::track::TrackParD* refLin) {
    if (frameArr[i]) { // update with cluster
      if (!prop->propagateToAlphaX(tr, refLin, frameArr[i]->alpha, frameArr[i]->x, false, mParams->maxSnp, mParams->maxStep, 1, mParams->corrType)) {
        return 2;
      }
      const auto& cluster = frameArr[i]->cluster;
      meas.dy = cluster.getY() - tr.getY();
      meas.dz = cluster.getZ() - tr.getZ();
      meas.sig2y = cluster.getSigmaY2();
      meas.sig2z = cluster.getSigmaZ2();
      meas.z = tr.getZ();
      meas.phi = tr.getPhi();
      o2::math_utils::bringTo02Pid(meas.phi);
      chi2 += (float)tr.getPredictedChi2Quiet(cluster);
      if (!tr.update(cluster)) {
        return 2;
      }
      if (refLin) { // displace the reference to the last updated cluster
        refLin->setY(cluster.getY());
        refLin->setZ(cluster.getZ());
      }
      return 0;
    }
    return 1;
  };

  FrameInfoExt pvInfo;
  if (mUsePVConstraint) { // add PV as constraint RSTODO: this should be a per-track decision, not global, should not be datamember
    const int iPV = mT2PVMC[iTrk];
    if (iPV < 0) {
      return false;
    }
    const auto& pv = mPVMC[iPV];
    auto tmp = convertTrack<double>(itsTrack.getParamIn());
    if (!prop->propagateToDCA(pv, tmp, bz)) {
      return false;
    }
    pvInfo.alpha = (float)tmp.getAlpha();
    double ca{0}, sa{0};
    o2::math_utils::bringToPMPi(pvInfo.alpha);
    o2::math_utils::sincosd(pvInfo.alpha, sa, ca);
    pvInfo.x = tmp.getX();
    pvInfo.cluster = ClusterF(-1, 0.f, -pv.getX() * sa + pv.getY() * ca, pv.getZ(), 0.5 * (pv.getSigmaX2() + pv.getSigmaY2()), pv.getSigmaY2(), 0.);
    pvInfo.lr = -1;
    frameArr[0] = &pvInfo;
  }

  // collect all track clusters to array, placing them to layer+1 slot
  int nCl = itsTrack.getNClusters();
  for (int i = 0; i < nCl; i++) { // clusters are ordered from the outermost to the innermost
    const auto& curInfo = mITS->getPointsInfo()[itsClRefs[itsTrack.getClusterEntry(i)]];
    frameArr[1 + curInfo.lr] = &curInfo;
  }

  // start refit
  resTrack.points.clear();
  resTrack.info.clear();
  trFit.resetCovariance();
  trFit.setCov(trFit.getQ2Pt() * trFit.getQ2Pt() * trFit.getCov()[14], 14);
  float chi2{0};
  for (int i{7}; i >= 0; --i) {
    Measurement point;
    int res = accountCluster(i, trFit, chi2, point, refLin);
    if (res == 2) {
      return false;
    } else if (res == 0) {
      resTrack.points.push_back(point);
      resTrack.info.push_back(*frameArr[i]);
      resTrack.track = trFit; // put track to whatever the IU is
    }
  }
  // reverse inserted points so they are in the same order as the track
  std::reverse(resTrack.info.begin(), resTrack.info.end());
  std::reverse(resTrack.points.begin(), resTrack.points.end());
  resTrack.kfFit.chi2 = chi2;
  resTrack.kfFit.ndf = (int)resTrack.info.size() * 2 - 5;
  resTrack.kfFit.chi2Ndf = chi2 / (float)resTrack.kfFit.ndf;

  return true;
}

#if 0
void AlignmentSpec::prepareITSPoints()
{
  // prepare TF ITS data for processing: convert clusters and build overlap info. // RSTODO At the moment it is not ITS3/staggering compatible!
  const auto clusITS = mRecoData->getITSClusters();
  const auto clusITSROF = mRecoData->getITSClustersROFRecords();
  const auto patterns = mRecoData->getITSClustersPatterns();
  auto pattIt = patterns.begin();
  mITSPointsInfo.clear();
  mITSPointsInfo.reserve(clusITS.size());
  if (mParams->ITSOverlapMargin > 0) {
    mITSOvlClusRef.clear();
    mITSOvlClusRef.resize(clusITS.size(), -1);
    mITSOvlCandidateID.clear();
    mITSOvlCandidateID.reserve(clusITS.size());
  }
  auto geom = its::GeometryTGeo::Instance();
  std::vector<int> edgeClusters;
  int ROFCount = 0, curSensID = -1;
  struct ROFChipEntry {
    int rofCount = -1;
    int chipFirstEntry = -1;
  };
  std::array<ROFChipEntry, o2::itsmft::ChipMappingITS::getNChips()> chipROFStart{}; // fill only for clusters with overlaps

  for (const auto& rof : clusITSROF) {
    int maxic = rof.getFirstEntry() + rof.getNEntries();
    edgeClusters.clear();
    for (int ic = rof.getFirstEntry(); ic < maxic; ic++) {
      const auto& cls = clusITS[ic];
      const auto sensID = cls.getSensorID();
      const auto lay = geom->getLayer(sensID);
      float sigmaY2{0.}, sigmaZ2{0.};
      math_utils::Point3D<float> locXYZ;
      if (mIsITS3) {
        locXYZ = o2::its3::ioutils::extractClusterData(cls, pattIt, mIT3Dict, sigmaY2, sigmaZ2);
      } else {
        locXYZ = o2::its::ioutils::extractClusterData(cls, pattIt, mITSDict, sigmaY2, sigmaZ2);
      }
      sigmaY2 += mParams->extraClsErrYITS[lay] * mParams->extraClsErrYITS[lay];
      sigmaZ2 += mParams->extraClsErrZITS[lay] * mParams->extraClsErrZITS[lay];
      const auto gloXYZ = geom->getMatrixL2G(sensID) * locXYZ; // local --> global
      auto trkXYZ = geom->getMatrixT2L(sensID) ^ locXYZ; // local --> tracking
      // Tracking alpha angle: We want that each cluster rotates its tracking frame to the clusters phi
      // that way the track linearization around the measurement is less biases to the arc
      // this means automatically that the measurement on the arc is at 0 for the curved layers
      double alpha = geom->getSensorRefAlpha(sensID); // RSTODO: what is this for ITS3 IB?
      double x = geom->getSensorRefX(sensID);
      if (mIsITS3 && o2::its3::constants::detID::isDetITS3(sensID)) {
        trkXYZ.SetY(0.f);
        x = std::hypot(gloXYZ.x(), gloXYZ.y());  // alpha, x always have to be defined wrt to the global Z axis!
        trkXYZ.SetX(x);
        alpha = std::atan2(gloXYZ.y(), gloXYZ.x());
      }
      math_utils::bringToPMPid(alpha);
      o2::BaseCluster<float> clus(sensID, trkXYZ, sigmaY2, sigmaZ2, 0.f);
      auto& pointInfo = mITSPointsInfo.emplace_back(lay, x, alpha, clus, {});
      if (mParams->ITSOverlapMargin > 0 && (!mIsITS3 || lay > 2)) {
        // fill chips overlaps info for clusters whose center is within of the mParams->ITSOverlapMargin distance from the chip min or max row edge
        // but the pixel closest to this edge has distance of at least mParams->ITSOverlapEdgeRows from the edge
        int row = 0, col = 0; // effective row/col of the cluster center
        o2::itsmft::SegmentationAlpide::localToDetectorUnchecked(locXYZ.X(), locXYZ.Z(), row, col);
        int drow = row < o2::itsmft::SegmentationAlpide::NRows / 2 ? row : o2::itsmft::SegmentationAlpide::NRows - row - 1; // distance to the edge
        if (drow * o2::itsmft::SegmentationAlpide::PitchRow < mParams->ITSOverlapMargin) {                                   // rough check is passed, check if the edge cluster is indeed good
          pointInfo.cluster.setBit(row < o2::itsmft::SegmentationAlpide::NRows / 2 ? EdgeFlags::LowRow : EdgeFlags::HighRow);            // flag that this is an edge cluster and indicate the low/high row side
          // check if it is not too close to the edge (to be biased)
          if (mParams->ITSOverlapEdgeRows > 0) { // is there a restriction?
            auto pattID = cls.getPatternID();
            drow = cls.getRow();
            if (pattID != itsmft::CompCluster::InvalidPatternID) {
              if (!mITSDict->isGroup(pattID)) {
                const auto& patt = mITSDict->getPattern(pattID); // reference pixel is min row/col corner
                if (row > o2::itsmft::SegmentationAlpide::NRows / 2) {
                  drow = o2::itsmft::SegmentationAlpide::NRows - 1 - (drow + patt.getRowSpan() - 1);
                }
              } else { // group: reference pixel is the one containing the COG
                o2::itsmft::ClusterPattern patt(pattItCopy);
                drow = row < o2::itsmft::SegmentationAlpide::NRows / 2 ? drow - patt.getRowSpan() / 2 : o2::itsmft::SegmentationAlpide::NRows - 1 - (drow + patt.getRowSpan() / 2 - 1);
              }
            } else {
              o2::itsmft::ClusterPattern patt(pattItCopy); // reference pixel is min row/col corner
              if (row > o2::itsmft::SegmentationAlpide::NRows / 2) {
                drow = o2::itsmft::SegmentationAlpide::NRows - 1 - (drow + patt.getRowSpan() - 1);
              }
            }
            if (drow < mParams->ITSOverlapEdgeRows) { // too close to the edge, flag this
              pointInfo.cluster.setBit(EdgeFlags::Biased);
            }
          }
          if (!pointInfo.cluster.isBitSet(EdgeFlags::Biased)) {
            if (chipROFStart[sensID].rofCount != ROFCount) { // remember 1st entry
              chipROFStart[sensID].rofCount = ROFCount;
              chipROFStart[sensID].chipFirstEntry = edgeClusters.size(); // remember 1st entry of edge cluster for this chip
            }
            edgeClusters.push_back(ic);
          }
        }
      }
    } // clusters of ROF
    // relate edge clusters of ROF to each other
    int prevSensID = -1;
    for (auto ic : edgeClusters) {
      auto& cl = mITSPointsInfo[ic].cluster;
        int sensID = cl.getSensorID();
        auto ovl = mOverlaps[sensID];
        int ovlCount = 0;
        for (int ir = 0; ir < OVL::NSides; ir++) {
          if (ovl.rowSide[ir] == OVL::NONE) { // no overlap from this row side
            continue;
        }
        int chipOvl = ovl.rowSide[ir]; // look for overlaps with this chip
        // are there clusters with overlaps on chipOvl?
        if (chipROFStart[chipOvl].rofCount == ROFCount) {
          auto oClusID = edgeClusters[chipROFStart[chipOvl].chipFirstEntry];
          while (oClusID < int(mITSPointsInfo.size())) {
            auto oClus = mITSPointsInfo[oClusID].cluster;
            if (oClus.getSensorID() != sensID) {
              break; // no more clusters on the overlapping chip
            }
            if (oClus.isBitSet(ovl.rowSideOverlap[ir]) &&                       // make sure that the edge cluster is on the right side of the row
                !oClus.isBitSet(EdgeFlags::Biased) &&                           // not too close to the edge
                std::abs(oClus.getZ() - cl.getZ()) < mParams->ITSOverlapMaxDZ) { // apply fiducial cut on Z distance of 2 clusters
              // register overlaping cluster
              if (!ovlCount) { // 1st overlap
                mITSOvlClusRef[ic] = mITSOvlCandidateID.size();
              }
              mITSOvlCandidateID.push_back(oClusID);
              ovlCount++;
            }
            oClusID++;
          }
        }
      }
      cl.setCount(std::min(127, ovlCount));
    }
    ROFCount++;
  } // loop over ROFs
}

#endif

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
    Label lbl(0, frame.cluster.getSensorID(), true);
    if (mChip2Hiearchy.find(lbl) == mChip2Hiearchy.end()) {
      return true; // sensor not in hierarchy, skip
    }
    const auto* tileVol = mChip2Hiearchy.at(lbl);

    // derivative in TRK frame (3x6: rows = dy, dz, dsnp)
    Matrix36 der = getRigidBodyBaseDerivatives(makeDerivativeContext(frame, wTrk));

    // TRK -> tile LOC
    const double posTrk[3] = {frame.x, 0., 0.};
    double posLoc[3];
    tileVol->getT2L().LocalToMaster(posTrk, posLoc);
    Matrix66 jacL2T;
    tileVol->computeJacobianL2T(posLoc, jacL2T);
    der *= jacL2T;

    // tile LOC -> halfBarrel LOC (same chain as GBL hierarchy walk)
    der *= tileVol->getJL2P();

    // apply: delta_res = der * delta_a_halfBarrel
    Eigen::Vector3d shift = der * mRigidBodyParams[sensorID];
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

  if (mOutOpt[o2::alignrs::OutputOpt::MisRes]) {
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

void AlignmentSpec::endOfStream(EndOfStreamContext& /*ec*/)
{
  mDBGOut->Close();
  mDBGOut.reset();
}

void AlignmentSpec::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    LOG(info) << "its cluster dictionary updated";
    mITSDict = (const o2::itsmft::TopologyDictionary*)obj;
    return;
  }
  if (matcher == ConcreteDataMatcher("IT3", "CLUSDICT", 0)) {
    LOG(info) << "it3 cluster dictionary updated";
    mIT3Dict = (const o2::its3::TopologyDictionary*)obj;
    return;
  }
}

DataProcessorSpec getAlignmentSpec(GTrackID::mask_t srcTracks, GTrackID::mask_t srcClusters, bool useMC, bool withITS, o2::alignrs::OutputEnum out)
{
  auto dataRequest = std::make_shared<DataRequest>();
  std::shared_ptr<o2::base::GRPGeomRequest> ggRequest{nullptr};
  if (!out[o2::alignrs::OutputOpt::MilleRes]) {
    dataRequest->requestTracks(srcTracks, useMC);
    if (!withITS) {
      dataRequest->requestIT3Clusters(useMC);
    } else {
      dataRequest->requestClusters(srcClusters, useMC);
    }
    dataRequest->requestPrimaryVertices(useMC);
    ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
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
    ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                           false,                             // GRPECS=true
                                                           false,                             // GRPLHCIF
                                                           false,                             // GRPMagField
                                                           false,                             // askMatLUT
                                                           o2::base::GRPGeomRequest::Aligned, // geometry
                                                           dataRequest->inputs);
  }

  Options opts{
    {"nthreads", VariantType::Int, 1, {"number of threads"}},
  };

  return DataProcessorSpec{
    .name = "its3-alignment",
    .inputs = dataRequest->inputs,
    .outputs = {},
    .algorithm = AlgorithmSpec{adaptFromTask<AlignmentSpec>(dataRequest, ggRequest, srcTracks, useMC, withITS, out)},
    .options = opts};
}
} // namespace o2::alignrs
