// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "ITSMFTReconstruction/ChipMappingITS.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DataFormatsTRD/Constants.h"
#include "TRDBase/Geometry.h"
#include "TOFBase/Geo.h"
#include <TVector3.h>
#include <stdexcept>
#endif

using DetID = o2::detectors::DetID;

void algDump(const std::string& geom = "", const std::string& outname = "algdump.root", int NX = 5, int NZ = 3)
{
  if (NX <= 0 || NZ <= 0) {
    throw std::invalid_argument("algDump: NX and NZ must be positive");
  }
  o2::base::GeometryManager::loadGeometry(geom.c_str());
  o2::utils::TreeStreamRedirector outstream(outname.c_str(), "recreate");
  TGeoHMatrix* matAlg = nullptr;
  TGeoHMatrix matOrig;
  TVector3 pos0, pos;
  DetID det;
  auto procSens = [](TGeoHMatrix& mat, double x = 0., double z = 0.) {
    double loc[3] = {x, 0., z}, glo[3];
    mat.LocalToMaster(loc, glo);
    return TVector3(glo[0], glo[1], glo[2]);
  };

  auto store = [&outstream, &det, &pos, &pos0](int lr, int sid, int sidLr, const char* outmod = "") {
    TString str = Form("gm%s", outmod);
    outstream << str.Data()
              << "det=" << det.getID() << "lr=" << lr << "sid=" << sid << "sidlr=" << sidLr << "pos0=" << pos0 << "pos=" << pos << "\n";
    //    printf("xx %d %d %d %f %f %f\n", det.getID(), lr, sid, pos0[0], pos0[1], pos0[1]);
  };

  det = DetID("ITS");
  o2::itsmft::ChipMappingITS mpits;
  for (int ic = 0; ic < mpits.getNChips(); ic++) {
    int lr = mpits.getLayer(ic);
    int ic0 = ic - mpits.getFirstChipsOnLayer(lr);
    matAlg = o2::base::GeometryManager::getMatrix(det, ic);
    o2::base::GeometryManager::getOriginalMatrix(det, ic, matOrig);
    pos0 = procSens(matOrig);
    pos = procSens(*matAlg);
    store(lr, ic, ic0);
    // Sample cell centers so all mesh points stay inside the active matrix.
    for (int ix = 0; ix < NX; ix++) {
      const double x = ((ix + 0.5) / NX - 0.5) * o2::itsmft::SegmentationAlpide::ActiveMatrixSizeRows;
      for (int iz = 0; iz < NZ; iz++) {
        const double z = ((iz + 0.5) / NZ - 0.5) * o2::itsmft::SegmentationAlpide::ActiveMatrixSizeCols;
        pos0 = procSens(matOrig, x, z);
        pos = procSens(*matAlg, x, z);
        store(lr, ic, ic0, "e");
      }
    }
  }

  det = DetID("TRD");
  for (int ilr = 0; ilr < o2::trd::constants::NLAYER; ilr++) {                                 // layer
    for (int ich = 0; ich < o2::trd::constants::NSTACK * o2::trd::constants::NSECTOR; ich++) { // chamber
      int isector = ich / o2::trd::constants::NSTACK;
      int istack = ich % o2::trd::constants::NSTACK;
      uint16_t sid = o2::trd::Geometry::getDetector(ilr, istack, isector);
      const char* symname = Form("TRD/sm%02d/st%d/pl%d", isector, istack, ilr);
      if (!gGeoManager->GetAlignableEntry(symname)) {
        continue;
      }
      matAlg = o2::base::GeometryManager::getMatrix(det, sid);
      o2::base::GeometryManager::getOriginalMatrix(det, sid, matOrig);
      pos0 = procSens(matOrig);
      pos = procSens(*matAlg);
      store(ilr, sid, ich);
    }
  }

  det = DetID("TOF");
  int cnt = -1;
  for (int isc = 0; isc < 18; isc++) {
    for (int istr = 1; istr <= o2::tof::Geo::NSTRIPXSECTOR; istr++) { // strip
      const char* symname = Form("TOF/sm%02d/strip%02d", isc, istr);
      cnt++;
      if (!gGeoManager->GetAlignableEntry(symname)) {
        continue;
      }
      matAlg = o2::base::GeometryManager::getMatrix(det, cnt);
      o2::base::GeometryManager::getOriginalMatrix(det, cnt, matOrig);
      pos0 = procSens(matOrig);
      pos = procSens(*matAlg);
      store(0, cnt, cnt);
    }
  }
  outstream.Close();
}
