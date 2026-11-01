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

// Per IB sensor: id layer halfBarrel stave chipInStave refAlpha refX centre(x y z) dz_glo/dz_loc, from a geometry file.
// Input table (sens.txt) of macro/milleSplitModels.py.
// root -l -b -q 'sensDump.C+("o2sim_geometry-aligned.root","sens.txt")'
#include "ITSBase/GeometryTGeo.h"
#include "DetectorsBase/GeometryManager.h"
#include "MathUtils/Utils.h"
#include <fstream>
void sensDump(const char* geo, const char* out)
{
  o2::base::GeometryManager::loadGeometry(geo);
  auto g = o2::its::GeometryTGeo::Instance();
  g->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G, o2::math_utils::TransformType::T2G));
  std::ofstream f(out);
  for (int i = 0; i < 432; i++) {
    const auto& m = g->getMatrixL2G(i);
    auto glo = m(o2::math_utils::Point3D<float>(0, 0, 0));
    auto gp = m(o2::math_utils::Point3D<float>(0, 0, 1.0));
    f << i << " " << g->getLayer(i) << " " << g->getHalfBarrel(i) << " " << g->getStave(i) << " " << g->getChipIdInStave(i) << " "
      << g->getSensorRefAlpha(i) << " " << g->getSensorRefX(i) << " " << glo.X() << " " << glo.Y() << " " << glo.Z() << " " << (gp.Z() - glo.Z()) << "\n";
  }
}
