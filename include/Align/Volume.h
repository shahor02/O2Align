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

#ifndef O2_ALIGN_VOLUME_H
#define O2_ALIGN_VOLUME_H

#include <memory>
#include <utility>
#include <vector>
#include <ostream>
#include <string>
#include <map>
#include <algorithm>

#include <Eigen/Dense>

#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>

#include "Align/Label.h"
#include "Align/DOFSet.h"

namespace o2::alignrs
{
using Matrix36 = Eigen::Matrix<double, 3, 6>;
using Matrix66 = Eigen::Matrix<double, 6, 6>;

class Volume
{
 public:
  using Ptr = std::unique_ptr<Volume>;
  using SensorMapping = std::map<Label, Volume*>;

  Volume(const Volume&) = delete;
  Volume(Volume&&) = delete;
  Volume& operator=(const Volume&) = delete;
  Volume& operator=(Volume&&) = delete;
  Volume(const char* symName, uint32_t label, uint32_t det, bool sens);
  Volume(const char* symName, Label label);
  virtual ~Volume() = default;

  static void applyDOFConfig(Volume* root, const std::string& jsonPath);  
  static void writeMillepedeResults(Volume* root, const std::string& milleResPath, const std::string& outJsonPath, const std::string& injectedJsonPath = "");

  
  void finalise(uint8_t level = 0);

  // steering file output
  void writeRigidBodyConstraints(std::ostream& os) const;
  void writeParameters(std::ostream& os) const;
  void writeTree(std::ostream& os, int indent = 0) const;

  // tree-like
  auto getLevel() const noexcept { return mLevel; }
  bool isRoot() const noexcept { return mParent == nullptr; }
  bool isLeaf() const noexcept { return mChildren.empty(); }
  
  template <class T = Volume>
    requires std::derived_from<T, Volume>
  Volume* addChild(const char* symName, uint32_t label, uint32_t det, bool sens)
  {
    auto c = std::make_unique<T>(symName, label, det, sens);
    return setParent(std::move(c));
  }

  template <class T = Volume>
    requires std::derived_from<T, Volume>
  Volume* addChild(const char* symName, Label lbl)
  {
    auto c = std::make_unique<T>(symName, lbl);
    return setParent(std::move(c));
  }

  // bfs traversal
  void traverse(const std::function<void(Volume*)>& visitor)
  {
    visitor(this);
    for (auto& c : mChildren) {
      c->traverse(visitor);
    }
  }

  std::string getSymName() const noexcept { return mSymName; }
  Label getLabel() const noexcept { return mLabel; }
  Volume* getParent() const { return mParent; }
  int getNChildren() const noexcept { return static_cast<int>(mChildren.size()); }

  // DOF management
  void setRigidBody(std::unique_ptr<DOFSet> rb) { mRigidBody = std::move(rb); }
  void setCalib(std::unique_ptr<DOFSet> cal) { mCalib = std::move(cal); }
  DOFSet* getRigidBody() const { return mRigidBody.get(); }
  DOFSet* getCalib() const { return mCalib.get(); }
  void setPseudo(bool p) noexcept { mIsPseudo = p; }
  bool isPseudo() const noexcept { return mIsPseudo; }
  void setSensorId(int id) noexcept { mSensorId = id; }
  int getSensorId() const noexcept { return mSensorId; }
  // true if this volume participates in the hierarchy (has DOFs or is pseudo)
  bool isActive() const noexcept { return mRigidBody != nullptr || mIsPseudo; }

  // transformation matrices
  virtual void defineMatrixL2G() {}
  virtual void defineMatrixT2L() {}
  virtual void computeJacobianL2T(const double* pos, Matrix66& jac) const {};
  const TGeoHMatrix& getL2P() const { return mL2P; }
  const TGeoHMatrix& getT2L() const { return mT2L; }
  const Matrix66& getJL2P() const { return mJL2P; }
  const Matrix66& getJP2L() const { return mJP2L; }

 protected:
  /// matrices
  Volume* mParent{nullptr}; // parent
  TGeoPNEntry* mPNE{nullptr};        // physical entry
  TGeoPhysicalNode* mPN{nullptr};    // physical node
  TGeoHMatrix mL2G;                  // (LOC) -> (GLO)
  TGeoHMatrix mL2P;                  // (LOC) -> (PAR)
  Matrix66 mJL2P;                    // jac (LOC) -> (PAR)
  Matrix66 mJP2L;                    // jac (PAR) -> (LOC)
  TGeoHMatrix mT2L;                  // (TRK) -> (LOC)

 private:
  std::string mSymName;
  Label mLabel;
  uint8_t mLevel{0};
  bool mIsPseudo{false};
  int mSensorId{-1}; // RS check if needed
  std::unique_ptr<DOFSet> mRigidBody;
  std::unique_ptr<DOFSet> mCalib;

  Volume* setParent(Ptr c)
  {
    c->mParent = this;
    mChildren.push_back(std::move(c));
    return mChildren.back().get();
  }
  std::vector<Ptr> mChildren; // children

  void init();
};

} // namespace o2::alignment

#endif
