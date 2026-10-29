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

// Helpers shared between Detector::MP2AlignParams and the weak-modes diagnostics to turn CCDB
// o2::detectors::AlignParam objects into the cumulative global delta they represent relative to
// the ideal geometry, see the derivation in Detector.cxx.
#ifndef O2ALIGN_ALIGNGEOMUTILS_H
#define O2ALIGN_ALIGNGEOMUTILS_H

#include <string>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <istream>
#include <nlohmann/json.hpp>

#include <TGeoMatrix.h>

#include "DetectorsCommonDataFormats/AlignParam.h"

namespace o2::alignrs::geomutils
{

/// Is this JSON object key a comment: "#<anything>" or "_comment<anything>"
inline bool isJsonCommentKey(const std::string& key)
{
  return key.starts_with('#') || key.starts_with("_comment");
}

/// Parse a JSON config file, dropping every object member whose key is a comment (see isJsonCommentKey)
/// at any depth. An object left empty by that (e.g. an array element with a comment only) stays in
/// its array as {}, which the consumers skip.
inline nlohmann::json parseJsonWithComments(std::istream& in)
{
  return nlohmann::json::parse(in, [](int, nlohmann::json::parse_event_t event, nlohmann::json& parsed) {
    return !(event == nlohmann::json::parse_event_t::key && isJsonCommentKey(parsed.get<std::string>()));
  });
}

/// Read the list of rules of a DOF-config-style JSON file (dofConfigJson, misAlgJson): the file is either
/// an object with a "rules" array or a bare array. Comment members are dropped (parseJsonWithComments)
/// and so are the elements left empty by that. An object-level "defaults" is inserted as the first rule
/// with match ".*" if allowDefaults, else it is fatal. what names the file kind in the messages.
nlohmann::json readJsonRules(const std::string& path, const char* what, bool allowDefaults);

/// A DOF-config-style symbolic-volume pattern: std::regex_match against either the full name or
/// the name with an implicit arbitrary prefix (so a suffix pattern can match at any hierarchy depth).
class VolumePattern
{
 public:
  explicit VolumePattern(const std::string& pattern);
  bool matches(const std::string& symName) const;

 private:
  std::regex mDirect;
  std::regex mAnyPrefix;
};

/// matrices keyed by the TGeo path of the alignable volume they belong to
using PathMatrices = std::unordered_map<std::string, TGeoHMatrix>;

/// TGeo path of an alignable volume, used to navigate the TGeo hierarchy by path prefixes
std::string getAlignablePath(const std::string& symName);

/// TGeo level of a path, used to order volumes parents first
int getPathDepth(const std::string& path);

/// L2G matrix of an alignable volume in the loaded (pre-aligned) geometry
const TGeoHMatrix& getGeometryMatrix(const std::string& symName);

/// Matrix of the TGeo volume `path` (if includeSelf) or else of its nearest TGeo ancestor present
/// in the map, the identity if there is none.
TGeoHMatrix getOnPath(const PathMatrices& matrices, std::string path, bool includeSelf);

/// global delta D^p equivalent to an initial AlignParam of either convention
TGeoHMatrix getInitialGlobalDelta(const o2::detectors::AlignParam& par);

/// cumulative global delta P^p of the initial alignment for every volume it contains, the objects
/// being accumulated parents first, like GeometryManager::applyAlignment applies them. An empty
/// map (identity) for every volume if initial is null.
PathMatrices collectInitialCumulativeDeltas(const std::vector<o2::detectors::AlignParam>* initial);

/// alignable entry of the loaded geometry
struct Alignable {
  std::string symName;
  std::string path;
};

/// all alignable entries of the loaded geometry whose TGeo path is one of topPaths or lies below
/// one of them, ordered parents first
std::vector<Alignable> collectAlignablesBelow(const std::unordered_set<std::string>& topPaths);

/// indices of the alignable children (nearest alignable TGeo descendants) of every alignable having
/// some, keyed by the path of the parent
std::unordered_map<std::string, std::vector<size_t>> collectAlignableChildren(const std::vector<Alignable>& alignables);

/// selection of the parents for transferChildrenCommonMotion
struct ChildrenTransferConfig {
  std::vector<VolumePattern> parents; // symbolic names of the parents to process, all if empty
  int minChildren = 2;                // skip the parents with fewer alignable children
};

/// A posteriori regularization of an alignment: to every selected parent (deepest first, so that the
/// common part cascades up) transfer the rigid motion robustly (L1, the multidimensional counterpart of
/// the median) describing, in the parent frame, the displacements of the bounding box corners of its TGeo
/// alignable children, keeping every child, hence every sensor, at its final place. A common motion of
/// the children which no rigid motion of the parent produces (e.g. all staves turned around their own
/// axes) stays in the children.
/// \param cumulative cumulative global deltas P (G = P * G^id) of the volumes having an alignment
///                   object, modified in place: the parents and their children which had none get one
/// \param alignables  the volumes among which parents and children are looked for (parents first)
/// \param ideal       ideal L2G matrix G^id of every alignable, by path
void transferChildrenCommonMotion(PathMatrices& cumulative, const std::vector<Alignable>& alignables, const PathMatrices& ideal, const ChildrenTransferConfig& cfg);

} // namespace o2::alignrs::geomutils

#endif
