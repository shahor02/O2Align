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

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include <pwd.h>
#include <unistd.h>

#include <nlohmann/json.hpp>

#include "Framework/Logger.h"
#include "O2Align/TimeSlotsSet.h"

namespace o2::alignrs
{

namespace
{
// expand a leading "~" or "~user" and make the path absolute (relative paths are resolved wrt the current directory)
std::filesystem::path expandPath(const std::string& fileName)
{
  std::string path = fileName;
  if (!path.empty() && path[0] == '~') {
    auto slash = path.find('/');
    std::string user = path.substr(1, slash == std::string::npos ? std::string::npos : slash - 1);
    std::string home;
    if (user.empty()) { // "~" or "~/..." : current user
      if (const char* env = std::getenv("HOME"); env && env[0]) {
        home = env;
      } else if (const auto* pw = getpwuid(getuid()); pw && pw->pw_dir) {
        home = pw->pw_dir;
      }
    } else if (const auto* pw = getpwnam(user.c_str()); pw && pw->pw_dir) { // "~user/..."
      home = pw->pw_dir;
    }
    if (home.empty()) {
      LOGP(fatal, "Cannot resolve home directory for path {}", fileName);
    }
    path = home + (slash == std::string::npos ? "" : path.substr(slash));
  }
  std::error_code ec;
  auto abs = std::filesystem::weakly_canonical(std::filesystem::absolute(path, ec), ec);
  return ec ? std::filesystem::path{path} : abs; // fall back to the expanded path if resolution failed
}
} // namespace

int TimeSlotsSet::readSlotsFromFile(const std::string& fileName)
{
  using json = nlohmann::json;
  slots.clear();
  const auto path = expandPath(fileName);
  std::ifstream inpf(path);
  if (!inpf.good()) {
    LOGP(fatal, "Cannot open time slots file {} (resolved to {})", fileName, path.native());
  }
  auto data = json::parse(inpf);
  if (!data.is_array()) {
    LOGP(fatal, "Time slots file {} must contain an array of slots", path.native());
  }
  slots.reserve(data.size());
  for (const auto& entry : data) {
    auto& slot = slots.emplace_back();
    slot.runNumber = entry.value("run", -1);
    slot.intervalID = entry.value("ID", -1);
    slot.timeStampS = entry.value("tsS", 0L);
    slot.timeStampE = entry.value("tsE", 0L);
    if (slot.timeStampE <= slot.timeStampS) {
      LOGP(fatal, "Time slot {} of run {} has invalid time range {}:{}", slot.intervalID, slot.runNumber, slot.timeStampS, slot.timeStampE);
    }
    if (slots.size() > 1 && slots[slots.size() - 2].timeStampE > slot.timeStampS) {
      LOGP(fatal, "Time slot {} starting at {} overlaps with the preceding slot {} ending at {}",
           slot.intervalID, slot.timeStampS, slots[slots.size() - 2].intervalID, slots[slots.size() - 2].timeStampE);
    }
  }
  LOGP(info, "Read {} time slots from {}", slots.size(), path.native());
  return static_cast<int>(slots.size());
}


const TimeSlot& TimeSlotsSet::getSlot(long timestamp) const
{
  static const TimeSlot emptySlot;
  if (slots.empty()) {
    return emptySlot;
  }
  // 1st slot starting after the timestamp: the timestamp can be covered only by the preceding one
  auto it = std::upper_bound(slots.begin(), slots.end(), timestamp, [](long ts, const TimeSlot& slot) {
    return ts < slot.timeStampS;
  });
  if (it != slots.begin()) {
    const auto& slot = *(--it);
    if (timestamp < slot.timeStampE) {
      return slot;
    }
  }
  LOGP(warn, "Did not find time slot for timestamp {}", timestamp);
  return emptySlot;
}

int TimeSlotsSet::getSlotID(long timestamp) const
{
  if (slots.empty()) {
    return 0;
  }
  return getSlot(timestamp).intervalID;
}

} // namespace o2::alignrs
