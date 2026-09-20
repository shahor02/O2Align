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

#ifndef O2_ALIGN_TIME_SLOTS_SET_H
#define O2_ALIGN_TIME_SLOTS_SET_H

#include <string>
#include <vector>

#include <Rtypes.h>

namespace o2::alignrs
{

struct TimeSlot {
  int runNumber{-1};
  int intervalID{-1};
  long timeStampS{0};
  long timeStampE{0};
  ClassDefNV(TimeSlot, 1);
};

struct TimeSlotsSet {
  std::vector<TimeSlot> slots; // slots ordered by timeStampS, with non-overlapping [timeStampS, timeStampE) ranges

  /// Fill `slots` from a json file with an array of {"run", "ID", "tsS", "tsE"} records.
  /// Returns the number of slots read.
  int readSlotsFromFile(const std::string& fileName);

  /// Return the slot with timeStampS <= timestamp < timeStampE,
  /// an empty slot if no slots are defined or no matching slot is found
  const TimeSlot& getSlot(long timestamp) const;

  /// Return the intervalID of the slot with timeStampS <= timestamp < timeStampE,
  /// 0 if no slots are defined, -1 if the timestamp is not covered by any slot.
  int getSlotID(long timestamp) const;

  ClassDefNV(TimeSlotsSet, 1);
};

} // namespace o2::alignrs

#endif
