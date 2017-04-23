/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2017 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include "kahypar/definitions.h"


namespace kahypar {
class MultiplicativePenalty {
 public:
  static HypernodeWeight penalty(const HypernodeWeight weight_u, const HypernodeWeight weight_v) {
    return weight_u * weight_v;
  }

  MultiplicativePenalty(const MultiplicativePenalty&) = delete;
  MultiplicativePenalty& operator= (const MultiplicativePenalty&) = delete;

  MultiplicativePenalty(MultiplicativePenalty&&) = delete;
  MultiplicativePenalty& operator= (MultiplicativePenalty&&) = delete;

 protected:
  ~MultiplicativePenalty() = default;
};

class NoWeightPenalty {
 public:
  static HypernodeWeight penalty(const HypernodeWeight, const HypernodeWeight) {
    return 1;
  }

  NoWeightPenalty(const NoWeightPenalty&) = delete;
  NoWeightPenalty& operator= (const NoWeightPenalty&) = delete;

  NoWeightPenalty(NoWeightPenalty&&) = delete;
  NoWeightPenalty& operator= (NoWeightPenalty&&) = delete;

 protected:
  ~NoWeightPenalty() = default;
};
}  // namespace kahypar