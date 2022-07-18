/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include "kahypar/meta/policy_registry.h"
#include "kahypar/macros.h"
#include "kahypar/meta/typelist.h"
#include "kahypar/partition/context.h"

#include <limits>
#include <stack>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <cmath>

namespace kahypar {
class FlowAcceptancePolicy : public meta::PolicyBase {
  protected:
    static constexpr bool debug = true;
    FlowAcceptancePolicy() : meta::PolicyBase(),
      _step0_smallest_block_weight(0),
      _step0_heaviest_block_weight(0),
      _total_num_steps(0),
      _ideal_block_weight(0) { }

    HypernodeWeight _step0_smallest_block_weight;
    HypernodeWeight _step0_heaviest_block_weight;
    HypernodeID _total_num_steps;
    HypernodeID _ideal_block_weight;
  
  public:
    virtual void init(HypernodeWeight step0_smallest_block_weight, HypernodeWeight step0_heaviest_block_weight, uint32_t total_num_steps, const Hypergraph& hg, const Context& context) {
      _step0_heaviest_block_weight = step0_heaviest_block_weight;
      _step0_smallest_block_weight = step0_smallest_block_weight;
      _total_num_steps = total_num_steps;
      _ideal_block_weight = hg.totalWeight() / context.partition.k;
    }
  
};

class BalanceApproachingAcceptancePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;
  //friend Base;

  public:
    BalanceApproachingAcceptancePolicy() : Base() { }

    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      HypernodeID current_step = hypergraph.currentNumNodes() - context.partition.k;
      double initialPureUpper = static_cast<double>(_ideal_block_weight + blockWeightDelta(context, 0));
      double initialImbalance = static_cast<double>(_step0_heaviest_block_weight);
      uint32_t current_pseudo_step = static_cast<uint32_t>(static_cast<double>(current_step) * (1 / (1 - context.local_search.fm.balance_convergence_time)));
      current_pseudo_step = std::min(current_pseudo_step, _total_num_steps);
      double exponent = 1 - ((static_cast<double>(current_pseudo_step)) / static_cast<double>((_total_num_steps)));
      double modifier = std::pow(initialImbalance / initialPureUpper, exponent);
      return static_cast<HypernodeWeight>((_ideal_block_weight + blockWeightDelta(context, current_step)) * modifier);
    }

    HypernodeWeight currentLowerBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      HypernodeID current_step = hypergraph.currentNumNodes() - context.partition.k;
      double initialPureLower = static_cast<double>(_ideal_block_weight - blockWeightDelta(context, 0));
      double initialImbalance = static_cast<double>(_step0_smallest_block_weight);
      uint32_t current_pseudo_step = static_cast<uint32_t>(static_cast<double>(current_step) * (1 / (1 - context.local_search.fm.balance_convergence_time)));
      current_pseudo_step = std::min(current_pseudo_step, _total_num_steps);
      double exponent = 1 - ((static_cast<double>(current_pseudo_step)) / static_cast<double>((_total_num_steps)));
      double modifier = std::pow(initialImbalance / initialPureLower, exponent);
      return static_cast<HypernodeWeight>((_ideal_block_weight - blockWeightDelta(context, current_step)) * modifier);
    }

  private:
    HypernodeWeight blockWeightDelta(const Context& context, HypernodeID current_step) {
      uint32_t current_pseudo_step = std::min(_total_num_steps, current_step + static_cast<uint32_t>(context.local_search.fm.balance_convergence_time * static_cast<double>(_total_num_steps)));
      uint32_t step_diff = _total_num_steps - current_pseudo_step;
      return static_cast<HypernodeWeight>(static_cast<double>(_ideal_block_weight)
        * std::pow((static_cast<double>(step_diff) / static_cast<double>(_total_num_steps)) + 1, context.local_search.fm.balance_convergence_speed)
        * context.partition.epsilon);
    }
};


class ImbalanceHoldingAcceptancePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;

  public:
    ImbalanceHoldingAcceptancePolicy() : Base() { }

    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      HypernodeID current_step = hypergraph.currentNumNodes() - context.partition.k;
      uint32_t imbalance_step_offset = static_cast<uint32_t>((context.local_search.fm.balance_convergence_time) * static_cast<double>(_total_num_steps));
      uint32_t balancing_start_step = _total_num_steps - 2 * imbalance_step_offset;
      if (current_step < balancing_start_step) {
        return _step0_heaviest_block_weight;
      }
      uint32_t goal_optimizing_step = _total_num_steps - imbalance_step_offset;
      if (current_step >= goal_optimizing_step) {
        return static_cast<HypernodeWeight>(static_cast<double>(_ideal_block_weight) * (1 + context.partition.epsilon));
      }
      //f(x) = mx + c
      double y0 = _step0_heaviest_block_weight;
      double x0 = balancing_start_step;
      double y1 = static_cast<double>(_ideal_block_weight) * (1 + context.partition.epsilon);
      double x1 = static_cast<double>(goal_optimizing_step);
      double x = static_cast<double>(current_step);
      double m = (y0 - y1) / (x0 - x1);
      double c = y0 - m * x0;
      return m * x + c;
    }
    
    HypernodeWeight currentLowerBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      return 2 * _ideal_block_weight - currentUpperBlockWeightBound(hypergraph, context);
    }

};

class StaircaseAcceptancePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;
  private:
    BalanceApproachingAcceptancePolicy _balance_approaching_policy;
    HypernodeWeight _final_upper_bound;

  public:
    StaircaseAcceptancePolicy() : Base(),
    _balance_approaching_policy(),
    _final_upper_bound(0) { }

    void init(HypernodeWeight step0_smallest_block_weight, HypernodeWeight step0_heaviest_block_weight, uint32_t total_num_steps, const Hypergraph& hg, const Context& context) override {
      _step0_heaviest_block_weight = step0_heaviest_block_weight;
      _step0_smallest_block_weight = step0_smallest_block_weight;
      _total_num_steps = total_num_steps;
      _ideal_block_weight = hg.totalWeight() / context.partition.k;
      _final_upper_bound = _ideal_block_weight * (context.partition.epsilon + 1);
      _balance_approaching_policy.init(step0_smallest_block_weight, step0_heaviest_block_weight, total_num_steps, hg, context);
    }

    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      return std::max(round(_balance_approaching_policy.currentUpperBlockWeightBound(hypergraph, context), context.local_search.flow.rounding_zeta), _final_upper_bound);
    }
    
    HypernodeWeight currentLowerBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      return 2 * _ideal_block_weight - currentUpperBlockWeightBound(hypergraph, context);
    }
  
  private:
    template <typename NumberType, typename RoundingType>
    static NumberType round(NumberType raw, RoundingType rounding_beta) {
      return static_cast<NumberType>(static_cast<int>(raw) / static_cast<int>(rounding_beta)) * rounding_beta;
    }
};

using FlowAcceptancePolicyClasses = meta::Typelist< BalanceApproachingAcceptancePolicy, 
                                                ImbalanceHoldingAcceptancePolicy, 
                                                StaircaseAcceptancePolicy>;
} // namespace kahypar

