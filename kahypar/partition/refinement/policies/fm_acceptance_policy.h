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
#include "kahypar/partition/metrics.h"

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
      _ideal_block_weight(0),
      _final_block_weight(0) { }

    HypernodeWeight _step0_smallest_block_weight;
    HypernodeWeight _step0_heaviest_block_weight;
    HypernodeID _total_num_steps;
    HypernodeWeight _ideal_block_weight;
    HypernodeWeight _final_block_weight;
  
  public:
    virtual void init(HypernodeWeight step0_smallest_block_weight, HypernodeWeight step0_heaviest_block_weight, uint32_t total_num_steps, const Hypergraph& hg, const Context& context) {
      _step0_heaviest_block_weight = step0_heaviest_block_weight;
      _step0_smallest_block_weight = step0_smallest_block_weight;
      _total_num_steps = total_num_steps;
      double raw_ideal_weight = static_cast<double>(hg.totalWeight()) / static_cast<double>(context.partition.k);
      _ideal_block_weight = static_cast<HypernodeWeight>(raw_ideal_weight);
      _final_block_weight = raw_ideal_weight * (1.0 + context.partition.epsilon);
    }
  
};

class BalanceApproachingAcceptancePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;
  
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
      //use shifted and stretched arctan:
      // f(x) = a * arctan(b * (x - c)) + d
      // a and d are parameters deriving from the initial hypergraph imbalance at step0 (first refinement step, right after finish of contraction) and the final bound
      // b and c are hyperparameters tuneable via .ini files
      // x is the current step  in refinement phase
      double x = static_cast<double>(hypergraph.currentNumNodes() - context.partition.k);
      double a = static_cast<double>(_step0_heaviest_block_weight - _final_block_weight) / 2.0;
      double d = _final_block_weight + a;
      double b = -0.0005 * context.local_search.fm.balance_convergence_speed; //must be negative
      //for large large x >> c it should hold: f(x) == _final_block_weight
      double c = context.local_search.fm.balance_convergence_time * hypergraph.initialNumNodes(); //must be positive
      double tan = arctan(b * (x - c));
      double value = a * tan + d;
      ASSERT(value >= _final_block_weight, V(value));
      return std::max(value, -1.0);
    }
    
    HypernodeWeight currentLowerBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      return 2 * _ideal_block_weight - currentUpperBlockWeightBound(hypergraph, context);
    }

  private:
    double arctan(double x) {
      if (x > 705) { //otherwise exp(x) overflows
        return 1.0;
      }
      return (exp(x) - exp(-x)) / (exp(x) + exp(-x));
    }

};

class StaircaseAcceptancePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;
  private:
    BalanceApproachingAcceptancePolicy _balance_approaching_policy;

  public:
    StaircaseAcceptancePolicy() : Base(),
    _balance_approaching_policy() { }

    void init(HypernodeWeight step0_smallest_block_weight, HypernodeWeight step0_heaviest_block_weight, uint32_t total_num_steps, const Hypergraph& hg, const Context& context) override {
      _step0_heaviest_block_weight = step0_heaviest_block_weight;
      _step0_smallest_block_weight = step0_smallest_block_weight;
      _total_num_steps = total_num_steps;
      double raw_ideal_weight = static_cast<double>(hg.totalWeight()) / static_cast<double>(context.partition.k);
      _ideal_block_weight = static_cast<HypernodeWeight>(raw_ideal_weight);
      _final_block_weight = raw_ideal_weight * (1.0 + context.partition.epsilon);
      _balance_approaching_policy.init(step0_smallest_block_weight, step0_heaviest_block_weight, total_num_steps, hg, context);
    }

    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      HypernodeID current_step = hypergraph.currentNumNodes() - context.partition.k;
      if (current_step < 50) { //TODO(fritsch) magic number
        return _balance_approaching_policy.currentUpperBlockWeightBound(hypergraph, context);
      }
      return std::max(round(_balance_approaching_policy.currentUpperBlockWeightBound(hypergraph, context), context.local_search.flow.rounding_zeta), _final_block_weight);
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



class HeaviestNodePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;

  public:
    HeaviestNodePolicy() : Base() { }

    //TODO(fritsch) heaviest node weight is already maintained in coarsener base and calculated again
    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      return std::max((1.0 + context.partition.epsilon) * static_cast<double>(hypergraph.totalWeight()) / static_cast<double>(context.partition.k),
                                                    hypergraph.weightOfHeaviestNode() + static_cast<double>(hypergraph.totalWeight()) / static_cast<double>(context.partition.k));
    }

    HypernodeWeight currentLowerBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      return 2 * _ideal_block_weight - currentUpperBlockWeightBound(hypergraph, context);
    }
};

using FlowAcceptancePolicyClasses = meta::Typelist< BalanceApproachingAcceptancePolicy, 
                                                ImbalanceHoldingAcceptancePolicy, 
                                                StaircaseAcceptancePolicy,
                                                HeaviestNodePolicy>;
} // namespace kahypar

