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
    static constexpr bool debug = false;
    FlowAcceptancePolicy() : meta::PolicyBase(),
      _step0_smallest_block_weight(0),
      _step0_heaviest_block_weight(0),
      _step0_num_nodes(0),
      _total_num_steps(0),
      _ideal_block_weight(0),
      _final_upper_bound(0),
      _final_lower_bound(0) { }

    HypernodeWeight _step0_smallest_block_weight;
    HypernodeWeight _step0_heaviest_block_weight;
    HypernodeID _step0_num_nodes;
    HypernodeID _total_num_steps;
    HypernodeWeight _ideal_block_weight;
    HypernodeWeight _final_upper_bound;
    HypernodeWeight _final_lower_bound;
  
  public:
    virtual void init(HypernodeWeight step0_smallest_block_weight, HypernodeWeight step0_heaviest_block_weight, HypernodeID step0_num_nodes, const Hypergraph& hg, const Context& context) {
      _step0_smallest_block_weight = step0_smallest_block_weight;
      _step0_heaviest_block_weight = step0_heaviest_block_weight;
      _step0_num_nodes = step0_num_nodes;
      _total_num_steps = hg.initialNumNodes() - step0_num_nodes;
      double raw_ideal_weight = static_cast<double>(hg.totalWeight()) / static_cast<double>(context.partition.k);
      _ideal_block_weight = static_cast<HypernodeWeight>(raw_ideal_weight);
      _final_upper_bound = raw_ideal_weight * (1.0 + context.partition.epsilon);
    }
  
};

class BalanceApproachingAcceptancePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;
  
  public:
    BalanceApproachingAcceptancePolicy() : Base() { }
    
    /*
    Concept:
    Interpolate a curve that starts at y_0 = step0_heaviest_block_weight (heaviest block at start of refining)
    and ends at y_final = _final_upper_bound
    Use a polynomial f of degree at least 3 defined by balance convergence speed to interpolate
    x_0 = 0;
    x = _hg.currentNumNodes() - 0step_num_nodes;
    x_final = _hg.NumNodes() - 0step_num_nodes;
    x_convergence = x_final * (1 - balance_convergence_time)
    f(x_0) = step0_heaviest_block_weight
    f(x_final) = f(x_convergence) = final_upper_bound;
    f(x) = a * (x - b)^(convergence_speed) + c 
    and the other way round fpr lower bound
    */


    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      if (_step0_heaviest_block_weight <= _final_upper_bound) {
        return _final_upper_bound;
      }
      double x = static_cast<double>(hypergraph.currentNumNodes() - _step0_num_nodes);
      double b = static_cast<double>(hypergraph.initialNumNodes() - _step0_num_nodes) * (1 - context.local_search.fm.balance_convergence_time);
      if (x >= b) {
        return _final_upper_bound;
      }
      double c = static_cast<double>(_final_upper_bound);
      ASSERT(b != 0, "Divide by zero error");
      double a = (_step0_heaviest_block_weight - c) / std::pow(-b, context.local_search.fm.balance_convergence_speed);
      double value = std::pow(x - b, context.local_search.fm.balance_convergence_speed);
      double bound = a * value + c;
      ASSERT(bound >= _final_upper_bound);
      return bound;
    }
};


class ImbalanceHoldingAcceptancePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;

  public:
    ImbalanceHoldingAcceptancePolicy() : Base() { }

    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      if (_step0_heaviest_block_weight <= _final_upper_bound) {
        return _final_upper_bound;
      }
      //use shifted and stretched tanh:
      // f(x) = a * tanh(b * (x - c)) + d
      // a and d are parameters deriving from the initial hypergraph imbalance at step0 (first refinement step, right after finish of contraction) and the final bound
      // b and c are hyperparameters tuneable via .ini files
      // x is the current step  in refinement phase
      double x = static_cast<double>(hypergraph.currentNumNodes() - _step0_num_nodes);
      double b = -0.0005 * context.local_search.fm.balance_convergence_speed; //must be negative
      double c = context.local_search.fm.balance_convergence_time * hypergraph.initialNumNodes(); //must be positive
      double a = static_cast<double>(_step0_heaviest_block_weight - _final_upper_bound) / 2.0;
      double d = _final_upper_bound + a;
      //for large large x >> c it should hold: f(x) == _final_upper_bound
      double tan = tanh(b * (x - c));
      double value = a * tan + d;
      ASSERT(value >= _final_upper_bound, V(value));
      return value;
    }

  private:
    double tanh(double x) {
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
    std::vector<HypernodeWeight> _imbalance_steps;

  public:
    StaircaseAcceptancePolicy() : Base(),
    _balance_approaching_policy(),
    _imbalance_steps() { }

     void init(HypernodeWeight step0_smallest_block_weight, HypernodeWeight step0_heaviest_block_weight, HypernodeID step0_num_nodes, const Hypergraph& hg, const Context& context) override {
      Base::init(step0_smallest_block_weight, step0_heaviest_block_weight, step0_num_nodes, hg, context);
      _balance_approaching_policy.init(step0_smallest_block_weight, step0_heaviest_block_weight, step0_num_nodes, hg, context);
      _imbalance_steps.resize(context.local_search.fm.num_staircase_steps, 0.0);
      double rounding_zeta = static_cast<double>(_step0_heaviest_block_weight - _final_upper_bound) / static_cast<double>(context.local_search.fm.num_staircase_steps);
      size_t i = 0;
      _imbalance_steps[i++] = _final_upper_bound;
      for (; i < context.local_search.fm.num_staircase_steps - 1; i++) {
        _imbalance_steps[i] = static_cast<HypernodeWeight>(i * rounding_zeta + _final_upper_bound);
      }
      _imbalance_steps[i] = _step0_heaviest_block_weight;
      DBG << V(rounding_zeta) << V(_step0_heaviest_block_weight) << V(_final_upper_bound);
    }

    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      if (_step0_heaviest_block_weight <= _final_upper_bound) {
        return _final_upper_bound;
      }
      return closestEntry(_balance_approaching_policy.currentUpperBlockWeightBound(hypergraph, context), _imbalance_steps);
    }
  
  private:
    template <typename NumberType>
    static NumberType closestEntry(NumberType value, std::vector<NumberType> entries) {
      if (entries.empty()) {
        LOG << "Illegal call entries are empty" << V(value);
      }
      NumberType closest_entry = entries[0];
      NumberType min_diff = std::abs(value - closest_entry);
      for (NumberType entry : entries) {
        NumberType diff = std::abs(value - entry);
        if (diff < min_diff) {
          min_diff = diff;
          closest_entry = entry;
          //this can be implemented more efficiently using the information that the here used entries are sorted. we could break in this case
        }
      }
      return closest_entry;
    }
};



class HeaviestNodePolicy : public FlowAcceptancePolicy {
  using Base = FlowAcceptancePolicy;

  public:
    HeaviestNodePolicy() : Base() { }

    //TODO(fritsch) heaviest node weight is already maintained in coarsener base and calculated again
    HypernodeWeight currentUpperBlockWeightBound(Hypergraph& hypergraph, const Context& context) {
      return std::max(static_cast<double>(_final_upper_bound), std::max((1.0 + context.partition.epsilon) * static_cast<double>(hypergraph.totalWeight()) / static_cast<double>(context.partition.k),
                                                    hypergraph.weightOfHeaviestNode() + static_cast<double>(hypergraph.totalWeight()) / static_cast<double>(context.partition.k)));
    }
};

using FlowAcceptancePolicyClasses = meta::Typelist< BalanceApproachingAcceptancePolicy, 
                                                ImbalanceHoldingAcceptancePolicy, 
                                                StaircaseAcceptancePolicy,
                                                HeaviestNodePolicy>;
} // namespace kahypar

