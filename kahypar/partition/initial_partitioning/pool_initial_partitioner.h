/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2015 Tobias Heuer <tobias.heuer@gmx.net>
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

#include <limits>
#include <string>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/partition/factories.h"
#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/partition/partitioner.h"
#include "kahypar/utils/randomize.h"

namespace kahypar {
class PoolInitialPartitioner : public IInitialPartitioner,
                               private InitialPartitionerBase<PoolInitialPartitioner>{
  using Base = InitialPartitionerBase<PoolInitialPartitioner>;
  friend Base;

 private:
  static constexpr bool debug = false;

  static constexpr HyperedgeWeight kInvalidCut = std::numeric_limits<HyperedgeWeight>::max();
  static constexpr double kInvalidImbalance = std::numeric_limits<double>::max();

  class PartitioningResult {
 public:
    PartitioningResult(InitialPartitionerAlgorithm algo, Objective objective,
                       HyperedgeWeight quality, double imbalance) :
      algo(algo),
      objective(objective),
      quality(quality),
      imbalance(imbalance) { }

    void print_result(const std::string& desc) const {
      if (objective == Objective::cut) {
        LOG << desc << "=" << "[ Cut=" << quality << "- Imbalance=" << imbalance << "- Algorithm="
            << algo << "]";
      } else {
        LOG << desc << "=" << "[ Km1=" << quality << "- Imbalance=" << imbalance << "- Algorithm="
            << algo << "]";
      }
    }

    InitialPartitionerAlgorithm algo;
    Objective objective;
    HyperedgeWeight quality;
    double imbalance;
  };

 public:
  PoolInitialPartitioner(Hypergraph& hypergraph, Context& context) :
    Base(hypergraph, context),
    _partitioner_pool() {
    // mix3 => pool_type = 1011110110111_{2} = 6071_{10}
    // Set bits in pool_type decides which partitioner is executed
    _partitioner_pool.push_back(InitialPartitionerAlgorithm::bin_packing);  // 13th bit set to 1
    _partitioner_pool.push_back(InitialPartitionerAlgorithm::greedy_global);  // 12th bit set to 1
    _partitioner_pool.push_back(InitialPartitionerAlgorithm::greedy_round);  // 11th bit set to 1
    _partitioner_pool.push_back(
      InitialPartitionerAlgorithm::greedy_sequential);  // 10th bit set to 1
    _partitioner_pool.push_back(
      InitialPartitionerAlgorithm::greedy_global_maxpin);  // 9th bit set to 1
    _partitioner_pool.push_back(
      InitialPartitionerAlgorithm::greedy_round_maxpin);  // 8th bit set to 1
    _partitioner_pool.push_back(
      InitialPartitionerAlgorithm::greedy_sequential_maxpin);  // 7th bit set to 1
    _partitioner_pool.push_back(
      InitialPartitionerAlgorithm::greedy_global_maxnet);  // 6th bit set to 1
    _partitioner_pool.push_back(
      InitialPartitionerAlgorithm::greedy_round_maxnet);  // 5th bit set to 1
    _partitioner_pool.push_back(
      InitialPartitionerAlgorithm::greedy_sequential_maxnet);  // 4th bit set to 1
    _partitioner_pool.push_back(InitialPartitionerAlgorithm::lp);  // 3th bit set to 1
    _partitioner_pool.push_back(InitialPartitionerAlgorithm::bfs);  // 2th bit set to 1
    _partitioner_pool.push_back(InitialPartitionerAlgorithm::random);  // 1th bit set to 1
  }

  ~PoolInitialPartitioner() override = default;

  PoolInitialPartitioner(const PoolInitialPartitioner&) = delete;
  PoolInitialPartitioner& operator= (const PoolInitialPartitioner&) = delete;

  PoolInitialPartitioner(PoolInitialPartitioner&&) = delete;
  PoolInitialPartitioner& operator= (PoolInitialPartitioner&&) = delete;

 private:
  void partitionImpl() override final {
    Base::multipleRunsInitialPartitioning();
  }

  void initialPartition() {
    Objective obj = _context.partition.objective;
    PartitioningResult best_cut(InitialPartitionerAlgorithm::pool, obj, kInvalidCut, kInvalidImbalance);
    PartitioningResult min_cut(InitialPartitionerAlgorithm::pool, obj, kInvalidCut, 0.0);
    PartitioningResult max_cut(InitialPartitionerAlgorithm::pool, obj, -1, 0.0);
    PartitioningResult min_imbalance(InitialPartitionerAlgorithm::pool, obj, kInvalidCut,
                                     kInvalidImbalance);
    PartitioningResult max_imbalance(InitialPartitionerAlgorithm::pool, obj, kInvalidCut, -0.1);

    std::vector<PartitionID> best_partition(_hg.initialNumNodes());
    double original_epsilon = _context.partition.epsilon;
    modifyEpsilon();
    _context.stats.add(StatTag::InitialPartitioning, "initial_epsilon", _context.partition.epsilon);
    unsigned int n = _partitioner_pool.size() - 1;
    for (unsigned int i = 0; i <= n; ++i) {
      // If the (n-i)th bit of pool_type is set we execute the corresponding
      // initial partitioner (see constructor)
      if (!((_context.initial_partitioning.pool_type >> (n - i)) & 1)) {
        continue;
      }
      InitialPartitionerAlgorithm algo = _partitioner_pool[i];
      if (algo == InitialPartitionerAlgorithm::greedy_round_maxpin ||
          algo == InitialPartitionerAlgorithm::greedy_global_maxpin ||
          algo == InitialPartitionerAlgorithm::greedy_sequential_maxpin) {
        DBG << "skipping maxpin";
        continue;
      }

      std::unique_ptr<IInitialPartitioner> partitioner(
        InitialPartitioningFactory::getInstance().createObject(algo, _hg, _context));
      partitioner->partition();
      HyperedgeWeight current_quality = obj == Objective::cut ?
                                        metrics::hyperedgeCut(_hg) : metrics::km1(_hg);
      double current_imbalance = metrics::imbalance(_hg, _context);
      DBG << algo << V(obj) << V(current_quality) << V(current_imbalance);

      const bool equal_metric = current_quality == best_cut.quality;
      const bool improved_metric = current_quality < best_cut.quality;
      const bool improved_imbalance = current_imbalance < best_cut.imbalance;
      const bool is_feasible_partition = current_imbalance <= _context.partition.epsilon;
      const bool is_best_cut_feasible_paritition = best_cut.imbalance <= _context.partition.epsilon;
      
      bool approve_new_partition = false;
      if (_context.initial_partitioning.modified_epsilon != ModifiedEpsilon::same &&
          _context.initial_partitioning.modified_epsilon != ModifiedEpsilon::UNDEFINED) {
        approve_new_partition = (improved_metric && is_feasible_partition) ||
                                (is_feasible_partition && !is_best_cut_feasible_paritition);
      } else {
        approve_new_partition = (improved_metric && (is_feasible_partition || improved_imbalance)) ||
                                (equal_metric && improved_imbalance) ||
                                (is_feasible_partition && !is_best_cut_feasible_paritition);
      }
      if (approve_new_partition) {
        for (const HypernodeID& hn : _hg.nodes()) {
          best_partition[hn] = _hg.partID(hn);
        }
        applyPartitioningResults(best_cut, current_quality, current_imbalance, algo);
      }
      if (current_quality < min_cut.quality) {
        applyPartitioningResults(min_cut, current_quality, current_imbalance, algo);
      }
      if (current_quality > max_cut.quality) {
        applyPartitioningResults(max_cut, current_quality, current_imbalance, algo);
      }
      if (current_imbalance < min_imbalance.imbalance) {
        applyPartitioningResults(min_imbalance, current_quality, current_imbalance, algo);
      }
      if (current_imbalance > max_imbalance.imbalance) {
        applyPartitioningResults(max_imbalance, current_quality, current_imbalance, algo);
      }
    }

    if (_context.initial_partitioning.verbose_output) {
      min_cut.print_result("Minimum Quality  ");
      max_cut.print_result("Maximum Quality  ");
      min_imbalance.print_result("Minimum Imbalance");
      max_imbalance.print_result("Maximum Imbalance");
      best_cut.print_result("==> Best Quality ");
    }

    const PartitionID unassigned_part = _context.initial_partitioning.unassigned_part;
    _context.initial_partitioning.unassigned_part = -1;
    Base::resetPartitioning();
    _context.initial_partitioning.unassigned_part = unassigned_part;
    for (const HypernodeID& hn : _hg.nodes()) {
      if (!_hg.isFixedVertex(hn)) {
        _hg.setNodePart(hn, best_partition[hn]);
      }
    }

    _hg.initializeNumCutHyperedges();

    ASSERT([&]() {
        for (const HypernodeID& hn : _hg.nodes()) {
          if (_hg.partID(hn) == -1) {
            return false;
          }
        }
        return true;
      } (), "There are unassigned hypernodes!");

    // Pool Partitioner executes each initial partitioner nruns times.
    // To prevent pool partitioner to execute himself nruns times, we
    // set the nruns parameter to 1.
    _context.initial_partitioning.nruns = 1;
    _context.partition.epsilon = original_epsilon;
  }

  void modifyEpsilon() {
    switch(_context.initial_partitioning.modified_epsilon) {
      case ModifiedEpsilon::custom: 
        _context.partition.epsilon = _context.initial_partitioning.custom_epsilon;
        break;
      case ModifiedEpsilon::same:
        break;
      case ModifiedEpsilon::node_avg: {
        //(1+e')avg = max((1+e)avg, avg+max) 
        double avg = ceil(static_cast<double>(_hg.totalWeight()) / _context.partition.k);
        double onePlusModifiedEpsilon = std::max((1.0 + _context.partition.epsilon) * avg, avg + _hg.weightOfHeaviestNode()) / avg;
        _context.partition.epsilon = onePlusModifiedEpsilon - 1;
        DBG0 << "use modified epsilon e' =" << (onePlusModifiedEpsilon - 1.0);
        break;
      }
      case ModifiedEpsilon::lpt: {
        //apply lpt, choose heaviest block weight as epsilon
        //ds::BinaryMinMaxHeap<PartitionID, HypernodeWeight> block_queue(_context.partition.k);
        std::vector<HypernodeWeight> block_queue(_context.partition.k, 0);
        std::vector<HypernodeID> descending_nodes(bin_packing::nodesInDescendingWeightOrder(_hg));
        /*for (PartitionID k = 0; k < _context.partition.k; k++) {
          block_queue.push(k, 0);
        }*/
        for (const HypernodeID& hn : descending_nodes) {
          ASSERT(_hg.nodeIsEnabled(hn));
          PartitionID min_block = 0;
          for (PartitionID k = 1; k < _context.partition.k; k++) {
            if (block_queue[k] < block_queue[min_block]) {
              min_block = k;
            }
          } 
          block_queue[min_block] += _hg.nodeWeight(hn);
        }
        HypernodeWeight heaviest_block_weight = block_queue[0];
        for (PartitionID k = 1; k < _context.partition.k; k++) {
          if (block_queue[k] > heaviest_block_weight) {
            heaviest_block_weight = block_queue[k];
          }
        }
        double avg = ceil(static_cast<double>(_hg.totalWeight()) / _context.partition.k);
        _context.partition.epsilon = (static_cast<double>(heaviest_block_weight) / avg) - 1.0;
        break;
      }
      case ModifiedEpsilon::UNDEFINED:
        LOG << "Undefined modified epsilon, using the same as before";
        break;
      // omit default case to trigger compiler warning for missing cases
      }
  }

  void applyPartitioningResults(PartitioningResult& result, const HyperedgeWeight quality,
                                const double imbalance,
                                const InitialPartitionerAlgorithm algo) const {
    result.quality = quality;
    result.imbalance = imbalance;
    result.algo = algo;
  }

  using Base::_hg;
  using Base::_context;
  std::vector<InitialPartitionerAlgorithm> _partitioner_pool;
};
}  // namespace kahypar
