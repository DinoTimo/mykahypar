#pragma once

#include <limits>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/meta/template_parameter_to_string.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/datastructure/binary_heap.h"
#include "kahypar/partition/refinement/i_refiner.h"
#include "kahypar/partition/refinement/move.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/utils/randomize.h"

namespace kahypar {

class Rebalancer {
  static constexpr bool debug = false;
  constexpr static PartitionID _invalid_part = -1;
  const static Gain _invalid_gain = std::numeric_limits<Gain>::min();
  using Queue = ds::BinaryMinMaxHeap<HypernodeWeight, HypernodeWeight>;
  private:

  public:
    Rebalancer(Hypergraph& hg, const Context& context) :
    _hg(hg),
    _context(context),
    k(context.partition.k),
    _queues(),
    _queue_weights(k, 0),
    _current_upper_bound(0),
    _adjacency_bitmap(k * k, false), //too large
    _node_to_part_map(hg.initialNumNodes(), _invalid_part) { }

    ~Rebalancer() = default;

    Rebalancer(const Rebalancer&) = delete;
    Rebalancer& operator= (const Rebalancer&) = delete;

    Rebalancer(Rebalancer&&) = default;
    Rebalancer& operator= (Rebalancer&&) = default;

    void initialize() {
      for (PartitionID part = 0; part < k; part++) {
        _queues.emplace_back(_hg.initialNumNodes());
      }
    }

    void rebalance(HypernodeWeight heaviest_node_weight, IRefiner& refiner, Metrics& current_metrics, std::vector<HypernodeID>& refinement_nodes, HypernodeWeight current_upper_bound) {
      _current_upper_bound = current_upper_bound;
      reset();
      // ------------------------------
      // calculate connected blocks
      // ------------------------------
      for (HyperedgeID edge : _hg.edges()) {
        //maybe optimize for low current num nodes, since its expected to be fully connected
        for (PartitionID block1 : _hg.connectivitySet(edge)) {
          for (PartitionID block2 : _hg.connectivitySet(edge)) {
            if (block1 == block2) continue;
            setBlocksAdjacent(block1, block2);
          }
        }
      }
      // ------------------------------
      // fill binary heaps
      // ------------------------------
      for (HypernodeID node : _hg.nodes()) {
        ASSERT(_hg.nodeIsEnabled(node));
        const PartitionID from_part = _hg.partID(node);
        if (_hg.partWeight(from_part) <= _current_upper_bound) continue;
        const std::pair<bool, std::pair<PartitionID, Gain>> relativeGainMove = highestGainMoveToNotOverloadedBlock(node);
        const PartitionID to_part = relativeGainMove.second.first;
        const Gain relativeGain = relativeGainMove.second.second;
        if (!relativeGainMove.first) {
          DBG << "no move found for node " << node << " from part " << from_part;
          ASSERT(relativeGain == _invalid_gain && to_part == _invalid_part);
          continue;
        } else {
          DBG0 << "found move for node " << node << " from part " << from_part << " to part " << to_part << " with gain " << relativeGain;
        }
        //try to insert move into
        if (_queue_weights.at(from_part) < excessWeight(from_part)) {
          insertIntoQ(node, to_part, relativeGain);
        } else if (relativeGain > _queues.at(from_part).minKey()) {
          insertIntoQ(node, to_part, relativeGain);
          if (_queue_weights.at(from_part) > excessWeight(from_part) + heaviest_node_weight) {
            HypernodeID removedNode = _queues.at(from_part).min();
            ASSERT(_node_to_part_map.at(removedNode) != _invalid_part);
            _queues.at(from_part).popMin();
            _node_to_part_map.at(removedNode) = _invalid_part;
            _queue_weights.at(from_part) -= _hg.nodeWeight(removedNode);
            ASSERT(_node_to_part_map.at(removedNode) != from_part);
            ASSERT(!_queues.at(from_part).contains(removedNode));
          }
        }
      }
      // ------------------------------
      // empty queues
      // ------------------------------
      std::vector<Move> moves;
      emptyQueues(moves);
      // ------------------------------
      // rollback moves
      // ------------------------------
      UncontractionGainChanges emptyGains;
      for (const Move& move : moves) {
        _hg.changeNodePart(move.hn, move.to, move.from);
      }
      if (!moves.empty()) {
        refiner.performMovesAndUpdateCache(moves, refinement_nodes, emptyGains);
        current_metrics.heaviest_block_weight = metrics::heaviest_block_weight(_hg);
        current_metrics.km1 = metrics::km1(_hg);
      } else {
        DBG << "No moves done";
      }
    }


  private:
    void emptyQueues(std::vector<Move>& moves) {
      std::vector<PartitionID> overloaded_blocks;
      for (PartitionID block = 0; block < k; block++) {
        if (_hg.partWeight(block) > _current_upper_bound) {
          overloaded_blocks.push_back(block);
        }
      }
      DBG << "Overloaded blocks:" << joinVector(overloaded_blocks, "{", ",", "}");
      ASSERT(!overloaded_blocks.empty());
      bool movedAnythingThisIteration = true;
      while (movedAnythingThisIteration) {
        updateOverloadedBlocks(overloaded_blocks);
        movedAnythingThisIteration = false;
        //block order policy
        //for now implement round robin in arbitrary order
        for (PartitionID block : overloaded_blocks) {
          movedAnythingThisIteration = movedAnythingThisIteration || tryToMoveOutOfBlock(block, moves);
        }
      }
    }

    void updateOverloadedBlocks(std::vector<PartitionID>& overloaded_blocks) const {
      overloaded_blocks.erase(std::remove_if(overloaded_blocks.begin(), overloaded_blocks.end(), [&](PartitionID& block) {
        if (_hg.partWeight(block) <= _current_upper_bound) {
          DBG << V(block) << "has been balanced";
          return true;
        }
        return false;
      }), overloaded_blocks.end());
    }

    bool tryToMoveOutOfBlock(const PartitionID& block, std::vector<Move>& moves) {
      if (_queues.at(block).empty()) {
        DBG << "No moves available for" << V(block) << " but still has weight " << _hg.partWeight(block);
        return false;
      }
      if (excessWeight(block) <= 0) {
        DBG << V(block) << " has been balanced and now weighs " << _hg.partWeight(block);
        return false;
      }
      const HypernodeID node = _queues.at(block).max();
      const Gain relativeGain = _queues.at(block).maxKey();
      const PartitionID to_part = _node_to_part_map.at(node);
      ASSERT(block == _hg.partID(node), V(block) << V(_hg.partID(node)) << V(node));
      ASSERT(to_part != _invalid_part);
      _queues.at(block).popMax();
      _node_to_part_map.at(node) = _invalid_part;
      ASSERT(!_queues.at(block).contains(node));
      _queue_weights.at(block) -= _hg.nodeWeight(node); 
      if (_hg.partWeight(to_part) + _hg.nodeWeight(node) > _current_upper_bound
      || gainChangedFor(node, to_part, relativeGain))  {
        if (_hg.partWeight(to_part) + _hg.nodeWeight(node) > _current_upper_bound) {
          DBG << V(to_part) << " would become overloaded if " << V(node) << "[" << _hg.nodeWeight(node) << "] was moved to it";
        } else {
          DBG << "gain changed for node " << node << " from part " << block << " to part " << to_part;
        }
        if (_hg.isBorderNode(node) || _queue_weights.at(block) < excessWeight(block)) {
          std::pair<bool, std::pair<PartitionID, Gain>> newMove = highestGainMoveToNotOverloadedBlock(node);
          PartitionID to_part = newMove.second.first;
          Gain newRelativeGain = newMove.second.second;
          if (!newMove.first) {
            ASSERT(newRelativeGain == _invalid_gain && to_part == _invalid_part);
          } else {
            insertIntoQ(node, to_part, newRelativeGain);
          }
          //Try to insert all neighbours from former block
          for (const HyperedgeID& edge : _hg.incidentEdges(node)) {
            for (const HypernodeID& neighbour : _hg.pins(edge)) {
              if (!_hg.active(neighbour) || neighbour == node || _hg.partID(neighbour) != block) {
                tryToInsertIntoCorrectQ(node);
              }
            }
          }
        }
      } else {
        DBG0 << "hypernode " << node << "[" << _hg.nodeWeight(node) << "] is moved from part " << block << "[" << _hg.partWeight(block) << "] to part " << to_part << "[" << _hg.partWeight(to_part) << "]";
        ASSERT(_hg.nodeIsEnabled(node));
        _hg.changeNodePart(node, block, to_part);
        moves.push_back(Move{node, block, to_part});
        ASSERT(!_queues.at(block).contains(node));
      }
      return true;
    }

    inline bool gainChangedFor(const HypernodeID node, const PartitionID to_part, const Gain relativeGain) const {
      Gain actualGain = gainInducedByHypergraph(node, to_part);
      Gain actualRelativeGain = (actualGain >= 0) ? actualGain * _hg.nodeWeight(node) : actualGain / _hg.nodeWeight(node);
      return actualRelativeGain != relativeGain;
    }

    inline bool tryToInsertIntoCorrectQ(HypernodeID node) {
      const std::pair<bool, std::pair<PartitionID, Gain>> relativeGainMove = highestGainMoveToNotOverloadedBlock(node);
      if (relativeGainMove.first) {
        const PartitionID to_part = relativeGainMove.second.first;
        const Gain relativeGain = relativeGainMove.second.second;
        return tryToInsertIntoQ(node, to_part, relativeGain);
      }
      return false;
     }

    inline bool tryToInsertIntoQ(HypernodeID node, PartitionID to_part, Gain relativeGain) {
      PartitionID from_part = _hg.partID(node);
      ASSERT(to_part != from_part);
      if (!_queues.at(from_part).contains(node)) {
        return insertIntoQ(node, to_part, relativeGain);
      }
      return false;
    }

    inline bool insertIntoQ(HypernodeID node, PartitionID to_part, Gain relativeGain) {
      if (_node_to_part_map.at(node) != _invalid_part) {
        return false;
      }
      PartitionID from_part = _hg.partID(node);
      ASSERT(to_part != from_part);
      ASSERT([&]() {
        int q_index = 0;
        for (const Queue& q : _queues) {
          ASSERT(!q.contains(node), V(node) << ", " << V(to_part) << ", " << V(q_index) << ", " << V(from_part));
          q_index++;
        }
        return true;
      } (), "Queues invalid");
      _node_to_part_map.at(node) = to_part;
      _queues.at(from_part).push(node, relativeGain);
      _queue_weights.at(from_part) += _hg.nodeWeight(node);
      return true;
    }

    inline void reset() {
      std::fill(_queue_weights.begin(), _queue_weights.end(), 0);
      std::fill(_adjacency_bitmap.begin(), _adjacency_bitmap.end(), false);
      std::fill(_node_to_part_map.begin(), _node_to_part_map.end(), _invalid_part);
      for (Queue&  q : _queues) {
        q.clear();
      }
    }

    std::pair<bool, std::pair<HypernodeID, Gain>> highestGainMoveToNotOverloadedBlock(HypernodeID node) const {
      return highestGainMoveToNotOverloadedBlock(node, false);
    }

    void setBlocksAdjacent(PartitionID block1, PartitionID block2) {
      _adjacency_bitmap.at(block1 * k + block2) = true;
      _adjacency_bitmap.at(block2 * k + block1) = true;
    }

    bool blockAdjacency(const PartitionID block1, const PartitionID block2) const {
      ASSERT(_adjacency_bitmap.at(block2 * k + block1) == _adjacency_bitmap.at(block1 * k + block2));
      return _adjacency_bitmap.at(block2 * k + block1);
    }

    std::pair<bool, std::pair<HypernodeID, Gain>> highestGainMoveToNotOverloadedBlock(HypernodeID node, bool ignore_neighbourhood) const {
      const PartitionID from_part = _hg.partID(node);
      Gain max_gain = _invalid_gain;
      PartitionID max_to_part = _invalid_part;
      bool foundAnyMove = false;
      for (PartitionID to_part = 0; to_part < k; to_part++) {
        if (from_part == to_part
        || (!blockAdjacency(from_part, to_part) && !ignore_neighbourhood)
        || excessWeight(from_part) <= 0
        || _hg.partWeight(to_part) + _hg.nodeWeight(node) > _current_upper_bound) {
          continue;
        }
        Gain newGain = gainInducedByHypergraph(node, to_part);
        if (!foundAnyMove) {
          max_gain = newGain;
          max_to_part = to_part;
          foundAnyMove = true;
        } else if (newGain == max_gain) {
          bool coin = Randomize::instance().flipCoin();
          if (coin) {
            max_gain = newGain;
            max_to_part = to_part;
          }
        } else if (newGain > max_gain) {
          max_gain = newGain;
          max_to_part = to_part;
        }
      }
      if (!foundAnyMove && !ignore_neighbourhood) {
        return highestGainMoveToNotOverloadedBlock(node, true);
      } 
      ASSERT((max_to_part == _invalid_part) == (max_gain == _invalid_gain));
      ASSERT((max_to_part < k || max_to_part == _invalid_part) && max_to_part != from_part);
      ASSERT(max_to_part == _invalid_part || _hg.partWeight(max_to_part) + _hg.nodeWeight(node) <= _current_upper_bound);
      if (foundAnyMove) {
        ASSERT(max_gain != _invalid_gain);
        max_gain = (max_gain >= 0) ? max_gain * _hg.nodeWeight(node) : max_gain / _hg.nodeWeight(node);
      }
      return std::make_pair(foundAnyMove, std::make_pair(max_to_part, max_gain));
    }

    Gain gainInducedByHypergraph(const HypernodeID hn, const PartitionID target_part) const {
      ASSERT(target_part != _hg.partID(hn), V(hn) << V(target_part));
      Gain gain = 0;
      for (const HyperedgeID& he : _hg.incidentEdges(hn)) {
        ASSERT(_hg.edgeSize(he) > 1, V(he));
        gain += gainInducedByHyperedge(hn, he, target_part);
      }
      return gain;
    }

    Gain gainInducedByHyperedge(const HypernodeID hn, const HyperedgeID he,
                              const PartitionID target_part) const {
      const HypernodeID pins_in_source_part = _hg.pinCountInPart(he, _hg.partID(hn));
      const HypernodeID pins_in_target_part = _hg.pinCountInPart(he, target_part);
      const HyperedgeWeight he_weight = _hg.edgeWeight(he);
      Gain gain = pins_in_source_part == 1 ? he_weight : 0;
      gain -= pins_in_target_part == 0 ? he_weight : 0;
      return gain;
    }

    HypernodeWeight excessWeight(PartitionID block) const {
      return _hg.partWeight(block) - _current_upper_bound;
    }

    Hypergraph& _hg;
    const Context& _context;
    PartitionID k;
    std::vector<Queue> _queues;
    std::vector<HypernodeWeight> _queue_weights;
    HypernodeWeight _current_upper_bound;
    std::vector<bool> _adjacency_bitmap;
    std::vector<PartitionID> _node_to_part_map;
};
}