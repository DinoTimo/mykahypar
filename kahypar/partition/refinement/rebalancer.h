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

namespace kahypar {

class Rebalancer {
  static constexpr bool debug = true;
  const static PartitionID _invalid_part = std::numeric_limits<PartitionID>::max();
  const static Gain _invalid_gain = std::numeric_limits<Gain>::min();
  using Queue = ds::BinaryMinMaxHeap<HypernodeWeight, HypernodeWeight>;
  private:
    struct NodeMove
    {
      HypernodeID node;
      PartitionID to_part;

      bool operator==(NodeMove& other) {
        return other.node == node && other.to_part == to_part;
      }
    };
    


  public:
    Rebalancer(Hypergraph& hg, const Context& context) :
    _hg(hg),
    _context(context),
    k(context.partition.k),
    _queues(),
    _queue_weights(k, 0),
    _heaviest_node_weight(0), //TODO(fritsch)
    _current_upper_bound(0),
    _adjacency_bitmap(k * (k - 1), false),
    _was_inserted_into_q_bitmap(hg.initialNumNodes()) { }

    ~Rebalancer() = default;

    Rebalancer(const Rebalancer&) = delete;
    Rebalancer& operator= (const Rebalancer&) = delete;

    Rebalancer(Rebalancer&&) = default;
    Rebalancer& operator= (Rebalancer&&) = default;

    void initialize() {
      for (PartitionID part = 0; part < k; part++) {
        _queues.emplace_back(_hg.initialNumNodes()); //this is too large TODO(fritsch)
      }
    }

    void rebalance(HypernodeWeight heaviest_node_weight, IRefiner& refiner, Metrics& current_metrics) {
      _current_upper_bound = std::max((1.0 + _context.partition.epsilon) * static_cast<double>(_hg.totalWeight()) / static_cast<double>(k),
                                                    heaviest_node_weight + static_cast<double>(_hg.totalWeight()) / static_cast<double>(k));
      reset();
      // ------------------------------
      // calculate connected blocks
      // ------------------------------
      for (HyperedgeID edge : _hg.edges()) {
        //maybe optimize for low current num nodes, since its expected to be fully connected
        for (PartitionID block1 : _hg.connectivitySet(edge)) {
          for (PartitionID block2 : _hg.connectivitySet(edge)) {
            if (block1 == block2) continue;
            _adjacency_bitmap[block1 * (k - 1) + block2] = true; //TODO(fritsch) optimize vector size in half and remove redundancy
            _adjacency_bitmap[block2 * (k - 1) + block1] = true;
          }
        }
      }
      std::fill(_adjacency_bitmap.begin(), _adjacency_bitmap.end(), true); //TODO(fritsch) this is temporary and must be removed.
      // ------------------------------
      // fill binary heaps
      // ------------------------------
      for (HypernodeID node : _hg.nodes()) {
        ASSERT(_hg.nodeIsEnabled(node));
        const PartitionID from_part = _hg.partID(node);
        if (_hg.partWeight(from_part) <= _current_upper_bound) continue;
        const std::pair<PartitionID, Gain> relativeGainMove = highestGainMoveToNotOverloadedBlock(node);
        const PartitionID to_part = relativeGainMove.first;
        const Gain relativeGain = relativeGainMove.second;
        if (relativeGain == _invalid_gain || to_part == _invalid_part) {
          DBG << "no move found for node " << node << " from part " << from_part;
          ASSERT(relativeGain == _invalid_gain && to_part == _invalid_part);
          continue;
        }
        //try to insert move into
        if (_queue_weights[from_part] < excessWeight(from_part)) {
          insertIntoQ(node, to_part, relativeGain);
        } else if (relativeGain > _queues[from_part].minKey()) {
          insertIntoQ(node, to_part, relativeGain);
          if (_queue_weights[from_part] > excessWeight(from_part) + _heaviest_node_weight) {
            NodeMove removedMove = nodeMoveFromInt(_queues[from_part].min());
            _queues[from_part].popMin();
            _queue_weights[from_part] -= _hg.nodeWeight(removedMove.node);
            ASSERT(removedMove.to_part != from_part);
            ASSERT(!_queues[from_part].contains(nodeMoveToInt(removedMove)));
          }
        }
      }
      // ------------------------------
      // empty binary heaps
      // ------------------------------
      bool movedAnythingThisIteration = false;
      do {
        movedAnythingThisIteration = false;
        //block order policy
        //for now implement round robin in arbitrary order
        for (PartitionID block = 0; block < k; block++) {
          if (_hg.partWeight(block) <= _current_upper_bound || _queues[block].empty()) {
            continue;
          }
          NodeMove move = nodeMoveFromInt(_queues[block].max());
          ASSERT(block == _hg.partID(move.node), V(block) << V(_hg.partID(move.node)));
          if (_hg.partWeight(move.to_part) + _hg.nodeWeight(move.node) > _current_upper_bound
          || gainChangedFor(move.node, move.to_part))  {
            _queues[block].popMax();
            ASSERT(!_queues[block].contains(nodeMoveToInt(move)));
            _queue_weights[block] -= _hg.nodeWeight(move.node);
            
            if (_hg.isBorderNode(move.node)) {
              std::pair<PartitionID, Gain> newMove = highestGainMoveToNotOverloadedBlock(move.node);
              PartitionID to_part = newMove.first;
              Gain newRelativeGain = newMove.second;
              if (newRelativeGain == _invalid_gain || to_part == _invalid_part) {
                ASSERT(newRelativeGain == _invalid_gain && to_part == _invalid_part);
              } else {
                insertIntoQ(move.node, to_part, newRelativeGain);
              }
              //Try to insert all neighbours from former block
              for (const HyperedgeID& edge : _hg.incidentEdges(move.node)) {
                for (const HypernodeID& neighbour : _hg.pins(edge)) {
                  if (!_hg.active(neighbour) || neighbour == move.node || _hg.partID(neighbour) != block) {
                    tryToInsertIntoCorrectQ(move.node);
                  }
                }
              }
            }

          } else {
            _queues[block].popMax();
            ASSERT(!_queues[block].contains(nodeMoveToInt(move)));
            _queue_weights[block] -= _hg.nodeWeight(move.node);  
            DBG << "hypernode " << move.node << "[" << _hg.nodeWeight(move.node) << "] is moved from part " << block << "[" << _hg.partWeight(block) << "] to part " << move.to_part<< "[" << _hg.partWeight(move.to_part) << "]";
            ASSERT(_hg.nodeIsEnabled(move.node));
            refiner.moveNodeExternallyAndKeepInternalCacheCorrect(move.node, block, move.to_part);
            // mark other gains as invalid
            ASSERT(!_queues[block].contains(nodeMoveToInt(NodeMove{move.node, block})));
            movedAnythingThisIteration = true;
          }
        }
      } while(movedAnythingThisIteration);
      current_metrics.heaviest_block_weight = metrics::heaviest_block_weight(_hg);
      current_metrics.standard_deviation = metrics::standard_deviation(_hg);
      //TODO(fritsch) update and dont fully compute again
      current_metrics.km1 = metrics::km1(_hg);
    }


  private:
    inline bool gainChangedFor(HypernodeID node, PartitionID to_part) {
      NodeMove move{node, to_part};
      HypernodeID moveID = nodeMoveToInt(move);
      ASSERT(_queues[_hg.partID(node)].contains(moveID));
      Gain actualGain = gainInducedByHypergraph(node, to_part);
      return actualGain == _queues[_hg.partID(node)].getKey(moveID);
    }

    inline bool tryToInsertIntoCorrectQ(HypernodeID node) {
      const std::pair<PartitionID, Gain> relativeGainMove = highestGainMoveToNotOverloadedBlock(node);
      const PartitionID to_part = relativeGainMove.first;
      const Gain relativeGain = relativeGainMove.second;
      return insertIntoQ(node, to_part, relativeGain);
    }

    inline bool insertIntoQ(HypernodeID node, PartitionID to_part, Gain relativeGain) {
      if (_was_inserted_into_q_bitmap[node]) {
        return false;
      }
      PartitionID from_part = _hg.partID(node);
      NodeMove move {node, to_part};
      ASSERT(!_queues[from_part].contains(nodeMoveToInt(move)));
      ASSERT(!_queues[to_part].contains(nodeMoveToInt(move)));
      ASSERT(to_part != from_part);
      _queues[from_part].push(nodeMoveToInt(move), relativeGain);
      _queue_weights[from_part] += _hg.nodeWeight(node);
      return true;
    }

    inline void reset() {
      std::fill(_queue_weights.begin(), _queue_weights.end(), 0);
      std::fill(_adjacency_bitmap.begin(), _adjacency_bitmap.end(), false);
      std::fill(_was_inserted_into_q_bitmap.begin(), _was_inserted_into_q_bitmap.end(), false);
      for (Queue&  q : _queues) {
        q.clear();
      }
    }

    std::pair<HypernodeID, Gain> highestGainMoveToNotOverloadedBlock(HypernodeID node) const {
      const PartitionID from_part = _hg.partID(node);
      Gain max_gain = _invalid_gain;
      PartitionID max_to_part = _invalid_part;
      for (PartitionID to_part = 0; to_part < k; to_part++) {
        if (from_part == to_part
        || !_adjacency_bitmap[from_part * (k - 1) + to_part]
        || excessWeight(from_part) <= 0
        || _hg.partWeight(to_part) + _hg.nodeWeight(node) > _current_upper_bound) {
          continue;
        }
        Gain newGain = gainInducedByHypergraph(node, to_part);
        if (newGain > max_gain) {
          max_gain = newGain;
          max_to_part = to_part;
        }
        //TODO(fritsch) tiebreaker if newGain == max_gain?
      }
      ASSERT((max_to_part == _invalid_part) == (max_gain == _invalid_gain));
      ASSERT((max_to_part < k || max_to_part == _invalid_part) && max_to_part != from_part);
      ASSERT(max_to_part == _invalid_part || _hg.partWeight(max_to_part) + _hg.nodeWeight(node) <= _current_upper_bound);
      if (max_gain != _invalid_gain) {
        max_gain = (max_gain >= 0) ? max_gain * _hg.nodeWeight(node) : max_gain / _hg.nodeWeight(node);
      }
      return std::make_pair(max_to_part, max_gain);
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

    NodeMove nodeMoveFromInt(uint32_t move) const {
      ASSERT(_hg.initialNumNodes() - k > 0);
      PartitionID to_part = move / _hg.initialNumNodes();
      HypernodeID node = move % _hg.initialNumNodes();
      NodeMove nodeMove {node, static_cast<PartitionID>(to_part)};
      ASSERT(nodeMoveToInt(nodeMove) == move);
      ASSERT(_hg.nodeIsEnabled(node));
      ASSERT(to_part < k && to_part >= 0);
      return nodeMove;
    }

    uint32_t nodeMoveToInt(NodeMove nodeMove) const {
      ASSERT(_hg.initialNumNodes() - k > 0);
      ASSERT(nodeMove.to_part < k && nodeMove.to_part >= 0, V(nodeMove.to_part));
      uint32_t move = nodeMove.node + nodeMove.to_part * _hg.initialNumNodes();
      ASSERT(move % _hg.initialNumNodes() == nodeMove.node);
      ASSERT(move / _hg.initialNumNodes() - nodeMove.to_part == 0);
      ASSERT(_hg.nodeIsEnabled(nodeMove.node));
      return move;
    }

    Hypergraph& _hg;
    const Context& _context;
    PartitionID k;
    std::vector<Queue> _queues;
    std::vector<HypernodeWeight> _queue_weights;
    HypernodeWeight _heaviest_node_weight;
    HypernodeWeight _current_upper_bound;
    std::vector<bool> _adjacency_bitmap;
    std::vector<bool> _was_inserted_into_q_bitmap;
};
}