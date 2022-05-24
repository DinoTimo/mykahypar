#pragma once

#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"
#include "kahypar/datastructure/kway_priority_queue.h"

namespace kahypar {
class DirectKInitialPartitioner : public IInitialPartitioner,
                                  private InitialPartitionerBase<DirectKInitialPartitioner> {
  using Base = InitialPartitionerBase<DirectKInitialPartitioner>;
  friend Base;

  public:
    DirectKInitialPartitioner(Hypergraph& hypergraph, Context& context) :
      Base(hypergraph, context),
      _descending_nodes(bin_packing::nodesInDescendingWeightOrder(hypergraph)),
      _block_queue(context.partition.k) { }

    ~DirectKInitialPartitioner() override = default;

    DirectKInitialPartitioner(const DirectKInitialPartitioner&) = delete;
    DirectKInitialPartitioner& operator= (const DirectKInitialPartitioner&) = delete;

    DirectKInitialPartitioner(DirectKInitialPartitioner&&) = delete;
    DirectKInitialPartitioner& operator= (DirectKInitialPartitioner&&) = delete;

  private:
    void partitionImpl() override final {
      for (PartitionID k = 0; k < _context.partition.k; k++) {
        _block_queue.push(k, 0);
      }
      for (const auto hn : _hg.nodes()) {
        if (_hg.nodeIsEnabled(hn)) {
          PartitionID lightestBlock = _block_queue.top();
          _hg.setNodePart(hn, lightestBlock);
          _block_queue.increaseKey(lightestBlock, _hg.nodeWeight(hn));
        }
      }
    }
    std::vector<HypernodeID> _descending_nodes;
    ds::BinaryMinHeap<PartitionID, HypernodeWeight> _block_queue;
};
}