#pragma once

#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"

namespace kahypar {
class DirectKInitialPartitioner : public IInitialPartitioner,
                                  private InitialPartitionerBase<DirectKInitialPartitioner> {
  using Base = InitialPartitionerBase<DirectKInitialPartitioner>;
  friend Base;

  public:
    DirectKInitialPartitioner(Hypergraph& hypergraph, Context& context) :
      Base(hypergraph, context) { }

    ~DirectKInitialPartitioner() override = default;

    DirectKInitialPartitioner(const DirectKInitialPartitioner&) = delete;
    DirectKInitialPartitioner& operator= (const DirectKInitialPartitioner&) = delete;

    DirectKInitialPartitioner(DirectKInitialPartitioner&&) = delete;
    DirectKInitialPartitioner& operator= (DirectKInitialPartitioner&&) = delete;

  private:
    void partitionImpl() override final {
      uint16_t current_k = 0;
      for (const auto hn : _hg.nodes()) {
        if (_hg.nodeIsEnabled(hn)) {
          _hg.setNodePart(hn, current_k++);
        }
      }
      if (_context.partition.k != current_k) {
        LOG << "THIS SHOULD NOT HAVE BEEN USED";
        std::exit(1);
      }
    }
};
}