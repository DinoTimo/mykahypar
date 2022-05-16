#pragma once

#include "kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "kahypar/partition/initial_partitioning/initial_partitioner_base.h"

namespace kahypar {
/* ! This class should only ever used, when the graph is contracted down to k vertices during coarsening.
   ! This class can be used with "i-algo=direct_k"
*/
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
        LOG << "!!! Direct k initial partitioning was used, although there were not k nodes !!!";
        std::exit(1);
      }
    }
};
}