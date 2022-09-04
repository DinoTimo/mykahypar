
#include "gmock/gmock.h"
#include "gtest/gtest_prod.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/refinement/flow/graph_flow_solver.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/partition/refinement/rebalancer.h"
#include "kahypar/partition/refinement/rebalancing_refiner.h"
#include "kahypar/partition/refinement/policies/fm_acceptance_policy.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"
#include "kahypar/partition/refinement/policies/fm_improvement_policy.h"

using ::testing::Test;
using ::testing::Eq;


namespace kahypar {
const static size_t k = 4;

class ARebalancer : public Test {
  using ThisRefiner = RebalancingKwayKMinusOneRefiner<AdvancedRandomWalkModelStopsSearch,
                                  MultilevelFlowExecution,
                                  BalanceApproachingAcceptancePolicy>;
  public:
    ARebalancer() :
    _hypergraph(new Hypergraph(7, 4, HyperedgeIndexVector { 0, 2, 6, 9,  /*sentinel*/ 12 },
                            HyperedgeVector { 0, 2, 0, 1, 3, 4, 3, 4, 6, 2, 5, 6 }, k)),
    _context(),
    _refiner(nullptr) {
      _context.partition.k = k;
      _context.local_search.fm.rebalancing_order_policy = RebalancerType::round_robin;
      _context.local_search.algorithm = RefinementAlgorithm::rebalancing_kway_fm_km1;
      _hypergraph->setNodePart(0, 0);
      _hypergraph->setNodePart(1, 1);
      _hypergraph->setNodePart(2, 2);
      _hypergraph->setNodePart(3, 3);
      _hypergraph->setNodePart(4, 3);
      _hypergraph->setNodePart(5, 3);
      _hypergraph->setNodePart(6, 3);
      _refiner = std::make_unique<ThisRefiner>(*_hypergraph, _context);
      _refiner->initialize(-1);
    }
    
    std::unique_ptr<Hypergraph> _hypergraph;
    Context _context;
    std::unique_ptr<ThisRefiner> _refiner;
};

TEST_F(ARebalancer, UnbalancedHypergraphRebalancingTest) { 

  Metrics best_metrics;
  best_metrics.heaviest_block_weight = metrics::heaviest_block_weight(*_hypergraph);
  best_metrics.km1 = metrics::km1(*_hypergraph);
  std::vector<HypernodeID> refinementNodes;

  HypernodeWeight previous_imbalance = best_metrics.heaviest_block_weight;
  HypernodeWeight target_weight = 2;
  _refiner->performRebalancing(best_metrics, refinementNodes, target_weight);
  LOG << "done rebalancing";
  for (size_t block = 0; block < k; block++) {
    ASSERT(_hypergraph->partWeight(block) <= target_weight, V(target_weight) << ", " << V(_hypergraph->partWeight(block)));
  }
  ASSERT(previous_imbalance > best_metrics.heaviest_block_weight, V(previous_imbalance) << ", " << V(best_metrics.heaviest_block_weight));
}

}