
#include <vector>
#include <random>

#include "gmock/gmock.h"
#include "gtest/gtest_prod.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/refinement/flow/graph_flow_solver.h"
#include "kahypar/datastructure/hypergraph.h"

using ::testing::Test;
using ::testing::Eq;


namespace kahypar {

class AFlowSolver : public Test {
  private:
  public:
    AFlowSolver() : 
    flowSolver() {     }
    
    FlowSolver<HypernodeWeight, PartitionID> flowSolver;
};

TEST_F(AFlowSolver, SimpleFlowSolverTest) { 
  PartitionID source = 0;
  PartitionID sink = 2;
  std::vector<HypernodeWeight> capacity_matrix({-1, 1, 1,
                                                1, 1, 1,
                                                1, 1, -1});

  std::vector<HypernodeWeight> flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({2, 1, 1,
                                                       -1, 1, 1,
                                                       -1, -1, 0}));
  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({2, 1, 1,
                                                       -1, 1, 1,
                                                       -1, -1, 0}));


  capacity_matrix = std::vector<HypernodeWeight>({-1, 1, 1,
                                                  1, 0, 1,
                                                  1, 1, -1});

  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({1, 0, 1,
                                                       0, 0, 0,
                                                       -1, 0, 0}));
  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({2, 1, 1,
                                                       -1, 1, 1,
                                                       -1, -1, 0}));
  const uint16_t num_runs = 50;
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<HypernodeWeight> rand;
  for (uint16_t run = 0; run < num_runs; run++) {
    capacity_matrix = std::vector<HypernodeWeight>({rand(rng), 1, 1,
                                                    1, rand(rng), 1,
                                                    1, 1, rand(rng)});
    flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, false);
    ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({2, 1, 1,
                                                        -1, 1, 1,
                                                        -1, -1, 0}));
  }
 
}

TEST_F(AFlowSolver, MoreComplexFlowSolverTest) { 
  std::vector<HypernodeWeight> capacity_matrix({0, 5, 5, 0,
                                                0, 3, 0, 6,
                                                0, 7, 5, 5,
                                                0, 0, 0, 0
  });
  PartitionID source = 0;
  PartitionID sink = 3;

  std::vector<HypernodeWeight> flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({10, 5, 5, 0,
                                                       -5,  5, 0, 5,
                                                       -5,  0, 5, 5,
                                                       0,  -5, -5, 0 
  }));

  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({8, 3, 5, 0,
                                                       -3, 3, 0, 3,
                                                       -5, 0, 5, 5,
                                                       0, -3, -5, 0 
  }));  

  capacity_matrix = std::vector<HypernodeWeight>({0, 4, 5, 0,
                                                  0, 6, 0, 6,
                                                  0, 5, 5, 3,
                                                  0, 0, 0, 0});
  
  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, true);
  
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({9, 4, 5, 0,
                                                       -4, 6, -2, 6,
                                                       -5, 2, 5, 3,
                                                       0, -6, -3, 0}));
  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({9, 4, 5, 0,
                                                       -4, 6, -2, 6,
                                                       -5, 2, 5, 3,
                                                       0, -6, -3, 0}));

  sink = 4;

  capacity_matrix = std::vector<HypernodeWeight>({0, 6, 4, 0, 0,
                                                  0, 4, 4, 3, 0,
                                                  0, 0, 5, 2, 3,
                                                  0, 0, 0, 5, 4,
                                                  0, 0, 0, 0, 0});
  
  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({7, 3, 4, 0, 0,
                                                       -3, 3, 0, 3, 0,
                                                       -4, 0, 4, 1, 3,
                                                       0, -3, -1, 4, 4, 
                                                       0, 0, -3, -4, 0
                                                        }));

  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({7, 3, 4, 0, 0,
                                                       -3, 3, 0, 3, 0,
                                                       -4, 0, 4, 1, 3,
                                                       0, -3, -1, 4, 4, 
                                                       0, 0, -3, -4, 0
                                                        }));


  capacity_matrix = std::vector<HypernodeWeight>({0, 6, 4, 0, 0,
                                                  0, 4, 4, 1, 0,
                                                  0, 0, 5, 3, 3,
                                                  0, 0, 0, 5, 4,
                                                  0, 0, 0, 0, 0});
  
  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({7, 3, 4, 0, 0,
                                                       -3, 3, 2, 1, 0,
                                                       -4, -2, 6, 3, 3,
                                                       0, -1, -3, 4, 4, 
                                                       0, 0, -3, -4, 0}));

  flow_matrix = flowSolver.solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flow_matrix, std::vector<HypernodeWeight>({6, 2, 4, 0, 0,
                                                       -2, 2, 1, 1, 0,
                                                       -4, -1, 5, 2, 3,
                                                       0, -1, -2, 3, 3, 
                                                       0, 0, -3, -3, 0}));

}
}