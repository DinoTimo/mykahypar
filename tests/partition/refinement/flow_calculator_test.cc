
#include <vector>
#include "gmock/gmock.h"
#include "gtest/gtest_prod.h"

#include "kahypar/definitions.h"
#include "kahypar/partition/refinement/flow/graph_flow_solver.h"
#include "kahypar/datastructure/hypergraph.h"
#include "kahypar/partition/refinement/kway_fm_cut_refiner.h"
#include "kahypar/partition/refinement/policies/fm_stop_policy.h"

using ::testing::Test;
using ::testing::Eq;


namespace kahypar {

class AFlowSolver : public Test {
  private:
  public:
    AFlowSolver() : 
    capacity_matrix(),
    flowSolver(nullptr) { 
      flowSolver = std::make_shared<FlowSolver<HypernodeWeight, PartitionID>>();
    }

    std::vector<HypernodeWeight> capacity_matrix;
    std::shared_ptr<FlowSolver<HypernodeWeight, PartitionID>> flowSolver;
};

TEST_F(AFlowSolver, SimpleFlowSolverTest) { 
  capacity_matrix.clear();
  PartitionID source = 0;
  PartitionID sink = 2;
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  capacity_matrix.push_back(1);
  flowSolver->solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flowSolver->flow(source), 2);
  flowSolver->solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flowSolver->flow(source), 2);
}

TEST_F(AFlowSolver, MoreComplexFlowSolverTest) { 
  capacity_matrix.clear();
  PartitionID source = 0;
  PartitionID sink = 3;

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(5);
  capacity_matrix.push_back(5);
  capacity_matrix.push_back(0);

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(3);
  capacity_matrix.push_back(0);
  capacity_matrix.push_back(6);

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(7);
  capacity_matrix.push_back(5);
  capacity_matrix.push_back(5);

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(0);
  capacity_matrix.push_back(0);
  capacity_matrix.push_back(0);

  flowSolver->solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flowSolver->flow(source), 10);
  flowSolver->solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flowSolver->flow(source), 8);

  capacity_matrix.clear();

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(4);
  capacity_matrix.push_back(5);
  capacity_matrix.push_back(0);

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(6);
  capacity_matrix.push_back(0);
  capacity_matrix.push_back(6);

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(5);
  capacity_matrix.push_back(5);
  capacity_matrix.push_back(3);

  capacity_matrix.push_back(0);
  capacity_matrix.push_back(0);
  capacity_matrix.push_back(0);
  capacity_matrix.push_back(0);

  flowSolver->solveFlow(capacity_matrix, source, sink, true);
  ASSERT_EQ(flowSolver->flow(source), 9);

  flowSolver->solveFlow(capacity_matrix, source, sink, false);
  ASSERT_EQ(flowSolver->flow(source), 9);

}

}