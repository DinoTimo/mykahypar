
#pragma once

#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/partition/refinement/fm_refiner_base.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/refinement/flow/graph_flow_solver.h"
#include "kahypar/partition/refinement/flow/flow_balancing_refiner.h"
#include "kahypar/partition/refinement/flow/quotient_graph_block_scheduler.h"
#include "kahypar/partition/refinement/matrices/matrix.h"
#include "kahypar/partition/refinement/matrices/matrix_solver.h"
#include "kahypar/partition/refinement/flow/policies/flow_execution_policy.h"

namespace kahypar {

//This class does not use the type 'Derived' and only passes it to the FMRefinerBase class
template <typename RollbackElement = Mandatory, typename FlowExecutionPolicy = Mandatory, typename Derived = Mandatory>
class FlowBalancingRefiner : protected FMRefinerBase<RollbackElement, Derived> {

  private:
    static constexpr bool debug = false;
    using Base = FMRefinerBase<RollbackElement, Derived>;

  public:
    FlowBalancingRefiner(const FlowBalancingRefiner&) = delete;
    FlowBalancingRefiner& operator= (const FlowBalancingRefiner&) = delete;

    FlowBalancingRefiner(FlowBalancingRefiner&&) = delete;
    FlowBalancingRefiner& operator= (FlowBalancingRefiner&&) = delete;

    ~FlowBalancingRefiner() = default;

  protected:
    FlowBalancingRefiner(Hypergraph& hypergraph, const Context& context) : FMRefinerBase<RollbackElement, Derived>(hypergraph, context),
      _flow_solver(),
      _step0_smallest_block_weight(0),
      _step0_heaviest_block_weight(0),
      _initial_imbalance_set(false),
      _total_num_steps(0),
      _current_step(0),
      _num_flow_nodes(context.partition.k + 2),
      _flow_matrix(_num_flow_nodes * _num_flow_nodes, 0),
      _capacity_matrix(_num_flow_nodes * _num_flow_nodes, 0),
      _quotient_edge_capacities(context.partition.k * context.partition.k, 0),
      _vertex_block_pair_bitvector(hypergraph.initialNumNodes() * context.partition.k, false),
      _adjacency_bitmap(context.partition.k * (context.partition.k - 1), 0),
      _degree_vector(context.partition.k, 0),
      _laplace_matrix(context.partition.k * context.partition.k, 0),
      _block_weight_diff_vector(context.partition.k, 0),
      _matrix_solver(),
      _flow_vector(),
      _flow_execution_policy() { }

    HypernodeWeight idealBlockWeight() {
      return _hg.totalWeight() / _context.partition.k;
    }

    void initQuotientEdgeCapacities() {
      std::fill(_vertex_block_pair_bitvector.begin(), _vertex_block_pair_bitvector.end(), false);
      std::fill(_quotient_edge_capacities.begin(), _quotient_edge_capacities.end(), 0);
      for (HyperedgeID edge : _hg.edges()) {
        if (_hg.connectivitySet(edge).size() <= 1) {
          continue;
        }
        for (HypernodeID pin : _hg.pins(edge)) {
          for (PartitionID block : _hg.connectivitySet(edge)) {
            if (!_vertex_block_pair_bitvector[pin * _context.partition.k + block]) {
              _quotient_edge_capacities[_hg.partID(pin) * _context.partition.k + block] += _hg.nodeWeight(pin);
              _vertex_block_pair_bitvector[pin * _context.partition.k + block] = true;
            }
          }
        }
      }
    }
    
    void calculateCapacityMatrix() {
      if (_context.local_search.fm.flow_model == BalancingFlowModel::finite_edges) {
        initQuotientEdgeCapacities();
      }
      PartitionID source = _context.partition.k;
      PartitionID sink = source + 1;
      std::fill(_capacity_matrix.begin(), _capacity_matrix.end(), 0);
      QuotientGraphBlockScheduler scheduler(_hg, _context);
      scheduler.buildQuotientGraph();
      for (const std::pair<PartitionID, PartitionID> edge : scheduler.quotientGraphEdges()) {
        PartitionID heavierBlockID = (_hg.partWeight(edge.first) > _hg.partWeight(edge.second)) ? edge.first : edge.second;
        PartitionID lighterBlockID = (_hg.partWeight(edge.first) > _hg.partWeight(edge.second)) ? edge.second : edge.first;
        if (_context.local_search.fm.flow_model == BalancingFlowModel::finite_edges) {
          HypernodeWeight heavierBlockWeight = _hg.partWeight(heavierBlockID);
          HypernodeWeight lighterBlockWeight = _hg.partWeight(lighterBlockID);
          HypernodeWeight idealWeight = idealBlockWeight();
          HypernodeWeight overload = heavierBlockWeight - idealWeight;
          HypernodeWeight underload = idealWeight - lighterBlockWeight;
          if (overload < 0 || underload < 0) {
            _capacity_matrix[heavierBlockID * _num_flow_nodes + lighterBlockID] = std::min(heavierBlockWeight - lighterBlockWeight, calculateQuotientEdgeCapacity(heavierBlockID, lighterBlockID));
            continue;
          }
          HypernodeWeight maxEdgeCapacity = std::min(overload, underload);
          _capacity_matrix[heavierBlockID * _num_flow_nodes + lighterBlockID] = std::min(maxEdgeCapacity, calculateQuotientEdgeCapacity(heavierBlockID, lighterBlockID));
        } else if(_context.local_search.fm.flow_model == BalancingFlowModel::infinity_edges) {
          _capacity_matrix[heavierBlockID * _num_flow_nodes + lighterBlockID] = std::numeric_limits<HypernodeWeight>::max();
        }
      }

      for (int i = 0; i < source; i++) {
        for (int j = i + 1; j < source; j++) {
          HypernodeWeight & capacity = _capacity_matrix[i * _num_flow_nodes + j]; 
          capacity = capacity > 0 ? capacity : _hg.totalWeight() / (source);
          capacity = _capacity_matrix[j * _num_flow_nodes + i]; 
          capacity = capacity > 0 ? capacity : _hg.totalWeight() / (source);
        }
      }
      
      // Add edges between source and overloaded blocks and sink and underloaded blocks
      for (PartitionID blockNode = 0; blockNode < _context.partition.k; blockNode++) {
        if (isOverloadedBlock(blockNode)) {
          _capacity_matrix[source * _num_flow_nodes + blockNode] = _hg.partWeight(blockNode) - idealBlockWeight(); 
        }
        if (isUnderloadedBlock(blockNode)) {
          _capacity_matrix[blockNode * _num_flow_nodes + sink] = idealBlockWeight() - _hg.partWeight(blockNode);
        }
      }
    }

    HypernodeWeight calculateQuotientNodeCapacity(PartitionID quotientNode) {
      return _hg.partWeight(quotientNode);
    }

    HypernodeWeight calculateQuotientEdgeCapacity(PartitionID first, PartitionID second) {
      return _quotient_edge_capacities[first * _context.partition.k + second];
    }
    
    bool isOverloadedBlock(PartitionID block) {
      return _hg.partWeight(block) > idealBlockWeight();
    }

    bool isUnderloadedBlock(PartitionID block) {
      return _hg.partWeight(block) < idealBlockWeight();
    }

    bool moveFeasibilityByFlow(PartitionID from, PartitionID to, HypernodeID node) {
      /*
      return (_hg.nodeWeight(node) <= _flow_matrix[from * _num_flow_nodes + to] * 2) 
      || ( (2 * (_hg.partWeight(from) - _hg.partWeight(to)) >= _hg.nodeWeight(node))
          && isOverloadedBlock(from)
          && isUnderloadedBlock(to));*/
      return _hg.nodeWeight(node) <= 2 * (_flow_vector[from] - _flow_vector[to]);
    }

    std::vector<HypernodeWeight> solveFlow(HypernodeID source, HypernodeID sink, bool useNodeCapacities) {
      return _flow_solver.solveFlow(_capacity_matrix, source, sink, useNodeCapacities);
    }

    void rollbackFlow(int last_index, const int min_cut_index) {
      DBG << "min_cut_index=" << min_cut_index;
      DBG << "last_index=" << last_index;
      while (last_index != min_cut_index) {
        const HypernodeID hn = _performed_moves[last_index].hn;
        const PartitionID from_part = _performed_moves[last_index].to_part;
        const PartitionID to_part = _performed_moves[last_index].from_part;
        updateFlow(hn, from_part, to_part);
        --last_index;
      }
    }

    void updateFlow(HypernodeID hn, PartitionID from_part, PartitionID to_part) {
      _flow_matrix[from_part * _num_flow_nodes + to_part] -= _hg.nodeWeight(hn);
      _flow_matrix[to_part * _num_flow_nodes + from_part] += _hg.nodeWeight(hn);
      _flow_vector[from_part] -= _hg.nodeWeight(hn);
      _flow_vector[to_part] += _hg.nodeWeight(hn);
    }

    // ------------------------------------------------------------------------------------------------------------------------
    // ABOVE ACTUAL FLOW CALCULATION ON QUTIOENT GRAPH
    // BELOW LAPLACE MATRIX BALANICNG FLOW CALCULATION
    // ------------------------------------------------------------------------------------------------------------------------
    
    void calculateLaplaceMatrix() {
      calculateAdjacency();
      std::vector<double> b(_context.partition.k, 0);
      for (int i = 0; i < _context.partition.k; i++) {
        for (int j = 0; j < _context.partition.k; j++) {
          _laplace_matrix[i * _context.partition.k + j] = (i == j) ? _degree_vector[i] : _adjacency_bitmap[i * (_context.partition.k - 1) + j];
        }  
      }
      for (int i = 0; i <_context.partition.k; i++) {
        _block_weight_diff_vector[i] = _hg.partWeight(i) - idealBlockWeight();
      }
    }

    void calculateAdjacency() {
      std::fill(_adjacency_bitmap.begin(), _adjacency_bitmap.end(), false);
      std::fill(_degree_vector.begin(), _degree_vector.end(), 0);
      for (HyperedgeID edge : _hg.edges()) {
        if (_hg.connectivitySet(edge).size() <= 1) {
          continue;
        }
        for (PartitionID block1 : _hg.connectivitySet(edge)) {
          for (PartitionID block2 : _hg.connectivitySet(edge)) {
            if (block1 == block2) continue;
            int index = block1 * (_context.partition.k - 1) + block2;
            if (_adjacency_bitmap[index]) continue;
            _degree_vector[block1] += 1;
            _degree_vector[block2]++;
            _adjacency_bitmap[index] = true;
          }
        }
        
      }
    }

    void solveBalancingEquations() {
      //solve Mx = b, where b is the diff from ideal weight and M the laplacian matrix
      _flow_vector = _matrix_solver.solve(
        matrices::full_square_matrix<double>(std::vector<double>(_laplace_matrix.begin(), _laplace_matrix.end())),
        std::vector<double>(_block_weight_diff_vector.begin(), _block_weight_diff_vector.end()),
        _context.partition.k);
      std::cout << "flow vector" << std::endl;
      std::cout << joinVector(_flow_vector, "[" , ",", "]") << std::endl;
    }

    FlowSolver<HypernodeWeight, PartitionID> _flow_solver;
    HypernodeWeight _step0_smallest_block_weight;
    HypernodeWeight _step0_heaviest_block_weight;
    bool _initial_imbalance_set;
    uint32_t _total_num_steps;
    uint32_t _current_step;
    PartitionID _num_flow_nodes;
    std::vector<HypernodeWeight> _flow_matrix;
    std::vector<HypernodeWeight> _capacity_matrix;
    std::vector<HypernodeWeight> _quotient_edge_capacities;
    std::vector<bool> _vertex_block_pair_bitvector;
    std::vector<bool> _adjacency_bitmap;
    std::vector<HypernodeWeight> _degree_vector;
    std::vector<HypernodeWeight> _laplace_matrix;
    std::vector<HypernodeWeight> _block_weight_diff_vector;
    matrices::LU_Decomp_matrix_solver _matrix_solver;
    std::vector<double> _flow_vector;
    FlowExecutionPolicy _flow_execution_policy;
    
    using Base::_context;
    using Base::_hg;
    using Base::_performed_moves;
};
}