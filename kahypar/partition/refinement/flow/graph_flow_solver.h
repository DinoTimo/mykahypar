

#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include <memory>
#include <limits>
#include "gtest/gtest_prod.h"

template <typename Capacity = int, typename Node = int>
class FlowSolver {
  using Flow = Capacity;
  struct Edge {
    Node start;
    Node end;
    
    bool operator==(const Edge& e) const {
      return start == e.start && end == e.end;
    }

    bool operator!=(const Edge& e) const {
      return !operator==(e);
    }
  };

  private:
    FRIEND_TEST(AFlowSolver, SimpleFlowSolverTest);  //TODO this does not seem to work right now
    FRIEND_TEST(AFlowSolver, MoreComplexFlowSolverTest);
    static constexpr bool debug = false;
    constexpr static Node _invalid_node = -1;
    const Edge _invalid_edge = Edge { _invalid_node, _invalid_node};
    std::vector<Capacity> _input_capacity_adjacency_matrix;
    std::vector<Flow> _output_flow_adjacency_matrix;
    std::vector<Edge> _pred;
    std::vector<Node> _queue;
    Node _source;
    Node _sink;
    int _num_nodes;
    Flow _total_flow;

    void setFlow(const Edge e, const Flow flow) {
      setFlow(e.start, e.end, flow);
    }

    void setFlow(const Node v, const Flow flow) {
      setFlow(v, v, flow);
    }

    void setFlow(const Node start, const Node end, const Flow flow) {
      _output_flow_adjacency_matrix[start * _num_nodes + end] = flow;
    }

    Edge reverseEdge(const Edge e) const {
      return Edge {e.end, e.start};
    }

    std::vector<Edge> getOutgoingEdgesOf(const Node start) const {
      std::vector<Edge> edges;
      for (Node end = 0; end < _num_nodes; end++) {
        if (end == start) {
          continue;
        }
        edges.push_back(Edge {start, end});
      }
      return edges;
    }

    Capacity capacity(const Edge e) const {
      return capacity(e.start, e.end);
    }

    Capacity capacity(const Node v) const {
      if (v == _source || v == _sink) {
        return std::numeric_limits<Capacity>::max();
      }
      return capacity(v, v);
    }

    Capacity capacity(const Node u, const Node v) const {
      return _input_capacity_adjacency_matrix[u * _num_nodes + v];
    }

        Flow flow(const Edge e) const {
      return flow(e.start, e.end);
    }
    
 

    Flow flow(const Node start, const Node end) const {
      return _output_flow_adjacency_matrix[start * _num_nodes + end];
    }

    void init(const std::vector<Capacity> input, const Node source, const Node sink) {
      _input_capacity_adjacency_matrix = move(input);
      _output_flow_adjacency_matrix = std::vector<Capacity>(input.size(), 0);
      _num_nodes = (int) sqrt(input.size());
      _pred = std::vector<Edge>(_num_nodes * (_num_nodes - 1), _invalid_edge);
      if(_num_nodes * _num_nodes  != input.size()) {
        LOG << "invalid input size of " << std::to_string(input.size()) <<". This is now undefined behaviour";
      }
      if(source < 0) {
        LOG << "invalid source: " << std::to_string(source) <<". This is now undefined behaviour";
      }
      if(sink >= _num_nodes) {
        LOG << "invalid sink: " << std::to_string(sink) <<". This is now undefined behaviour";
      }
      _queue = std::vector<Node>();
      //maybe reserve some space for the _queue? log n?
      _total_flow = 0;
      _source = source;
      _sink = sink;
    } 

  public:
    FlowSolver() : 
    _input_capacity_adjacency_matrix(),
    _output_flow_adjacency_matrix(),
    _pred(),
    _queue(),
    _source(0),
    _sink(0),
    _num_nodes(0),
    _total_flow(0) { }

    FlowSolver(const FlowSolver&) = delete;
    FlowSolver(FlowSolver&&) = delete;
    FlowSolver& operator= (const FlowSolver&) = delete;
    FlowSolver& operator= (FlowSolver&&) = delete;

    //this should be private and shared to test classes via Friend_test TODO
    Flow flow(const Node v) const {
      Flow total = 0;
      for (const Edge e : getOutgoingEdgesOf(v)) {
        total += flow(e);
      }
      return total;
    }

    std::vector<Capacity> solveFlow(const std::vector<Capacity> input, const Node source, const Node sink, const bool respect_Node_Capacity) {
      init(input, source, sink);
      do {
        _queue.insert(_queue.begin(), _source);
        fill(_pred.begin(), _pred.end(), _invalid_edge);
        while(!_queue.empty()) {
          Node cur = _queue.back();
          _queue.pop_back();
          std::vector<Edge> curEdges = getOutgoingEdgesOf(cur);
          for (Edge e : curEdges) {
            if (_pred[e.end] == _invalid_edge && e.end != _source &&
              capacity(e) > flow(e) &&
              (!respect_Node_Capacity || capacity(e.start) > flow(e.start))) { 
              _pred[e.end] = e;
              _queue.insert(_queue.begin(), e.end);
            }
          }
        } 
        
        if (_pred[_sink] != _invalid_edge) {
          size_t pathLen = 0;
          DBG << "path found";
          Flow cur_flow = std::numeric_limits<Flow>::max();
          for (Edge e = _pred[_sink]; e != _invalid_edge; e = _pred[e.start]) {
            cur_flow = std::min(cur_flow, capacity(e) - flow(e));
            if (respect_Node_Capacity) {
              cur_flow = std::min(cur_flow, capacity(e.start) - flow(e.start));
            }
            pathLen++;
          }
          DBG << "found flow of size " << std::to_string(cur_flow) << " and of length " << std::to_string(pathLen);
          for (Edge e = _pred[_sink]; e != _invalid_edge; e = _pred[e.start]) {
            setFlow(e, flow(e) + cur_flow);
            setFlow(reverseEdge(e), flow(e) - cur_flow);
          }
          DBG << "flow updated\n";
          _total_flow += cur_flow;
        }
      } while (_pred[_sink] != _invalid_edge);
      DBG << "total flow is = " << std::to_string(_total_flow);
      DBG << "input flow at target is = " << std::to_string(flow(sink)) << " and at source is = " << std::to_string(flow(source));
      return _output_flow_adjacency_matrix;
    }
};