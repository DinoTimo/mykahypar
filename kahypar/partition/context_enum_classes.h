/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
******************************************************************************/

#pragma once

#include <iostream>
#include <string>

namespace kahypar {
enum class ContextType : bool {
  main,
  initial_partitioning
};

enum class Mode : uint8_t {
  recursive_bisection,
  direct_kway,
  UNDEFINED
};

enum class InitialPartitioningTechnique : uint8_t {
  multilevel,
  flat,
  UNDEFINED
};

enum class RatingFunction : uint8_t {
  heavy_edge,
  edge_frequency,
  UNDEFINED
};

enum class CommunityPolicy : uint8_t {
  use_communities,
  ignore_communities,
  UNDEFINED
};

enum class IgnoreMaxNodeWeight : uint8_t {
  respect_max_node_weight,
  ignore_max_node_weight,
  UNDEFINED
};

enum class ContractCommunities : uint8_t {
  contract_communities,
  dont_contract_communities,
  UNDEFINED
};

enum class HeavyNodePenaltyPolicy : uint8_t {
  no_penalty,
  multiplicative_penalty,
  edge_frequency_penalty,
  UNDEFINED
};

enum class AcceptancePolicy : uint8_t {
  best,
  best_prefer_unmatched,
  UNDEFINED
};

enum class RatingPartitionPolicy : uint8_t {
  normal,
  evolutionary
};

enum class FixVertexContractionAcceptancePolicy : uint8_t {
  free_vertex_only,
  fixed_vertex_allowed,
  equivalent_vertices,
  UNDEFINED
};

enum class CoarseningAlgorithm : uint8_t {
  heavy_full,
  heavy_lazy,
  ml_style,
  do_nothing,
  UNDEFINED
};

enum class RefinementAlgorithm : uint8_t {
  twoway_fm,
  kway_fm,
  kway_fm_km1,
  balance_approaching_kway_fm_km1,
  imbalance_holding_kway_fm_km1,
  twoway_fm_hyperflow_cutter,
  twoway_hyperflow_cutter,
  kway_hyperflow_cutter,
  kway_fm_hyperflow_cutter,
  kway_fm_hyperflow_cutter_km1,
  do_nothing,
  UNDEFINED
};

enum class InitialPartitionerAlgorithm : uint8_t {
  greedy_sequential,
  greedy_global,
  greedy_round,
  greedy_sequential_maxpin,
  greedy_global_maxpin,
  greedy_round_maxpin,
  greedy_sequential_maxnet,
  greedy_global_maxnet,
  greedy_round_maxnet,
  bfs,
  random,
  lp,
  bin_packing,
  pool,
  direct_k,
  UNDEFINED
};

enum class BinPackingAlgorithm : uint8_t {
  worst_fit,
  first_fit,
  UNDEFINED
};

enum class LouvainEdgeWeight : uint8_t {
  hybrid,
  uniform,
  non_uniform,
  degree,
  UNDEFINED
};

enum class RefinementStoppingRule : uint8_t {
  simple,
  adaptive_opt,
  UNDEFINED
};

enum class BalancingFlowModel : uint8_t {
  laplace_matrix,
  quotient_flow,
  UNDEFINED
};

enum class Objective : uint8_t {
  cut,
  km1,
  UNDEFINED
};
enum class EvoReplaceStrategy : uint8_t {
  worst,
  diverse,
  strong_diverse
};


enum class EvoCombineStrategy : uint8_t {
  basic,
  edge_frequency,
  UNDEFINED
};
enum class EvoMutateStrategy : uint8_t {
  new_initial_partitioning_vcycle,
  vcycle,
  UNDEFINED
};

enum class EvoDecision : uint8_t {
  normal,
  mutation,
  combine
};

enum class FileLogLevel : uint8_t {
  no_file_logging,
  write_imbalance_km1,
  write_imbalance_km1_target
};

enum class FlowExecutionMode : uint8_t {
  constant,
  multilevel,
  exponential,
  UNDEFINED
};

enum class AcceptanceRule : uint8_t {
  balance_approaching,
  imbalance_holding,
  staircase,
  UNDEFINED
};

enum class FlowHypergraphSizeConstraint : uint8_t {
  part_weight_fraction,
  max_part_weight_fraction,
  scaled_max_part_weight_fraction_minus_opposite_side
};

static std::ostream& operator<< (std::ostream& os, const EvoReplaceStrategy& replace) {
  switch (replace) {
    case EvoReplaceStrategy::worst: return os << "worst";
    case EvoReplaceStrategy::diverse: return os << "diverse";
    case EvoReplaceStrategy::strong_diverse: return os << "strong_diverse";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(replace);
}

static std::ostream& operator<< (std::ostream& os, const EvoCombineStrategy& combine) {
  switch (combine) {
    case EvoCombineStrategy::basic: return os << "basic";
    case EvoCombineStrategy::edge_frequency: return os << "edge_frequency";
    case EvoCombineStrategy::UNDEFINED: return os << "-";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(combine);
}

static std::ostream& operator<< (std::ostream& os, const EvoMutateStrategy& mutation) {
  switch (mutation) {
    case EvoMutateStrategy::new_initial_partitioning_vcycle:
      return os << "new_initial_partitioning_vcycle";
    case EvoMutateStrategy::vcycle: return os << "vcycle";
    case EvoMutateStrategy::UNDEFINED:  return os << "-";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(mutation);
}


static std::ostream& operator<< (std::ostream& os, const EvoDecision& decision) {
  switch (decision) {
    case EvoDecision::normal:  return os << "normal";
    case EvoDecision::mutation:  return os << "mutation";
    case EvoDecision::combine:  return os << "combine";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(decision);
}

static std::ostream& operator<< (std::ostream& os, const FileLogLevel& log_level) {
  switch (log_level) {
    case FileLogLevel::no_file_logging:  return os << "no file logging";
    case FileLogLevel::write_imbalance_km1:  return os << "only write imbalance and km1 goal to files";
    case FileLogLevel::write_imbalance_km1_target :  return os << "write imbalance, target imbalance and km1 goal to files";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(log_level);
}

static std::ostream& operator<< (std::ostream& os, const RatingPartitionPolicy& policy) {
  switch (policy) {
    case RatingPartitionPolicy::normal: return os << "normal";
    case RatingPartitionPolicy::evolutionary: return os << "evolutionary";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(policy);
}

static std::ostream& operator<< (std::ostream& os, const Mode& mode) {
  switch (mode) {
    case Mode::recursive_bisection: return os << "recursive";
    case Mode::direct_kway: return os << "direct";
    case Mode::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(mode);
}

static std::ostream& operator<< (std::ostream& os, const ContextType& type) {
  if (type == ContextType::main) {
    return os << "main";
  } else {
    return os << "ip";
  }
  return os << static_cast<uint8_t>(type);
}

static std::ostream& operator<< (std::ostream& os, const CommunityPolicy& comm_policy) {
  switch (comm_policy) {
    case CommunityPolicy::use_communities: return os << "true";
    case CommunityPolicy::ignore_communities: return os << "false";
    case CommunityPolicy::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(comm_policy);
}

static std::ostream& operator<< (std::ostream& os, const HeavyNodePenaltyPolicy& heavy_hn_policy) {
  switch (heavy_hn_policy) {
    case HeavyNodePenaltyPolicy::multiplicative_penalty: return os << "multiplicative";
    case HeavyNodePenaltyPolicy::no_penalty: return os << "no_penalty";
    case HeavyNodePenaltyPolicy::edge_frequency_penalty: return os << "edge_frequency_penalty";
    case HeavyNodePenaltyPolicy::UNDEFINED: return os << "UNDEFINED";
  }
  return os << static_cast<uint8_t>(heavy_hn_policy);
}

static std::ostream& operator<< (std::ostream& os, const AcceptancePolicy& acceptance_policy) {
  switch (acceptance_policy) {
    case AcceptancePolicy::best: return os << "best";
    case AcceptancePolicy::best_prefer_unmatched: return os << "best_prefer_unmatched";
    case AcceptancePolicy::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(acceptance_policy);
}

static std::ostream& operator<< (std::ostream& os, const FixVertexContractionAcceptancePolicy& acceptance_policy) {
  switch (acceptance_policy) {
    case FixVertexContractionAcceptancePolicy::free_vertex_only: return os << "free_vertex_only";
    case FixVertexContractionAcceptancePolicy::fixed_vertex_allowed: return os << "fixed_vertex_allowed";
    case FixVertexContractionAcceptancePolicy::equivalent_vertices: return os << "equivalent_vertices";
    case FixVertexContractionAcceptancePolicy::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(acceptance_policy);
}

static std::ostream& operator<< (std::ostream& os, const RatingFunction& func) {
  switch (func) {
    case RatingFunction::heavy_edge: return os << "heavy_edge";
    case RatingFunction::edge_frequency: return os << "edge_frequency";
    case RatingFunction::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(func);
}

static std::ostream& operator<< (std::ostream& os, const IgnoreMaxNodeWeight& ignore) {
  switch (ignore) {
    case IgnoreMaxNodeWeight::ignore_max_node_weight: return os << "ignore_max_node_weight";
    case IgnoreMaxNodeWeight::respect_max_node_weight: return os << "respect_max_node_weight";
    case IgnoreMaxNodeWeight::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(ignore);
}

static std::ostream& operator<< (std::ostream& os, const ContractCommunities& contract) {
  switch (contract) {
    case ContractCommunities::contract_communities: return os << "contract_communities";
    case ContractCommunities::dont_contract_communities: return os << "dont_contract_communities";
    case ContractCommunities::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(contract);
}

static std::ostream& operator<< (std::ostream& os, const Objective& objective) {
  switch (objective) {
    case Objective::cut: return os << "cut";
    case Objective::km1: return os << "km1";
    case Objective::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(objective);
}

static std::ostream& operator<< (std::ostream& os, const InitialPartitioningTechnique& technique) {
  switch (technique) {
    case InitialPartitioningTechnique::flat: return os << "flat";
    case InitialPartitioningTechnique::multilevel: return os << "multilevel";
    case InitialPartitioningTechnique::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(technique);
}

static std::ostream& operator<< (std::ostream& os, const CoarseningAlgorithm& algo) {
  switch (algo) {
    case CoarseningAlgorithm::heavy_full: return os << "heavy_full";
    case CoarseningAlgorithm::heavy_lazy: return os << "heavy_lazy";
    case CoarseningAlgorithm::ml_style: return os << "ml_style";
    case CoarseningAlgorithm::do_nothing: return os << "do_nothing";
    case CoarseningAlgorithm::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(algo);
}

static std::ostream& operator<< (std::ostream& os, const RefinementAlgorithm& algo) {
  switch (algo) {
    case RefinementAlgorithm::twoway_fm: return os << "twoway_fm";
    case RefinementAlgorithm::kway_fm: return os << "kway_fm";
    case RefinementAlgorithm::kway_fm_km1: return os << "kway_fm_km1";
    case RefinementAlgorithm::balance_approaching_kway_fm_km1: return os << "balance_approaching_kway_fm_km1";
    case RefinementAlgorithm::imbalance_holding_kway_fm_km1: return os << "imbalance_holding_kway_fm_km1";
    case RefinementAlgorithm::twoway_hyperflow_cutter: return os << "twoway_hyperflow_cutter";
    case RefinementAlgorithm::twoway_fm_hyperflow_cutter: return os << "twoway_fm_hyperflow_cutter";
    case RefinementAlgorithm::kway_hyperflow_cutter: return os << "kway_hyperflow_cutter";
    case RefinementAlgorithm::kway_fm_hyperflow_cutter: return os << "kway_fm_hyperflow_cutter";
    case RefinementAlgorithm::kway_fm_hyperflow_cutter_km1: return os << "kway_fm_hyperflow_cutter_km1";
    case RefinementAlgorithm::do_nothing: return os << "do_nothing";
    case RefinementAlgorithm::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(algo);
}

static std::ostream& operator<< (std::ostream& os, const InitialPartitionerAlgorithm& algo) {
  switch (algo) {
    case InitialPartitionerAlgorithm::greedy_sequential: return os << "greedy_sequential";
    case InitialPartitionerAlgorithm::greedy_global: return os << "greedy_global";
    case InitialPartitionerAlgorithm::greedy_round: return os << "greedy_round";
    case InitialPartitionerAlgorithm::greedy_sequential_maxpin: return os << "greedy_maxpin";
    case InitialPartitionerAlgorithm::greedy_global_maxpin: return os << "greedy_global_maxpin";
    case InitialPartitionerAlgorithm::greedy_round_maxpin: return os << "greedy_round_maxpin";
    case InitialPartitionerAlgorithm::greedy_sequential_maxnet: return os << "greedy_maxnet";
    case InitialPartitionerAlgorithm::greedy_global_maxnet: return os << "greedy_global_maxnet";
    case InitialPartitionerAlgorithm::greedy_round_maxnet: return os << "greedy_round_maxnet";
    case InitialPartitionerAlgorithm::bfs: return os << "bfs";
    case InitialPartitionerAlgorithm::random: return os << "random";
    case InitialPartitionerAlgorithm::lp: return os << "lp";
    case InitialPartitionerAlgorithm::bin_packing: return os << "bin_packing";
    case InitialPartitionerAlgorithm::direct_k: return os << "direct_k";
    case InitialPartitionerAlgorithm::pool: return os << "pool";
    case InitialPartitionerAlgorithm::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(algo);
}

static std::ostream& operator<< (std::ostream& os, const LouvainEdgeWeight& weight) {
  switch (weight) {
    case LouvainEdgeWeight::hybrid: return os << "hybrid";
    case LouvainEdgeWeight::uniform: return os << "uniform";
    case LouvainEdgeWeight::non_uniform: return os << "non_uniform";
    case LouvainEdgeWeight::degree: return os << "degree";
    case LouvainEdgeWeight::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(weight);
}

static std::ostream& operator<< (std::ostream& os, const RefinementStoppingRule& rule) {
  switch (rule) {
    case RefinementStoppingRule::simple: return os << "simple";
    case RefinementStoppingRule::adaptive_opt: return os << "adaptive_opt";
    case RefinementStoppingRule::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(rule);
}

static std::ostream& operator<< (std::ostream& os, const BalancingFlowModel& model) {
  switch (model) {
    case BalancingFlowModel::quotient_flow : return os << "quotient_flow";
    case BalancingFlowModel::laplace_matrix : return os << "laplace_matrix";
    case BalancingFlowModel::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(model);
}

static std::ostream& operator<< (std::ostream& os, const FlowExecutionMode& mode) {
  switch (mode) {
    case FlowExecutionMode::constant: return os << "constant";
    case FlowExecutionMode::multilevel: return os << "multilevel";
    case FlowExecutionMode::exponential: return os << "exponential";
    case FlowExecutionMode::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(mode);
}

static std::ostream& operator<< (std::ostream& os, const AcceptanceRule& rule) {
  switch (rule) {
    case AcceptanceRule::balance_approaching: return os << "balance_approaching";
    case AcceptanceRule::imbalance_holding: return os << "imbalance_holding";
    case AcceptanceRule::staircase: return os << "staircase";
    case AcceptanceRule::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(rule);
}


std::ostream& operator<< (std::ostream& os, const BinPackingAlgorithm& bp_algo) {
  switch (bp_algo) {
    case BinPackingAlgorithm::worst_fit: return os << "worst_fit";
    case BinPackingAlgorithm::first_fit: return os << "first_fit";
    case BinPackingAlgorithm::UNDEFINED: return os << "UNDEFINED";
      // omit default case to trigger compiler warning for missing cases
  }
  return os << static_cast<uint8_t>(bp_algo);
}

static EvoMutateStrategy mutateStrategyFromString(const std::string& strat) {
  if (strat == "new-initial-partitioning-vcycle") {
    return EvoMutateStrategy::new_initial_partitioning_vcycle;
  } else if (strat == "vcycle") {
    return EvoMutateStrategy::vcycle;
  }
  LOG << "No valid mutate strategy. ";
  exit(0);
}
static EvoCombineStrategy combineStrategyFromString(const std::string& strat) {
  if (strat == "basic") {
    return EvoCombineStrategy::basic;
  } else if (strat == "edge-frequency") {
    return EvoCombineStrategy::edge_frequency;
  }
  LOG << "No valid combine strategy. ";
  exit(0);
}
static EvoReplaceStrategy replaceStrategyFromString(const std::string& strat) {
  if (strat == "worst") {
    return EvoReplaceStrategy::worst;
  } else if (strat == "diverse") {
    return EvoReplaceStrategy::diverse;
  } else if (strat == "strong-diverse") {
    return EvoReplaceStrategy::strong_diverse;
  }
  LOG << "No valid replace strategy. ";
  exit(0);
}

static AcceptancePolicy acceptanceCriterionFromString(const std::string& crit) {
  if (crit == "best") {
    return AcceptancePolicy::best;
  } else if (crit == "best_prefer_unmatched") {
    return AcceptancePolicy::best_prefer_unmatched;
  }
  LOG << "No valid acceptance criterion for rating.";
  exit(0);
}

static IgnoreMaxNodeWeight ignoreMaxNodeWeightFromString(const std::string& crit) {
  if (crit == "true") {
    return IgnoreMaxNodeWeight::ignore_max_node_weight;
  } else if (crit == "false") {
    return IgnoreMaxNodeWeight::respect_max_node_weight;
  }
  LOG << "No valid boolean if max node weight should be ignored while rating.";
  exit(0);
}

static ContractCommunities contractCommunitiesFromString(const std::string& crit) {
  if (crit == "true") {
    return ContractCommunities::contract_communities;
  } else if (crit == "false") {
    return ContractCommunities::dont_contract_communities;
  }
  LOG << "No valid boolean if communities should be contracted while rating.";
  exit(0);
}

static RatingPartitionPolicy ratingPartitionPolicyFromString(const std::string& partition) {
  if (partition == "normal") {
    return RatingPartitionPolicy::normal;
  } else if (partition == "evolutionary") {
    return RatingPartitionPolicy::evolutionary;
  }
  LOG << "No valid partition policy for rating.";
  exit(0);
  return RatingPartitionPolicy::normal;
}

static FixVertexContractionAcceptancePolicy fixedVertexAcceptanceCriterionFromString(const std::string& crit) {
  if (crit == "free_vertex_only") {
    return FixVertexContractionAcceptancePolicy::free_vertex_only;
  } else if (crit == "fixed_vertex_allowed") {
    return FixVertexContractionAcceptancePolicy::fixed_vertex_allowed;
  } else if (crit == "equivalent_vertices") {
    return FixVertexContractionAcceptancePolicy::equivalent_vertices;
  }
  LOG << "No valid fixed vertex acceptance criterion for rating.";
  exit(0);
}

static HeavyNodePenaltyPolicy heavyNodePenaltyFromString(const std::string& penalty) {
  if (penalty == "multiplicative") {
    return HeavyNodePenaltyPolicy::multiplicative_penalty;
  } else if (penalty == "no_penalty") {
    return HeavyNodePenaltyPolicy::no_penalty;
  } else if (penalty == "edge_frequency_penalty") {
    return HeavyNodePenaltyPolicy::edge_frequency_penalty;
    // omit default case to trigger compiler warning for missing cases
  }
  LOG << "No valid edge penalty policy for rating.";
  exit(0);
  return HeavyNodePenaltyPolicy::multiplicative_penalty;
}

static RatingFunction ratingFunctionFromString(const std::string& function) {
  if (function == "heavy_edge") {
    return RatingFunction::heavy_edge;
  } else if (function == "edge_frequency") {
    return RatingFunction::edge_frequency;
  }
  LOG << "No valid rating function for rating.";
  exit(0);
  return RatingFunction::heavy_edge;
}

static RefinementStoppingRule stoppingRuleFromString(const std::string& rule) {
  if (rule == "simple") {
    return RefinementStoppingRule::simple;
  } else if (rule == "adaptive_opt") {
    return RefinementStoppingRule::adaptive_opt;
  }
  LOG << "No valid stopping rule for FM.";
  exit(0);
  return RefinementStoppingRule::simple;
}

static BalancingFlowModel balancingFlowModelFromString(const std::string& model) {
  if (model == "laplace_matrix") {
    return BalancingFlowModel::laplace_matrix;
  } else if (model == "quotient_flow") {
    return BalancingFlowModel::quotient_flow;
  }
  LOG << "No valid flow model for flow balancing.";
  exit(0);
  return BalancingFlowModel::UNDEFINED;
}

static CoarseningAlgorithm coarseningAlgorithmFromString(const std::string& type) {
  if (type == "heavy_full") {
    return CoarseningAlgorithm::heavy_full;
  } else if (type == "heavy_lazy") {
    return CoarseningAlgorithm::heavy_lazy;
  } else if (type == "ml_style") {
    return CoarseningAlgorithm::ml_style;
  } else if (type == "do_nothing") {
    return CoarseningAlgorithm::do_nothing;
  }
  LOG << "Illegal option:" << type;
  exit(0);
  return CoarseningAlgorithm::heavy_lazy;
}

static RefinementAlgorithm refinementAlgorithmFromString(const std::string& type) {
  if (type == "twoway_fm") {
    return RefinementAlgorithm::twoway_fm;
  } else if (type == "kway_fm") {
    return RefinementAlgorithm::kway_fm;
  } else if (type == "kway_fm_km1") {
    return RefinementAlgorithm::kway_fm_km1;
  } else if (type == "twoway_hyperflow_cutter") {
    return RefinementAlgorithm::twoway_hyperflow_cutter;
  } else if (type == "kway_hyperflow_cutter") {
    return RefinementAlgorithm::kway_hyperflow_cutter;
  } else if (type == "kway_fm_hyperflow_cutter") {
    return RefinementAlgorithm::kway_fm_hyperflow_cutter;
  } else if (type == "twoway_fm_hyperflow_cutter") {
    return RefinementAlgorithm::twoway_fm_hyperflow_cutter;
  } else if (type == "kway_fm_hyperflow_cutter_km1") {
    return RefinementAlgorithm::kway_fm_hyperflow_cutter_km1;
  } else if (type == "balance_approaching_kway_fm_km1") {
    return RefinementAlgorithm::balance_approaching_kway_fm_km1;
  } else if (type == "imbalance_holding_kway_fm_km1") {
    return RefinementAlgorithm::imbalance_holding_kway_fm_km1;
  } else if (type == "do_nothing") {
    return RefinementAlgorithm::do_nothing;
  }
  LOG << "Illegal option:" << type;
  exit(0);
  return RefinementAlgorithm::kway_fm;
}

static InitialPartitionerAlgorithm initialPartitioningAlgorithmFromString(const std::string& mode) {
  if (mode == "greedy_sequential") {
    return InitialPartitionerAlgorithm::greedy_sequential;
  } else if (mode == "greedy_global") {
    return InitialPartitionerAlgorithm::greedy_global;
  } else if (mode == "greedy_round") {
    return InitialPartitionerAlgorithm::greedy_round;
  } else if (mode == "greedy_sequential_maxpin") {
    return InitialPartitionerAlgorithm::greedy_sequential_maxpin;
  } else if (mode == "greedy_global_maxpin") {
    return InitialPartitionerAlgorithm::greedy_global_maxpin;
  } else if (mode == "greedy_round_maxpin") {
    return InitialPartitionerAlgorithm::greedy_round_maxpin;
  } else if (mode == "greedy_sequential_maxnet") {
    return InitialPartitionerAlgorithm::greedy_sequential_maxnet;
  } else if (mode == "greedy_global_maxnet") {
    return InitialPartitionerAlgorithm::greedy_global_maxnet;
  } else if (mode == "greedy_round_maxnet") {
    return InitialPartitionerAlgorithm::greedy_round_maxnet;
  } else if (mode == "lp") {
    return InitialPartitionerAlgorithm::lp;
  } else if (mode == "bfs") {
    return InitialPartitionerAlgorithm::bfs;
  } else if (mode == "random") {
    return InitialPartitionerAlgorithm::random;
  } else if (mode == "pool") {
    return InitialPartitionerAlgorithm::pool;
  } else if (mode == "bin_packing") {
    return InitialPartitionerAlgorithm::bin_packing;
  } else if (mode == "direct_k") {
    return InitialPartitionerAlgorithm::direct_k;
  }
  LOG << "Illegal option:" << mode;
  exit(0);
  return InitialPartitionerAlgorithm::greedy_global;
}

static InitialPartitioningTechnique initialPartitioningTechniqueFromString(const std::string& technique) {
  if (technique == "flat") {
    return InitialPartitioningTechnique::flat;
  } else if (technique == "multi") {
    return InitialPartitioningTechnique::multilevel;
  }
  LOG << "Illegal option:" << technique;
  exit(0);
  return InitialPartitioningTechnique::multilevel;
}

static LouvainEdgeWeight edgeWeightFromString(const std::string& type) {
  if (type == "hybrid") {
    return LouvainEdgeWeight::hybrid;
  } else if (type == "uniform") {
    return LouvainEdgeWeight::uniform;
  } else if (type == "non_uniform") {
    return LouvainEdgeWeight::non_uniform;
  } else if (type == "degree") {
    return LouvainEdgeWeight::degree;
  }
  LOG << "Illegal option:" << type;
  exit(0);
  return LouvainEdgeWeight::uniform;
}

static FileLogLevel fileLogLevelFromString(std::string log_level) {
  if (log_level == "km1-imbalance-target") {
    return FileLogLevel::write_imbalance_km1_target;
  } else if (log_level == "km1-imbalance") {
    return FileLogLevel::write_imbalance_km1;
  } if (log_level == "none") {
    return FileLogLevel::no_file_logging;
  }
  LOG << "Illegal option:" << log_level;
  exit(0);
  return FileLogLevel::no_file_logging;
}

static Mode modeFromString(const std::string& mode) {
  if (mode == "recursive") {
    return Mode::recursive_bisection;
  } else if (mode == "direct") {
    return Mode::direct_kway;
  }
  LOG << "Illegal option:" << mode;
  exit(0);
  return Mode::direct_kway;
}

static FlowExecutionMode flowExecutionPolicyFromString(const std::string& mode) {
  if (mode == "constant") {
    return FlowExecutionMode::constant;
  } else if (mode == "multilevel") {
    return FlowExecutionMode::multilevel;
  } else if (mode == "exponential") {
    return FlowExecutionMode::exponential;
  }
  LOG << "No valid flow execution mode.";
  exit(0);
  return FlowExecutionMode::exponential;
}

static AcceptanceRule flowAcceptancePolicyFromString(const std::string& mode) {
  if (mode == "balance_approaching") {
    return AcceptanceRule::balance_approaching;
  } else if (mode == "imbalance_holding") {
    return AcceptanceRule::imbalance_holding;
  } else if (mode == "staircase") {
    return AcceptanceRule::staircase;
  }
  LOG << "No valid acceptance policy.";
  exit(0);
  return AcceptanceRule::UNDEFINED;
}

static BinPackingAlgorithm binPackingAlgorithmFromString(const std::string& type) {
  if (type == "worst_fit") {
    return BinPackingAlgorithm::worst_fit;
  } else if (type == "first_fit") {
    return BinPackingAlgorithm::first_fit;
  }
  LOG << "Illegal option:" << type;
  exit(0);
  return BinPackingAlgorithm::worst_fit;
}
}  // namespace kahypar
