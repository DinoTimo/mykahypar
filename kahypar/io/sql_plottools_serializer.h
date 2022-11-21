/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014-2016 Sebastian Schlag <sebastian.schlag@kit.edu>
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
#include <array>
#include <chrono>
#include <fstream>
#include <limits>
#include <string>

#include "kahypar/definitions.h"
#include "kahypar/git_revision.h"
#include "kahypar/partition/context.h"
#include "kahypar/partition/evolutionary/individual.h"
#include "kahypar/partition/metrics.h"
#include "kahypar/partition/partitioner.h"
#include "kahypar/partition/bin_packing/bin_packing_utils.h"

namespace kahypar {
namespace io {
namespace serializer {
static inline void serialize(const Context& context, const Hypergraph& hypergraph,
                             const std::chrono::duration<double>& elapsed_seconds,
                             const size_t iteration = 0, bool interrupted = false) {
  if (!context.partition.sp_process_output) {
    return;
  }

  std::ostringstream oss;
  oss << "RESULT"
      << " graph=" << context.partition.graph_filename.substr(context.partition.graph_filename.find_last_of('/') + 1)
      << " numHNs=" << hypergraph.initialNumNodes()
      << " numHEs=" << hypergraph.initialNumEdges()
      << " numComponents=" << metrics::amountConnectedComponents(hypergraph)
      << " " << hypergraph.typeAsString();
  oss << " mode=" << context.partition.mode
      << " objective=" << context.partition.objective
      << " k=" << context.partition.k
      << " epsilon=" << context.partition.epsilon
      << " seed=" << context.partition.seed
      << " total_graph_weight=" << hypergraph.totalWeight()
      << " he_size_threshold=" << context.partition.hyperedge_size_threshold;

  oss << " pre_enable_deduplication=" << std::boolalpha
      << context.preprocessing.enable_deduplication
      << " pre_enable_min_hash_sparsifier=" << std::boolalpha
      << context.preprocessing.enable_min_hash_sparsifier
      << " pre_min_hash_max_hyperedge_size="
      << context.preprocessing.min_hash_sparsifier.max_hyperedge_size
      << " pre_min_hash_max_cluster_size="
      << context.preprocessing.min_hash_sparsifier.max_cluster_size
      << " pre_min_hash_min_cluster_size="
      << context.preprocessing.min_hash_sparsifier.min_cluster_size
      << " pre_min_hash_num_hash_functions="
      << context.preprocessing.min_hash_sparsifier.num_hash_functions
      << " pre_min_hash_combined_num_hash_functions="
      << context.preprocessing.min_hash_sparsifier.combined_num_hash_functions
      << " pre_min_sparsifier_is_active="
      << context.preprocessing.min_hash_sparsifier.is_active
      << " pre_min_sparsifier_activation_median_he_size="
      << context.preprocessing.min_hash_sparsifier.min_median_he_size
      << " enable_community_detection=" << std::boolalpha
      << context.preprocessing.enable_community_detection
      << " enable_louvain_in_initial_partitioning=" << std::boolalpha
      << context.preprocessing.community_detection.enable_in_initial_partitioning
      << " max_louvain_pass_iterations="
      << context.preprocessing.community_detection.max_pass_iterations
      << " min_louvain_eps_improvement="
      << context.preprocessing.community_detection.min_eps_improvement
      << " louvain_edge_weight=" << context.preprocessing.community_detection.edge_weight
      << " reuse_community_structure=" << std::boolalpha
      << context.preprocessing.community_detection.reuse_communities
      << " contract_communities="
      << context.coarsening.rating.contract_communities
      << " ignore_max_node_weight="
      << context.coarsening.rating.ignore_max_node_weight
      << " coarsening_algo=" << context.coarsening.algorithm
      << " coarsening_max_allowed_weight_multiplier="
      << context.coarsening.max_allowed_weight_multiplier
      << " coarsening_contraction_limit_multiplier="
      << context.coarsening.contraction_limit_multiplier
      << " coarsening_hypernode_weight_fraction=" << context.coarsening.hypernode_weight_fraction
      << " coarsening_max_allowed_node_weight=" << context.coarsening.max_allowed_node_weight
      << " coarsening_contraction_limit=" << context.coarsening.contraction_limit
      << " coarsening_rating_function=" << context.coarsening.rating.rating_function
      << " coarsening_rating_use_communities="
      << context.coarsening.rating.community_policy
      << " coarsening_rating_heavy_node_penalty="
      << context.coarsening.rating.heavy_node_penalty_policy
      << " coarsening_rating_acceptance_policy="
      << context.coarsening.rating.acceptance_policy
      << " coarsening_rating_fixed_vertex_acceptance_policy="
      << context.coarsening.rating.fixed_vertex_acceptance_policy;
  if (!context.partition_evolutionary &&
      !context.partition.time_limited_repeated_partitioning) {
    oss << " " << context.stats.serialize().str();
  }
  oss << " git=" << STR(KaHyPar_BUILD_VERSION)
      << std::endl;

  std::cout << oss.str() << std::endl;
}

static inline void serializeEvolutionary(const Context& context, const Hypergraph& hg) {
  std::ostringstream oss;
  if (context.partition.quiet_mode) {
    return;
  }
  EvoCombineStrategy combine_strat = EvoCombineStrategy::UNDEFINED;
  EvoMutateStrategy mutate_strat = EvoMutateStrategy::UNDEFINED;
  switch (context.evolutionary.action.decision()) {
    case EvoDecision::combine:
      combine_strat = context.evolutionary.combine_strategy;
      break;
    case EvoDecision::mutation:
      mutate_strat = context.evolutionary.mutate_strategy;
      break;
    case EvoDecision::normal:
      break;
    default:
      LOG << "Trying to print a nonintentional action:" << context.evolutionary.action.decision();
  }

  std::string graph_name = context.partition.graph_filename;
  std::string truncated_graph_name = graph_name.substr(graph_name.find_last_of("/") + 1);
  oss << "RESULT "
      << "connectivity=" << metrics::km1(hg)
      << " action=" << context.evolutionary.action.decision()
      << " time-total=" << Timer::instance().evolutionaryResult().total_evolutionary
      << " iteration=" << context.evolutionary.iteration
      << " replace-strategy=" << context.evolutionary.replace_strategy
      << " combine-strategy=" << combine_strat
      << " mutate-strategy=" << mutate_strat
      << " population-size=" << context.evolutionary.population_size
      << " mutation-chance=" << context.evolutionary.mutation_chance
      << " diversify-interval=" << context.evolutionary.diversify_interval
      << " dynamic-pop-size=" << context.evolutionary.dynamic_population_size
      << " dynamic-pop-percentile=" << context.evolutionary.dynamic_population_amount_of_time
      << " seed=" << context.partition.seed
      << " graph-name=" << truncated_graph_name
      << " SOED=" << metrics::soed(hg)
      << " cut=" << metrics::hyperedgeCut(hg)
      << " absorption=" << metrics::absorption(hg)
      << " imbalance=" << metrics::imbalance(hg, context)
      << " k=" << context.partition.k
      << std::endl;

  std::cout << oss.str() << std::endl;
}
}  // namespace serializer
}  // namespace io
}  // namespace kahypar
