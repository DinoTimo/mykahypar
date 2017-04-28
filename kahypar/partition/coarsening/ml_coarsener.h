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

#include <limits>
#include <string>
#include <vector>


#include "kahypar/datastructure/sparse_map.h"
#include "kahypar/definitions.h"
#include "kahypar/macros.h"
#include "kahypar/partition/coarsening/heavy_edge_rater.h"
#include "kahypar/partition/coarsening/policies/rating_acceptance_policy.h"
#include "kahypar/partition/coarsening/policies/rating_heavy_node_penalty_policy.h"
#include "kahypar/partition/coarsening/policies/rating_tie_breaking_policy.h"
#include "kahypar/partition/preprocessing/louvain.h"
#include "kahypar/partition/preprocessing/modularity.h"
#include "kahypar/utils/float_compare.h"

namespace kahypar {
class MLCoarsener final : public ICoarsener,
                          private VertexPairCoarsenerBase<>{
 private:
  static constexpr bool debug = false;

  static constexpr HypernodeID kInvalidTarget = std::numeric_limits<HypernodeID>::max();

  using Base = VertexPairCoarsenerBase;
  using Rater = HeavyEdgeRater<RatingType, RandomRatingWins,
                               BestRatingPreferringUnmatched<RandomRatingWins>,
                               NoWeightPenalty>;
  using Rating = typename Rater::Rating;

 public:
  MLCoarsener(Hypergraph& hypergraph, const Context& context,
              const HypernodeWeight weight_of_heaviest_node) :
    Base(hypergraph, context, weight_of_heaviest_node),
    _rater(_hg, _context) { }

  ~MLCoarsener() override = default;

  MLCoarsener(const MLCoarsener&) = delete;
  MLCoarsener& operator= (const MLCoarsener&) = delete;

  MLCoarsener(MLCoarsener&&) = delete;
  MLCoarsener& operator= (MLCoarsener&&) = delete;

 private:
  void coarsenImpl(const HypernodeID limit) override final {
    int pass_nr = 0;
    std::vector<HypernodeID> current_hns;

    while (_hg.currentNumNodes() > limit) {
      DBG << V(pass_nr);
      DBG << V(_hg.currentNumNodes());
      DBG << V(_hg.currentNumEdges());

      _rater.resetMatches();
      current_hns.clear();

      const HypernodeID num_hns_before_pass = _hg.currentNumNodes();
      for (const HypernodeID& hn : _hg.nodes()) {
        current_hns.push_back(hn);
      }
      Randomize::instance().shuffleVector(current_hns, current_hns.size());
      // std::sort(std::begin(current_hns), std::end(current_hns),
      //           [&](const HypernodeID l, const HypernodeID r) {
      //             return _hg.nodeDegree(l) < _hg.nodeDegree(r);
      //           });

      for (const HypernodeID& hn : current_hns) {
        if (_hg.nodeIsEnabled(hn)) {
          const Rating rating = _rater.rate(hn);

          if (rating.target != kInvalidTarget) {
            _rater.markAsMatched(hn);
            _rater.markAsMatched(rating.target);
            // if (_hg.nodeDegree(hn) > _hg.nodeDegree(rating.target)) {

            performContraction(hn, rating.target);
            // } else {
            //   contract(rating.target, hn);
            // }
          }

          if (_hg.currentNumNodes() <= limit) {
            break;
          }
        }
      }

      if (num_hns_before_pass == _hg.currentNumNodes()) {
        break;
      }
      ++pass_nr;
    }
    Stats::instance().addToTotal(_context, "hns_after_coarsening", _hg.currentNumNodes());
  }

  bool uncoarsenImpl(IRefiner& refiner) override final {
    return doUncoarsen(refiner);
  }

  std::string policyStringImpl() const override final {
    return std::string("");
  }

  using Base::_pq;
  using Base::_hg;
  using Base::_context;
  using Base::_history;
  Rater _rater;
};
}  // namespace kahypar
