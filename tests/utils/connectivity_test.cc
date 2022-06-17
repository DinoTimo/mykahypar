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

#include "gmock/gmock.h"
#include "kahypar/io/hypergraph_io.h"
#include "kahypar/partition/metrics.h"

using ::testing::Eq;

namespace kahypar {
namespace metrics {
TEST(IsGraphConnected, WorksAsExpected) {
  Hypergraph hg(kahypar::io::createHypergraphFromFile("test_instances/ISPD98_ibm01.hgr", 4));
  ASSERT_THAT(metrics::amountConnectedComponents(hg), Eq(1));
  
  Hypergraph hg2(kahypar::io::createHypergraphFromFile("test_instances//Oregon-1.mtx.hgr", 4));
  ASSERT_THAT(metrics::amountConnectedComponents(hg2), Eq(319));
  
 }

TEST(VectorElementFinding, WorksAsExpected) {
  std::vector<int> vec({4, 5, 6, 20, 4, 3});
  ASSERT_THAT(metrics::internal::findElem(vec, 4), Eq(0));
  ASSERT_THAT(metrics::internal::findElem(vec, 5), Eq(1));
  ASSERT_THAT(metrics::internal::findElem(vec, 6), Eq(2));
  ASSERT_THAT(metrics::internal::findElem(vec, 20), Eq(3));
  ASSERT_THAT(metrics::internal::findElem(vec, 12), Eq(-1));

}

}  // namespace connectivity
}  // namespace kahypar
