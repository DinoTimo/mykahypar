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

#include <vector>

#include "gmock/gmock.h"
#include "kahypar/macros.h"
#include "kahypar/partition/refinement/matrices/matrix.h"
#include "kahypar/partition/refinement/matrices/matrix_solver.h"
#include "kahypar/utils/float_compare.h"
#include "kahypar/utils/logger.h"

using ::testing::Test;

namespace kahypar {
namespace matrices {

TEST(ModifyingMatrix, WorksAsExpected) {
  size_t size = 5;
  full_square_matrix<int> matrix(size);
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      ASSERT(matrix.at(i, j) == 0);
    } 
  }

  matrix.at(1, 2) = 5;
  ASSERT(matrix.at(1, 2) == 5);
  matrix.at(1, 2) = 3;
  ASSERT(matrix.at(1, 2) == 3);

  matrix.at(0, 0) = 5;
  ASSERT(matrix.at(0, 0) == 5);
  matrix.at(0, 0) = 3;
  ASSERT(matrix.at(0, 0) == 3);

  matrix.at(2, 2) = 12;
  ASSERT(matrix.at(2, 2) == 12);
}

bool almostEquals(double x, double y) {
  return std::abs(x - y) < 0.0001;
}

bool almostEquals(const full_square_matrix<double> matrix, const std::vector<double> vector, const size_t size) {
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      ASSERT(almostEquals(matrix.read(i, j), vector.at(i * size + j)));
    }
  }
  return true;
}

bool almostEquals(const std::vector<double> vector2, const std::vector<double> vector) {
  ASSERT(vector2.size() == vector.size());
  for (size_t i = 0; i < vector2.size(); i++) {
    ASSERT(almostEquals(vector2.at(i), vector.at(i)));
  }
  return true;
}

TEST(ForwardSubstitution, WorksAsExpected) {
  size_t size = 3;
  full_square_matrix<double> matrix(std::vector<double>
  { 3, 0, 0,
    2, 4, 0,
    1, 7, 9 }
  );
  std::vector<double> b{1, 2, 3};
  LU_Decomp_matrix_solver solver;
  std::vector<double> x = solver.forwardSubstitute(matrix, b, size);
  ASSERT(almostEquals(x, std::vector<double>{ 1.0 / 3.0, 1.0 / 3.0, 1.0 / 27.0}));
}

TEST(BackwardSubstitution, WorksAsExpected) {
  size_t size = 3;
  full_square_matrix<double> matrix(std::vector<double>
  { 3, 2, 8,
    0, 1, 8,
    0, 0, 9 }
  );
  std::vector<double> b{3, 2, 1};
  LU_Decomp_matrix_solver solver;
  std::vector<double> x = solver.backwardSubstitute(matrix, b, size);
  ASSERT(almostEquals(x, std::vector<double>{ -1.0 / 27.0,  10.0 / 9.0, 1.0 / 9.0}));
}


TEST(LU_Decomposition_test, WorksAsExpected) {
  size_t size = 3;
  full_square_matrix<double> matrix(std::vector<double> 
  { 15, 25, 5,
    6, 28, 11,
    3, 23, 12} 
  );
  LU_Decomp_matrix_solver solver;
  const std::pair<full_square_matrix<double>, full_square_matrix<double>> pair(solver.LU_decompose(matrix, size));
  const full_square_matrix<double> L = pair.first;
  const full_square_matrix<double> U = pair.second;
  ASSERT(almostEquals(L, std::vector<double>{ 1, 0, 0,
                                              0.4, 1, 0,
                                              0.2, 1, 1}, size));

  ASSERT(almostEquals(U, std::vector<double>{15, 25, 5,
                                              0, 18, 9,
                                              0, 0, 2}, size));
}

TEST(FullEquationSolvingTest, WorksAsExpected) {
  size_t size = 3;
  full_square_matrix<double> matrix(std::vector<double> 
  { 15, 25, 5,
    6, 28, 11,
    3, 23, 12} 
  );
  std::vector<double> b{80, 95, 85};
  LU_Decomp_matrix_solver solver;
  std::vector<double> x(solver.solve(matrix, b, size));
  ASSERT(almostEquals(x, std::vector<double>{1, 2, 3}));  
}
}
}