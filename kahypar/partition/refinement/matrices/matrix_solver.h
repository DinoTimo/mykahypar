/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
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

#include <vector>

#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/macros.h"
#include "kahypar/partition/refinement/matrices/matrix.h"

namespace kahypar {
namespace matrices {



template <typename Derived = Mandatory>
class matrix_solver {
  
  public:
    //template <typename MatrixType>
    std::vector<double> solve(full_square_matrix<double> A, std::vector<double> b, size_t N) {
      return static_cast<Derived*>(this)->solveImpl(A, b, N);
    }

};

/**
 * This Decomposition class only does A = LU with no pivot matrix, as it is expected for the main diagonal
 * of A to be strictly non 0.
 * 
 */
class LU_Decomp_matrix_solver final : public matrix_solver<LU_Decomp_matrix_solver> {

  public:
    //template <typename MatrixType>
    std::vector<double> solveImpl(full_square_matrix<double> A, std::vector<double> b, size_t N) {
      std::pair<full_square_matrix<double>, full_square_matrix<double>> LU_pair(LU_decompose(A, N)); // TODO use proper matrix types
      full_square_matrix<double> L(LU_pair.first);
      full_square_matrix<double> U(LU_pair.second);
      std::vector<double> y = forwardSubstitute(L, b, N);
      std::vector<double> x = backwardSubstitute(U, y, N);
      return x;
    }

    //Use Doolittle algorithm, no pivoting, not in-place, return pair of matrices, not one matrix where both upper/lower matrix are implicitly stored inside
    //template <typename MatrixType>
    std::pair<full_square_matrix<double>, full_square_matrix<double>> LU_decompose(const full_square_matrix<double> A, const size_t N) {
      // TODO use lower matrix
      full_square_matrix<double> L(N, 1);
      // TODO use upper matrix
      full_square_matrix<double> U(A);
      for (size_t i = 0; i < N - 1; i++) {
        for (size_t k = i + 1; k < N; k++) {
          ASSERT(U.at(i, i) > 0.0001 || U.at(i, i) < -0.0001);
          L.at(k, i) = U.at(k, i) / U.at(i, i);
          for (size_t j = i; j < N; j++ ) {
            U.at(k, j) -= L.at(k, i) * U.at(i, j);
          }
        }
      }
      static std::pair<full_square_matrix<double>, full_square_matrix<double>> pair(L, U);
      return pair;
    }

    //TODO use lower matrix type and remove template argument
    //template <typename MatrixType>
    std::vector<double> forwardSubstitute(full_square_matrix<double> A, std::vector<double> b, size_t N) {
      ASSERT(b.size() == N);
      std::vector<double> x(N, 0);
      for (size_t i = 0; i < N; i++) {
        x[i] = b[i];
        for (size_t j = 0; j < i; j++) {
          x[i] -= A.at(i, j) * x[j];
        }
        x[i] /= A.at(i, i); 
      }
      return x;
    }
    
    //TODO use upper matrix type and remove template argument
    //template <typename MatrixType>
    std::vector<double> backwardSubstitute(full_square_matrix<double> A, std::vector<double> b, size_t n) {
      //int is used instead of size_t, since int is unsigned and the loop going to zero.
      int N = static_cast<int>(n);
      ASSERT(b.size() == n);
      std::vector<double> x(n, 0);
      for (int i = N - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < N; j++) {
          x[i] -= A.at(i, j) * x[j];
        }
        x[i] /= A.at(i, i); 
      }
      return x;
    }
    
};

}
}