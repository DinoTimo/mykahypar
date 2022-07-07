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

namespace kahypar {
namespace matrices {

template <typename Entry = Mandatory, typename Derived = Mandatory>
class matrix {

  public:
    matrix() = default;

    matrix(const matrix&) = delete;
    matrix& operator= (const matrix&) = delete;

    matrix(matrix&&) = delete;
    matrix& operator= (matrix&&) = delete;

    ~matrix() = default;
    
    Entry& at(size_t i, size_t j) {
      return static_cast<Derived*>(this)->atImpl(i, j);
    }

};

template<typename Entry = Mandatory>
class full_square_matrix final : public matrix<Entry, full_square_matrix<Entry>> {
  using Base = matrix<Entry, full_square_matrix<Entry>>;

  private:
    std::vector<Entry> elements;
    size_t size;

  public:
    full_square_matrix& operator= (full_square_matrix& other) {
      this->elements = other.elements;
      this->size = other.size;
    }

    full_square_matrix(size_t N) : Base(),
     elements(N*N, 0),
     size(N) { }

    full_square_matrix(size_t N, Entry diag_value) : Base(),
     elements(N*N, 0),
     size(N) { 
      for (size_t i = 0; i < N; i++) {
          atImpl(i, i) = diag_value; 
      }
     }



    full_square_matrix(const full_square_matrix<Entry>& other) : Base(),
     elements(other.elements),
     size(other.size) { }

    full_square_matrix(const std::vector<Entry>& vector) : Base(),
     elements(vector),
     size(static_cast<size_t>(std::sqrt(static_cast<double>(vector.size())))) {
      ASSERT(size * size == elements.size());
    }

    Entry& atImpl(const size_t i, const size_t j) {
      ASSERT((0 <= i) && (i < size));
      ASSERT((0 <= j) && (j < size));
      size_t index = i * size + j;
      ASSERT((index >= 0) && (index < elements.size()));
      return elements[index];
    }

    Entry read(const size_t i, const size_t j) const {
      ASSERT((0 <= i) && (i < size));
      ASSERT((0 <= j) && (j < size));
      size_t index = i * size + j;
      ASSERT((index >= 0) && (index < elements.size()));
      return elements[index];
    }

    void print() const {
      for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
          std::cout << read(i, j) << ", ";
        } 
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

};



}
} //namespace kahypar