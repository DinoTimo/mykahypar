/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2015-2017 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "kahypar/meta/pack_expand.h"
#include "kahypar/definitions.h"

namespace kahypar {
class Logger {
 public:
  explicit Logger(const bool newline) :
    _newline(newline),
    _oss() { }

  template <typename Arg, typename ... Args>
  Logger(const bool newline, Arg&& arg, Args&& ... args) :
    _newline(newline),
    _oss() {
    _oss << "[" << std::forward<Arg>(arg);
    meta::expandPack((_oss << ":" << std::forward<Args>(args)) ...);
    _oss << "]: ";
  }

  template <typename T>
  Logger& operator<< (const T& output) {
    _oss << output << ' ';
    return *this;
  }

  template <typename T>
  Logger& operator<< (const T* output) {
    _oss << output << ' ';
    return *this;
  }


  Logger& operator<< (decltype(std::left)& output) {
    _oss << output;
    return *this;
  }

  Logger& operator<< (const decltype(std::setw(1))& output) {
    _oss << output;
    return *this;
  }

  ~Logger() {
    std::cout << _oss.str();
    if (_newline) {
      std::cout << std::endl;
    } else {
      std::cout << ' ';
    }
  }

 private:
  bool _newline;
  std::ostringstream _oss;
};

class LoggerVoidify {
 public:
  void operator& (Logger&) { }
};


inline std::string joinVector(const std::vector<std::string> vec, const std::string prefix, const std::string delim, const std::string postfix) {
  std::string str;
  if (vec.empty()) {
    return "";
  }
  str.append(prefix);
  unsigned long int index = 0;
  while (index < vec.size() - 1) {
    str.append(vec[index]);
    str.append(delim);
    index++;
  }
  str.append(vec[index]);
  str.append(postfix);
  return str;
}


template<typename Content>
inline std::string joinVector(const std::vector<Content> vec, const std::string prefix, const std::string delim, const std::string postfix) {
  std::vector<std::string> stringVec;
  std::for_each(vec.begin(), vec.end(), [&stringVec](Content elem) {
    stringVec.push_back(std::to_string(elem));
  });
  return joinVector(stringVec, prefix, delim, postfix);
}

inline void writeToFile(const std::string str, const std::string fileName) {
  std::ofstream file;
  file.open(fileName);
  file << str;
  file.close();
}

template <typename Content>
inline void writeVectorToFile(const std::vector<Content> vec, const std::string fileName) {
  writeToFile(joinVector(vec, "", "\n", ""), fileName);
}

inline void addToFile(const std::string newContent, const std::string filePath) {
  std::ofstream file;
  file.open(filePath, std::ios::app);
  file << newContent;
  file.close();
}

template <typename Number, typename Number2>
inline void writeVectorsToCSV(
  const std::vector<Number>& num_nodes,
  const std::vector<Number2>& km1s,
  const std::vector<double>& imbalances, 
  const std::vector<Number>& rebalance_steps,
  const std::vector<double>& target_imbalances,
  std::string file) {
  if (! (num_nodes.size() == km1s.size() && km1s.size() == imbalances.size() && imbalances.size() == rebalance_steps.size() && rebalance_steps.size() == target_imbalances.size()) ) {
    std::cerr << "Incorrect sizes of vectors to be plotted" << std::endl;
    return;
  }
  if (file.empty()) {
    std::cerr << "Invalid file";
    return;
  }
  addToFile("num_nodes,km1,imbalance,rebalance_step,target_imbalance\n", file);
  for (size_t i = 0; i < num_nodes.size(); i++) {
    std::stringstream s;
    s << num_nodes[i] << ",";
    s << km1s[i] << ",";
    s << imbalances[i] << ",";
    s << rebalance_steps[i] << ",";
    s << target_imbalances[i] << "\n";
    addToFile(s.str(), file);
  }
}

}  // namespace kahypar
