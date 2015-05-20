/***************************************************************************
 *  Copyright (C) 2015 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_CORE_STATICDISPATCHER_H_
#define SRC_LIB_CORE_STATICDISPATCHER_H_

#include "lib/core/Mandatory.h"
#include "lib/core/Parameters.h"
#include "lib/core/Typelist.h"

namespace core {
template <
  class Executor = Mandatory,
  class TypesLhs = Mandatory,
  class TypesRhs = TypesLhs,
  typename ResultType = void
  >
class StaticDispatcher {
  using Head = typename TypesLhs::Head;
  using Tail = typename TypesLhs::Tail;

 public:
  template <typename BaseTypeLhs,
            typename BaseTypeRhs>
  static ResultType go(BaseTypeLhs& lhs, BaseTypeRhs& rhs, Executor exec,
                       const Parameters& parameters) {
    if (Head* p1 = dynamic_cast<Head*>(&lhs)) {
      return StaticDispatcher<Executor, TypesLhs, TypesRhs,
                              ResultType>::dispatchRhs(*p1, rhs, exec, parameters);
    } else {
      return StaticDispatcher<Executor, Tail, TypesRhs,
                              ResultType>::go(lhs, rhs, exec, parameters);
    }
  }

  template <class SomeLhs,
            typename BaseTypeRhs>
  static ResultType dispatchRhs(SomeLhs& lhs, BaseTypeRhs& rhs, Executor exec,
                                const Parameters& parameters) {
    using Head = typename TypesRhs::Head;
    using Tail = typename TypesRhs::Tail;
    if (Head* p2 = dynamic_cast<Head*>(&rhs)) {
      return exec.fire(lhs, *p2, parameters);
    } else {
      return StaticDispatcher<Executor, TypesLhs, Tail,
                              ResultType>::dispatchRhs(lhs, rhs, exec, parameters);
    }
  }
};

template <
  class Executor,
  class TypesRhs,
  typename ResultType>
class StaticDispatcher<Executor, NullType, TypesRhs, ResultType>{
 public:
  template <typename BaseTypeLhs,
            typename BaseTypeRhs>
  static ResultType go(BaseTypeLhs& lhs, BaseTypeRhs& rhs, Executor exec,
                       const Parameters& parameters) {
    return exec.onError(lhs, rhs, parameters);
  }
};

template <
  class Executor,
  class TypesLhs,
  typename ResultType>
class StaticDispatcher<Executor, TypesLhs, NullType, ResultType>{
 public:
  template <typename BaseTypeLhs,
            typename BaseTypeRhs>
  static ResultType dispatchRhs(BaseTypeLhs& lhs, BaseTypeRhs& rhs, Executor exec,
                                const Parameters& parameters) {
    return exec.onError(lhs, rhs, parameters);
  }
};
}  // namespace core

#endif  // SRC_LIB_CORE_STATICDISPATCHER_H_
