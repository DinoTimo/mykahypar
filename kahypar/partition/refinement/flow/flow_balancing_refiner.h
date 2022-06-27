
#pragma once

#include "kahypar/definitions.h"
#include "kahypar/partition/context.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/partition/refinement/fm_refiner_base.h"

namespace kahypar {



//This class does not use the type 'Derived' and only directly gives it to the FMRefinerBase class
template <typename RollbackElement = Mandatory, typename Derived = Mandatory>
class FlowBalancingRefiner : public FMRefinerBase<RollbackElement, Derived> {

  private:
    static constexpr bool debug = false;

  public:
    FlowBalancingRefiner(const FlowBalancingRefiner&) = delete;
    FlowBalancingRefiner& operator= (const FlowBalancingRefiner&) = delete;

    FlowBalancingRefiner(FlowBalancingRefiner&&) = delete;
    FlowBalancingRefiner& operator= (FlowBalancingRefiner&&) = delete;

    ~FlowBalancingRefiner() = default;

  protected:
    FlowBalancingRefiner(Hypergraph& hypergraph, const Context& context) : FMRefinerBase<RollbackElement, Derived>(hypergraph, context) { }



};
}