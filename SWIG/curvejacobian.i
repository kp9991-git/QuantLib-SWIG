/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Kyrylo Protsenko

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_curve_jacobian_i
#define quantlib_curve_jacobian_i

%include common.i
%include types.i
%include vectors.i
%include linearalgebra.i
%include termstructures.i

namespace std {
    %template(YieldTermStructureVector) vector<ext::shared_ptr<YieldTermStructure> >;
    %template(ArrayVector) vector<Array>;
}

%{
#include <ql/termstructures/curvejacobiangraph.hpp>

// The C++ graph is registered with statically-typed curves, since both
// bootstrap access and derived-curve inspection need the concrete type.
// The wrappers only ever see ext::shared_ptr<YieldTermStructure>, so
// supported exported curves register a downcasting adder here and add()
// tries them in turn.
class CurveJacobianGraphProxy;

typedef bool (*CurveJacobianAdder)(CurveJacobianGraphProxy&,
                                   const ext::shared_ptr<YieldTermStructure>&);

inline std::vector<CurveJacobianAdder>& curveJacobianAdders() {
    static std::vector<CurveJacobianAdder> adders;
    return adders;
}

/* Adds the registration order and the shared ownership of the curves
   to QuantLib::CurveJacobianGraph, so that results can be returned as
   vectors aligned with the registered curves and the curves are kept
   alive for as long as the graph is.
*/
class CurveJacobianGraphProxy {
  public:
    void add(const ext::shared_ptr<YieldTermStructure>& curve) {
        QL_REQUIRE(curve, "null curve");
        for (auto adder : curveJacobianAdders())
            if (adder(*this, curve))
                return;
        QL_FAIL("the given curve type is not supported by "
                "CurveJacobianGraph.add(); supported curves must provide "
                "bootstrap Jacobians or expose their underlying curve to C++");
    }

    void validateDependencies(bool requireAnalyticMetadata = false) const {
        graph_.validateDependencies(requireAnalyticMetadata);
    }

    template <class Curve>
    void addCurve(const ext::shared_ptr<Curve>& curve) {
        graph_.add(curve);
        if constexpr (QuantLib::detail::supportsCurveJacobianNode<Curve>) {
            for (auto& c : curves_) {
                if (c.get() == curve.get())
                    return;  // already registered; graph_.add() replaced it
            }
            curves_.push_back(curve);
        } else {
            for (auto& c : inspectedCurves_)
                if (c.get() == curve.get())
                    return;
            inspectedCurves_.push_back(curve);
        }
    }

    const std::vector<ext::shared_ptr<YieldTermStructure> >& curves() const {
        return curves_;
    }

    Matrix crossJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                         const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                         std::vector<bool>* analyticRows = nullptr) const {
        QL_REQUIRE(of && withRespectTo, "null curve");
        return graph_.crossJacobian(*of, *withRespectTo, analyticRows);
    }

    Matrix nodeQuoteJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                             const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                             std::vector<bool>* analyticRows = nullptr) const {
        QL_REQUIRE(of && withRespectTo, "null curve");
        return graph_.nodeQuoteJacobian(*of, *withRespectTo, analyticRows);
    }

    std::map<const YieldTermStructure*, Array>
    parRiskMap(const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
               const std::vector<Array>& nodeRisk,
               std::vector<bool>* analyticRows = nullptr) const {
        QL_REQUIRE(curves.size() == nodeRisk.size(),
                   "the number of curves (" << curves.size() <<
                   ") does not match the number of node-risk vectors (" <<
                   nodeRisk.size() << ")");
        std::map<const YieldTermStructure*, Array> input;
        for (Size i = 0; i < curves.size(); ++i) {
            QL_REQUIRE(curves[i], "null curve");
            auto [entry, inserted] = input.emplace(curves[i].get(), nodeRisk[i]);
            if (!inserted) {
                QL_REQUIRE(entry->second.size() == nodeRisk[i].size(),
                           "node-risk vectors for the same curve have different sizes");
                entry->second += nodeRisk[i];
            }
        }
        return graph_.parRisk(input, analyticRows);
    }

  private:
    QuantLib::CurveJacobianGraph graph_;
    std::vector<ext::shared_ptr<YieldTermStructure> > curves_;
    std::vector<ext::shared_ptr<YieldTermStructure> > inspectedCurves_;
};

template <class Curve>
bool addCurveToJacobianGraph(CurveJacobianGraphProxy& graph,
                             const ext::shared_ptr<YieldTermStructure>& curve) {
    if (auto c = ext::dynamic_pointer_cast<Curve>(curve)) {
        graph.addCurve(c);
        return true;
    }
    return false;
}

template <class Curve>
bool registerCurveJacobianAdder() {
    curveJacobianAdders().push_back(&addCurveToJacobianGraph<Curve>);
    return true;
}
%}

/* Registers a piecewise curve type with CurveJacobianGraph::add().
   Invoked by the export_piecewise_curve macros; the generated variable
   has external linkage so that its initializer is guaranteed to run
   when the module is loaded.
*/
%define export_curve_to_jacobian_graph(Name)
%{
bool Name ## _registeredWithCurveJacobianGraph = registerCurveJacobianAdder<Name>();
%}
%enddef

/* Registers a supported concrete derived curve under a stable token name.
   Such curves are inspected by add() but do not appear in curves() or add a
   block to the bootstrap Jacobian system.
*/
%define export_derived_curve_to_jacobian_graph(Name,Curve)
%{
bool Name ## _registeredAsDerivedWithCurveJacobianGraph =
    registerCurveJacobianAdder<Curve>();
%}
%enddef

/* Jacobian accessors shared by all the piecewise curves.  Expanded
   inside the class body of the curve being exported.
*/
%define export_curve_jacobian_methods
%extend {
    /*! Jacobian of the implied quotes with respect to the curve nodes.

        Element (i,j) is the derivative of the i-th alive helper's
        implied quote with respect to the curve value data()[j+1] (the
        value at the reference date, data()[0], is not a free
        variable).  Rows follow the order of the helpers as stored in
        the curve; columns follow the curve nodes, so the matrix is
        square.

        When the curve is part of a multi-curve group or built with
        exogenous curves, this is the partial derivative with the other
        curves' nodes kept fixed.

        Rows that could not be computed analytically fall back to
        numerical differentiation; the flag vector, when given, reports
        which rows were analytical.
    */
    Matrix jacobian() {
        return self->jacobian();
    }
    Matrix jacobian(std::vector<bool>& analyticRows) {
        return self->jacobian(&analyticRows);
    }

    /*! Jacobian of the curve nodes with respect to the quotes.

        Element (j,i) is the derivative of the curve value data()[j+1]
        with respect to the i-th alive helper's quote.  For a
        stand-alone curve, this is the inverse of jacobian().

        When the curve is bootstrapped jointly with other curves as
        part of a multi-curve group (see MultiCurve), the returned
        derivatives account for the feedback through the whole group:
        the columns then span the quotes of all member curves --- this
        curve's quotes first, followed by those of the other members in
        the group's registration order.
    */
    Matrix inverseJacobian() {
        return self->inverseJacobian();
    }
    Matrix inverseJacobian(std::vector<bool>& analyticRows) {
        return self->inverseJacobian(&analyticRows);
    }
}
%enddef

%rename(CurveJacobianGraph) CurveJacobianGraphProxy;
class CurveJacobianGraphProxy {
  public:
    //! cross-curve Jacobians of a set of bootstrapped curves
    /*! Given a set of bootstrapped curves depending on each other (for
        instance, a projection curve built with an exogenous discount
        curve, or co-dependent curves bootstrapped jointly through
        MultiCurve), this class provides the Jacobians of the helper
        quotes and curve nodes across the whole set, so that a
        sensitivity with respect to the nodes of any curve can be
        propagated back through the curve dependencies and expressed as
        a sensitivity with respect to the quoted instruments of any
        curve.

        Rows and columns follow the order of the alive helpers and of
        the curve nodes (excluding the value at the reference date,
        which is not a free variable) as stored in each curve.
    */
    CurveJacobianGraphProxy();

    //! registers a curve; adding the same curve again replaces its entry
    void add(const ext::shared_ptr<YieldTermStructure>& curve);

    //! validates all dependencies reported by registered helpers
    void validateDependencies(bool requireAnalyticMetadata = false) const;

    //! the registered curves, in registration order
    std::vector<ext::shared_ptr<YieldTermStructure> > curves() const;

    %extend {
        void add(const Handle<YieldTermStructure>& curve) {
            QL_REQUIRE(!curve.empty(), "empty curve handle");
            self->add(curve.currentLink());
        }

        /*! Jacobian of the implied quotes of the first curve's helpers
            with respect to the nodes of the second curve, all other
            curves' nodes being kept fixed.  When the two curves
            coincide, this is the curve's own Jacobian.  Rows that could
            not be computed analytically fall back to numerical
            differentiation; the flag vector, when given, reports which
            rows were analytical.
        */
        Matrix crossJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                             const ext::shared_ptr<YieldTermStructure>& withRespectTo) {
            return self->crossJacobian(of, withRespectTo);
        }
        Matrix crossJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                             const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                             std::vector<bool>& analyticRows) {
            return self->crossJacobian(of, withRespectTo, &analyticRows);
        }

        /*! Jacobian of the nodes of the first curve with respect to the
            quotes of the second curve's helpers, obtained by solving
            the differentiated bootstrap conditions of all registered
            curves at once.  Co-dependent (cyclical) sets of curves are
            handled naturally.
        */
        Matrix nodeQuoteJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                                 const ext::shared_ptr<YieldTermStructure>& withRespectTo) {
            return self->nodeQuoteJacobian(of, withRespectTo);
        }
        Matrix nodeQuoteJacobian(const ext::shared_ptr<YieldTermStructure>& of,
                                 const ext::shared_ptr<YieldTermStructure>& withRespectTo,
                                 std::vector<bool>& analyticRows) {
            return self->nodeQuoteJacobian(of, withRespectTo, &analyticRows);
        }

        /*! Propagates sensitivities with respect to curve nodes (for
            instance, obtained by repricing a trade under bumps of the
            node values) into sensitivities with respect to the helper
            quotes of every registered curve.  The given curves must be
            registered and each node-risk vector must have one entry per
            node of the corresponding curve.  Repeated curve entries are
            summed.  The results are returned in the order given by curves().
        */
        std::vector<Array>
        parRisk(const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                const std::vector<Array>& nodeRisk) {
            std::map<const YieldTermStructure*, Array> risk =
                self->parRiskMap(curves, nodeRisk);
            std::vector<Array> result;
            result.reserve(self->curves().size());
            for (const auto& c : self->curves())
                result.push_back(risk.at(c.get()));
            return result;
        }

        std::vector<Array>
        parRisk(const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                const std::vector<Array>& nodeRisk,
                std::vector<bool>& analyticRows) {
            std::map<const YieldTermStructure*, Array> risk =
                self->parRiskMap(curves, nodeRisk, &analyticRows);
            std::vector<Array> result;
            result.reserve(self->curves().size());
            for (const auto& c : self->curves())
                result.push_back(risk.at(c.get()));
            return result;
        }

        //! par risk on one registered curve; see the overload above
        Array parRisk(const std::vector<ext::shared_ptr<YieldTermStructure> >& curves,
                      const std::vector<Array>& nodeRisk,
                      const ext::shared_ptr<YieldTermStructure>& onCurve) {
            QL_REQUIRE(onCurve, "null curve");
            std::map<const YieldTermStructure*, Array> risk =
                self->parRiskMap(curves, nodeRisk);
            auto i = risk.find(onCurve.get());
            QL_REQUIRE(i != risk.end(),
                       "the given curve was not added to the graph");
            return i->second;
        }
    }
};

#endif
