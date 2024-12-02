// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include "constraint.hpp"

namespace crab {

enum class refinement_type_t {
    NUM,    // numeric refinement
    PACKET, // packet refinement
    ANY     // represents any refinement
};

class refinement_t {
    refinement_type_t _type;
    expression_t _value;
    constraint_t meta_begin_constraint; // meta + sth <= begin
    constraint_t begin_end_constraint; // begin + sth <= end
    bool has_constraints = false;

  public:
    refinement_t() : _type(refinement_type_t::ANY) {}
    refinement_t(refinement_type_t type, expression_t value,
                 constraint_t meta_begin_constraint, constraint_t begin_end_constraint) :
        _type(type), _value(value), meta_begin_constraint(meta_begin_constraint),
        begin_end_constraint(begin_end_constraint), has_constraints(true) {}
    refinement_t(refinement_type_t type, expression_t value) :
        _type(type), _value(value), meta_begin_constraint(constraint_t::get_true()),
        begin_end_constraint(constraint_t::get_true()) {}
 
    [[nodiscard]] interval_t get_interval_value() const;
    [[nodiscard]] refinement_type_t get_type() const { return _type; }
    [[nodiscard]] expression_t get_value() const { return _value; }
    [[nodiscard]] expression_t get_equivalent_expression() const;
    refinement_t operator+(int n) const;
    refinement_t operator+(interval_t) const;
    refinement_t operator+(const refinement_t &other) const;
    refinement_t operator-(const refinement_t &other) const;
    refinement_t operator|(const refinement_t &other) const;
    refinement_t widen(const refinement_t &other) const;
    bool operator<=(const refinement_t &other) const;
    constraint_t assume_le(const refinement_t &other) const;
    constraint_t assume_gt(const refinement_t &other) const;
    bool operator==(const refinement_t &other) const;
    void write(std::ostream &o) const;
    bool same_type(const refinement_t &other) const;
    bool safe_access(const expression_t&, const expression_t&, bool) const;
    bool check_consistent(const constraint_t&) const;
    constraint_t construct_meta_end_constraint() const;
    bool is_numeric_refinement() const;
    bool is_packet_refinement() const;
    interval_t simplify_for_subtraction(const symbol_t&, const symbol_t&) const;
    std::map<symbol_t, mock_interval_t> get_slack_intervals() const;
    void add_constraint(constraint_t);
    friend std::ostream &operator<<(std::ostream &, const refinement_t&);

    static refinement_t begin_with_constraints() {
        constraint_t meta_begin_constraint = constraint_t::meta_begin_init_constraint();
        constraint_t begin_end_constraint = constraint_t::begin_end_init_constraint();
        return refinement_t(refinement_type_t::PACKET, expression_t::begin(),
                            meta_begin_constraint, begin_end_constraint);
    }

    static refinement_t begin(std::shared_ptr<slacks_t> slacks = nullptr) {
        return refinement_t(refinement_type_t::PACKET, expression_t::begin(slacks));
    }

    static refinement_t end(std::shared_ptr<slacks_t> slacks) {
        return refinement_t(refinement_type_t::PACKET, expression_t::end(slacks));
    }

    static refinement_t meta(std::shared_ptr<slacks_t> slacks) {
        return refinement_t(refinement_type_t::PACKET, expression_t::meta(slacks));
    }

    static refinement_t numeric_refinement(expression_t value) {
        return refinement_t(refinement_type_t::NUM, value);
    }

    static refinement_t numeric_refinement(interval_t i, std::shared_ptr<slacks_t> slacks) {
        symbol_t s = symbol_t::make();
        (*slacks)[s] = i;
        return numeric_refinement(expression_t(s, slacks));
    }

    static refinement_t numeric_refinement_top(std::shared_ptr<slacks_t> slacks = nullptr) {
        // expression_t::top() represents a top interval
        return numeric_refinement(expression_t::top(slacks));
    }
};


} // namespace crab
