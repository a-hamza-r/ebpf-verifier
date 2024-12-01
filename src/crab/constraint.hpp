// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include "expression.hpp"

namespace crab {

// represents a constraint of the form lhs <= rhs, where rhs is always 0
class constraint_t {
    expression_t _lhs;
    expression_t _rhs = expression_t(0);

  private:
    // true constraint is 0 <= 0, false constraint is 1 <= 0
    constraint_t(bool is_true) : _lhs(is_true ? 0 : 1), _rhs(0) {}

  public:
    constraint_t() : _lhs(1), _rhs(0) {} // false constraint
    constraint_t(expression_t lhs, expression_t rhs)
        : _lhs(lhs - rhs), _rhs(expression_t(0)) {}
    constraint_t(expression_t lhs) : _lhs(lhs), _rhs(expression_t(0)) {}
    static constraint_t get_false() { return constraint_t(false); }
    static constraint_t get_true() { return constraint_t(true); }
    static constraint_t meta_begin_init_constraint() {
        // construct meta <= begin
        return constraint_t(expression_t::meta(), expression_t::begin());
    }
    static constraint_t begin_end_init_constraint() {
        // construct begin <= end
        return constraint_t(expression_t::begin(), expression_t::end());
    }

    bool is_meta_begin_constraint() const;
    bool is_begin_end_constraint() const;
    [[nodiscard]] expression_t get_lhs() const { return _lhs; }
    std::map<symbol_t, mock_interval_t> get_slack_intervals() const;
    bool contains(const symbol_t&) const;
    bool implies(const constraint_t&) const;
    constraint_t operator+(constraint_t) const;
    bool is_unsat(constraint_t) const;
    bool is_bottom() const;
    bool is_top() const;
    interval_t compute_subtraction(const symbol_t&, const symbol_t&) const;
    [[nodiscard]] constraint_t negate() const;
    constraint_t operator|(const constraint_t&) const;
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream &, const constraint_t&);
};


} // namespace crab
