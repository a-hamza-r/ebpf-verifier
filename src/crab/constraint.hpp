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
    constraint_t() : constraint_t(true) {} // default is true constraint
    constraint_t(expression_t lhs, expression_t rhs)
        : _lhs(lhs - rhs), _rhs(expression_t(0)) {}
    constraint_t(expression_t lhs) : _lhs(lhs), _rhs(expression_t(0)) {}
    static constraint_t false_constraint() { return constraint_t(false); }
    static constraint_t true_constraint() { return constraint_t(true); }
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
    bool contains_single_pkt_symbol() const;
    constraint_t substitute_for_pkt_symbols() const;
    [[nodiscard]] expression_t get_lhs() const { return _lhs; }
    std::vector<symbol_t> get_slacks() const;
    bool contains(const symbol_t&) const;
    bool implies(const constraint_t&, std::shared_ptr<slacks_t>) const;
    constraint_t operator+(constraint_t) const;
    bool check_eq(const constraint_t&, std::shared_ptr<slacks_t>) const;
    bool is_inconsistent(constraint_t, std::shared_ptr<slacks_t>) const;
    bool is_sat(std::shared_ptr<slacks_t>) const;
    bool is_unsat(std::shared_ptr<slacks_t>) const;
    interval_t compute_subtraction(const symbol_t&, const symbol_t&,
                                   std::shared_ptr<slacks_t>) const;
    [[nodiscard]] constraint_t negate() const;
    constraint_t join(const constraint_t&, std::shared_ptr<slacks_t>) const;
    constraint_t widen(const constraint_t&, std::shared_ptr<slacks_t>) const;
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream &, const constraint_t&);
};


} // namespace crab
