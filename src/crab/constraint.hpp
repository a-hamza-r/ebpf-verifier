// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include "expression.hpp"
#include "common.hpp"

namespace crab {

using slacks_t = std::map<symbol_t, mock_interval_t>;

// represents a constraint of the form lhs <= rhs
class constraint_t {
  private:
    expression_t _expression_lhs;
    expression_t _expression_rhs;

  public:
    constraint_t(expression_t lhs, expression_t rhs)
        : _expression_lhs(lhs), _expression_rhs(rhs) {}

    bool is_bottom(std::shared_ptr<slacks_t>);

    expression_t get_lhs() const { return _expression_lhs; }
    expression_t get_rhs() const { return _expression_rhs; }

    constraint_t operator+(const constraint_t&) const;    
    constraint_t operator+(int) const;
    constraint_t operator+(interval_t) const;
    constraint_t operator-(const constraint_t&) const;
    constraint_t operator|(const constraint_t&) const;
    constraint_t operator<=(const constraint_t&) const;
    constraint_t operator>(const constraint_t&) const;
    bool is_equality() const;
    bool contains(const symbol_t&) const;
    void simplify(std::shared_ptr<slacks_t>);
    [[nodiscard]] constraint_t negate() const;
    void write(std::ostream&) const;
};

std::ostream& operator<<(std::ostream &, const constraint_t&);

} // namespace crab
