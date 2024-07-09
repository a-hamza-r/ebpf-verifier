// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include "expression.hpp"

namespace crab {

// represents a constraint of the form lhs <= rhs
class constraint_t {
  private:
    expression_t _lhs;
    expression_t _rhs;

  public:
    constraint_t(expression_t lhs, expression_t rhs)
        : _lhs(lhs), _rhs(rhs) {}

    bool is_bottom();

    expression_t get_lhs() const { return _lhs; }
    expression_t get_rhs() const { return _rhs; }

    constraint_t operator+(const constraint_t&) const;    
    constraint_t operator+(int) const;
    constraint_t operator+(interval_t) const;
    constraint_t operator-(const constraint_t&) const;
    constraint_t operator|(const constraint_t&) const;
    constraint_t operator<=(const constraint_t&) const;
    constraint_t operator>(const constraint_t&) const;
    bool is_equality() const;
    bool contains(const symbol_t&) const;
    void normalize();
    [[nodiscard]] constraint_t negate() const;
    void write(std::ostream&) const;
};

std::ostream& operator<<(std::ostream &, const constraint_t&);

} // namespace crab
