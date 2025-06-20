// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "constraint.hpp"

namespace crab {

std::ostream &operator<<(std::ostream &o, const constraint_t &c) {
    c.write(o);
    return o;
}

// if a certain symbol is present in the constraint
bool constraint_t::contains(const symbol_t &s) const {
    return _lhs.contains(s);
}

// get slack intervals for all symbols in the constraint
std::map<symbol_t, interval_t> constraint_t::get_slack_intervals() const {
    return _lhs.get_slack_intervals();
}

interval_t constraint_t::compute_subtraction(const symbol_t &a, const symbol_t &b) const {
    expression_t lhs = _lhs.get_equivalent_expression();
    int coeff_a = static_cast<int>(lhs.get_coefficient(a));
    int coeff_b = static_cast<int>(lhs.get_coefficient(b));
    interval_t constant = lhs.get_constant_term();
    if (coeff_a == 1 && coeff_b == -1) {
        // a - b + [lb, ub] <= 0 -> a - b = [-oo, -lb]
        return interval_t{bound_t::minus_infinity(), -constant.lb()};
    } else if (coeff_a == -1 && coeff_b == 1) {
        // b - a + [lb, ub] <= 0 -> a - b = [lb, +oo]
        return interval_t{constant.lb(), bound_t::plus_infinity()};
    } else {
        return interval_t::top();
    }
}

constraint_t constraint_t::negate() const {
    // neg(x <= y) -> x > y -> y < x -> y <= x - 1
    return constraint_t(_rhs, _lhs - 1);
}

constraint_t constraint_t::operator|(const constraint_t &other) const {
    if (other.implies(*this)) {
        // other => this
        return *this;
    } else if (implies(other)) {
        // this => other
        return other;
    } else {
        if (is_meta_begin_constraint()) {
            // meta <= begin
            return constraint_t::meta_begin_init_constraint();
        }
        else if (is_begin_end_constraint()) {
            // begin <= end
            return constraint_t::begin_end_init_constraint();
        }
        else {
            return constraint_t::get_true();
        }
    }
}

constraint_t constraint_t::widen(const constraint_t &other) const {
    if (this->implies(other) && other.implies(*this)) {
        // this == other
        return *this;
    } else {
        if (is_meta_begin_constraint()) {
            // meta <= begin
            return constraint_t::meta_begin_init_constraint();
        }
        else if (is_begin_end_constraint()) {
            // begin <= end
            return constraint_t::begin_end_init_constraint();
        }
        else {
            return constraint_t::get_true();
        }
    }
}

// inclusion operator
bool constraint_t::operator<=(const constraint_t &other) const {
    return this->implies(other);
}

bool constraint_t::is_meta_begin_constraint() const {
    return _lhs.contains(symbol_t::meta()) && _lhs.contains(symbol_t::begin());
}

bool constraint_t::is_begin_end_constraint() const {
    return _lhs.contains(symbol_t::begin()) && _lhs.contains(symbol_t::end());
}

bool constraint_t::is_bottom() const {
    // lhs > rhs is bottom
    return _lhs > _rhs;
}

bool constraint_t::is_top() const {
    // lhs <= rhs is top
    // in expression_t, we need to support <= operator as inclusion operator,
    // hence we use check_le instead of operator<=
    return _lhs.check_le(_rhs);
}

constraint_t constraint_t::operator+(constraint_t c2) const {
    return constraint_t(_lhs + c2.get_lhs());
}

bool constraint_t::is_unsat(constraint_t other) const {
    // check if the other constraint is unsat with this constraint
    return operator+(other).is_bottom();
}

bool constraint_t::implies(const constraint_t& other) const {
    // c1 => c2 iff c1 and not c2 is unsat
    constraint_t negated = other.negate();
    return is_unsat(std::move(negated));
}

void constraint_t::write(std::ostream &o) const {
    o << _lhs << " <= " << _rhs;
}

} // namespace crab
