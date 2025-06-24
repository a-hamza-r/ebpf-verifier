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
std::vector<symbol_t> constraint_t::get_slacks() const {
    return _lhs.get_slacks();
}

// compute the subtraction of two symbols present in the constraint
interval_t constraint_t::compute_subtraction(const symbol_t &a, const symbol_t &b,
                                             std::shared_ptr<slacks_t> slacks) const {
    if (!contains(a) || !contains(b)) {
        // if either symbol is not present, return top interval
        return interval_t::top();
    }
    expression_t lhs = _lhs.get_equivalent_expression(slacks);
    if (lhs.get_num_terms() == 2) {
        int coeff_a = static_cast<int>(lhs.get_coefficient(a));
        int coeff_b = static_cast<int>(lhs.get_coefficient(b));
        interval_t constant = lhs.get_constant_term();
        if (coeff_a == 1 && coeff_b == -1) {
            // a - b + [lb, ub] <= 0
            // -> a - b + k <= 0, k in [lb, ub]
            // -> a - b <= -k, k in [lb, ub]
            // -> a - b <= -ub, a - b <= -lb
            // -> a - b <= -lb is a weaker constraint
            // -> a - b = [-oo, -lb]
            return interval_t{bound_t::minus_infinity(), -constant.lb()};
        } else if (coeff_a == -1 && coeff_b == 1) {
            // b - a + [lb, ub] <= 0
            // -> b - a + k <= 0, k in [lb, ub]
            // -> b - a <= -k, k in [lb, ub]
            // -> b - a <= -ub, b - a <= -lb
            // -> lb <= a - b, ub <= a - b
            // -> a - b = [lb, +oo]
            return interval_t{constant.lb(), bound_t::plus_infinity()};
        }
    }
    return interval_t::top();
}

constraint_t constraint_t::negate() const {
    // neg(x <= y) -> x > y -> y < x -> y <= x - 1
    return constraint_t(_rhs, _lhs - 1);
}

constraint_t constraint_t::join(const constraint_t &other, std::shared_ptr<slacks_t> slacks) const {
    if (other.implies(*this, slacks)) {
        // other => this
        return *this;
    } else if (implies(other, slacks)) {
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
            return constraint_t::true_constraint();
        }
    }
}

constraint_t constraint_t::widen(const constraint_t &other, std::shared_ptr<slacks_t> slacks) const {
    if (this->implies(other, slacks) && other.implies(*this, slacks)) {
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
            return constraint_t::true_constraint();
        }
    }
}

bool constraint_t::is_meta_begin_constraint() const {
    return _lhs.contains(symbol_t::meta()) && _lhs.contains(symbol_t::begin());
}

bool constraint_t::is_begin_end_constraint() const {
    return _lhs.contains(symbol_t::begin()) && _lhs.contains(symbol_t::end());
}

bool constraint_t::check_eq(const constraint_t &other, std::shared_ptr<slacks_t> slacks) const {
    // check if this constraint is equal to the other constraint
    return _lhs.check_eq(other.get_lhs(), slacks);
}

bool constraint_t::is_bottom(std::shared_ptr<slacks_t> slacks) const {
    if (this->check_eq(constraint_t::false_constraint(), slacks)) {
        return true; // false constraint is bottom
    }
    // lhs > rhs is bottom
    return _lhs.check_gt(_rhs, slacks);
}

constraint_t constraint_t::operator+(constraint_t c2) const {
    return constraint_t(_lhs + c2.get_lhs());
}

bool constraint_t::is_unsat(constraint_t other, std::shared_ptr<slacks_t> slacks) const {
    // check if the other constraint is unsat with this constraint
    return operator+(other).is_bottom(slacks);
}

bool constraint_t::implies(const constraint_t& other, std::shared_ptr<slacks_t> slacks) const {
    if (this->check_eq(constraint_t::true_constraint(), slacks)
        && other.check_eq(constraint_t::false_constraint(), slacks)) {
        return false;
    }
    if (other.check_eq(constraint_t::true_constraint(), slacks)) {
        return true;
    }
    // c1 => c2 iff c1 and not c2 is unsat
    return is_unsat(std::move(other.negate()), slacks);
}

void constraint_t::write(std::ostream &o) const {
    o << _lhs << " <= " << _rhs;
}

} // namespace crab
