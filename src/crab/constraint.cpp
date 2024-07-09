
#include "constraint.hpp"

namespace crab {

std::ostream &operator<<(std::ostream &o, const constraint_t &c) {
    c.write(o);
    return o;
}

bool constraint_t::contains(const symbol_t &s) const {
    return _lhs.contains(s) || _rhs.contains(s);
}

// move all terms to the left side
void constraint_t::normalize() {
    _lhs = _lhs + _rhs.negate();
    _rhs = expression_t(0);
}

bool constraint_t::is_bottom() {
    normalize();
    symbol_terms_t terms = _lhs.get_symbol_terms();
    if (terms.size() == 1 && terms.begin()->first == symbol_t::begin()) {
        // we can sometimes simplify if we have a single symbol begin, since it is always offset 0
        _lhs -= terms.begin()->first;
    }

    return _lhs.is_greater_than(expression_t(0));
}

constraint_t constraint_t::negate() const {
    // neg(x <= y) -> x > y -> y < x -> y <= x - 1
    return constraint_t(_rhs, _lhs + (-1));
}

void constraint_t::write(std::ostream &o) const {
    o << _lhs << " <= " << _rhs;
}

} // namespace crab
