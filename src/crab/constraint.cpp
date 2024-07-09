
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

std::vector<std::pair<symbol_t, interval_t>> constraint_t::get_slack_intervals() const {
    std::vector<std::pair<symbol_t, interval_t>> result = _lhs.get_slack_intervals();
    auto rhs_intervals = _rhs.get_slack_intervals();
    result.insert(result.end(), rhs_intervals.begin(), rhs_intervals.end());
    return result;
}

constraint_t constraint_t::negate() const {
    // neg(x <= y) -> x > y -> y < x -> y <= x - 1
    return constraint_t(_rhs, _lhs + (-1));
}

void constraint_t::write(std::ostream &o) const {
    o << _lhs << " <= " << _rhs;
}

} // namespace crab
