
#include "constraint.hpp"

namespace crab {

std::ostream &operator<<(std::ostream &o, const constraint_t &c) {
    c.write(o);
    return o;
}

bool constraint_t::contains(const symbol_t &s) const {
    return _lhs.contains(s) || _rhs.contains(s);
}

void constraint_t::simplify(std::shared_ptr<slacks_t> slacks) {
    auto neg = _rhs.negate();
    _lhs = _lhs + _rhs.negate();
    _rhs = expression_t(0);
    for (auto &e : _lhs.get_symbol_terms()) {
        if (e.second == 0) {
            _lhs -= e.first;
        }
        else if (e.first.is_slack()) {
            auto val = slacks->find(e.first);
            if (val != slacks->end()) {
                _lhs -= e.first;
                _lhs = _lhs + val->second.to_interval() * interval_t{(int)e.second};
            }
        }
    }
}

bool constraint_t::is_bottom(std::shared_ptr<slacks_t> slacks) {
    simplify(slacks);
    symbol_terms_t terms = _lhs.get_symbol_terms();
    if (terms.size() == 1 && terms.begin()->first == symbol_t::begin()) {
        // we can sometimes simplify if we have a single symbol begin, since it is always offset 0
        _lhs -= terms.begin()->first;
    }
    if (!_lhs.is_constant()) {
        return false;
    }
    auto constant = _lhs.get_constant_term();
    //auto interval_0 = _rhs.get_constant_term();
    return (bound_t{number_t{0}} < constant.lb());
    //return (constant.ub() > interval_0.lb()) && ((constant & interval_0) == interval_t::bottom());
}

constraint_t constraint_t::negate() const {
    // neg(x <= y) -> x > y -> y < x -> y <= x - 1
    return constraint_t(_rhs, _lhs + (-1));
}

void constraint_t::write(std::ostream &o) const {
    o << _lhs << " <= " << _rhs;
}

} // namespace crab
