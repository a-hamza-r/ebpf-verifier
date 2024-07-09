
#include "refinement.hpp"

namespace crab {

constexpr int PROPAGATE_INEQUALITIES = 2;

static void propagate_inequalities(std::vector<constraint_t>& constraints_to_check) {
    std::vector<constraint_t> new_constraints;
    for (int i = 0; i < (int)constraints_to_check.size()-1; i++) {
        for (int j = i + 1; j < (int)constraints_to_check.size(); j++) {
            auto c1 = constraints_to_check[i];
            auto c2 = constraints_to_check[j];
            expression_t lhs1 = c1.get_lhs();
            expression_t rhs1 = c1.get_rhs();
            expression_t lhs2 = c2.get_lhs();
            expression_t rhs2 = c2.get_rhs();
            if (rhs1.is_singleton()) {
                symbol_t y = rhs1.get_singleton();
                if (lhs2.contains(y)) {
                    expression_t new_lhs2 = lhs2.substitute(y, lhs1);
                    new_constraints.push_back(constraint_t(new_lhs2, rhs2));
                }
            }
            if (rhs2.is_singleton()) {
                symbol_t y = rhs2.get_singleton();
                if (lhs1.contains(y)) {
                    expression_t new_lhs1 = lhs1.substitute(y, lhs2);
                    new_constraints.push_back(constraint_t(new_lhs1, rhs1));
                }
            }
        }
    }
    constraints_to_check.insert(constraints_to_check.end(),
            new_constraints.begin(), new_constraints.end());
}

std::ostream& operator<<(std::ostream& o, const refinement_t& r) {
    r.write(o);
    return o;
}

bool refinement_t::is_bottom() {
    std::vector<constraint_t> constraints_to_check = _constraints;
    int i = 0;
    while (i <= PROPAGATE_INEQUALITIES) {
        for (constraint_t c : constraints_to_check) {
            if (c.is_bottom()) {
                return true;
            }
        }
        if (i < PROPAGATE_INEQUALITIES) {
            propagate_inequalities(constraints_to_check);
        }
        i++;
    }
    return false;
}

bool refinement_t::has_value(refinement_t&& other) const {
    return (_type == other._type && _value == other._value);
}

void refinement_t::write(std::ostream& o) const {
    symbol_t nu = symbol_t::nu();
    symbol_t i = symbol_t::pkt_symbol();
    o << "{" << nu << " : ";
    if (_type == data_type_t::NUM) {
        o << "num | " << nu << " = ";
    }
    else if (_type == data_type_t::PACKET) {
        o << "pkt<" << i << "> | " << i << " = ";
    }
    else {
        o << "_";
    }
    o << _value;
    if (_constraints.size() > 0) {
        o << " & ";
        for (size_t i = 0; i < _constraints.size(); i++) {
            o << _constraints[i];
            if (i < _constraints.size() - 1) {
                o << " & ";
            }
        }
    }
    o << "}";
}

void refinement_t::simplify() {
    std::set<int> to_remove;
    std::vector<crab::constraint_t> constraints;
    for (int i = 0; i < (int)_constraints.size()-1; i++) {
        for (int j = i+1; j < (int)_constraints.size(); j++) {
            auto c = _constraints[i];
            auto c1 = _constraints[j];
            if (c.get_rhs() == c1.get_rhs()) {
                if (c.get_lhs() < c1.get_lhs()) {
                    to_remove.insert(i);
                }
                else if (c1.get_lhs() < c.get_lhs()) {
                    to_remove.insert(j);
                }
            }
        }
    }
    for (size_t i = 0; i < _constraints.size(); i++) {
        if (to_remove.find(i) == to_remove.end()) {
            constraints.push_back(_constraints[i]);
        }
    }
    _constraints = constraints;
}

refinement_t refinement_t::operator+(interval_t i) const {
    expression_t added_value = _value + i;
    return refinement_t(_type, added_value, _constraints);
}

refinement_t refinement_t::operator+(int n) const {
    return operator+(interval_t{n});
}

refinement_t refinement_t::operator+(const refinement_t &other) const {
    expression_t new_value = _value + other._value;
    std::vector<crab::constraint_t> new_constraints;
    new_constraints.insert(new_constraints.end(), _constraints.begin(), _constraints.end());
    new_constraints.insert(new_constraints.end(),
            other._constraints.begin(), other._constraints.end());
    data_type_t new_type = (_type == other._type) ? _type
        : (_type == data_type_t::PACKET || other._type == data_type_t::PACKET) ? data_type_t::PACKET
        : data_type_t::ANY;
    return refinement_t(new_type, new_value, new_constraints);
}

interval_t refinement_t::simplify_for_subtraction(const symbol_t& dst, const symbol_t& src,
        const std::vector<constraint_t>& constraints) const {
    bound_t max_packet_size = bound_t{number_t{MAX_PACKET_SIZE}};
    bound_t max_meta_size = bound_t{number_t{4098}};
    bound_t zero = bound_t{number_t{0}};
    symbol_t begin = symbol_t::begin();
    symbol_t end = symbol_t::end();
    symbol_t meta = symbol_t::meta();

    if (src != begin && src != meta && src != end && dst != begin && dst != meta && dst != end) {
        // We currently only support subtraction between begin, meta, and end
        return interval_t::top();
    }

    // we can do better in case both end and meta are involved, but for now we should be okay
    if (dst == end && src == meta) {
        return interval_t{zero, max_meta_size + max_packet_size};
    }
    else if (dst == meta && src == end) {
        return interval_t{-(max_meta_size + max_packet_size), zero};
    }
    // we only have cases (begin, meta), and (begin, end) left
    bound_t max_size = (dst == meta || src == meta) ? max_meta_size : max_packet_size;

    interval_t result = interval_t::top();
    for (constraint_t c : constraints) {
        c.simplify();
        expression_t lhs = c.get_lhs();
        if (lhs.contains(dst) && lhs.contains(src)) {
            symbol_terms_t terms = lhs.get_symbol_terms();
            if (terms.size() == 2) {
                if (terms[dst] == 1 && terms[src] == -1) {
                    // dst - src + [lb, ub] <= 0 -> dst - src <= -lb
                    interval_t lhs_interval = lhs.get_constant_term();
                    result = result & interval_t{-max_size, -lhs_interval.lb()};
                }
                else if (terms[dst] == -1 && terms[src] == 1) {
                    // src - dst + [lb, ub] <= 0 -> dst - src >= lb
                    interval_t lhs_interval = lhs.get_constant_term();
                    result = result & interval_t{lhs_interval.lb(), max_size};
                }
            }
        }
    }
    // If no constraints were found, we can still infer some bounds
    if (result == interval_t::top()) {
        if (src == end)
            result = interval_t{-max_packet_size, zero};
        else if (dst == end)
            result = interval_t{zero, max_packet_size};
        else if (src == meta)
            result = interval_t{zero, max_meta_size};
        else if (dst == meta)
            result = interval_t{-max_meta_size, zero};
    }
    return result;
}

refinement_t refinement_t::operator-(const refinement_t &other) const {
    expression_t new_value = _value + other._value.negate();
    std::vector<crab::constraint_t> new_constraints;
    new_constraints.insert(new_constraints.end(), _constraints.begin(), _constraints.end());
    new_constraints.insert(new_constraints.end(),
            other._constraints.begin(), other._constraints.end());
    data_type_t new_type = (_type == other._type) ? data_type_t::NUM
        : (_type == data_type_t::PACKET || other._type == data_type_t::PACKET) ? data_type_t::PACKET
        : data_type_t::ANY;
    return refinement_t(new_type, new_value, new_constraints);
}

bool refinement_t::same_type(const refinement_t &other) const {
    return _type == other._type;
}

refinement_t refinement_t::operator|(const refinement_t &other) const {
    assert(same_type(other));
    auto joined_value = _value | other._value;
    std::vector<crab::constraint_t> new_constraints;
    std::set<int> to_keep, to_keep1;
    for (size_t i = 0; i < _constraints.size(); i++) {
        for (size_t j = 0; j < other._constraints.size(); j++) {
            auto c = _constraints[i];
            auto c1 = other._constraints[j];
            auto lhs = c.get_lhs();
            auto lhs1 = c1.get_lhs();
            if (lhs.is_equal(lhs1)) {
                to_keep.insert(i);
            }
            else if (lhs.is_less_than(lhs1)) {
                to_keep.insert(i);
            }
            else if (lhs1.is_less_than(lhs)) {
                to_keep1.insert(j);
            }
        }
    }
    for (size_t i = 0; i < _constraints.size(); i++) {
        if (to_keep.find(i) != to_keep.end()) {
            new_constraints.push_back(_constraints[i]);
        }
    }
    for (size_t i = 0; i < other._constraints.size(); i++) {
        if (to_keep1.find(i) != to_keep1.end()) {
            new_constraints.push_back(other._constraints[i]);
        }
    }
    // remove duplicates
    for (int i = 0; i < (int)new_constraints.size()-1; i++) {
        for (int j = i+1; j < (int)new_constraints.size(); j++) {
            auto c = new_constraints[i];
            auto c1 = new_constraints[j];
            c.simplify();   c1.simplify();
            if (c.get_lhs().is_equal(c1.get_lhs())) {
                new_constraints.erase(new_constraints.begin() + j);
            }
        }
    }
    return refinement_t(_type, joined_value, new_constraints);
}

constraint_t refinement_t::operator<=(const refinement_t &other) const {
    assert(same_type(other));
    return constraint_t(_value, other._value);
}

constraint_t refinement_t::operator>(const refinement_t &other) const {
    assert(same_type(other));
    return constraint_t(other._value, _value + (-1));
}

void refinement_t::add_constraint(const constraint_t& c) {
    _constraints.push_back(c);
}

bool refinement_t::is_safe_with(refinement_t begin, bool is_comparison_check) const {
    refinement_t check_lb = begin;
    auto lb = constraint_t(expression_t::meta(), _value);
    constraint_t neg_lb = lb.negate();
    check_lb.add_constraint(neg_lb);
    bool lb_satisfied = check_lb.is_bottom();

    refinement_t check_ub = std::move(begin);
    auto ub = is_comparison_check ? constraint_t(_value, expression_t(interval_t{MAX_PACKET_SIZE}))
        : constraint_t(_value, expression_t::end());
    constraint_t neg_ub = ub.negate();
    check_ub.add_constraint(neg_ub);
    bool ub_satisfied = check_ub.is_bottom();

    return lb_satisfied && ub_satisfied;
}

} // namespace crab
