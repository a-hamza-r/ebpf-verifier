// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "refinement.hpp"

namespace crab {

std::ostream& operator<<(std::ostream& o, const refinement_t& r) {
    r.write(o);
    return o;
}

bool refinement_t::is_numeric_refinement() const {
    return _type == refinement_type_t::NUM;
}

bool refinement_t::is_packet_refinement() const {
    return _type == refinement_type_t::PACKET;
}

interval_t refinement_t::get_interval_value() const {
    expression_t e = _value.get_equivalent_expression();
    if (e.only_has_interval()) {
        return e.get_constant_term();
    }
    return interval_t::bottom();
}

std::map<symbol_t, mock_interval_t> refinement_t::get_slack_intervals() const {
    std::map<symbol_t, mock_interval_t> slack_intervals = _value.get_slack_intervals();
    if (has_constraints) {
        slack_intervals.merge(meta_begin_constraint.get_slack_intervals());
        slack_intervals.merge(begin_end_constraint.get_slack_intervals());
    }
    return slack_intervals;
}

void refinement_t::write(std::ostream& o) const {
    std::map<symbol_t, mock_interval_t> slack_intervals = get_slack_intervals();
    bool has_extra_info = has_constraints || !slack_intervals.empty();
    if (has_extra_info) o << "{";
    if (is_numeric_refinement()) {
        o << "num<" << _value << ">";
    }
    else if (is_packet_refinement()) {
        o << "pkt<" << _value << ">";
    }
    else {
        o << "_";
    }
    if (has_extra_info) o << " | ";
    if (has_constraints) {
        o << meta_begin_constraint << " & " << begin_end_constraint;
    }
    size_t i = 0;
    size_t n = slack_intervals.size();
    if (n > 0 && has_constraints) {
        o << " & ";
    }
    for (auto [s, mock_interv] : slack_intervals) {
        o << s << " in " << mock_interv.to_interval();
        if (i < n - 1) {
            o << " & ";
        }
        i++;
    }
    if (has_extra_info) o << "}";
}

refinement_t refinement_t::operator+(interval_t i) const {
    expression_t added_value = _value + i;
    return refinement_t(_type, added_value);
}

refinement_t refinement_t::operator+(int n) const {
    return operator+(interval_t{n});
}

refinement_t refinement_t::operator+(const refinement_t &other) const {
    expression_t new_value = _value + other._value;
    refinement_type_t new_type = (_type == other._type) ? _type
        : (_type == refinement_type_t::PACKET || other._type == refinement_type_t::PACKET)
          ? refinement_type_t::PACKET : refinement_type_t::ANY;
    return refinement_t(new_type, new_value);
}

// use the constraints to compute a value for the subtraction, if possible
interval_t refinement_t::simplify_for_subtraction(const symbol_t& dst, const symbol_t& src) const {
    bound_t max_packet_size = bound_t{number_t{MAX_PACKET_SIZE}};
    bound_t max_meta_size = bound_t{number_t{MAX_META_SIZE}};
    bound_t zero = bound_t{number_t{0}};
    symbol_t begin = symbol_t::begin();
    symbol_t end = symbol_t::end();
    symbol_t meta = symbol_t::meta();

    if (dst == meta && src == begin) {
        // compute meta - begin
        interval_t i = meta_begin_constraint.compute_subtraction(meta, begin);
        return i & interval_t{-max_meta_size, zero};
    }
    else if (dst == begin && src == meta) {
        // compute begin - meta
        interval_t i = meta_begin_constraint.compute_subtraction(begin, meta);
        return i & interval_t{zero, max_meta_size};
    }
    else if (dst == end && src == begin) {
        // compute end - begin
        interval_t i = begin_end_constraint.compute_subtraction(end, begin);
        return i & interval_t{zero, max_packet_size};
    }
    else if (dst == begin && src == end) {
        // compute begin - end
        interval_t i = begin_end_constraint.compute_subtraction(begin, end);
        return i & interval_t{-max_packet_size, zero};
    }
    else if (dst == end && src == meta) {
        // compute end - meta
        constraint_t meta_end_constraint = construct_meta_end_constraint();
        interval_t i = meta_end_constraint.compute_subtraction(end, meta);
        return i & interval_t{zero, max_meta_size + max_packet_size};
    }
    else if (dst == meta && src == end) {
        // compute meta - end
        constraint_t meta_end_constraint = construct_meta_end_constraint();
        interval_t i = meta_end_constraint.compute_subtraction(meta, end);
        return i & interval_t{-max_meta_size - max_packet_size, zero};
    }
    else {
        // We currently only support subtraction between begin, meta, and end
        return interval_t::top();
    }
}

static constraint_t solve_constraints(constraint_t c1, constraint_t c2) {
    // We are adding c2 to a system already containing c1
    if (c1.is_unsat(c2)) {
        // c1 and c2 are unsatisfiable
        // e.g., c1 := begin + 14 <= end and c2 := begin + 14 > end
        return constraint_t::get_false();
    }
    if (c1.implies(c2)) {
        // c1 implies c2, so c2 is a weaker constraint. We keep c1.
        // e.g., c1 := begin + 14 <= end and c2 := begin + 12 <= end
        return c1;
    }
    if (c2.implies(c1)) {
        // c2 implies c1, so c1 is a weaker constraint. We keep c2.
        // e.g., c1 := begin + 14 <= end and c2 := begin + 18 <= end
        return c2;
    }
    // some heuristic to resolve certain cases like the following:
    // c1 := begin + 34 <= end and c2 := begin + a_0 + 18 <= end, a_0 in [0, 60]
    // it will not fall in any category above
    // We can prefer to keep c2 in this case as it might contain more information than c1
    // as it contains a slack variable a_0
    expression_t c2_lhs = c2.get_lhs();
    expression_t c1_lhs = c1.get_lhs();
    auto c2_slacks = c2_lhs.get_slack_intervals();
    auto c1_slacks = c1_lhs.get_slack_intervals();
    if (!c2_slacks.empty() && c1_slacks.empty()) {
        // very specific heuristic but could be beneficial if expressed properly
        // TODO: fix later
        expression_t c2_eq = c2_lhs.get_equivalent_expression();
        interval_t c2_interval = c2_eq.get_constant_term();
        auto symbol_terms = c2_eq.get_symbol_terms();
        if (c2_interval.singleton()) {
            return c1;
        }
        auto lb = interval_t{c2_interval.lb()};
        auto ub = interval_t{c2_interval.ub()};
        expression_t lb_c2 = expression_t{symbol_terms, lb};
        expression_t ub_c2 = expression_t{symbol_terms, ub};
        if (c1.implies(constraint_t(lb_c2)) && constraint_t(ub_c2).implies(c1)) {
            return c2;
        }
        return c1;
    }
    // c1 and c2 are not comparable, i.e., c2 is providing information not affecting c1
    // e.g., c1 := begin + 14 <= end and c2 := begin + 18 > end
    return c1;
}

void refinement_t::add_constraint(constraint_t c) {
    if (c.is_meta_begin_constraint()) {
        meta_begin_constraint = solve_constraints(meta_begin_constraint, std::move(c));
    }
    else if (c.is_begin_end_constraint()) {
        begin_end_constraint = solve_constraints(begin_end_constraint, std::move(c));
    }
}

refinement_t refinement_t::operator-(const refinement_t &other) const {
    expression_t new_value = _value - other._value;
    refinement_type_t new_type = (_type == other._type) ? refinement_type_t::NUM
        : (_type == refinement_type_t::PACKET || other._type == refinement_type_t::PACKET)
          ? refinement_type_t::PACKET : refinement_type_t::ANY;
    return refinement_t(new_type, new_value);
}

bool refinement_t::same_type(const refinement_t &other) const {
    return _type == other._type && _type != refinement_type_t::ANY;
}


// inclusion operator
bool refinement_t::operator<=(const refinement_t &other) const {
    if (!has_constraints) {
        return _value <= other._value;
    }
    return _value <= other._value &&
           meta_begin_constraint <= other.meta_begin_constraint &&
           begin_end_constraint <= other.begin_end_constraint;
}

refinement_t refinement_t::widen(const refinement_t &other) const {
    if (!has_constraints) {
        return refinement_t(_type, _value.widen(other._value));
    }
    return refinement_t(_type, _value.widen(other._value),
                        meta_begin_constraint.widen(other.meta_begin_constraint),
                        begin_end_constraint.widen(other.begin_end_constraint));
}

refinement_t refinement_t::operator|(const refinement_t &other) const {
    if (!has_constraints) {
        return refinement_t(_type, _value | other._value);
    }
    return refinement_t(_type, _value | other._value,
                        meta_begin_constraint | other.meta_begin_constraint,
                        begin_end_constraint | other.begin_end_constraint);
}

// assume that the current refinement is less than or equal to the other refinement
constraint_t refinement_t::assume_le(const refinement_t &other) const {
    return constraint_t(_value, other._value);
}

// assume that the current refinement is greater than the other refinement
constraint_t refinement_t::assume_gt(const refinement_t &other) const {
    return assume_le(other).negate();
}

constraint_t refinement_t::construct_meta_end_constraint() const {
    return meta_begin_constraint + begin_end_constraint;
}

bool refinement_t::check_consistent(const constraint_t& c) const {
    if (c.is_bottom()) {
        return false;
    }
    else if (c.is_top()) {
        return true;
    }
    else if (c.is_begin_end_constraint()) {
        return begin_end_constraint.implies(c);
    }
    else if (c.is_meta_begin_constraint()) {
        return meta_begin_constraint.implies(c);
    }
    else {
        return construct_meta_end_constraint().implies(c);
    }
}

bool refinement_t::safe_access(const expression_t& access_lb, const expression_t& access_ub,
                               bool is_comparison_check) const {
    // access is at least meta
    constraint_t lb = constraint_t(expression_t::meta(), access_lb);
    bool lb_satisfied = check_consistent(lb);

    // access is at most end or MAX_PACKET_SIZE, depending on the check
    constraint_t ub = is_comparison_check ?
        constraint_t(access_ub, expression_t(interval_t{MAX_PACKET_SIZE}))
        : constraint_t(access_ub, expression_t::end());
    bool ub_satisfied = check_consistent(ub);

    return lb_satisfied && ub_satisfied;
}

bool refinement_t::operator==(const refinement_t &other) const {
    // TODO: we need to compare constraints as well
    return _type == other._type && _value == other._value;
}

} // namespace crab
