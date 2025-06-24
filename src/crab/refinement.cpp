// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "refinement.hpp"

namespace crab {

bool refinement_t::is_numeric_refinement() const {
    return _type == refinement_type_t::NUM;
}

bool refinement_t::is_packet_refinement() const {
    return _type == refinement_type_t::PACKET;
}

// compute an interval value for the refinement
interval_t refinement_t::get_interval_value(std::shared_ptr<slacks_t> slacks) const {
    if (is_packet_refinement()) {
        // cannot compute an interval value for packet refinements
        return interval_t::bottom();
    }
    expression_t e = _value.get_equivalent_expression(slacks);
    if (e.is_constant()) {
        return e.get_constant_term();
    }
    return interval_t::bottom();
}

std::vector<symbol_t> refinement_t::get_all_slacks() const {
    std::vector<symbol_t> slacks = _value.get_slacks();
    if (_contains_pkt_constraints) {
        auto meta_begin_slacks = _meta_begin_constraint.get_slacks();
        auto begin_end_slacks = _begin_end_constraint.get_slacks();
        slacks.insert(slacks.end(), meta_begin_slacks.begin(), meta_begin_slacks.end());
        slacks.insert(slacks.end(), begin_end_slacks.begin(), begin_end_slacks.end());
    }
    return slacks;
}

void refinement_t::write(std::ostream& o, std::shared_ptr<slacks_t> slack_intervs) const {
    std::vector<symbol_t> slack_vars = get_all_slacks();
    bool has_extra_info = _contains_pkt_constraints || !slack_vars.empty();
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
    if (_contains_pkt_constraints) {
        o << _meta_begin_constraint << " & " << _begin_end_constraint;
    }
    size_t i = 0;
    size_t n = slack_vars.size();
    if (n > 0 && _contains_pkt_constraints) {
        o << " & ";
    }
    for (const auto s : slack_vars) {
        auto interv = slack_intervs->at(s);
        o << s << " in " << interv;
        if (i < n - 1) {
            o << " & ";
        }
        i++;
    }
    if (has_extra_info) o << "}";
}

refinement_t refinement_t::operator+(interval_t i) const {
    expression_t added_value = _value + i;
    return refinement_t(_type, added_value, _meta_begin_constraint, _begin_end_constraint,
                        _contains_pkt_constraints);
}

refinement_t refinement_t::operator+(int n) const {
    return operator+(interval_t{n});
}

refinement_t refinement_t::add(const refinement_t &other, std::shared_ptr<slacks_t> slacks) const {
    expression_t new_value = _value + other._value;
    refinement_type_t new_type = (_type == other._type) ? _type
        : (_type == refinement_type_t::PACKET || other._type == refinement_type_t::PACKET)
          ? refinement_type_t::PACKET : refinement_type_t::ANY;
    if (_contains_pkt_constraints || other._contains_pkt_constraints) {
        auto new_rf = refinement_t(new_type, new_value,
                            _meta_begin_constraint, _begin_end_constraint, true);
        new_rf.add_constraint(other._meta_begin_constraint, slacks);
        new_rf.add_constraint(other._begin_end_constraint, slacks);
        return new_rf;
    }
    return refinement_t(new_type, new_value);
}

// use the constraints to compute a value for the subtraction, if possible
interval_t refinement_t::simplify_for_subtraction(const symbol_t& dst, const symbol_t& src,
                                                  std::shared_ptr<slacks_t> slacks) const {
    bound_t max_packet_size = bound_t{MAX_PACKET_SIZE};
    bound_t max_meta_size = bound_t{MAX_META_SIZE};
    bound_t zero = bound_t{number_t{0}};
    symbol_t begin = symbol_t::begin();
    symbol_t end = symbol_t::end();
    symbol_t meta = symbol_t::meta();

    if (dst == meta && src == begin) {
        // compute meta - begin
        interval_t i = _meta_begin_constraint.compute_subtraction(meta, begin, slacks);
        return i & interval_t{-max_meta_size, zero};
    }
    else if (dst == begin && src == meta) {
        // compute begin - meta
        interval_t i = _meta_begin_constraint.compute_subtraction(begin, meta, slacks);
        return i & interval_t{zero, max_meta_size};
    }
    else if (dst == end && src == begin) {
        // compute end - begin
        interval_t i = _begin_end_constraint.compute_subtraction(end, begin, slacks);
        return i & interval_t{zero, max_packet_size};
    }
    else if (dst == begin && src == end) {
        // compute begin - end
        interval_t i = _begin_end_constraint.compute_subtraction(begin, end, slacks);
        return i & interval_t{-max_packet_size, zero};
    }
    else if (dst == end && src == meta) {
        // compute end - meta
        constraint_t meta_end_constraint = construct_meta_end_constraint();
        interval_t i = meta_end_constraint.compute_subtraction(end, meta, slacks);
        return i & interval_t{zero, max_meta_size + max_packet_size};
    }
    else if (dst == meta && src == end) {
        // compute meta - end
        constraint_t meta_end_constraint = construct_meta_end_constraint();
        interval_t i = meta_end_constraint.compute_subtraction(meta, end, slacks);
        return i & interval_t{-max_meta_size - max_packet_size, zero};
    }
    else {
        // We currently only support subtraction between begin, meta, and end
        return interval_t::top();
    }
}

static constraint_t solve_constraints(constraint_t c1, constraint_t c2, 
                                      std::shared_ptr<slacks_t> slacks) {
    // We are adding c2 to a system already containing c1
    if (c1.is_unsat(c2, slacks)) {
        // c1 and c2 are unsatisfiable
        // e.g., c1 := begin + 14 <= end and c2 := begin + 14 > end
        return constraint_t::false_constraint();
    }
    if (c1.implies(c2, slacks)) {
        // c1 implies c2, so c2 is a weaker constraint. We keep c1.
        // e.g., c1 := begin + 14 <= end and c2 := begin + 12 <= end
        return c1;
    }
    if (c2.implies(c1, slacks)) {
        // c2 implies c1, so c1 is a weaker constraint. We keep c2.
        // e.g., c1 := begin + 14 <= end and c2 := begin + 18 <= end
        return c2;
    }
    // some heuristic to resolve certain cases like the following:
    // c1 := begin + 34 <= end and c2 := begin + a_0 + 18 <= end, a_0 in [0, 60]
    // it will not fall in any category above
    // We can prefer to keep c2 in this case as it might contain more information than c1
    // as it contains a slack variable a_0, which might be needed in some context.
    expression_t c1_lhs = c1.get_lhs();
    expression_t c2_lhs = c2.get_lhs();
    if (c2_lhs.get_num_terms() > c1_lhs.get_num_terms()) {
        // very specific heuristic but could be beneficial if expressed properly
        expression_t c2_eq = c2_lhs.get_equivalent_expression(slacks);
        interval_t c2_interval = c2_eq.get_constant_term();
        auto symbol_terms = c2_eq.get_symbol_terms();
        if (c2_interval.singleton()) {
            return c1;
        }
        auto lb = interval_t{c2_interval.lb()};
        auto ub = interval_t{c2_interval.ub()};
        expression_t lb_c2 = expression_t{symbol_terms, lb};
        expression_t ub_c2 = expression_t{symbol_terms, ub};
        if (c1.implies(constraint_t(lb_c2), slacks) && constraint_t(ub_c2).implies(c1, slacks)) {
            return c2;
        }
    }
    // c1 and c2 are not comparable, i.e., c2 is providing information not affecting c1
    // e.g., c1 := begin + 14 <= end and c2 := begin + 18 > end,
    // or the previous heuristic did not apply;
    // We can lose some information here, but this keeps the constraints simple
    return c1;
}

void refinement_t::add_constraint(constraint_t c, std::shared_ptr<slacks_t> slacks) {
    if (c.is_meta_begin_constraint()) {
        _meta_begin_constraint = solve_constraints(_meta_begin_constraint, std::move(c), slacks);
    }
    else if (c.is_begin_end_constraint()) {
        _begin_end_constraint = solve_constraints(_begin_end_constraint, std::move(c), slacks);
    }
}

refinement_t refinement_t::subtract(const refinement_t &other,
                                    std::shared_ptr<slacks_t> slacks) const {
    expression_t new_value = _value - other._value;
    refinement_type_t new_type = (_type == other._type) ? refinement_type_t::NUM
        : (_type == refinement_type_t::PACKET) ? refinement_type_t::PACKET
        : refinement_type_t::ANY;
    if (_contains_pkt_constraints || other._contains_pkt_constraints) {
        auto new_rf = refinement_t(new_type, new_value,
                            _meta_begin_constraint, _begin_end_constraint, true);
        new_rf.add_constraint(other._meta_begin_constraint, slacks);
        new_rf.add_constraint(other._begin_end_constraint, slacks);
        return new_rf;
    }
    return refinement_t(new_type, new_value);
}

bool refinement_t::same_type(const refinement_t &other) const {
    return _type == other._type && _type != refinement_type_t::ANY;
}


// inclusion operator
bool refinement_t::inclusion(const refinement_t &other, std::shared_ptr<slacks_t> slacks) const {
    return _value.inclusion(other._value, slacks) &&
              _meta_begin_constraint.implies(other._meta_begin_constraint, slacks) &&
              _begin_end_constraint.implies(other._begin_end_constraint, slacks);
}

refinement_t refinement_t::widen(const refinement_t &other, std::shared_ptr<slacks_t> slacks) const {
    return refinement_t(_type, _value.widen(other._value, slacks),
                        _meta_begin_constraint.widen(other._meta_begin_constraint, slacks),
                        _begin_end_constraint.widen(other._begin_end_constraint, slacks),
                        _contains_pkt_constraints || other._contains_pkt_constraints);
}

refinement_t refinement_t::join(const refinement_t &other, std::shared_ptr<slacks_t> slacks) const {
    return refinement_t(_type, _value.join(other._value, slacks),
                        _meta_begin_constraint.join(other._meta_begin_constraint, slacks),
                        _begin_end_constraint.join(other._begin_end_constraint, slacks),
                        _contains_pkt_constraints || other._contains_pkt_constraints);
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
    return _meta_begin_constraint + _begin_end_constraint;
}

bool refinement_t::check_consistent(const constraint_t& c, std::shared_ptr<slacks_t> slacks) const {
    if (c.is_bottom(slacks)) {
        return false;
    }
    else if (c.is_begin_end_constraint()) {
        return _begin_end_constraint.implies(c, slacks);
    }
    else if (c.is_meta_begin_constraint()) {
        return _meta_begin_constraint.implies(c, slacks);
    }
    else {
        return construct_meta_end_constraint().implies(c, slacks);
    }
}

bool refinement_t::safe_access(const expression_t& access_lb, const expression_t& access_ub,
                               bool is_comparison_check, std::shared_ptr<slacks_t> slacks) const {
    // access is at least meta
    constraint_t lb = constraint_t(expression_t::meta(), access_lb);
    bool lb_satisfied = check_consistent(lb, slacks);

    // access is at most end or MAX_PACKET_SIZE, depending on the check
    constraint_t ub = is_comparison_check ?
        constraint_t(access_ub, expression_t(interval_t{MAX_PACKET_SIZE}))
        : constraint_t(access_ub, expression_t::end());
    bool ub_satisfied = check_consistent(ub, slacks);

    return lb_satisfied && ub_satisfied;
}

bool refinement_t::check_eq(const refinement_t &other, std::shared_ptr<slacks_t> slacks) const {
    // Checking strict equality between constraints might be too strict
    // TODO: revisit this
    bool meta_begin_eq = _meta_begin_constraint.implies(other._meta_begin_constraint, slacks) &&
                         other._meta_begin_constraint.implies(_meta_begin_constraint, slacks);
    bool begin_end_eq = _begin_end_constraint.implies(other._begin_end_constraint, slacks) &&
                        other._begin_end_constraint.implies(_begin_end_constraint, slacks);
    return _type == other._type && _value.check_eq(other._value, slacks) 
           && meta_begin_eq && begin_end_eq;
}

} // namespace crab
