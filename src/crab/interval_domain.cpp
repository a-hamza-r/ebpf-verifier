// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/interval_domain.hpp"
#include "boost/endian/conversion.hpp"

namespace crab {

bool interval_domain_t::is_bottom() const {
    if (m_is_bottom) return true;
    return (m_signed.is_bottom() || m_unsigned.is_bottom());
}

bool interval_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_signed.is_top() && m_unsigned.is_top());
}

interval_domain_t interval_domain_t::bottom() {
    interval_domain_t interval;
    interval.set_to_bottom();
    return interval;
}

void interval_domain_t::set_to_bottom() {
    m_is_bottom = true;
    m_signed.set_to_bottom();
    m_unsigned.set_to_bottom();
}

void interval_domain_t::set_to_top() {
    m_is_bottom = false;
    m_signed.set_to_top();
    m_unsigned.set_to_top();
}

void interval_domain_t::set_registers_to_bottom() {
    m_signed.set_registers_to_bottom();
    m_unsigned.set_registers_to_bottom();
}

void interval_domain_t::set_registers_to_top() {
    m_signed.set_registers_to_top();
    m_unsigned.set_registers_to_top();
}

std::optional<interval_cells_t> interval_domain_t::find_in_stack_signed(uint64_t key) const {
    return m_signed.find_in_stack(key);
}

std::optional<interval_cells_t> interval_domain_t::find_in_stack_unsigned(uint64_t key) const {
    return m_unsigned.find_in_stack(key);
}

void interval_domain_t::adjust_bb_for_types(location_t loc) {
    m_signed.adjust_bb_for_types(loc);
    m_unsigned.adjust_bb_for_types(loc);
}

std::vector<uint64_t> interval_domain_t::get_stack_keys() const {
    // likewise for both domains, so just use signed
    return m_signed.get_stack_keys();
}

bool interval_domain_t::all_numeric_in_stack(uint64_t start_loc, int width) const {
    // likewise for both domains, so just use signed
    return m_signed.all_numeric_in_stack(start_loc, width);
}

std::vector<uint64_t> interval_domain_t::find_overlapping_cells_in_stack(uint64_t start_loc,
        int width) const {
    // likewise for both domains, so just use signed
    return m_signed.find_overlapping_cells_in_stack(start_loc, width);
}

void interval_domain_t::remove_overlap_in_stack(const std::vector<uint64_t>& overlap,
        uint64_t start_loc, int width) {
    m_signed.remove_overlap_in_stack(overlap, start_loc, width);
    m_unsigned.remove_overlap_in_stack(overlap, start_loc, width);
}

void interval_domain_t::fill_values_in_stack(const std::vector<uint64_t>& overlap,
        uint64_t start_loc, int width) {
    m_signed.fill_values_in_stack(overlap, start_loc, width);
    m_unsigned.fill_values_in_stack(overlap, start_loc, width);
}

std::optional<mock_interval_t> interval_domain_t::find_interval_value(register_t reg) const {
    // in this case, it does not matter which domain we look into, hence we look into the signed one
    return m_signed.find_interval_value(reg);
}

std::optional<mock_interval_t> interval_domain_t::find_signed_interval_value(register_t reg) const {
    return m_signed.find_interval_value(reg);
}

std::optional<mock_interval_t> interval_domain_t::find_unsigned_interval_value(register_t reg) const {
    return m_unsigned.find_interval_value(reg);
}

std::optional<mock_interval_t> interval_domain_t::find_signed_interval_at_loc(
        const reg_with_loc_t reg) const {
    return m_signed.find_interval_at_loc(reg);
}

std::optional<mock_interval_t> interval_domain_t::find_unsigned_interval_at_loc(
        const reg_with_loc_t reg) const {
    return m_unsigned.find_interval_at_loc(reg);
}

void interval_domain_t::insert_in_registers(register_t reg, location_t loc, interval_t interval) {
    insert_in_registers_unsigned(reg, loc, interval);
    insert_in_registers_signed(reg, loc, interval);
}

void interval_domain_t::insert_in_registers_signed(register_t reg, location_t loc,
        interval_t interval) {
    m_signed.insert_in_registers(reg, loc, interval);
    auto v = m_unsigned.find_interval_at_loc(reg_with_loc_t{reg, loc});
    if (!v) {
        m_unsigned.insert_in_registers(reg, loc, interval_t::top());
    }
}

void interval_domain_t::insert_in_registers_unsigned(register_t reg, location_t loc,
        interval_t interval) {
    m_unsigned.insert_in_registers(reg, loc, interval);
    auto v = m_signed.find_interval_at_loc(reg_with_loc_t{reg, loc});
    if (!v) {
        m_signed.insert_in_registers(reg, loc, interval_t::top());
    }
}

void interval_domain_t::store_in_stack(uint64_t key, mock_interval_t interval, int width) {
    store_in_stack_signed(key, interval, width);
    store_in_stack_unsigned(key, interval, width);
}

void interval_domain_t::store_in_stack_signed(uint64_t key, mock_interval_t interval, int width) {
    m_signed.store_in_stack(key, interval, width);
}

void interval_domain_t::store_in_stack_unsigned(uint64_t key, mock_interval_t interval, int width) {
    m_unsigned.store_in_stack(key, interval, width);
}

bool interval_domain_t::operator<=(const interval_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return true;
}

void interval_domain_t::operator|=(const interval_domain_t& abs) {
    interval_domain_t tmp{abs};
    operator|=(std::move(tmp));
}

void interval_domain_t::operator|=(interval_domain_t&& abs) {
    if (is_bottom()) {
        *this = abs;
        return;
    }
    *this = *this | std::move(abs);
}

interval_domain_t interval_domain_t::operator|(const interval_domain_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return interval_domain_t(m_signed | other.m_signed, m_unsigned | other.m_unsigned);
}

interval_domain_t interval_domain_t::operator|(interval_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return interval_domain_t(
            m_signed | std::move(other.m_signed), m_unsigned | std::move(other.m_unsigned));
}

interval_domain_t interval_domain_t::operator&(const interval_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

interval_domain_t interval_domain_t::widen(const interval_domain_t& abs, bool to_constants) {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

interval_domain_t interval_domain_t::narrow(const interval_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

crab::bound_t interval_domain_t::get_loop_count_upper_bound() {
    /* WARNING: The operation is not implemented yet.*/
    return crab::bound_t{crab::number_t{0}};
}

void interval_domain_t::initialize_loop_counter(const label_t& label) {
    /* WARNING: The operation is not implemented yet.*/
}

string_invariant interval_domain_t::to_set() {
    return string_invariant{};
}

interval_domain_t interval_domain_t::setup_entry() {
    auto&& _signed = signed_interval_domain_t::setup_entry();
    auto&& _unsigned = unsigned_interval_domain_t::setup_entry();
    interval_domain_t interval(std::move(_signed), std::move(_unsigned));
    return interval;
}

void interval_domain_t::operator()(const Un& u, location_t loc) {
    m_signed(u, loc);
    m_unsigned(u, loc);
}

void interval_domain_t::operator()(const LoadMapFd& u, location_t loc) {
    operator-=(register_t{u.dst.v});
}

void interval_domain_t::scratch_caller_saved_registers() {
    for (uint8_t i = R1_ARG; i <= R5_ARG; i++) {
        operator-=(register_t{i});
    }
}

void interval_domain_t::do_call(const Call& u, const stack_cells_t& store_in_stack,
        location_t loc) {
    for (const auto& kv : store_in_stack) {
        auto offset = kv.first;
        auto width = kv.second;
        auto overlapping_cells = find_overlapping_cells_in_stack(offset, width);
        if (overlapping_cells.empty()) {
            m_signed.store_in_stack(offset, interval_t::top(), width);
            m_unsigned.store_in_stack(offset, interval_t::top(), width);
        }
        else {
            fill_values_in_stack(overlapping_cells, offset, width);
        }
        //remove_overlap_in_stack(overlapping_cells, offset, width);
        //m_signed.store_in_stack(offset, interval_t::top(), width);
        //m_unsigned.store_in_stack(offset, interval_t::top(), width);
    }
    auto r0 = register_t{R0_RETURN_VALUE};
    // TODO: Check if packet_reallocate() function call needs handling separately
    if (u.is_map_lookup) {
        operator-=(r0);
    }
    else {
        insert_in_registers(r0, loc, interval_t::top());
    }
    scratch_caller_saved_registers();
}

void interval_domain_t::operator()(const Packet& u, location_t loc) {
    auto r0 = register_t{R0_RETURN_VALUE};
    insert_in_registers(r0, loc, interval_t::top());
    scratch_caller_saved_registers();
}

// Given left and right values, get the left and right intervals
static void get_unsigned_intervals(bool is64, const interval_t& dst_signed,
        const interval_t& dst_unsigned, const interval_t& src_unsigned,
        interval_t& left_interval, interval_t& right_interval) {

    // Get intervals as 32-bit or 64-bit as appropriate.
    left_interval = dst_unsigned;
    right_interval = src_unsigned;
    if (!is64) {
        for (interval_t* interval : {&left_interval, &right_interval}) {
            if (!(*interval <= interval_t::unsigned_int(false))) {
                *interval = interval->truncate_to_uint(false);
            }
        }
    }

    if (left_interval.is_top()) {
        left_interval = dst_signed;
        if (left_interval.is_top()) {
            left_interval = interval_t::unsigned_int(is64);
        }
        else {
            // make left interval as union of two intervals:
                // [0, left_interval.ub()] truncated to uint
                // [left_interval.lb(), -1] truncated to uint, as negative_int <=> unsigned_high
            left_interval = interval_t{number_t{0}, left_interval.ub()}.truncate_to_uint(true) |
                    interval_t{left_interval.lb(), number_t{-1}}.truncate_to_uint(true);
        }
    }

    for (interval_t* interval : {&left_interval, &right_interval}) {
        if (!(*interval <= interval_t::unsigned_int(true))) {
            *interval = interval->truncate_to_uint(true);
        }
    }
}

void interval_domain_t::assume_unsigned_lt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    if (right_interval <= interval_t::nonnegative_int(is64)) {
        // Both left_interval and right_interval fit in [0, INT_MAX],
        // and can be treated as both signed and unsigned values
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, std::move(positive), std::move(positive), true, true, false, false);
    }
    else if (left_interval <= interval_t::unsigned_int(is64) &&
            right_interval <= interval_t::unsigned_int(is64)) {
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), true, true, false, false);
    }
    else if (left_interval <= interval_t::unsigned_int(is64)) {
        // interval can only be represented as uvalue
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::top(),
                false, true, false, false);
        insert_in_registers_signed(left, loc, interval_t::top());
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_signed(std::get<Reg>(right).v, loc, interval_t::top());
        }
    }
    // possibly redundant case, since when left interval is negative, it is converted to
    // unsigned high representation, while right interval likely is not negative
    /*
    else if (left_signed <= interval_t::negative_int(is64) &&
            right_signed <= interval_t::negative_int(is64)) {
        // right_signed and left_signed fit in [INT_MIN, -1], and can be treated as
        // both signed and unsigned, since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::negative_int(is64),
                interval_t::unsigned_high(is64), true, true);
    }
    */
    else {
        // interval can only be represented as uvalue
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), false, true, false, false);
        insert_in_registers_signed(left, loc, interval_t::top());
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_signed(std::get<Reg>(right).v, loc, interval_t::top());
        }
    }
}

void interval_domain_t::assume_unsigned_gt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    if (left_interval <= interval_t::unsigned_int(is64) &&
            right_interval <= interval_t::unsigned_int(is64)) {
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), true, true, false, false);
    }
    // possibly redundant analysis, see unsigned_lt
    /*
    else if (right_signed <= interval_t::negative_int(is64)
            && left_signed <= interval_t::negative_int(is64)) {
        // Both left_signed and right_signed fit in [INT_MIN, -1], and can be treated as both
        // signed and unsigned values since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::negative_int(is64),
                interval_t::unsigned_high(is64), true, true);
    }
    */
    else {
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), false, true, false, false);
        insert_in_registers_signed(left, loc, interval_t::top());
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers(std::get<Reg>(right).v, loc, interval_t::top());
        }
    }
}

void interval_domain_t::assume_unsigned_cst(Condition::Op op, bool is64,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto left_interval = interval_t::bottom();
    auto right_interval = interval_t::bottom();
    get_unsigned_intervals(is64, left_signed, left_unsigned, right_unsigned,
            left_interval, right_interval);

    // Handle uvalue != right.
    if (op == Condition::Op::NE) {
        if (auto rn = right_interval.singleton()) {
            if (rn == left_interval.truncate_to_uint(is64).lb().number()) {
                // "NE lower bound" is equivalent to "GT lower bound".
                op = Condition::Op::GT;
                right_interval = interval_t{left_interval.lb()};
            } else if (rn == left_interval.ub().number()) {
                // "NE upper bound" is equivalent to "LT upper bound".
                op = Condition::Op::LT;
                right_interval = interval_t{left_interval.ub()};
            } else {
                return;
            }
        } else {
            return;
        }
    }
    const bool is_lt = op == Condition::Op::LT || op == Condition::Op::LE;
    bool strict = op == Condition::Op::LT || op == Condition::Op::GT;

    auto llb = left_interval.lb();
    auto lub = left_interval.ub();
    auto rlb = right_interval.lb();
    auto rub = right_interval.ub();
    if (!is_lt && (strict ? (lub <= rlb) : (lub < rlb))) {
        // Left unsigned interval is lower than right unsigned interval.
        set_registers_to_bottom();
        return;
    } else if (is_lt && (strict ? (llb >= rub) : (llb > rub))) {
        // Left unsigned interval is higher than right unsigned interval.
        set_registers_to_bottom();
        return;
    }
    if (is_lt && (strict ? (lub < rlb) : (lub <= rlb))) {
        // Left unsigned interval is lower than right unsigned interval.
        // TODO: verify if setting to top is the correct equivalent of returning linear cst true
        // set_registers_to_top();
        return;
    } else if (!is_lt && (strict ? (llb > rub) : (llb >= rub))) {
        // Left unsigned interval is higher than right unsigned interval.
        // set_registers_to_top();
        return;
    }

    if (is_lt)
        assume_unsigned_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc);
    else
        assume_unsigned_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc);
}

// Given left and right values, get the left and right intervals
static void get_signed_intervals(bool is64, const interval_t& dst_signed,
        const interval_t& dst_unsigned, const interval_t& src_signed,
        interval_t& left_interval, interval_t& right_interval) {

    // Get intervals as 32-bit or 64-bit as appropriate.
    left_interval = dst_signed;
    right_interval = src_signed;
    if (!is64) {
        for (interval_t* interval : {&left_interval, &right_interval}) {
            if (!(*interval <= interval_t::signed_int(false))) {
                *interval = interval->truncate_to_sint(false);
            }
        }
    }

    if (left_interval.is_top()) {
        left_interval = dst_unsigned;
        if (left_interval.is_top()) {
            left_interval = interval_t::signed_int(is64);
        }
        else {
            auto low = (left_interval & interval_t::unsigned_high(is64)).truncate_to_sint(true);
            auto high = (left_interval & interval_t::nonnegative_int(is64)).truncate_to_sint(true);
            left_interval = low | high;
        }
    }

    for (interval_t* interval : {&left_interval, &right_interval}) {
        if (!(*interval <= interval_t::signed_int(true))) {
            *interval = interval->truncate_to_sint(true);
        }
    }
}

void interval_domain_t::update_lt(bool is64, bool strict, interval_t&& left_interval,
        interval_t&& right_interval, const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc,
        interval_t&& restrict_signed, interval_t&& restrict_unsigned,
        bool update_signed, bool update_unsigned, bool mk_equal_unsigned, bool is_signed) {

    auto rlb = right_interval.lb();
    auto rub = right_interval.ub();
    auto llb = left_interval.lb();
    auto lub = left_interval.ub();
    auto llbs = left_signed.lb();
    auto lubs = left_signed.ub();
    auto rlbs = right_signed.lb();
    auto rubs = right_signed.ub();
    auto llbu = left_unsigned.lb();
    auto lubu = left_unsigned.ub();
    auto rlbu = right_unsigned.lb();
    auto rubu = right_unsigned.ub();

    bool holds_reg = std::holds_alternative<Reg>(right);

    if (strict ? llb < rlb : llb <= rlb && lub >= rlb) {
        if (update_signed) {
            auto to_insert_signed = interval_t{llbs, strict ? rlbs - number_t{1} : rlbs};
            to_insert_signed = to_insert_signed & restrict_signed;
            insert_in_registers_signed(left, loc, to_insert_signed);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, to_insert_signed);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{llbu, strict ? rlbu - number_t{1} : rlbu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(left, loc, to_insert_unsigned);
        }
    }
    else if (left_interval <= right_interval && strict ? lub < rub : lub <= rub && holds_reg) {
        auto right_reg = std::get<Reg>(right).v;
        if (update_signed) {
            auto to_insert_signed = interval_t{strict ? lubs + number_t{1} : lubs, rubs};
            to_insert_signed = to_insert_signed & restrict_signed;
            insert_in_registers_signed(right_reg, loc, to_insert_signed);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(right_reg, loc, to_insert_signed);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{strict ? lubu + number_t{1} : lubu, rubu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(right_reg, loc, to_insert_unsigned);
        }
    }
    else if (lub >= rub && strict ? llb < rub : llb <= rub) {
        if (update_signed) {
            auto to_insert_left_signed = interval_t{llbs, strict ? rubs - number_t{1} : rubs};
            to_insert_left_signed = to_insert_left_signed & restrict_signed;
            insert_in_registers_signed(left, loc, to_insert_left_signed);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, to_insert_left_signed);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_left_unsigned = interval_t{llbu, strict ? rubu - number_t{1} : rubu};
            to_insert_left_unsigned = to_insert_left_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(left, loc, to_insert_left_unsigned);
        }

        // this is only one way to resolve this scenario, i.e. set right to singleton value (rub)
        // and set left to the rest of the interval < (or <=) of right
        // a more sound analysis is needed
        if (holds_reg) {
            auto right_reg = std::get<Reg>(right).v;
            if (update_signed) {
                auto to_insert_right_signed = interval_t{rubs} & restrict_signed;
                insert_in_registers_signed(right_reg, loc, to_insert_right_signed);
                if (mk_equal_unsigned && !is_signed) {
                    insert_in_registers_unsigned(right_reg, loc, to_insert_right_signed);
                }
            }
            if (update_unsigned && !is_signed) {
                insert_in_registers_unsigned(right_reg, loc, interval_t{rubu} & restrict_unsigned);
            }
        }
    }
    else {
        // TODO: verify if any legitimate case can fall into here
        set_registers_to_bottom();
    }
    if (is_signed) {
        insert_in_registers_unsigned(left, loc, left_unsigned);
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_unsigned(std::get<Reg>(right).v, loc, right_unsigned);
        }
    }
}

void interval_domain_t::assume_signed_lt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    if (right_interval <= interval_t::negative_int(is64)) {
        // right_interval fits in [INT_MIN, -1], and can be treated as both signed and unsigned
        // since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        // likewise for left_interval, as it is not > right_interval, and truncated to signed int
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::unsigned_high(is64), true, true, false, true);
    }
    else if (left_interval <= interval_t::nonnegative_int(is64) &&
            right_interval <= interval_t::nonnegative_int(is64)) {
        // Both left_interval and right_interval fit in [0, INT_MAX],
        // and can be treated as both signed and unsigned values
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, std::move(positive), std::move(positive), true, false, true, true);
    }
    else  {
        // left_interval and right_interval can be treated as signed values only
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::top(), true, false, false, false);
        insert_in_registers_unsigned(left, loc, left_unsigned);
        if (std::holds_alternative<Reg>(right)) {
            auto right_reg = std::get<Reg>(right).v;
            insert_in_registers_unsigned(right_reg, loc, right_unsigned);
        }
    }
}

void interval_domain_t::update_gt(bool is64, bool strict, interval_t&& left_interval,
        interval_t&& right_interval, const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc, interval_t&& restrict_signed,
        interval_t&& restrict_unsigned, bool update_signed, bool update_unsigned,
        bool mk_equal_unsigned, bool is_signed) {

    auto rlb = right_interval.lb();
    auto rub = right_interval.ub();
    auto llb = left_interval.lb();
    auto lub = left_interval.ub();
    auto llbs = left_signed.lb();
    auto lubs = left_signed.ub();
    auto rlbs = right_signed.lb();
    auto rubs = right_signed.ub();
    auto llbu = left_unsigned.lb();
    auto lubu = left_unsigned.ub();
    auto rlbu = right_unsigned.lb();
    auto rubu = right_unsigned.ub();

    bool holds_reg = std::holds_alternative<Reg>(right);

    if (strict ? lub > rub : lub >= rub && llb <= rub) {
        if (update_signed) {
            auto to_insert_signed = interval_t{strict ? rubs + number_t{1} : rubs, lubs};
            to_insert_signed = to_insert_signed & restrict_signed;
            insert_in_registers_signed(left, loc, to_insert_signed);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, to_insert_signed);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{strict ? rubu + number_t{1} : rubu, lubu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(left, loc, to_insert_unsigned);
        }
    }
    else if (left_interval <= right_interval && strict ? llb > rlb : llb >= rlb && holds_reg) {
        auto right_reg = std::get<Reg>(right).v;
        if (update_signed) {
            auto to_insert_signed = interval_t{rlbs, strict ? llbs - number_t{1} : llbs};
            to_insert_signed = to_insert_signed & restrict_signed;
            insert_in_registers_signed(right_reg, loc, to_insert_signed);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(right_reg, loc, to_insert_signed);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{rlbu, strict ? llbu - number_t{1} : llbu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(right_reg, loc, to_insert_unsigned);
        }
    }
    else if (llb <= rlb && strict ? lub > rlb : lub >= rlb) {
        if (update_signed) {
            auto to_insert_signed_left = interval_t{strict ? rlbs + number_t{1} : rlbs, lubs};
            to_insert_signed_left = to_insert_signed_left & restrict_signed;
            insert_in_registers_signed(left, loc, to_insert_signed_left);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, to_insert_signed_left);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned_left = interval_t{strict ? rlbu + number_t{1} : rlbu, lubu};
            to_insert_unsigned_left = to_insert_unsigned_left & restrict_unsigned;
            insert_in_registers_unsigned(left, loc, to_insert_unsigned_left);
        }

        // this is only one way to resolve this scenario, i.e. set right to singleton value (rlb)
        // and set left to the rest of the interval > (or >=) of right
        // a more sound analysis is needed
        if (holds_reg) {
            auto right_reg = std::get<Reg>(right).v;
            if (update_signed) {
                auto to_insert_signed_right = interval_t{rlbs} & restrict_signed;
                insert_in_registers_signed(right_reg, loc, to_insert_signed_right);
                if (mk_equal_unsigned && !is_signed) {
                    insert_in_registers_unsigned(right_reg, loc, to_insert_signed_right);
                }
            }
            if (update_unsigned && !is_signed) {
                insert_in_registers_unsigned(right_reg, loc, interval_t{rlbu} & restrict_unsigned);
            }
        }
    }
    else {
        // TODO: verify if any legitimate case can fall into here
        set_registers_to_bottom();
    }
    if (is_signed) {
        insert_in_registers_unsigned(left, loc, left_unsigned);
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_unsigned(std::get<Reg>(right).v, loc, right_unsigned);
        }
    }
}

void interval_domain_t::assume_signed_gt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    if (right_interval <= interval_t::nonnegative_int(is64)) {
        // right_interval fits in [0, INT_MAX], and can be treated as both signed and unsigned
        // likewise fits in [0, UINT_MAX], as it is not < right_interval,
        // and truncated to signed int
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, std::move(positive), std::move(positive), true, false, true, true);
    }
    else if (right_interval <= interval_t::negative_int(is64)
            && left_interval <= interval_t::negative_int(is64)) {
        // Both left_interval and right_interval fit in [INT_MIN, -1], and can be treated as both
        // signed and unsigned values since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::unsigned_high(is64),
                true, true, false, true);
    }
    else {
        // left_interval and right_interval can be treated as signed values only
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::top(), true, false, false, true);
        insert_in_registers_unsigned(left, loc, left_unsigned);
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_unsigned(std::get<Reg>(right).v, loc, right_unsigned);
        }
    }
}

void interval_domain_t::assume_signed_cst(Condition::Op op, bool is64,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto left_interval = interval_t::bottom();
    auto right_interval = interval_t::bottom();
    get_signed_intervals(is64, left_signed, left_unsigned, right_signed,
            left_interval, right_interval);

    const bool is_lt = op == Condition::Op::SLT || op == Condition::Op::SLE;
    bool strict = op == Condition::Op::SLT || op == Condition::Op::SGT;

    auto llb = left_interval.lb();
    auto lub = left_interval.ub();
    auto rlb = right_interval.lb();
    auto rub = right_interval.ub();
    if (!is_lt && (strict ? (lub <= rlb) : (lub < rlb))) {
        // Left unsigned interval is lower than right unsigned interval.
        set_registers_to_bottom();
        return;
    } else if (is_lt && (strict ? (llb >= rub) : (llb > rub))) {
        // Left unsigned interval is higher than right unsigned interval.
        set_registers_to_bottom();
        return;
    }
    if (is_lt && (strict ? (lub < rlb) : (lub <= rlb))) {
        // Left unsigned interval is lower than right unsigned interval.
        // TODO: verify if setting to top is the correct equivalent of returning linear cst true
        // set_registers_to_top();
        return;
    } else if (!is_lt && (strict ? (llb > rub) : (llb >= rub))) {
        // Left unsigned interval is higher than right unsigned interval.
        // set_registers_to_top();
        return;
    }

    if (is_lt)
        assume_signed_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc);
    else
        assume_signed_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc);
}

void interval_domain_t::assume_cst(Condition::Op op, bool is64, register_t left,
        Value right, location_t loc) {
    using Op = Condition::Op;

    auto left_signed = find_signed_interval_value(left)->to_interval();
    auto left_unsigned = find_unsigned_interval_value(left)->to_interval();
    auto right_signed = interval_t::bottom();
    auto right_unsigned = interval_t::bottom();
    if (std::holds_alternative<Reg>(right)) {
        auto right_reg = register_t{std::get<Reg>(right).v};
        right_signed = find_signed_interval_value(right_reg)->to_interval();
        right_unsigned = find_unsigned_interval_value(right_reg)->to_interval();
    } else if (std::holds_alternative<Imm>(right)) {
        auto right_imm = std::get<Imm>(right).v;
        right_signed = interval_t{number_t{right_imm}};
        right_unsigned = interval_t(number_t{(uint64_t)right_imm});
    }

    switch (op) {
        case Op::EQ: {
            auto interval_signed = left_signed & right_signed;
            auto interval_unsigned = left_unsigned & right_unsigned;
            insert_in_registers_signed(left, loc, interval_signed);
            insert_in_registers_unsigned(left, loc, interval_unsigned);
            if (std::holds_alternative<Reg>(right)) {
                auto right_reg = std::get<Reg>(right).v;
                insert_in_registers_signed(right_reg, loc, interval_signed);
                insert_in_registers_unsigned(right_reg, loc, interval_unsigned);
            }
            break;
        }
        case Op::SGE:
        case Op::SLE:
        case Op::SGT:
        case Op::SLT: {
            assume_signed_cst(op, is64, left_signed, left_unsigned, right_signed,
                    right_unsigned, left, right, loc);
            break;
        }
        case Op::SET:
        case Op::NSET: {
            // TODO: implement SET and NSET
            break;
        }
        case Op::NE:
        case Op::GE:
        case Op::LE:
        case Op::GT:
        case Op::LT: {
            assume_unsigned_cst(op, is64, left_signed, left_unsigned, right_signed,
                    right_unsigned, left, right, loc);
            break;
        }
    }
}

void interval_domain_t::operator()(const Assume& s, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator-=(register_t reg) {
    m_signed.operator-=(reg);
    m_unsigned.operator-=(reg);
}

bool interval_domain_t::load_from_stack(register_t reg, interval_t load_at, int width,
        location_t loc) {
    uint64_t start_offset = 0;
    if (auto load_at_singleton = load_at.singleton()) {
        start_offset = (uint64_t)(*load_at_singleton);
        bool loaded_signed = m_signed.load_from_stack(reg, start_offset, loc);
        bool loaded_unsigned = m_unsigned.load_from_stack(reg, start_offset, loc);
        if (loaded_signed && loaded_unsigned) return true;
    }
    else {
        auto load_at_lb = load_at.lb();
        auto load_at_ub = load_at.ub();
        if (auto finite_size = load_at.finite_size()) {
            if (auto load_at_lb = load_at.lb().number()) {
                start_offset = (uint64_t)(*load_at_lb);
                width = (int)(*finite_size + number_t{width});
            }
        }
    }
    if (all_numeric_in_stack(start_offset, width)) {
        insert_in_registers(reg, loc, interval_t::top());
        return true;
    }
    return false;
}

void interval_domain_t::do_load(const Mem& b, const register_t& target_register,
        std::optional<ptr_or_mapfd_t> basereg_type, bool load_in_region, location_t loc) {

    if (!basereg_type) {
        operator-=(target_register);
        return;
    }

    // we check if we already loaded a pointer from ctx or stack in region domain,
        // we then do not store a number
    if (load_in_region) {
        operator-=(target_register);
        return;
    }
    int width = b.access.width;
    int offset = b.access.offset;
    auto basereg_ptr_or_mapfd_type = basereg_type.value();

    if (is_ctx_ptr(basereg_type)) {
        insert_in_registers(target_register, loc, interval_t::top());
        return;
    }
    if (is_packet_ptr(basereg_type) || is_shared_ptr(basereg_type)) {
        if (width == 1) {
            interval_t to_insert = interval_t(number_t{0}, number_t{UINT8_MAX});
            insert_in_registers(target_register, loc, to_insert);
        }
        else if (width == 2) {
            interval_t to_insert = interval_t(number_t{0}, number_t{UINT16_MAX});
            insert_in_registers(target_register, loc, to_insert);
        }
        else {
            insert_in_registers(target_register, loc, interval_t::top());
        }
        return;
    }

    if (is_stack_ptr(basereg_type)) {
        auto ptr_with_off = std::get<ptr_with_off_t>(basereg_ptr_or_mapfd_type);
        auto p_offset = ptr_with_off.get_offset();
        auto load_at = p_offset.to_interval() + interval_t(number_t{offset});
        if (load_from_stack(target_register, load_at, width, loc)) return;
    }
    operator-=(target_register);
}

void interval_domain_t::store_in_stack(const Mem& b, uint64_t store_at, int width) {
    m_signed.store_in_stack(b, store_at, width);
    m_unsigned.store_in_stack(b, store_at, width);
}

void interval_domain_t::do_mem_store(const Mem& b, std::optional<ptr_or_mapfd_t> basereg_type) {
    int offset = b.access.offset;
    int width = b.access.width;

    if (!is_stack_ptr(basereg_type)) {
        // we only store for stack pointers
        return;
    }

    auto basereg_ptr_with_off_type = std::get<ptr_with_off_t>(*basereg_type);
    auto offset_singleton = basereg_ptr_with_off_type.get_offset().to_interval().singleton();
    if (!offset_singleton) {
        m_errors.push_back("doing a store with unknown offset");
        return;
    }
    auto store_at = (uint64_t)(*offset_singleton + offset);
    auto overlapping_cells = find_overlapping_cells_in_stack(store_at, width);
    remove_overlap_in_stack(overlapping_cells, store_at, width);
    store_in_stack(b, store_at, width);
}

void interval_domain_t::check_valid_access(const ValidAccess& s, interval_t&& interval,
        int width, bool check_stack_all_numeric) {
    // access can be checked only in the signed domain
    m_signed.check_valid_access(s, std::move(interval), width, check_stack_all_numeric);
}

static void overflow_bounds(interval_t& interval, number_t span, int finite_width,
        bool issigned) {
    if (interval.ub() - interval.lb() >= span) {
        // Interval covers the full space.
        interval = interval_t::top();
        return;
    }
    if (interval.is_bottom()) {
        interval = interval_t::top();
        return;
    }
    number_t lb_value = interval.lb().number().value();
    number_t ub_value = interval.ub().number().value();

    // Compute the interval, taking overflow into account.
    // For a signed result, we need to ensure the signed and unsigned results match
    // so for a 32-bit operation, 0x80000000 should be a positive 64-bit number not
    // a sign extended negative one.
    number_t lb = lb_value.truncate_to_unsigned_finite_width(finite_width);
    number_t ub = ub_value.truncate_to_unsigned_finite_width(finite_width);
    if (issigned) {
        lb = lb.truncate_to_sint64();
        ub = ub.truncate_to_sint64();
    }
    if (lb > ub) {
        // Range wraps in the middle, so we cannot represent as an unsigned interval.
        interval = interval_t::top();
        return;
    }
    interval = crab::interval_t{lb, ub};
}

static void overflow_unsigned(interval_t& interval, int finite_width) {
    auto span{finite_width == 64   ? crab::z_number{std::numeric_limits<uint64_t>::max()}
        : finite_width == 32 ? crab::z_number{std::numeric_limits<uint32_t>::max()}
                                   : throw std::exception()};
    overflow_bounds(interval, span, finite_width, false);
}

static void overflow_signed(interval_t& interval, int finite_width) {
    auto span{finite_width == 64   ? crab::z_number{std::numeric_limits<int64_t>::max()}
        : finite_width == 32 ? crab::z_number{std::numeric_limits<int32_t>::max()}
                                   : throw std::exception()};
    overflow_bounds(interval, span, finite_width, true);
}

static void apply_unsigned(interval_t& dst_signed, interval_t& dst_unsigned,
        const interval_t& src, int finite_width, const interval_t& updated_interval, Bin::Op op) {
    switch (op) {
        case Bin::Op::UDIV: {
            // this hack is required because UDiv call returns [+oo, +oo], at least in case when 
            // updated_interval is top,
            // which causes the error: CRAB ERROR: Bound: undefined operation -oo + +oo
            // SplitDBM (possibly) uses a normalize() to avoid this issue,
            // but we don't have that here
            // TODO: fix this
            if (updated_interval.is_top()) {
                dst_unsigned = interval_t::top();
                dst_signed = interval_t::top();
                return;
            }
            dst_unsigned = updated_interval.UDiv(src);
            break;
        }
        case Bin::Op::UMOD: {
            dst_unsigned = updated_interval.URem(src);
            break;
        }
        case Bin::Op::AND: {
            dst_unsigned = updated_interval.And(src);
            break;
        }
        case Bin::Op::OR: {
            dst_unsigned = updated_interval.Or(src);
            break;
        }
        case Bin::Op::XOR: {
            dst_unsigned = updated_interval.Xor(src);
            break;
        }
        case Bin::Op::LSH: {
            dst_unsigned = updated_interval.Shl(src);
            break;
        }
        case Bin::Op::RSH: {
            dst_unsigned = updated_interval.LShr(src);
            break;
        }
        case Bin::Op::ARSH: {
            dst_unsigned = updated_interval.AShr(src);
            break;
        }
        default: {
            break;
        }
    }
    if (finite_width) {
        dst_signed = dst_unsigned;
        overflow_unsigned(dst_unsigned, finite_width);
        overflow_signed(dst_signed, finite_width);
    }
}

static void apply_signed(interval_t& dst_signed, interval_t& dst_unsigned, const interval_t& src,
        int finite_width, const interval_t& updated_interval, Bin::Op op) {
    switch (op) {
        case Bin::Op::ADD: {
            dst_signed = updated_interval + src;
            break;
        }
        case Bin::Op::SUB: {
            dst_signed = updated_interval - src;
            break;
        }
        case Bin::Op::MUL: {
            dst_signed = updated_interval * src;
            break;
        }
        default: {
            break;
        }
    }
    if (finite_width) {
        dst_unsigned = dst_signed;
        overflow_signed(dst_signed, finite_width);
        overflow_unsigned(dst_unsigned, finite_width);
    }
}

static void shl(register_t reg, int imm, int finite_width,
        interval_t& dst_signed, interval_t& dst_unsigned, location_t loc) {
    // The BPF ISA requires masking the imm.
    imm &= finite_width - 1;

    if (dst_unsigned.finite_size()) {
        number_t lb = dst_unsigned.lb().number().value();
        number_t ub = dst_unsigned.ub().number().value();
        uint64_t lb_n = lb.cast_to_uint64();
        uint64_t ub_n = ub.cast_to_uint64();
        uint64_t uint_max = (finite_width == 64) ? UINT64_MAX : UINT32_MAX;
        if ((lb_n >> (finite_width - imm)) != (ub_n >> (finite_width - imm))) {
            // The bits that will be shifted out to the left are different,
            // which means all combinations of remaining bits are possible.
            lb_n = 0;
            ub_n = (uint_max << imm) & uint_max;
        } else {
            // The bits that will be shifted out to the left are identical
            // for all values in the interval, so we can safely shift left
            // to get a new interval.
            lb_n = (lb_n << imm) & uint_max;
            ub_n = (ub_n << imm) & uint_max;
        }
        dst_unsigned = interval_t{number_t{lb_n}, number_t{ub_n}};
        dst_signed = interval_t::top();
        if ((int64_t)ub_n >= (int64_t)lb_n) {
            dst_signed = dst_unsigned;
        }
        return;
    }
    apply_unsigned(dst_signed, dst_unsigned, interval_t{number_t{imm}}, 64,
            dst_unsigned, Bin::Op::LSH);
}

static void lshr(register_t reg, int imm, int finite_width,
        interval_t& dst_signed, interval_t& dst_unsigned, location_t loc) {
    // The BPF ISA requires masking the imm.
    imm &= finite_width - 1;

    number_t lb_n{0};
    number_t ub_n{UINT64_MAX >> imm};
    if (dst_unsigned.finite_size()) {
        number_t lb = dst_unsigned.lb().number().value();
        number_t ub = dst_unsigned.ub().number().value();
        if (finite_width == 64) {
            lb_n = lb.cast_to_uint64() >> imm;
            ub_n = ub.cast_to_uint64() >> imm;
        } else {
            number_t lb_w = lb.cast_to_signed_finite_width(finite_width);
            number_t ub_w = ub.cast_to_signed_finite_width(finite_width);
            lb_n = lb_w.cast_to_uint32() >> imm;
            ub_n = ub_w.cast_to_uint32() >> imm;

            // The interval must be valid since a signed range crossing 0
            // was earlier converted to a full unsigned range.
            assert(lb_n <= ub_n);
        }
    }
    dst_unsigned = interval_t{lb_n, ub_n};
    if ((int64_t)ub_n >= (int64_t)lb_n) {
        dst_signed = dst_unsigned;
    } else {
        dst_signed = interval_t::top();
    }
    return;
}

static void ashr(register_t reg, interval_t src, int finite_width,
        interval_t& dst_signed, interval_t& dst_unsigned, location_t loc) {
    interval_t left_interval = interval_t::bottom();
    interval_t right_interval = interval_t::bottom();
    get_signed_intervals(finite_width == 64, dst_signed, dst_unsigned,
            src, left_interval, right_interval);
    if (auto sn = right_interval.singleton()) {
        // The BPF ISA requires masking the imm.
        int64_t imm = sn->cast_to_sint64() & (finite_width - 1);

        int64_t lb_n = INT64_MIN >> imm;
        int64_t ub_n = INT64_MAX >> imm;
        if (left_interval.finite_size()) {
            number_t lb = left_interval.lb().number().value();
            number_t ub = left_interval.ub().number().value();
            if (finite_width == 64) {
                lb_n = lb.cast_to_sint64() >> imm;
                ub_n = ub.cast_to_sint64() >> imm;
            } else {
                number_t lb_w = lb.cast_to_signed_finite_width(finite_width) >> (int)imm;
                number_t ub_w = ub.cast_to_signed_finite_width(finite_width) >> (int)imm;
                if (lb_w.cast_to_uint32() <= ub_w.cast_to_uint32()) {
                    lb_n = lb_w.cast_to_uint32();
                    ub_n = ub_w.cast_to_uint32();
                }
            }
        }
        dst_signed = interval_t{number_t{lb_n}, number_t{ub_n}};
        dst_unsigned = interval_t::top();
        if ((uint64_t)ub_n >= (uint64_t)lb_n) {
            dst_unsigned = dst_signed;
        }
        return;
    }
    dst_signed = interval_t::top();
    dst_unsigned = interval_t::top();
}

void interval_domain_t::do_bin(const Bin& bin,
        const std::optional<interval_t>& src_signed_interval_opt,
        const std::optional<interval_t>& src_unsigned_interval_opt,
        const std::optional<ptr_or_mapfd_t>& src_ptr_or_mapfd_opt,
        const std::optional<interval_t>& dst_signed_interval_opt,
        const std::optional<interval_t>& dst_unsigned_interval_opt,
        const std::optional<ptr_or_mapfd_t>& dst_ptr_or_mapfd_opt,
        const interval_t& subtracted, location_t loc) {
    
    using Op = Bin::Op;

    auto dst_register = register_t{bin.dst.v};
    auto finite_width = (bin.is64 ? 64 : 32);

    interval_t dst_signed = subtracted;
    interval_t dst_unsigned = subtracted;
    if (subtracted != interval_t::bottom()) {
        if (!(dst_signed <= interval_t::signed_int(bin.is64))) {
            dst_signed = dst_signed.truncate_to_sint(bin.is64);
        }
        insert_in_registers_signed(dst_register, loc, dst_signed);
        if (!(dst_unsigned <= interval_t::unsigned_int(bin.is64))) {
            dst_unsigned = dst_unsigned.truncate_to_uint(bin.is64);
        }
        insert_in_registers_unsigned(dst_register, loc, dst_unsigned);
        return;
    }
    if (dst_signed_interval_opt) dst_signed = std::move(*dst_signed_interval_opt);
    if (dst_unsigned_interval_opt) dst_unsigned = std::move(*dst_unsigned_interval_opt);

    if (std::holds_alternative<Imm>(bin.v)) {
        int64_t imm;
        if (bin.is64) {
            // Use the full signed value.
            imm = static_cast<int64_t>(std::get<Imm>(bin.v).v);
        } else {
            // Use only the low 32 bits of the value.
            imm = static_cast<int>(std::get<Imm>(bin.v).v);
            if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                apply_unsigned(dst_signed, dst_unsigned, interval_t{number_t{UINT32_MAX}},
                        64, dst_unsigned, Op::AND);
            }
        }
        auto imm_number = number_t{imm};
        auto imm_interval = interval_t{imm_number};
        auto imm_unsigned_interval = interval_t{number_t{imm_number.cast_to_uint64()}};
        auto imm_int_interval = interval_t{number_t{(int)imm}};
        switch (bin.op) {
            case Op::MOV: {
                // ra = imm
                dst_signed = imm_interval;
                overflow_unsigned(imm_interval, (bin.is64 ? 64 : 32));
                dst_unsigned = imm_interval;
                break;
            }
            case Op::ADD: {
                // ra += imm
                if (imm == 0) {
                    return;
                }
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    auto interval = dst_signed.is_top() ? dst_unsigned : dst_signed;
                    apply_signed(dst_signed, dst_unsigned, imm_int_interval, finite_width,
                            interval, Op::ADD);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::SUB: {
                // ra -= imm
                if (imm == 0) {
                    return;
                }
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    auto interval = dst_signed.is_top() ? dst_unsigned : dst_signed;
                    apply_signed(dst_signed, dst_unsigned, imm_int_interval, finite_width,
                            interval, Op::SUB);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::MUL: {
                // ra *= imm
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_signed(dst_signed, dst_unsigned, imm_interval, finite_width, dst_signed,
                            Op::MUL);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::UDIV: {
                // ra /= imm
                if (dst_unsigned_interval_opt && dst_signed_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, imm_interval, finite_width,
                            dst_unsigned, Op::UDIV);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::UMOD: {
                // ra %= imm
                if (dst_unsigned_interval_opt && dst_signed_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, imm_interval, finite_width,
                            dst_unsigned, Op::UMOD);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::AND: {
                // ra &= imm
                // might or might not be needed, but some case occurred where left was negative
                dst_unsigned = dst_unsigned.truncate_to_uint(finite_width);
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, imm_unsigned_interval, finite_width,
                            dst_unsigned, Op::AND);
                    if ((int32_t)imm > 0) {
                        // AND with immediate is only a 32-bit operation so svalue and uvalue
                        // are the same.
                        dst_signed = dst_signed & interval_t{number_t{0}, number_t{imm}};
                        dst_unsigned = dst_unsigned & interval_t{number_t{0}, number_t{imm}};
                    }
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::OR: {
                // ra |= imm
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, imm_unsigned_interval, finite_width,
                            dst_unsigned, Op::OR);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::XOR: {
                // ra ^= imm
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, imm_unsigned_interval, finite_width,
                            dst_unsigned, Op::XOR);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::LSH: {
                // ra <<= imm
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    shl(dst_register, (int32_t)imm, finite_width, dst_signed, dst_unsigned, loc);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::RSH: {
                // ra >>= imm
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    lshr(dst_register, (int32_t)imm, finite_width, dst_signed, dst_unsigned, loc);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::ARSH: {
                // ra >>>= imm
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    ashr(dst_register, interval_t{number_t{(int32_t)imm}}, finite_width,
                            dst_signed, dst_unsigned, loc);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            default: {
                break;
            }
        }
    }
    else {
        if (!src_signed_interval_opt || !src_unsigned_interval_opt) {
            operator-=(dst_register);
            return;
        }
        interval_t src_signed = *src_signed_interval_opt;
        interval_t src_unsigned = *src_unsigned_interval_opt;
        switch (bin.op) {
            case Op::MOV: {
                // ra = rb
                dst_signed = src_signed;
                dst_unsigned = src_unsigned;
                break;
            }
            case Op::ADD: {
                // ra += rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    auto interval = dst_signed.is_top() ? dst_unsigned : dst_signed;
                    apply_signed(dst_signed, dst_unsigned, src_signed, finite_width, interval,
                            Op::ADD);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::SUB: {
                // ra -= rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    auto interval = dst_signed.is_top() ? dst_unsigned : dst_signed;
                    apply_signed(dst_signed, dst_unsigned, src_signed, finite_width, interval,
                            Op::SUB);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::MUL: {
                // ra *= rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_signed(dst_signed, dst_unsigned, src_signed, finite_width, dst_signed,
                            Op::MUL);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::UDIV: {
                // ra /= rb
                if (dst_unsigned_interval_opt && dst_signed_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, src_unsigned, finite_width,
                            dst_unsigned, Op::UDIV);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::UMOD: {
                // ra %= rb
                if (dst_unsigned_interval_opt && dst_signed_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, src_unsigned, finite_width,
                            dst_unsigned, Op::UMOD);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::AND: {
                // ra &= rb
                // might or might not be needed
                dst_unsigned = dst_unsigned.truncate_to_uint(finite_width);
                src_unsigned = src_unsigned.truncate_to_uint(finite_width);
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, src_unsigned, finite_width,
                            dst_unsigned, Op::AND);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::OR: {
                // ra |= rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, src_unsigned, finite_width,
                            dst_unsigned, Op::OR);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::XOR: {
                // ra ^= rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    apply_unsigned(dst_signed, dst_unsigned, src_unsigned, finite_width,
                            dst_unsigned, Op::XOR);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::LSH: {
                // ra <<= rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    if (std::optional<number_t> sn = src_unsigned.singleton()) {
                        uint64_t imm = sn->cast_to_uint64() & (bin.is64 ? 63 : 31);
                        if (imm <= INT32_MAX) {
                            if (!bin.is64) {
                                // Use only the low 32 bits of the value.
                                dst_signed = dst_signed & interval_t{number_t{UINT32_MAX}};
                                dst_unsigned = dst_unsigned & interval_t{number_t{UINT32_MAX}};
                            }
                            shl(dst_register, (int32_t)imm, finite_width, dst_signed, dst_unsigned, loc);
                            break;
                        }
                    }
                    apply_unsigned(dst_signed, dst_unsigned, src_unsigned, 64,
                            dst_unsigned, Op::LSH);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::RSH: {
                // ra >>= rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    if (std::optional<number_t> sn = src_unsigned.singleton()) {
                        uint64_t imm = sn->cast_to_uint64() & (bin.is64 ? 63 : 31);
                        if (imm <= INT32_MAX) {
                            if (!bin.is64) {
                                // Use only the low 32 bits of the value.
                                dst_signed = dst_signed & interval_t{number_t{UINT32_MAX}};
                                dst_unsigned = dst_unsigned & interval_t{number_t{UINT32_MAX}};
                            }
                            lshr(dst_register, (int32_t)imm, finite_width, dst_signed, dst_unsigned, loc);
                            break;
                        }
                    }
                    dst_signed = interval_t::top();
                    dst_unsigned = interval_t::top();
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            case Op::ARSH: {
                // ra >>>= rb
                if (dst_signed_interval_opt && dst_unsigned_interval_opt) {
                    ashr(dst_register, src_signed, finite_width, dst_signed, dst_unsigned, loc);
                }
                else {
                    operator-=(dst_register);
                }
                break;
            }
            default: {
                break;
            }
        }
    }
    if (!dst_signed.is_bottom() && !dst_unsigned.is_bottom()) {
        if (!bin.is64) {
            apply_unsigned(dst_signed, dst_unsigned, interval_t{number_t{UINT32_MAX}},
                    finite_width, dst_unsigned, Op::AND);
        }
        insert_in_registers_signed(dst_register, loc, dst_signed);
        insert_in_registers_unsigned(dst_register, loc, dst_unsigned);
    }
    else {
        operator-=(dst_register);
    }
}

void interval_domain_t::operator()(const Undefined& u, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator()(const Bin& b, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator()(const Call&, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator()(const Exit&, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator()(const Jmp&, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator()(const Mem&, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator()(const Assert&, location_t loc) {
    // nothing to do here
}

void interval_domain_t::operator()(const basic_block_t& bb, int print) {
    // nothing to do here
}
void interval_domain_t::set_require_check(check_require_func_t f) {}

} // namespace crab
