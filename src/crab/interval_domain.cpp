// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/interval_domain.hpp"
#include "boost/endian/conversion.hpp"

namespace crab {

bool interval_domain_t::is_bottom() const {
    return (m_signed.is_bottom() || m_unsigned.is_bottom());
}

bool interval_domain_t::is_top() const {
    return (m_signed.is_top() && m_unsigned.is_top());
}

interval_domain_t interval_domain_t::bottom() {
    interval_domain_t interval;
    interval.set_to_bottom();
    return interval;
}

void interval_domain_t::set_to_bottom() {
    m_signed.set_to_bottom();
    m_unsigned.set_to_bottom();
}

void interval_domain_t::set_to_top() {
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

std::optional<signed_interval_stack_cell_t> interval_domain_t::find_in_stack_signed(uint64_t key) const {
    return m_signed.find_in_stack(key);
}

std::optional<unsigned_interval_stack_cell_t> interval_domain_t::find_in_stack_unsigned(uint64_t key) const {
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

std::optional<refinement_t> interval_domain_t::find_interval_value(register_t reg) const {
    // in this case, it does not matter which domain we look into, hence we look into the signed one
    return m_signed.find_interval_value(reg);
}

std::optional<refinement_t> interval_domain_t::find_signed_interval_value(register_t reg) const {
    return m_signed.find_interval_value(reg);
}

std::optional<refinement_t> interval_domain_t::find_unsigned_interval_value(register_t reg) const {
    return m_unsigned.find_interval_value(reg);
}

std::optional<refinement_t> interval_domain_t::find_signed_interval_at_loc(
        const register_location_t reg) const {
    return m_signed.find_interval_at_loc(reg);
}

std::optional<refinement_t> interval_domain_t::find_unsigned_interval_at_loc(
        const register_location_t reg) const {
    return m_unsigned.find_interval_at_loc(reg);
}

void interval_domain_t::insert_in_registers(register_t reg, location_t loc, refinement_t rf) {
    insert_in_registers_unsigned(reg, loc, rf);
    insert_in_registers_signed(reg, loc, rf);
}

void interval_domain_t::insert_in_registers_signed(register_t reg, location_t loc,
        interval_t interval) {
    refinement_t rf = refinement_t::numeric_refinement(interval, m_slacks);
    m_signed.insert_in_registers(reg, loc, rf);
}

void interval_domain_t::insert_in_registers_signed(register_t reg, location_t loc,
        refinement_t rf) {
    m_signed.insert_in_registers(reg, loc, rf);
}

void interval_domain_t::insert_in_registers_unsigned(register_t reg, location_t loc,
        interval_t interval) {
    refinement_t rf = refinement_t::numeric_refinement(interval, m_slacks);
    m_unsigned.insert_in_registers(reg, loc, rf);
}

void interval_domain_t::insert_in_registers_unsigned(register_t reg, location_t loc,
        refinement_t rf) {
    m_unsigned.insert_in_registers(reg, loc, rf);
}

void interval_domain_t::store_in_stack(uint64_t key, refinement_t rf, int width) {
    store_in_stack_signed(key, rf, width);
    store_in_stack_unsigned(key, rf, width);
}

void interval_domain_t::store_in_stack_signed(uint64_t key, refinement_t rf, int width) {
    m_signed.store_in_stack(key, rf, width);
}

void interval_domain_t::store_in_stack_unsigned(uint64_t key, refinement_t rf, int width) {
    m_unsigned.store_in_stack(key, rf, width);
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
    return interval_domain_t(m_signed | other.m_signed, m_unsigned | other.m_unsigned, m_slacks);
}

interval_domain_t interval_domain_t::operator|(interval_domain_t&& other) const {
    return interval_domain_t(m_signed | std::move(other.m_signed),
                             m_unsigned | std::move(other.m_unsigned), std::move(other.m_slacks));
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

crab::bound_t interval_domain_t::get_loop_count_upper_bound() const {
    /* WARNING: The operation is not implemented yet.*/
    return crab::bound_t{crab::number_t{0}};
}

void interval_domain_t::initialize_loop_counter(const label_t& label) {
    /* WARNING: The operation is not implemented yet.*/
}

string_invariant interval_domain_t::to_set() {
    return string_invariant{};
}

interval_domain_t interval_domain_t::setup_entry(std::shared_ptr<slacks_t> slacks) {
    return interval_domain_t{
        signed_interval_domain_t::setup_entry(slacks),
        unsigned_interval_domain_t::setup_entry(slacks),
        slacks
    };
}

void interval_domain_t::overflow_bounds(const register_t& lhs, number_t span, const int finite_width, location_t loc, bool is_signed) {
    auto rf_opt = is_signed ? m_signed.find_interval_value(lhs) : m_unsigned.find_interval_value(lhs);
    if (!rf_opt) return;
    interval_t interval = rf_opt->get_interval_value();
    // numeric_refinement_top() represents interval_t::top()
    refinement_t top_rf = refinement_t::numeric_refinement_top(m_slacks);
    if (interval.ub() - interval.lb() >= span) {
        // Interval covers the full space.
        // We do not forget the interval, as it will remove the information that it is a number.
        // We only set the interval to top.
        if (is_signed) {
            m_signed.insert_in_registers(lhs, loc, top_rf);
        } else {
            m_unsigned.insert_in_registers(lhs, loc, top_rf);
        }
        return;
    }
    if (interval.is_bottom()) {
        if (is_signed) {
            m_signed.insert_in_registers(lhs, loc, top_rf);
        } else {
            m_unsigned.insert_in_registers(lhs, loc, top_rf);
        }
        return;
    }
    number_t lb_value = interval.lb().number().value();
    number_t ub_value = interval.ub().number().value();

    // Compute the interval, taking overflow into account.
    // For a signed result, we need to ensure the signed and unsigned results match
    // so for a 32-bit operation, 0x80000000 should be a positive 64-bit number not
    // a sign extended negative one.
    number_t lb = lb_value.truncate_to_uint(finite_width);
    number_t ub = ub_value.truncate_to_uint(finite_width);
    if (is_signed) {
        lb = lb.truncate_to<int64_t>();
        ub = ub.truncate_to<int64_t>();
    }
    auto new_interval = interval_t{lb, ub};
    if (lb > ub) {
        // Range wraps in the middle, so we cannot represent as an unsigned interval.
        new_interval = interval_t::top();
    }
    if (is_signed) {
        m_signed.insert_in_registers(lhs, loc, new_interval);
    } else {
        m_unsigned.insert_in_registers(lhs, loc, new_interval);
    }
}

void interval_domain_t::overflow(const register_t& lhs, const int finite_width, location_t loc, bool is_signed) {
    const auto span{finite_width == 64   ? number_t{std::numeric_limits<uint64_t>::max()}
                    : finite_width == 32 ? number_t{std::numeric_limits<uint32_t>::max()}
                                         : throw std::exception()};
    overflow_bounds(lhs, span, finite_width, loc, is_signed);
}

// As defined in the BPF ISA specification, the immediate value of an unsigned modulo and division is treated
// differently depending on the width:
// * for 32 bit, as a 32-bit unsigned integer
// * for 64 bit, as a 32-bit (not 64 bit) signed integer
static number_t read_imm_for_udiv_or_umod(const number_t& imm, const int width) {
    assert(width == 32 || width == 64);
    if (width == 32) {
        return number_t{imm.cast_to<uint32_t>()};
    }
    return number_t{imm.cast_to<int32_t>()};
}

// As defined in the BPF ISA specification, the immediate value of a signed modulo and division is treated
// differently depending on the width:
// * for 32 bit, as a 32-bit signed integer
// * for 64 bit, as a 64-bit signed integer
static number_t read_imm_for_sdiv_or_smod(const number_t& imm, const int width) {
    assert(width == 32 || width == 64);
    if (width == 32) {
        return number_t{imm.cast_to<int32_t>()};
    }
    return number_t{imm.cast_to<int64_t>()};
}


void interval_domain_t::apply(const arith_binaryop_t& op, const register_t& x, const register_t& y, const register_t& z, const int finite_width, location_t loc, bool is_signed) {
    // performing arithmatic operation
    interval_t xi = interval_t::bottom();
    interval_t yi = interval_t::bottom();
    interval_t zi = interval_t::bottom();
    if (is_signed) {
        auto yi_opt = m_signed.find_interval_value(y);
        auto zi_opt = m_signed.find_interval_value(z);
        if (!yi_opt || !zi_opt) {
            std::cerr << "Error: registers not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
        zi = zi_opt->get_interval_value();
    } else {
        auto yi_opt = m_unsigned.find_interval_value(y);
        auto zi_opt = m_unsigned.find_interval_value(z);
        if (!yi_opt || !zi_opt) {
            std::cerr << "Error: registers not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
        zi = zi_opt->get_interval_value();
    }

    switch (op) {
        case arith_binaryop_t::ADD: xi = yi + zi; break;
        case arith_binaryop_t::SUB: xi = yi - zi; break;
        case arith_binaryop_t::MUL: xi = yi * zi; break;
        case arith_binaryop_t::SDIV: xi = yi.SDiv(zi); break;
        case arith_binaryop_t::UDIV: {
            // this hack is required because UDiv call returns [+oo, +oo], at least in case when 
            // updated_interval is top,
            // which causes the error: CRAB ERROR: Bound: undefined operation -oo + +oo
            // SplitDBM (possibly) uses a normalize() to avoid this issue,
            // but we don't have that here
            // TODO: fix this
            if (yi.is_top()) {
                if (is_signed) {
                    m_signed -= x;
                } else {
                    m_unsigned -= x;
                }
                return;
            }
            xi = yi.UDiv(zi);
            break;
        }
        case arith_binaryop_t::SREM: xi = yi.SRem(zi); break;
        case arith_binaryop_t::UREM: xi = yi.URem(zi); break;
        default: break;
    }
    if (is_signed) {
        m_signed.insert_in_registers(x, loc, xi);
    } else {
        m_unsigned.insert_in_registers(x, loc, xi);
    }
}


void interval_domain_t::apply(const arith_binaryop_t& op, const register_t& x, const register_t& y, const number_t& k, const int finite_width, location_t loc, bool is_signed) {
    // performing arithmatic operation
    interval_t xi = interval_t::bottom();
    interval_t yi = interval_t::bottom();
    if (is_signed) {
        auto yi_opt = m_signed.find_interval_value(y);
        if (!yi_opt) {
            std::cerr << "Error: register " << y << " not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
    } else {
        auto yi_opt = m_unsigned.find_interval_value(y);
        if (!yi_opt) {
            std::cerr << "Error: register " << y << " not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
    }

    switch (op) {
        case arith_binaryop_t::ADD: xi = yi + interval_t{k}; break;
        case arith_binaryop_t::SUB: xi = yi - interval_t{k}; break;
        case arith_binaryop_t::MUL: xi = yi * interval_t{k}; break;
        case arith_binaryop_t::SDIV: xi = yi.SDiv(interval_t{read_imm_for_sdiv_or_smod(k, finite_width)}); break;
        case arith_binaryop_t::UDIV: xi = yi.UDiv(interval_t{read_imm_for_udiv_or_umod(k, finite_width)}); break;
        case arith_binaryop_t::SREM: xi = yi.SRem(interval_t{read_imm_for_sdiv_or_smod(k, finite_width)}); break;
        case arith_binaryop_t::UREM: xi = yi.URem(interval_t{read_imm_for_udiv_or_umod(k, finite_width)}); break;
        default: break;
    }
    if (is_signed) {
        m_signed.insert_in_registers(x, loc, xi);
    } else {
        m_unsigned.insert_in_registers(x, loc, xi);
    }
}


void interval_domain_t::apply(const bitwise_binaryop_t& op, const register_t& x, const register_t& y, const number_t& k, const int finite_width, location_t loc, bool is_signed) {
    // performing bitwise operation 
    interval_t xi = interval_t::bottom();
    interval_t yi = interval_t::bottom();
    if (is_signed) {
        auto yi_opt = m_signed.find_interval_value(y);
        if (!yi_opt) {
            std::cerr << "Error: register " << y << " not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
    } else {
        auto yi_opt = m_unsigned.find_interval_value(y);
        if (!yi_opt) {
            std::cerr << "Error: register " << y << " not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
    }
    interval_t zi{number_t{k.cast_to<uint64_t>()}};

    switch (op) {
        case bitwise_binaryop_t::AND: xi = yi.And(zi); break;
        case bitwise_binaryop_t::OR: xi = yi.Or(zi); break;
        case bitwise_binaryop_t::XOR: xi = yi.Xor(zi); break;
        case bitwise_binaryop_t::SHL: xi = yi.Shl(zi); break;
        case bitwise_binaryop_t::LSHR: xi = yi.LShr(zi); break;
        case bitwise_binaryop_t::ASHR: xi = yi.AShr(zi); break;
        default: break;
    }
    if (is_signed) {
        m_signed.insert_in_registers(x, loc, xi);
    } else {
        m_unsigned.insert_in_registers(x, loc, xi);
    }
}


void interval_domain_t::apply(const bitwise_binaryop_t& op, const register_t& x, const register_t& y, const register_t& z, const int finite_width, location_t loc, bool is_signed) {
    // performing bitwise operation
    interval_t xi = interval_t::bottom();
    interval_t yi = interval_t::bottom();
    interval_t zi = interval_t::bottom();
    if (is_signed) {
        auto yi_opt = m_signed.find_interval_value(y);
        auto zi_opt = m_signed.find_interval_value(z);
        if (!yi_opt || !zi_opt) {
            std::cerr << "Error: registers not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
        zi = zi_opt->get_interval_value();
    } else {
        auto yi_opt = m_unsigned.find_interval_value(y);
        auto zi_opt = m_unsigned.find_interval_value(z);
        if (!yi_opt || !zi_opt) {
            std::cerr << "Error: registers not found in the interval environment\n";
            return;
        }
        yi = yi_opt->get_interval_value();
        zi = zi_opt->get_interval_value();
    }

    switch (op) {
        case bitwise_binaryop_t::AND: xi = yi.And(zi); break;
        case bitwise_binaryop_t::OR: xi = yi.Or(zi); break;
        case bitwise_binaryop_t::XOR: xi = yi.Xor(zi); break;
        case bitwise_binaryop_t::SHL: xi = yi.Shl(zi); break;
        case bitwise_binaryop_t::LSHR: xi = yi.LShr(zi); break;
        case bitwise_binaryop_t::ASHR: xi = yi.AShr(zi); break;
        default: break;
    }
    if (is_signed) {
        m_signed.insert_in_registers(x, loc, xi);
    } else {
        m_unsigned.insert_in_registers(x, loc, xi);
    }
}


void interval_domain_t::apply_signed(const binaryop_t& op, const register_t& result, const register_t& lhs, const number_t& k_rhs, const int finite_width, location_t loc) {
    apply(op, result, lhs, k_rhs, finite_width, loc, true);
    if (finite_width) {
        auto signed_result = m_signed.find_interval_value(result);
        if (signed_result) {
            m_unsigned.insert_in_registers(result, loc, *signed_result);
        }
        overflow(result, finite_width, loc, true);
        overflow(result, finite_width, loc, false);
    }
}

void interval_domain_t::apply_signed(const binaryop_t& op, const register_t& result, const register_t& lhs, const register_t& rhs, const int finite_width, location_t loc) {
    apply(op, result, lhs, rhs, finite_width, loc, true);
    if (finite_width) {
        auto signed_result = m_signed.find_interval_value(result);
        if (signed_result) {
            m_unsigned.insert_in_registers(result, loc, *signed_result);
        }
        overflow(result, finite_width, loc, true);
        overflow(result, finite_width, loc, false);
    }
}

void interval_domain_t::apply_unsigned(const binaryop_t& op, const register_t& result, const register_t& lhs, const number_t& k_rhs, const int finite_width, location_t loc) {
    apply(op, result, lhs, k_rhs, finite_width, loc, false);
    if (finite_width) {
        auto unsigned_result = m_unsigned.find_interval_value(result);
        if (unsigned_result) {
            m_signed.insert_in_registers(result, loc, *unsigned_result);
        }
        overflow(result, finite_width, loc, true);
        overflow(result, finite_width, loc, false);
    }
}

void interval_domain_t::apply_unsigned(const binaryop_t& op, const register_t& result, const register_t& lhs, const register_t& rhs, const int finite_width, location_t loc) {
    apply(op, result, lhs, rhs, finite_width, loc, false);
    if (finite_width) {
        auto unsigned_result = m_unsigned.find_interval_value(result);
        if (unsigned_result) {
            m_signed.insert_in_registers(result, loc, *unsigned_result);
        }
        overflow(result, finite_width, loc, true);
        overflow(result, finite_width, loc, false);
    }
}

void interval_domain_t::add(const register_t& lhs, const register_t& op2, location_t loc) {
    apply_signed(arith_binaryop_t::ADD, lhs, lhs, op2, 0, loc);
}

void interval_domain_t::add(const register_t& lhs, const number_t& op2, location_t loc) {
    apply_signed(arith_binaryop_t::ADD, lhs, lhs, op2, 0, loc);
}

void interval_domain_t::sub(const register_t& lhs, const register_t& op2, location_t loc) {
    apply_signed(arith_binaryop_t::SUB, lhs, lhs, op2, 0, loc);
}

void interval_domain_t::sub(const register_t& lhs, const number_t& op2, location_t loc) {
    apply_signed(arith_binaryop_t::SUB, lhs, lhs, op2, 0, loc);
}

// Add/subtract with overflow are both signed and unsigned. We can use either one of the two to compute the
// result before adjusting for overflow, though if one is top we want to use the other to retain precision.
void interval_domain_t::add_overflow(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    interval_t lhs_signed = m_signed.find_interval_value(lhs)->get_interval_value();
    if (!lhs_signed.is_top()) {
        apply_signed(arith_binaryop_t::ADD, lhs, lhs, op2, finite_width, loc);
    } else {
        apply_unsigned(arith_binaryop_t::ADD, lhs, lhs, op2, finite_width, loc);
    }
}

void interval_domain_t::add_overflow(const register_t& lhs, const number_t& op2, const int finite_width, location_t loc) {
    interval_t lhs_signed = m_signed.find_interval_value(lhs)->get_interval_value();
    if (!lhs_signed.is_top()) {
        apply_signed(arith_binaryop_t::ADD, lhs, lhs, op2, finite_width, loc);
    } else {
        apply_unsigned(arith_binaryop_t::ADD, lhs, lhs, op2, finite_width, loc);
    }
}

void interval_domain_t::sub_overflow(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    interval_t lhs_signed = m_signed.find_interval_value(lhs)->get_interval_value();
    if (!lhs_signed.is_top()) {
        apply_signed(arith_binaryop_t::SUB, lhs, lhs, op2, finite_width, loc);
    } else {
        apply_unsigned(arith_binaryop_t::SUB, lhs, lhs, op2, finite_width, loc);
    }
}

void interval_domain_t::sub_overflow(const register_t& lhs, const number_t& op2, const int finite_width, location_t loc) {
    interval_t lhs_signed = m_signed.find_interval_value(lhs)->get_interval_value();
    if (!lhs_signed.is_top()) {
        apply_signed(arith_binaryop_t::SUB, lhs, lhs, op2, finite_width, loc);
    } else {
        apply_unsigned(arith_binaryop_t::SUB, lhs, lhs, op2, finite_width, loc);
    }
}

void interval_domain_t::neg(const register_t& lhs, const int finite_width, location_t loc) {
    apply_signed(arith_binaryop_t::MUL, lhs, lhs, number_t{-1}, finite_width, loc);
}

void interval_domain_t::mul(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_signed(arith_binaryop_t::MUL, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::mul(const register_t& lhs, const number_t& op2, const int finite_width, location_t loc) {
    apply_signed(arith_binaryop_t::MUL, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::udiv(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_unsigned(arith_binaryop_t::UDIV, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::udiv(const register_t& lhs, const number_t& op2, const int finite_width, location_t loc) {
    apply_unsigned(arith_binaryop_t::UDIV, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::sdiv(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_signed(arith_binaryop_t::SDIV, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::sdiv(const register_t& lhs, const number_t& op2, const int finite_width, location_t loc) {
    apply_signed(arith_binaryop_t::SDIV, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::srem(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_signed(arith_binaryop_t::SREM, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::srem(const register_t& lhs, const number_t& op2, const int finite_width, location_t loc) {
    apply_signed(arith_binaryop_t::SREM, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::urem(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_unsigned(arith_binaryop_t::UREM, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::urem(const register_t& lhs, const number_t& op2, const int finite_width, location_t loc) {
    apply_unsigned(arith_binaryop_t::UREM, lhs, lhs, op2, finite_width, loc);
}


void interval_domain_t::bitwise_and(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_unsigned(bitwise_binaryop_t::AND, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::bitwise_and(const register_t& lhs, const number_t& op2, location_t loc) {
    // Use finite width 64 to make the svalue be set as well as the uvalue.
    apply_unsigned(bitwise_binaryop_t::AND, lhs, lhs, op2, 64, loc);
}

void interval_domain_t::bitwise_or(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_unsigned(bitwise_binaryop_t::OR, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::bitwise_or(const register_t& lhs, const number_t& op2, location_t loc) {
    apply_unsigned(bitwise_binaryop_t::OR, lhs, lhs, op2, 64, loc);
}

void interval_domain_t::bitwise_xor(const register_t& lhs, const register_t& op2, const int finite_width, location_t loc) {
    apply_unsigned(bitwise_binaryop_t::XOR, lhs, lhs, op2, finite_width, loc);
}

void interval_domain_t::bitwise_xor(const register_t& lhs, const number_t& op2, location_t loc) {
    apply_unsigned(bitwise_binaryop_t::XOR, lhs, lhs, op2, 64, loc);
}

void interval_domain_t::shl_overflow(const register_t& lhs, const register_t& op2, location_t loc) {
    apply_unsigned(bitwise_binaryop_t::SHL, lhs, lhs, op2, 64, loc);
}

void interval_domain_t::shl_overflow(const register_t& lhs, const number_t& op2, location_t loc) {
    apply_unsigned(bitwise_binaryop_t::SHL, lhs, lhs, op2, 64, loc);
}

void interval_domain_t::operator()(const Un& u, location_t loc) {
    if (u.op == Un::Op::NEG) {
        auto dst = register_t{u.dst.v};
        neg(dst, u.is64 ? 64 : 32, loc);
        return;
    }
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
    refinement_t top_rf = refinement_t::numeric_refinement_top(m_slacks);
    for (const auto& kv : store_in_stack) {
        auto offset = kv.first;
        auto width = kv.second;
        auto overlapping_cells = find_overlapping_cells_in_stack(offset, width);
        if (overlapping_cells.empty()) {
            m_signed.store_in_stack(offset, top_rf, width);
            m_unsigned.store_in_stack(offset, top_rf, width);
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
        insert_in_registers(r0, loc, top_rf);
    }
    scratch_caller_saved_registers();
}

void interval_domain_t::operator()(const Packet& u, location_t loc) {
    auto r0 = register_t{R0_RETURN_VALUE};
    insert_in_registers(r0, loc, refinement_t::numeric_refinement_top(m_slacks));
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
            if (!(*interval <= interval_t::unsigned_int(32))) {
                *interval = interval->truncate_to_uint(32);
            }
        }
    }

    if (left_interval.is_top()) {
        left_interval = dst_signed;
        if (left_interval.is_top()) {
            left_interval = interval_t::unsigned_int(64);
        }
        else {
            // make left interval as union of two intervals:
                // [0, left_interval.ub()] truncated to uint
                // [left_interval.lb(), -1] truncated to uint, as negative_int <=> unsigned_high
            left_interval = interval_t{number_t{0}, left_interval.ub()}.truncate_to_uint(64) |
                    interval_t{left_interval.lb(), number_t{-1}}.truncate_to_uint(64);
        }
    }

    for (interval_t* interval : {&left_interval, &right_interval}) {
        if (!(*interval <= interval_t::unsigned_int(64))) {
            *interval = interval->truncate_to_uint(64);
        }
    }
}

void interval_domain_t::assume_unsigned_lt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    refinement_t top_rf = refinement_t::numeric_refinement_top(m_slacks);
    if (right_interval <= interval_t::nonnegative(64)) {
        // Both left_interval and right_interval fit in [0, INT_MAX],
        // and can be treated as both signed and unsigned values
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, std::move(positive), std::move(positive), true, true, false, false);
    }
    else if (left_interval <= interval_t::unsigned_int(64) &&
            right_interval <= interval_t::unsigned_int(64)) {
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), true, true, false, false);
    }
    else if (left_interval <= interval_t::unsigned_int(64)) {
        // interval can only be represented as uvalue
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::top(),
                false, true, false, false);
        insert_in_registers_signed(left, loc, top_rf);
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_signed(std::get<Reg>(right).v, loc, top_rf);
        }
    }
    // possibly redundant case, since when left interval is negative, it is converted to
    // unsigned high representation, while right interval likely is not negative
    /*
    else if (left_signed <= interval_t::negative(64) &&
            right_signed <= interval_t::negative(64)) {
        // right_signed and left_signed fit in [INT_MIN, -1], and can be treated as
        // both signed and unsigned, since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::negative(64),
                interval_t::unsigned_high(64), true, true);
    }
    */
    else {
        // interval can only be represented as uvalue
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), false, true, false, false);
        insert_in_registers_signed(left, loc, top_rf);
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_signed(std::get<Reg>(right).v, loc, top_rf);
        }
    }
}

void interval_domain_t::assume_unsigned_gt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    refinement_t top_rf = refinement_t::numeric_refinement_top(m_slacks);
    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    if (left_interval <= interval_t::unsigned_int(64) &&
            right_interval <= interval_t::unsigned_int(64)) {
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), true, true, false, false);
    }
    // possibly redundant analysis, see unsigned_lt
    /*
    else if (right_signed <= interval_t::negative(64)
            && left_signed <= interval_t::negative(64)) {
        // Both left_signed and right_signed fit in [INT_MIN, -1], and can be treated as both
        // signed and unsigned values since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::negative(64),
                interval_t::unsigned_high(64), true, true);
    }
    */
    else {
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), std::move(positive), false, true, false, false);
        insert_in_registers_signed(left, loc, top_rf);
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers(std::get<Reg>(right).v, loc, top_rf);
        }
    }
}

void interval_domain_t::assume_unsigned_cst(Condition::Op op, bool is64,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    if (left_unsigned.is_top() && right_unsigned.is_top()) {
        // this is a drastic heuristic
        return;
    }

    auto left_interval = interval_t::bottom();
    auto right_interval = interval_t::bottom();
    get_unsigned_intervals(is64, left_signed, left_unsigned, right_unsigned,
            left_interval, right_interval);

    // Handle uvalue != right.
    if (op == Condition::Op::NE) {
        if (auto rn = right_interval.singleton()) {
            if (rn == left_interval.truncate_to_uint(64).lb().number()) {
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
            if (!(*interval <= interval_t::signed_int(32))) {
                *interval = interval->truncate_to_sint(32);
            }
        }
    }

    if (left_interval.is_top()) {
        left_interval = dst_unsigned;
        if (left_interval.is_top()) {
            left_interval = interval_t::signed_int(64);
        }
        else {
            auto low = (left_interval & interval_t::unsigned_high(64)).truncate_to_sint(64);
            auto high = (left_interval & interval_t::nonnegative(64)).truncate_to_sint(64);
            left_interval = low | high;
        }
    }

    for (interval_t* interval : {&left_interval, &right_interval}) {
        if (!(*interval <= interval_t::signed_int(64))) {
            *interval = interval->truncate_to_sint(64);
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
            refinement_t rf = refinement_t::numeric_refinement(to_insert_signed, m_slacks);
            insert_in_registers_signed(left, loc, rf);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, rf);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{llbu, strict ? rlbu - number_t{1} : rlbu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(left, loc,
                                refinement_t::numeric_refinement(to_insert_unsigned, m_slacks));
        }
    }
    else if (left_interval <= right_interval && strict ? lub < rub : lub <= rub && holds_reg) {
        auto right_reg = std::get<Reg>(right).v;
        if (update_signed) {
            auto to_insert_signed = interval_t{strict ? lubs + number_t{1} : lubs, rubs};
            to_insert_signed = to_insert_signed & restrict_signed;
            refinement_t rf = refinement_t::numeric_refinement(to_insert_signed, m_slacks);
            insert_in_registers_signed(right_reg, loc, rf);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(right_reg, loc, rf);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{strict ? lubu + number_t{1} : lubu, rubu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(right_reg, loc,
                            refinement_t::numeric_refinement(to_insert_unsigned, m_slacks));
        }
    }
    else if (lub >= rub && strict ? llb < rub : llb <= rub) {
        if (update_signed) {
            auto to_insert_left_signed = interval_t{llbs, strict ? rubs - number_t{1} : rubs};
            to_insert_left_signed = to_insert_left_signed & restrict_signed;
            refinement_t rf = refinement_t::numeric_refinement(to_insert_left_signed, m_slacks);
            insert_in_registers_signed(left, loc, rf);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, rf);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_left_unsigned = interval_t{llbu, strict ? rubu - number_t{1} : rubu};
            to_insert_left_unsigned = to_insert_left_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(left, loc,
                            refinement_t::numeric_refinement(to_insert_left_unsigned, m_slacks));
        }

        // this is only one way to resolve this scenario, i.e. set right to singleton value (rub)
        // and set left to the rest of the interval < (or <=) of right
        // a more sound analysis is needed
        if (holds_reg) {
            auto right_reg = std::get<Reg>(right).v;
            if (update_signed) {
                auto to_insert_right_signed = interval_t{rubs} & restrict_signed;
                refinement_t rf = refinement_t::numeric_refinement(to_insert_right_signed, m_slacks);
                insert_in_registers_signed(right_reg, loc, rf);
                if (mk_equal_unsigned && !is_signed) {
                    insert_in_registers_unsigned(right_reg, loc, rf);
                }
            }
            if (update_unsigned && !is_signed) {
                interval_t interval = interval_t{rubu} & restrict_unsigned;
                insert_in_registers_unsigned(right_reg, loc,
                                             refinement_t::numeric_refinement(interval, m_slacks));
            }
        }
    }
    else {
        // TODO: verify if any legitimate case can fall into here
        set_registers_to_bottom();
    }
    if (is_signed) {
        insert_in_registers_unsigned(left, loc,
                                     refinement_t::numeric_refinement(left_unsigned, m_slacks));
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_unsigned(std::get<Reg>(right).v, loc,
                                     refinement_t::numeric_refinement(right_unsigned, m_slacks));
        }
    }
}

void interval_domain_t::assume_signed_lt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    if (right_interval <= interval_t::negative(64)) {
        // right_interval fits in [INT_MIN, -1], and can be treated as both signed and unsigned
        // since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        // likewise for left_interval, as it is not > right_interval, and truncated to signed int
        update_lt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::unsigned_high(64), true, true, false, true);
    }
    else if (left_interval <= interval_t::nonnegative(64) &&
            right_interval <= interval_t::nonnegative(64)) {
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
        insert_in_registers_unsigned(left, loc,
                                     refinement_t::numeric_refinement(left_unsigned, m_slacks));
        if (std::holds_alternative<Reg>(right)) {
            auto right_reg = std::get<Reg>(right).v;
            insert_in_registers_unsigned(right_reg, loc,
                                     refinement_t::numeric_refinement(right_unsigned, m_slacks));
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
            refinement_t rf = refinement_t::numeric_refinement(to_insert_signed, m_slacks);
            insert_in_registers_signed(left, loc, rf);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, rf);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{strict ? rubu + number_t{1} : rubu, lubu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(left, loc,
                                    refinement_t::numeric_refinement(to_insert_unsigned, m_slacks));
        }
    }
    else if (left_interval <= right_interval && strict ? llb > rlb : llb >= rlb && holds_reg) {
        auto right_reg = std::get<Reg>(right).v;
        if (update_signed) {
            auto to_insert_signed = interval_t{rlbs, strict ? llbs - number_t{1} : llbs};
            to_insert_signed = to_insert_signed & restrict_signed;
            refinement_t rf = refinement_t::numeric_refinement(to_insert_signed, m_slacks);
            insert_in_registers_signed(right_reg, loc, rf);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(right_reg, loc, rf);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned = interval_t{rlbu, strict ? llbu - number_t{1} : llbu};
            to_insert_unsigned = to_insert_unsigned & restrict_unsigned;
            insert_in_registers_unsigned(right_reg, loc,
                                    refinement_t::numeric_refinement(to_insert_unsigned, m_slacks));
        }
    }
    else if (llb <= rlb && strict ? lub > rlb : lub >= rlb) {
        if (update_signed) {
            auto to_insert_signed_left = interval_t{strict ? rlbs + number_t{1} : rlbs, lubs};
            to_insert_signed_left = to_insert_signed_left & restrict_signed;
            refinement_t rf = refinement_t::numeric_refinement(to_insert_signed_left, m_slacks);
            insert_in_registers_signed(left, loc, rf);
            if (mk_equal_unsigned && !is_signed) {
                insert_in_registers_unsigned(left, loc, rf);
            }
        }

        if (update_unsigned && !is_signed) {
            auto to_insert_unsigned_left = interval_t{strict ? rlbu + number_t{1} : rlbu, lubu};
            to_insert_unsigned_left = to_insert_unsigned_left & restrict_unsigned;
            insert_in_registers_unsigned(left, loc,
                            refinement_t::numeric_refinement(to_insert_unsigned_left, m_slacks));
        }

        // this is only one way to resolve this scenario, i.e. set right to singleton value (rlb)
        // and set left to the rest of the interval > (or >=) of right
        // a more sound analysis is needed
        if (holds_reg) {
            auto right_reg = std::get<Reg>(right).v;
            if (update_signed) {
                auto to_insert_signed_right = interval_t{rlbs} & restrict_signed;
                refinement_t rf = refinement_t::numeric_refinement(to_insert_signed_right, m_slacks);
                insert_in_registers_signed(right_reg, loc, rf);
                if (mk_equal_unsigned && !is_signed) {
                    insert_in_registers_unsigned(right_reg, loc, rf);
                }
            }
            if (update_unsigned && !is_signed) {
                interval_t interval = interval_t{rlbu} & restrict_unsigned;
                insert_in_registers_unsigned(right_reg, loc,
                                             refinement_t::numeric_refinement(interval, m_slacks));
            }
        }
    }
    else {
        // TODO: verify if any legitimate case can fall into here
        set_registers_to_bottom();
    }
    if (is_signed) {
        insert_in_registers_unsigned(left, loc,
                                     refinement_t::numeric_refinement(left_unsigned, m_slacks));
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_unsigned(std::get<Reg>(right).v, loc,
                                     refinement_t::numeric_refinement(right_unsigned, m_slacks));
        }
    }
}

void interval_domain_t::assume_signed_gt(bool is64, bool strict,
        interval_t&& left_interval, interval_t&& right_interval,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    auto positive = interval_t{number_t{0}, bound_t::plus_infinity()};
    if (right_interval <= interval_t::nonnegative(64)) {
        // right_interval fits in [0, INT_MAX], and can be treated as both signed and unsigned
        // likewise fits in [0, UINT_MAX], as it is not < right_interval,
        // and truncated to signed int
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, std::move(positive), std::move(positive), true, false, true, true);
    }
    else if (right_interval <= interval_t::negative(64)
            && left_interval <= interval_t::negative(64)) {
        // Both left_interval and right_interval fit in [INT_MIN, -1], and can be treated as both
        // signed and unsigned values since [INT_MIN, -1] <=> [INT_MAX+1, UINT_MAX]
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::unsigned_high(64),
                true, true, false, true);
    }
    else {
        // left_interval and right_interval can be treated as signed values only
        update_gt(is64, strict, std::move(left_interval), std::move(right_interval),
                left_signed, left_unsigned, right_signed, right_unsigned,
                left, right, loc, interval_t::top(), interval_t::top(), true, false, false, true);
        insert_in_registers_unsigned(left, loc,
                                     refinement_t::numeric_refinement(left_unsigned, m_slacks));
        if (std::holds_alternative<Reg>(right)) {
            insert_in_registers_unsigned(std::get<Reg>(right).v, loc,
                                     refinement_t::numeric_refinement(right_unsigned, m_slacks));
        }
    }
}

void interval_domain_t::assume_signed_cst(Condition::Op op, bool is64,
        const interval_t& left_signed, const interval_t& left_unsigned,
        const interval_t& right_signed, const interval_t& right_unsigned,
        register_t left, Value right, location_t loc) {

    if (left_unsigned.is_top() && right_unsigned.is_top()) {
        // this is a drastic heuristic
        return;
    }

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

    auto left_signed = find_signed_interval_value(left)->get_interval_value();
    auto left_unsigned = find_unsigned_interval_value(left)->get_interval_value();
    auto right_signed = interval_t::bottom();
    auto right_unsigned = interval_t::bottom();
    if (std::holds_alternative<Reg>(right)) {
        auto right_reg = register_t{std::get<Reg>(right).v};
        right_signed = find_signed_interval_value(right_reg)->get_interval_value();
        right_unsigned = find_unsigned_interval_value(right_reg)->get_interval_value();
    } else if (std::holds_alternative<Imm>(right)) {
        auto right_imm = std::get<Imm>(right).v;
        right_signed = interval_t{number_t{right_imm}};
        right_unsigned = interval_t(number_t{(uint64_t)right_imm});
    }

    switch (op) {
        case Op::EQ: {
            auto interval_signed = left_signed & right_signed;
            auto interval_unsigned = left_unsigned & right_unsigned;
            refinement_t rf_signed = refinement_t::numeric_refinement(interval_signed, m_slacks);
            refinement_t rf_unsigned = refinement_t::numeric_refinement(interval_unsigned, m_slacks);
            insert_in_registers_signed(left, loc, rf_signed);
            insert_in_registers_unsigned(left, loc, rf_unsigned);
            if (std::holds_alternative<Reg>(right)) {
                auto right_reg = std::get<Reg>(right).v;
                insert_in_registers_signed(right_reg, loc, rf_signed);
                insert_in_registers_unsigned(right_reg, loc, rf_unsigned);
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
    uint start_offset = 0;
    if (auto load_at_singleton = load_at.singleton()) {
        start_offset = load_at_singleton->cast_to<uint>();
        bool loaded_signed = m_signed.load_from_stack(reg, start_offset, loc);
        bool loaded_unsigned = m_unsigned.load_from_stack(reg, start_offset, loc);
        if (loaded_signed && loaded_unsigned) return true;
    }
    else {
        auto load_at_lb = load_at.lb();
        auto load_at_ub = load_at.ub();
        if (auto finite_size = load_at.finite_size()) {
            if (auto load_at_lb = load_at.lb().number()) {
                start_offset = load_at_lb->cast_to<uint>();
                width = (*finite_size + number_t{width}).cast_to<int>();
            }
        }
    }
    auto overlapping_cells = find_overlapping_cells_in_stack(start_offset, width);
    if (overlapping_cells.size() == 1) {
        // only allow loading from a single cell
        if (all_numeric_in_stack(start_offset, width)) {
            insert_in_registers(reg, loc, refinement_t::numeric_refinement_top(m_slacks));
            return true;
        }
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

    refinement_t top_rf = refinement_t::numeric_refinement_top(m_slacks);
    if (is_ctx_ptr(basereg_type)) {
        insert_in_registers(target_register, loc, top_rf);
        return;
    }
    if (is_packet_ptr(basereg_type) || is_shared_ptr(basereg_type)) {
        if (width == 1) {
            interval_t to_insert = interval_t(number_t{0}, number_t{UINT8_MAX});
            expression_t e(to_insert, m_slacks);
            insert_in_registers(target_register, loc, refinement_t::numeric_refinement(e));
        }
        else if (width == 2) {
            interval_t to_insert = interval_t(number_t{0}, number_t{UINT16_MAX});
            expression_t e(to_insert, m_slacks);
            insert_in_registers(target_register, loc, refinement_t::numeric_refinement(e));
        }
        else {
            insert_in_registers(target_register, loc, refinement_t::numeric_refinement_top(m_slacks));
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
    auto store_at = (*offset_singleton + offset).cast_to<uint64_t>();
    auto overlapping_cells = find_overlapping_cells_in_stack(store_at, width);
    remove_overlap_in_stack(overlapping_cells, store_at, width);
    store_in_stack(b, store_at, width);
}

void interval_domain_t::check_valid_access(const ValidAccess& s, interval_t&& interval,
        int width, bool check_stack_all_numeric) {
    // access can be checked only in the signed domain
    m_signed.check_valid_access(s, std::move(interval), width, check_stack_all_numeric);
}


void interval_domain_t::shl(const register_t& reg, int imm, const int finite_width, location_t loc) {
    // The BPF ISA requires masking the imm.
    imm &= finite_width - 1;

    if (auto interval_opt = find_unsigned_interval_value(reg)) {
        interval_t interval = interval_opt->get_interval_value();
        if (interval.finite_size()) {
            const number_t lb = interval.lb().number().value();
            const number_t ub = interval.ub().number().value();
            uint64_t lb_n = lb.cast_to<uint64_t>();
            uint64_t ub_n = ub.cast_to<uint64_t>();
            uint64_t uint_max = (finite_width == 64) ? uint64_t{std::numeric_limits<uint64_t>::max()} :
                    uint64_t{std::numeric_limits<uint32_t>::max()};
            if (lb_n >> (finite_width - imm) != ub_n >> (finite_width - imm)) {
                // The bits that will be shifted out to the left are different,
                // which means all combinations of remaining bits are possible.
                lb_n = 0;
                ub_n = uint_max << imm & uint_max;
            } else {
                // The bits that will be shifted out to the left are identical
                // for all values in the interval, so we can safely shift left
                // to get a new interval.
                lb_n = lb_n << imm & uint_max;
                ub_n = ub_n << imm & uint_max;
            }
            insert_in_registers_unsigned(reg, loc, interval_t{lb_n, ub_n});
            if (to_signed(ub_n) >= to_signed(lb_n)) {
                insert_in_registers_signed(reg, loc, interval_t{lb_n, ub_n});
            } else {
                insert_in_registers_signed(reg, loc, refinement_t::numeric_refinement_top(m_slacks));
            }
            return;
        }
    }
    shl_overflow(reg, number_t{imm}, loc);
}

void interval_domain_t::lshr(const register_t& reg, int imm, const int finite_width, location_t loc) {
    // The BPF ISA requires masking the imm.
    imm &= finite_width - 1;

    if (auto interval_opt = find_unsigned_interval_value(reg)) {
        interval_t interval = interval_opt->get_interval_value();
        number_t lb_n{0};
        number_t ub_n{std::numeric_limits<uint64_t>::max() >> imm};
        if (interval.finite_size()) {
            number_t lb = interval.lb().number().value();
            number_t ub = interval.ub().number().value();
            if (finite_width == 64) {
                lb_n = lb.cast_to<uint64_t>() >> imm;
                ub_n = ub.cast_to<uint64_t>() >> imm;
            } else {
                number_t lb_w = lb.cast_to_sint(finite_width);
                number_t ub_w = ub.cast_to_sint(finite_width);
                lb_n = lb_w.cast_to<uint32_t>() >> imm;
                ub_n = ub_w.cast_to<uint32_t>() >> imm;

                // The interval must be valid since a signed range crossing 0
                // was earlier converted to a full unsigned range.
                assert(lb_n <= ub_n);
            }
        }
        insert_in_registers_unsigned(reg, loc, interval_t{lb_n, ub_n});
        if (ub_n.narrow<int64_t>() >= lb_n.narrow<int64_t>()) {
            insert_in_registers_signed(reg, loc, interval_t{lb_n, ub_n});
        } else {
            insert_in_registers_signed(reg, loc, refinement_t::numeric_refinement_top(m_slacks));
        }
        return;
    }
    insert_in_registers_unsigned(reg, loc, refinement_t::numeric_refinement_top(m_slacks));
    insert_in_registers_signed(reg, loc, refinement_t::numeric_refinement_top(m_slacks));
}

void interval_domain_t::do_bin(const Bin& bin, const std::optional<interval_t>& subtracted_opt,
                               location_t loc) {
    
    using Op = Bin::Op;

    auto dst_register = register_t{bin.dst.v};
    auto finite_width = (bin.is64 ? 64 : 32);

    if (subtracted_opt.has_value()) {
        interval_t dst_signed = *subtracted_opt;
        interval_t dst_unsigned = *subtracted_opt;
        if (!(dst_signed <= interval_t::signed_int(64))) {
            dst_signed = dst_signed.truncate_to_sint(64);
        }
        insert_in_registers_signed(dst_register, loc, dst_signed);
        if (!(dst_unsigned <= interval_t::unsigned_int(64))) {
            dst_unsigned = dst_unsigned.truncate_to_uint(64);
        }
        insert_in_registers_unsigned(dst_register, loc,
                                     refinement_t::numeric_refinement(dst_unsigned, m_slacks));
        return;
    }

    bool is_numeric_dst = find_signed_interval_value(dst_register).has_value();
    bool is_numeric_src = std::holds_alternative<Imm>(bin.v) ||
            find_signed_interval_value(register_t{std::get<Reg>(bin.v).v}).has_value();

    if (!is_numeric_dst && bin.op != Op::MOV) {
        operator-=(dst_register);
        return;
    }

    if (auto pimm = std::get_if<Imm>(&bin.v)) {
        int64_t imm;
        if (bin.is64) {
            // Use the full signed value.
            imm = to_signed(pimm->v);
        } else {
            // Use only the low 32 bits of the value.
            imm = gsl::narrow_cast<int32_t>(pimm->v);
            bitwise_and(dst_register, number_t{std::numeric_limits<uint32_t>::max()}, loc);
        }
        auto imm_interval = interval_t{number_t{imm}};
        switch (bin.op) {
            case Op::MOV: {
                // ra = imm
                refinement_t rf = refinement_t::numeric_refinement(imm_interval, m_slacks);
                m_signed.insert_in_registers(dst_register, loc, rf);
                m_unsigned.insert_in_registers(dst_register, loc, rf);
                overflow(dst_register, finite_width, loc, false);
                break;
            }
            case Op::MOVSX8:
            case Op::MOVSX16:
            case Op::MOVSX32: m_errors.push_back("MOVSX not implemented"); break;
            case Op::ADD: {
                // ra += imm
                if (imm == 0) {
                    return;
                }
                add_overflow(dst_register, number_t{gsl::narrow<int>(imm)}, finite_width, loc);
                break;
            }
            case Op::SUB: {
                // ra -= imm
                if (imm == 0) {
                    return;
                }
                add_overflow(dst_register, number_t{gsl::narrow<int>(-imm)}, finite_width, loc);
                break;
            }
            case Op::MUL: {
                // ra *= imm
                mul(dst_register, number_t{imm}, finite_width, loc);
                break;
            }
            case Op::UDIV: {
                // ra /= imm
                udiv(dst_register, number_t{imm}, finite_width, loc);
                break;
            }
            case Op::SDIV: {
                // ra s/= imm
                sdiv(dst_register, number_t{imm}, finite_width, loc);
                break;
            }
            case Op::UMOD: {
                // ra %= imm
                urem(dst_register, number_t{imm}, finite_width, loc);
                break;
            }
            case Op::SMOD: {
                // ra s%= imm
                srem(dst_register, number_t{imm}, finite_width, loc);
                break;
            }
            case Op::AND: {
                // ra &= imm
                bitwise_and(dst_register, number_t{imm}, loc);
                if (gsl::narrow<int32_t>(imm) > 0) {
                    // AND with immediate is only a 32-bit operation so svalue and uvalue
                    // are the same.
                    auto dst_signed = m_signed.find_interval_value(dst_register)->get_interval_value();
                    auto lb = dst_signed.lb().number().value();
                    auto ub = dst_signed.ub().number().value();
                    dst_signed = dst_signed & interval_t{number_t{0}, number_t{imm}};
                    refinement_t rf = refinement_t::numeric_refinement(dst_signed, m_slacks);
                    m_signed.insert_in_registers(dst_register, loc, rf);
                    m_unsigned.insert_in_registers(dst_register, loc, rf);
                }
                break;
            }
            case Op::OR: {
                // ra |= imm
                bitwise_or(dst_register, number_t{imm}, loc);
                break;
            }
            case Op::XOR: {
                // ra ^= imm
                bitwise_xor(dst_register, number_t{imm}, loc);
                break;
            }
            case Op::LSH: {
                // ra <<= imm
                shl(dst_register, gsl::narrow<int32_t>(imm), finite_width, loc);
                break;
            }
            case Op::RSH: {
                // ra >>= imm
                lshr(dst_register, gsl::narrow<int32_t>(imm), finite_width, loc);
                break;
            }
            case Op::ARSH: {
                // ra >>>= imm
                //ashr(dst_register, gsl::narrow<int32_t>(imm), finite_width, loc);
                m_signed.insert_in_registers(dst_register, loc,
                                             refinement_t::numeric_refinement_top(m_slacks));
                m_unsigned.insert_in_registers(dst_register, loc,
                                               refinement_t::numeric_refinement_top(m_slacks));
                break;
            }
            default: {
                break;
            }
        }
    }
    else {
        if (!is_numeric_src) {
            operator-=(dst_register);
            return;
        }
        register_t src_register = register_t{std::get<Reg>(bin.v).v};
        switch (bin.op) {
            case Op::MOVSX8:
            case Op::MOVSX16:
            case Op::MOVSX32:
                m_signed.insert_in_registers(dst_register, loc,
                                            refinement_t::numeric_refinement_top(m_slacks));
                m_unsigned.insert_in_registers(dst_register, loc,
                                            refinement_t::numeric_refinement_top(m_slacks));
                break;
            case Op::MOV: {
                // ra = rb
                auto src_signed_rf = m_signed.find_interval_value(src_register);
                auto src_unsigned_rf = m_unsigned.find_interval_value(src_register);
                m_signed.insert_in_registers(dst_register, loc, *src_signed_rf);
                m_unsigned.insert_in_registers(dst_register, loc, *src_unsigned_rf);
                break;
            }
            case Op::ADD: {
                // ra += rb
                add_overflow(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::SUB: {
                // ra -= rb
                if (src_register == dst_register) {
                    // an extra check only to pass a test
                    interval_t zero = interval_t{number_t{0}};
                    m_signed.insert_in_registers(dst_register, loc, zero);
                    m_unsigned.insert_in_registers(dst_register, loc, zero);
                    break;
                }
                apply_signed(arith_binaryop_t::SUB, dst_register, dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::MUL: {
                // ra *= rb
                mul(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::UDIV: {
                // ra /= rb
                udiv(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::SDIV: {
                // ra s/= rb
                sdiv(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::UMOD: {
                // ra %= rb
                urem(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::SMOD: {
                // ra s%= rb
                srem(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::AND: {
                // ra &= rb
                bitwise_and(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::OR: {
                // ra |= rb
                bitwise_or(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::XOR: {
                // ra ^= rb
                bitwise_xor(dst_register, src_register, finite_width, loc);
                break;
            }
            case Op::LSH: {
                // ra <<= rb
                if (auto src_unsigned_interval_opt = find_unsigned_interval_value(src_register)) {
                    auto src_unsigned = src_unsigned_interval_opt->get_interval_value();
                    if (std::optional<number_t> sn = src_unsigned.singleton()) {
                        uint64_t imm = sn->cast_to<int32_t>() & (bin.is64 ? 63 : 31);
                        if (imm <= std::numeric_limits<int32_t>::max()) {
                            if (!bin.is64) {
                                // Use only the low 32 bits of the value.
                                bitwise_and(dst_register, std::numeric_limits<uint32_t>::max(), loc);
                            }
                            shl(dst_register, gsl::narrow_cast<int32_t>(imm), finite_width, loc);
                            break;
                        }
                    }
                }
                shl_overflow(dst_register, register_t{src_register}, loc);
                break;
            }
            case Op::RSH: {
                // ra >>= rb
                if (auto src_unsigned_interval_opt = find_unsigned_interval_value(src_register)) {
                    auto src_unsigned = src_unsigned_interval_opt->get_interval_value();
                    if (std::optional<number_t> sn = src_unsigned.singleton()) {
                        uint64_t imm = sn->cast_to<uint64_t>() & (bin.is64 ? 63 : 31);
                        if (imm <= std::numeric_limits<int32_t>::max()) {
                            if (!bin.is64) {
                                // Use only the low 32 bits of the value.
                                bitwise_and(dst_register, std::numeric_limits<uint32_t>::max(), loc);
                            }
                            lshr(dst_register, gsl::narrow_cast<int32_t>(imm), finite_width, loc);
                            break;
                        }
                    }
                }
                m_signed.insert_in_registers(dst_register, loc,
                                             refinement_t::numeric_refinement_top(m_slacks));
                m_unsigned.insert_in_registers(dst_register, loc,
                                             refinement_t::numeric_refinement_top(m_slacks));
                break;
            }
            case Op::ARSH: {
                // ra >>>= rb
                //if (auto src_signed_interval_opt = find_signed_interval_value(src_register)) {
                //    ashr(dst_register, src_register, finite_width, loc);
                //    break;
                //}
                m_signed.insert_in_registers(dst_register, loc,
                                             refinement_t::numeric_refinement_top(m_slacks));
                m_unsigned.insert_in_registers(dst_register, loc,
                                             refinement_t::numeric_refinement_top(m_slacks));
                break;
            }
            default: {
                break;
            }
        }
    }
    if (!bin.is64) {
        bitwise_and(dst_register, std::numeric_limits<uint32_t>::max(), loc);
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

void interval_domain_t::operator()(const basic_block_t& bb) {
    // nothing to do here
}
void interval_domain_t::set_require_check(check_require_func_t f) {}

} // namespace crab
