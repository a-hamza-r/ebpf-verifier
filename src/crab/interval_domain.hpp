// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "signed_interval_domain.hpp"
#include "unsigned_interval_domain.hpp"

namespace crab {

using check_require_func_t = std::function<bool(crab::domains::NumAbsDomain&, const crab::linear_constraint_t&, std::string)>;
enum class arith_binaryop_t { ADD, SUB, MUL, SDIV, UDIV, SREM, UREM };
enum class bitwise_binaryop_t { AND, OR, XOR, SHL, LSHR, ASHR };
using binaryop_t = std::variant<arith_binaryop_t, bitwise_binaryop_t>;

class interval_domain_t final {
    std::shared_ptr<slacks_t> m_slacks = nullptr;
    signed_interval_domain_t m_signed;
    unsigned_interval_domain_t m_unsigned;
    std::vector<std::string> m_errors;

  public:

    interval_domain_t() : m_slacks(std::make_shared<slacks_t>()), m_signed(m_slacks),
        m_unsigned(m_slacks) {}
    interval_domain_t(std::shared_ptr<slacks_t> slacks) : m_slacks(std::move(slacks)),
        m_signed(m_slacks), m_unsigned(m_slacks) {}
    interval_domain_t(signed_interval_domain_t signed_domain,
            unsigned_interval_domain_t unsigned_domain, std::shared_ptr<slacks_t> slacks) :
        m_slacks(std::move(slacks)), m_signed(std::move(signed_domain)),
        m_unsigned(std::move(unsigned_domain)) {}

    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static interval_domain_t setup_entry(std::shared_ptr<slacks_t>);
    // bottom/top
    static interval_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    void set_registers_to_bottom();
    void set_registers_to_top();
    bool is_bottom() const;
    bool is_top() const;
    // inclusion
    bool operator<=(const interval_domain_t& other) const;
    // join
    void operator|=(const interval_domain_t& abs);
    void operator|=(interval_domain_t&& abs);
    interval_domain_t operator|(const interval_domain_t& other) const;
    interval_domain_t operator|(interval_domain_t&& abs) const;
    // meet
    interval_domain_t operator&(const interval_domain_t& other) const;
    // widening
    interval_domain_t widen(const interval_domain_t& other, bool);
    // narrowing
    interval_domain_t narrow(const interval_domain_t& other) const;
    void operator-=(register_t reg);


    void apply_signed(const binaryop_t&, const register_t&, const register_t&, const number_t&, const int, location_t);
    void apply_unsigned(const binaryop_t&, const register_t&, const register_t&, const number_t&, const int, location_t);
    void apply_signed(const binaryop_t&, const register_t&, const register_t&, const register_t&, const int, location_t);
    void apply_unsigned(const binaryop_t&, const register_t&, const register_t&, const register_t&, const int, location_t);

    void overflow(const register_t&, const int, location_t, bool);
    void overflow_bounds(const register_t&, number_t, const int, location_t, bool);

    void apply(const bitwise_binaryop_t&, const register_t&, const register_t&, const number_t&, const int, location_t, bool);
    void apply(const bitwise_binaryop_t&, const register_t&, const register_t&, const register_t&, const int, location_t, bool);
    void apply(const arith_binaryop_t&, const register_t&, const register_t&, const number_t&, const int, location_t, bool);
    void apply(const arith_binaryop_t&, const register_t&, const register_t&, const register_t&, const int, location_t, bool);

    void apply(const binaryop_t& op, const register_t& result, const register_t& lhs, const register_t& rhs, const int finite_width, location_t loc, bool is_signed) {
        std::visit([&](auto top) { apply(top, result, lhs, rhs, finite_width, loc, is_signed); }, op);
    }
    void apply(const binaryop_t& op, const register_t& result, const register_t& lhs, const number_t& rhs, const int finite_width, location_t loc, bool is_signed) {
        std::visit([&](auto top) { apply(top, result, lhs, rhs, finite_width, loc, is_signed); }, op);
    }


    void neg(const register_t&, const int, location_t);
    void add(const register_t&, const register_t&, location_t);
    void add(const register_t&, const number_t&, location_t);
    void sub(const register_t&, const register_t&, location_t);
    void sub(const register_t&, const number_t&, location_t);
    void add_overflow(const register_t&, const register_t&, const int, location_t);
    void add_overflow(const register_t&, const number_t&, const int, location_t);
    void sub_overflow(const register_t&, const register_t&, const int, location_t);
    void sub_overflow(const register_t&, const number_t&, const int, location_t);
    void mul(const register_t&, const register_t&, const int, location_t);
    void mul(const register_t&, const number_t&, const int, location_t);
    void udiv(const register_t&, const register_t&, const int, location_t);
    void udiv(const register_t&, const number_t&, const int, location_t);
    void sdiv(const register_t&, const register_t&, const int, location_t);
    void sdiv(const register_t&, const number_t&, const int, location_t);
    void srem(const register_t&, const register_t&, const int, location_t);
    void srem(const register_t&, const number_t&, const int, location_t);
    void urem(const register_t&, const register_t&, const int, location_t);
    void urem(const register_t&, const number_t&, const int, location_t);
    void bitwise_and(const register_t&, const register_t&, const int, location_t);
    void bitwise_and(const register_t&, const number_t&, location_t);
    void bitwise_or(const register_t&, const register_t&, const int, location_t);
    void bitwise_or(const register_t&, const number_t&, location_t);
    void bitwise_xor(const register_t&, const register_t&, const int, location_t);
    void bitwise_xor(const register_t&, const number_t&, location_t);
    void shl_overflow(const register_t&, const register_t&, location_t);
    void shl_overflow(const register_t&, const number_t&, location_t);

    void shl(const register_t&, int, const int, location_t);
    void lshr(const register_t&, int, const int, location_t);

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = location_t::top());
    void operator()(const Bin&, location_t loc = location_t::top());
    void operator()(const Un&, location_t loc = location_t::top());
    void operator()(const LoadMapFd&, location_t loc = location_t::top());
    void operator()(const LoadVariable&, location_t loc = location_t::top());
    void operator()(const Call&, location_t loc = location_t::top());
    void operator()(const Exit&, location_t loc = location_t::top());
    void operator()(const Jmp&, location_t loc = location_t::top());
    void operator()(const Mem&, location_t loc = location_t::top());
    void operator()(const Packet&, location_t loc = location_t::top());
    void operator()(const Assume&, location_t loc = location_t::top());
    void operator()(const Assert&, location_t loc = location_t::top());
    void operator()(const basic_block_t& bb);
    void write(std::ostream& os) const {}
    crab::bound_t get_loop_count_upper_bound() const;
    void initialize_loop_counter(const label_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f);

    void do_load(const Mem&, const register_t&, std::optional<ptr_or_mapfd_t>, bool, bool, location_t);
    void do_mem_store(const Mem&, std::optional<ptr_or_mapfd_t>);
    void do_call(const Call&, const stack_cells_t&, location_t);
    void do_bin(const Bin&, std::optional<interval_t>, location_t);
    void check_valid_access(const ValidAccess&, interval_t, interval_t, int, bool, location_t);
    void assume_cst(Condition::Op, bool, register_t, Value, location_t);
    void assume_signed_cst(Condition::Op, bool, const interval_t&, const interval_t&,
            const interval_t&, const interval_t&, register_t, Value, location_t);
    void assume_signed_lt(bool, bool, interval_t&&, interval_t&&, const interval_t&,
            const interval_t&, const interval_t&, const interval_t&, register_t, Value, location_t);
    void assume_signed_gt(bool, bool, interval_t&&, interval_t&&, const interval_t&,
            const interval_t&, const interval_t&, const interval_t&, register_t, Value, location_t);
    void assume_unsigned_cst(Condition::Op, bool, const interval_t&, const interval_t&,
            const interval_t&, const interval_t&, register_t, Value, location_t);
    void assume_unsigned_lt(bool, bool, interval_t&&, interval_t&&, const interval_t&,
            const interval_t&, const interval_t&, const interval_t&, register_t, Value, location_t);
    void assume_unsigned_gt(bool, bool, interval_t&&, interval_t&&, const interval_t&,
            const interval_t&, const interval_t&, const interval_t&, register_t, Value, location_t);
    void update_gt(bool, bool, interval_t&&, interval_t&&, const interval_t&, const interval_t&,
            const interval_t&, const interval_t&, register_t, Value, location_t,
            interval_t&&, interval_t&&, bool, bool, bool, bool);
    void update_lt(bool, bool, interval_t&&, interval_t&&, const interval_t&, const interval_t&,
            const interval_t&, const interval_t&, register_t, Value, location_t,
            interval_t&&, interval_t&&, bool, bool, bool, bool);
    std::optional<refinement_t> find_interval_value(register_t) const;
    std::optional<refinement_t> find_signed_interval_value(register_t) const;
    std::optional<refinement_t> find_unsigned_interval_value(register_t) const;
    std::optional<refinement_t> find_signed_interval_at_loc(const register_location_t reg) const;
    std::optional<refinement_t> find_unsigned_interval_at_loc(const register_location_t reg) const;
    std::optional<signed_interval_stack_cell_t> find_in_stack_signed(uint64_t) const;
    std::optional<unsigned_interval_stack_cell_t> find_in_stack_unsigned(uint64_t) const;
    void insert_in_registers(register_t, location_t, refinement_t);
    void insert_in_registers_signed(register_t, location_t, refinement_t);
    void insert_in_registers_signed(register_t, location_t, interval_t);
    void insert_in_registers_unsigned(register_t, location_t, refinement_t);
    void insert_in_registers_unsigned(register_t, location_t, interval_t);
    void store_in_stack(uint64_t, refinement_t, int);
    void store_in_stack_signed(uint64_t, refinement_t, int);
    void store_in_stack_unsigned(uint64_t, refinement_t, int);
    void adjust_bb_for_types(location_t);
    [[nodiscard]] std::vector<uint64_t> get_stack_keys() const;
    bool all_numeric_in_stack(uint64_t, int) const;
    std::vector<uint64_t> find_overlapping_cells_in_stack(uint64_t, int) const;
    void remove_overlap_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    void fill_values_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    void print_stack(std::ostream& o) const {};
    void print_state_init(std::ostream& o, label_t) const {}
    void print_ctx(std::ostream& o) const {};
    void print_annotated_bb(std::ostream& o, const basic_block_t& bb) const {};
    [[nodiscard]] std::vector<std::string>& get_errors() {
        operator+=(m_signed.get_errors());
        operator+=(m_unsigned.get_errors());
        return m_errors;
    }
    void reset_errors() { 
        m_errors.clear();
        m_signed.reset_errors();
        m_unsigned.reset_errors();
    }
    void operator+=(std::vector<std::string>& errs) {
        m_errors.insert(m_errors.end(), errs.begin(), errs.end());
    }

  private:
    void scratch_caller_saved_registers();
    bool load_from_stack(register_t, interval_t, int, location_t);
    void store_in_stack(const Mem&, uint64_t, int);
}; // end interval_domain_t

} // namespace crab
