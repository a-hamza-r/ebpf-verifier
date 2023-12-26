// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "crab/common.hpp"
#include "crab/signed_interval_domain.hpp"
#include "crab/unsigned_interval_domain.hpp"

namespace crab {

class interval_domain_t final {
    signed_interval_domain_t m_signed;
    unsigned_interval_domain_t m_unsigned;
    std::vector<std::string> m_errors;
    bool m_is_bottom = false;

  public:

    interval_domain_t() = default;
    interval_domain_t(interval_domain_t&& o) = default;
    interval_domain_t(const interval_domain_t& o) = default;
    explicit interval_domain_t(signed_interval_domain_t&& signed_domain,
            unsigned_interval_domain_t&& unsigned_domain, bool is_bottom = false) :
        m_signed(std::move(signed_domain)), m_unsigned(std::move(unsigned_domain)),
        m_is_bottom(is_bottom) {}
    interval_domain_t& operator=(interval_domain_t&& o) = default;
    interval_domain_t& operator=(const interval_domain_t& o) = default;
    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static interval_domain_t setup_entry();
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

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = boost::none);
    void operator()(const Bin&, location_t loc = boost::none);
    void operator()(const Un&, location_t loc = boost::none);
    void operator()(const LoadMapFd&, location_t loc = boost::none);
    void operator()(const Call&, location_t loc = boost::none);
    void operator()(const Exit&, location_t loc = boost::none);
    void operator()(const Jmp&, location_t loc = boost::none);
    void operator()(const Mem&, location_t loc = boost::none);
    void operator()(const Packet&, location_t loc = boost::none);
    void operator()(const Assume&, location_t loc = boost::none);
    void operator()(const Assert&, location_t loc = boost::none);
    void operator()(const basic_block_t& bb, int print = 0);
    void write(std::ostream& os) const {}
    crab::bound_t get_loop_count_upper_bound();
    void initialize_loop_counter(const label_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f);

    void do_load(const Mem&, const register_t&, std::optional<ptr_or_mapfd_t>, bool, location_t);
    void do_mem_store(const Mem&, std::optional<ptr_or_mapfd_t>);
    void do_call(const Call&, const stack_cells_t&, location_t);
    void do_bin(const Bin&, const std::optional<interval_t>&,
            const std::optional<interval_t>&, const std::optional<ptr_or_mapfd_t>&,
            const std::optional<interval_t>&, const std::optional<interval_t>&,
            const std::optional<ptr_or_mapfd_t>&, const interval_t&, location_t);
    void check_valid_access(const ValidAccess&, interval_t&&, int = -1, bool = false);
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
    std::optional<mock_interval_t> find_interval_value(register_t) const;
    std::optional<mock_interval_t> find_signed_interval_value(register_t) const;
    std::optional<mock_interval_t> find_unsigned_interval_value(register_t) const;
    std::optional<mock_interval_t> find_signed_interval_at_loc(const reg_with_loc_t reg) const;
    std::optional<mock_interval_t> find_unsigned_interval_at_loc(const reg_with_loc_t reg) const;
    std::optional<interval_cells_t> find_in_stack_signed(uint64_t) const;
    std::optional<interval_cells_t> find_in_stack_unsigned(uint64_t) const;
    void insert_in_registers(register_t, location_t, interval_t);
    void insert_in_registers_signed(register_t, location_t, interval_t);
    void insert_in_registers_unsigned(register_t, location_t, interval_t);
    void store_in_stack(uint64_t, mock_interval_t, int);
    void store_in_stack_signed(uint64_t, mock_interval_t, int);
    void store_in_stack_unsigned(uint64_t, mock_interval_t, int);
    void adjust_bb_for_types(location_t);
    std::vector<uint64_t> get_stack_keys() const;
    bool all_numeric_in_stack(uint64_t, int) const;
    std::vector<uint64_t> find_overlapping_cells_in_stack(uint64_t, int) const;
    void remove_overlap_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    void fill_values_in_stack(const std::vector<uint64_t>&, uint64_t, int);
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
