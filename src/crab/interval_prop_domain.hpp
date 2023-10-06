// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>
#include <map>

#include "crab/abstract_domain.hpp"
#include "crab/region_domain.hpp"
#include "crab/cfg.hpp"
#include "linear_constraint.hpp"
#include "string_constraints.hpp"
#include "crab/common.hpp"

namespace crab {

using live_registers_t = std::array<std::shared_ptr<reg_with_loc_t>, 11>;
using global_interval_env_t = std::unordered_map<reg_with_loc_t, mock_interval_t>;

class registers_cp_state_t {

    live_registers_t m_cur_def;
    std::shared_ptr<global_interval_env_t> m_interval_env;
    bool m_is_bottom = false;

    public:
        bool is_bottom() const;
        bool is_top() const;
        void set_to_bottom();
        void set_to_top();
        std::optional<mock_interval_t> find(reg_with_loc_t reg) const;
        std::optional<mock_interval_t> find(register_t key) const;
        void insert(register_t, const reg_with_loc_t&, interval_t);
        void operator-=(register_t);
        registers_cp_state_t operator|(const registers_cp_state_t& other) const;
        registers_cp_state_t(bool is_bottom = false) : m_interval_env(nullptr),
            m_is_bottom(is_bottom) {}
        explicit registers_cp_state_t(live_registers_t&& vars,
                std::shared_ptr<global_interval_env_t> interval_env, bool is_bottom = false)
            : m_cur_def(std::move(vars)), m_interval_env(interval_env), m_is_bottom(is_bottom) {}
        explicit registers_cp_state_t(std::shared_ptr<global_interval_env_t> interval_env,
                bool is_bottom = false)
            : m_interval_env(interval_env), m_is_bottom(is_bottom) {}
        void adjust_bb_for_registers(location_t);
};

using interval_cells_t = std::pair<mock_interval_t, int>;    // intervals with width
using interval_values_stack_t = std::map<uint64_t, interval_cells_t>;

class stack_cp_state_t {

    interval_values_stack_t m_interval_values;
    bool m_is_bottom = false;

    public:
        bool is_bottom() const;
        bool is_top() const;
        void set_to_bottom();
        void set_to_top();
        static stack_cp_state_t top();
        std::optional<interval_cells_t> find(uint64_t) const;
        void store(uint64_t, mock_interval_t, int);
        void operator-=(uint64_t);
        bool all_numeric(uint64_t, int) const;
        stack_cp_state_t operator|(const stack_cp_state_t& other) const;
        stack_cp_state_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        explicit stack_cp_state_t(interval_values_stack_t&& interval_values, bool is_bottom = false)
            : m_interval_values(std::move(interval_values)), m_is_bottom(is_bottom) {}
        std::vector<uint64_t> get_keys() const;
        size_t size() const;
        std::vector<uint64_t> find_overlapping_cells(uint64_t, int) const;
        void remove_overlap(const std::vector<uint64_t>&, uint64_t, int);
};

class interval_prop_domain_t final {
    registers_cp_state_t m_registers_interval_values;
    stack_cp_state_t m_stack_slots_interval_values;
    std::vector<std::string> m_errors;
    bool m_is_bottom = false;

  public:

    interval_prop_domain_t() = default;
    interval_prop_domain_t(interval_prop_domain_t&& o) = default;
    interval_prop_domain_t(const interval_prop_domain_t& o) = default;
    explicit interval_prop_domain_t(registers_cp_state_t&& consts_regs,
            stack_cp_state_t&& interval_stack_slots, bool is_bottom = false) :
        m_registers_interval_values(std::move(consts_regs)), m_stack_slots_interval_values(std::move(interval_stack_slots)),
        m_is_bottom(is_bottom) {}
    interval_prop_domain_t& operator=(interval_prop_domain_t&& o) = default;
    interval_prop_domain_t& operator=(const interval_prop_domain_t& o) = default;
    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static interval_prop_domain_t setup_entry();
    // bottom/top
    static interval_prop_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    bool is_bottom() const;
    bool is_top() const;
    // inclusion
    bool operator<=(const interval_prop_domain_t& other) const;
    // join
    void operator|=(const interval_prop_domain_t& abs);
    void operator|=(interval_prop_domain_t&& abs);
    interval_prop_domain_t operator|(const interval_prop_domain_t& other) const;
    interval_prop_domain_t operator|(interval_prop_domain_t&& abs) const;
    // meet
    interval_prop_domain_t operator&(const interval_prop_domain_t& other) const;
    // widening
    interval_prop_domain_t widen(const interval_prop_domain_t& other, bool);
    // narrowing
    interval_prop_domain_t narrow(const interval_prop_domain_t& other) const;
    //forget
    void operator-=(variable_t var);
    void operator-=(register_t reg) { m_registers_interval_values -= reg; }

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = boost::none, int print = 0);
    void operator()(const Bin&, location_t loc = boost::none, int print = 0);
    void operator()(const Un&, location_t loc = boost::none, int print = 0);
    void operator()(const LoadMapFd&, location_t loc = boost::none, int print = 0);
    void operator()(const Atomic&, location_t loc = boost::none, int print = 0) {}
    void operator()(const Call&, location_t loc = boost::none, int print = 0);
    void operator()(const Callx&, location_t loc = boost::none, int print = 0) {}
    void operator()(const Exit&, location_t loc = boost::none, int print = 0);
    void operator()(const Jmp&, location_t loc = boost::none, int print = 0);
    void operator()(const Mem&, location_t loc = boost::none, int print = 0);
    void operator()(const Packet&, location_t loc = boost::none, int print = 0);
    void operator()(const Assume&, location_t loc = boost::none, int print = 0);
    void operator()(const Assert&, location_t loc = boost::none, int print = 0);
    void operator()(const ValidAccess&, location_t loc = boost::none, int print = 0);
    void operator()(const Comparable&, location_t loc = boost::none, int print = 0);
    void operator()(const Addable&, location_t loc = boost::none, int print = 0);
    void operator()(const ValidStore&, location_t loc = boost::none, int print = 0);
    void operator()(const TypeConstraint&, location_t loc = boost::none, int print = 0);
    void operator()(const ValidSize&, location_t loc = boost::none, int print = 0);
    void operator()(const ValidMapKeyValue&, location_t loc = boost::none, int print = 0);
    void operator()(const ZeroCtxOffset&, location_t loc = boost::none, int print = 0);
    void operator()(const FuncConstraint&, location_t loc = boost::none, int print = 0) {}
    void operator()(const IncrementLoopCounter&, location_t loc = boost::none, int print = 0) {}
    void operator()(const basic_block_t& bb, int print = 0);
    void write(std::ostream& os) const;
    crab::bound_t get_loop_count_upper_bound();
    void initialize_loop_counter(const label_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f);

    void do_load(const Mem&, const Reg&, std::optional<ptr_or_mapfd_t>, bool, location_t);
    void do_mem_store(const Mem&, std::optional<ptr_or_mapfd_t>);
    void do_call(const Call&, const interval_values_stack_t&, location_t);
    void do_bin(const Bin&, const std::optional<interval_t>&, const std::optional<interval_t>&,
            const std::optional<ptr_or_mapfd_t>&, const std::optional<ptr_or_mapfd_t>&,
            const interval_t&, location_t);
    void check_valid_access(const ValidAccess&, interval_t&&, int, bool = false);
    void assume_unsigned_cst_interval(Condition::Op, bool, interval_t&&, interval_t&&,
            interval_t&&, interval_t&&, register_t, Value, location_t);
    void assume_signed_cst_interval(Condition::Op, bool, interval_t&&, interval_t&&,
            interval_t&&, interval_t&&, register_t, Value, location_t);
    void assume_cst(Condition::Op, bool, register_t, Value, location_t);
    void assume_lt(bool, interval_t&&, interval_t&&, interval_t&&, interval_t&&,
            register_t, Value, location_t, bool);
    void assume_gt(bool, interval_t&&, interval_t&&, interval_t&&, interval_t&&,
            register_t, Value, location_t);
    void assume_gt_and_lt(bool, bool, bool, interval_t&&, interval_t&&,
            interval_t&&, interval_t&&, register_t, Value, location_t, bool = true);
    std::optional<mock_interval_t> find_interval_value(register_t) const;
    std::optional<mock_interval_t> find_interval_at_loc(const reg_with_loc_t reg) const;
    std::optional<interval_cells_t> find_in_stack(uint64_t) const;
    void adjust_bb_for_types(location_t);
    std::vector<uint64_t> get_stack_keys() const;
    bool all_numeric_in_stack(uint64_t, int) const;
    std::vector<uint64_t> find_overlapping_cells_in_stack(uint64_t, int) const;
    void remove_overlap_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void reset_errors() { m_errors.clear(); }
}; // end interval_prop_domain_t

} // namespace crab
