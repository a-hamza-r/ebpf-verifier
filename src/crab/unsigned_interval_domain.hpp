// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "crab/common.hpp"

namespace crab {

using live_registers_t = std::array<std::shared_ptr<reg_with_loc_t>, NUM_REGISTERS>;
using global_interval_env_t = std::unordered_map<reg_with_loc_t, mock_interval_t>;

class registers_unsigned_state_t {

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
        void insert(register_t, const location_t&, interval_t);
        void operator-=(register_t);
        bool operator<=(const registers_unsigned_state_t& other) const;
        registers_unsigned_state_t operator|(const registers_unsigned_state_t& other) const;
        registers_unsigned_state_t(bool is_bottom = false) : m_interval_env(nullptr),
            m_is_bottom(is_bottom) {}
        explicit registers_unsigned_state_t(live_registers_t&& vars,
                std::shared_ptr<global_interval_env_t> interval_env, bool is_bottom = false)
            : m_cur_def(std::move(vars)), m_interval_env(interval_env), m_is_bottom(is_bottom) {}
        explicit registers_unsigned_state_t(std::shared_ptr<global_interval_env_t> interval_env,
                bool is_bottom = false)
            : m_interval_env(interval_env), m_is_bottom(is_bottom) {}
        void adjust_bb_for_registers(location_t);
};

using interval_cells_t = std::pair<mock_interval_t, int>;    // intervals with width
using interval_values_stack_t = std::map<uint64_t, interval_cells_t>;

class stack_slots_unsigned_state_t {

    interval_values_stack_t m_interval_values;
    bool m_is_bottom = false;

    public:
        bool is_bottom() const;
        bool is_top() const;
        void set_to_bottom();
        void set_to_top();
        static stack_slots_unsigned_state_t top();
        std::optional<interval_cells_t> find(uint64_t) const;
        void store(uint64_t, mock_interval_t, int);
        void operator-=(uint64_t);
        bool operator<=(const stack_slots_unsigned_state_t& other) const;
        stack_slots_unsigned_state_t operator|(const stack_slots_unsigned_state_t& other) const;
        stack_slots_unsigned_state_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        explicit stack_slots_unsigned_state_t(interval_values_stack_t&& interval_values, bool is_bottom = false)
            : m_interval_values(std::move(interval_values)), m_is_bottom(is_bottom) {}
        std::vector<uint64_t> get_keys() const;
        size_t size() const;
        void remove_overlap(const std::vector<uint64_t>&, uint64_t, int);
        void fill_values(const std::vector<uint64_t>&, uint64_t, int);
};

class unsigned_interval_domain_t final {
    registers_unsigned_state_t m_registers_values;
    stack_slots_unsigned_state_t m_stack_slots_values;
    std::vector<std::string> m_errors;
    bool m_is_bottom = false;

  public:

    unsigned_interval_domain_t() = default;
    unsigned_interval_domain_t(unsigned_interval_domain_t&& o) = default;
    unsigned_interval_domain_t(const unsigned_interval_domain_t& o) = default;
    explicit unsigned_interval_domain_t(registers_unsigned_state_t&& consts_regs,
            stack_slots_unsigned_state_t&& interval_stack_slots, bool is_bottom = false) :
        m_registers_values(std::move(consts_regs)), m_stack_slots_values(std::move(interval_stack_slots)),
        m_is_bottom(is_bottom) {}
    unsigned_interval_domain_t& operator=(unsigned_interval_domain_t&& o) = default;
    unsigned_interval_domain_t& operator=(const unsigned_interval_domain_t& o) = default;
    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static unsigned_interval_domain_t setup_entry();
    // bottom/top
    static unsigned_interval_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    void set_registers_to_bottom();
    void set_registers_to_top();
    bool is_bottom() const;
    bool is_top() const;
    // inclusion
    bool operator<=(const unsigned_interval_domain_t& other) const;
    // join
    void operator|=(const unsigned_interval_domain_t& abs);
    void operator|=(unsigned_interval_domain_t&& abs);
    unsigned_interval_domain_t operator|(const unsigned_interval_domain_t& other) const;
    unsigned_interval_domain_t operator|(unsigned_interval_domain_t&& abs) const;
    // meet
    unsigned_interval_domain_t operator&(const unsigned_interval_domain_t& other) const;
    // widening
    unsigned_interval_domain_t widen(const unsigned_interval_domain_t& other, bool);
    // narrowing
    unsigned_interval_domain_t narrow(const unsigned_interval_domain_t& other) const;
    //forget
    void operator-=(register_t reg) { m_registers_values -= reg; }

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
    void write(std::ostream& os) const;
    crab::bound_t get_loop_count_upper_bound();
    void initialize_loop_counter(const label_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f);

    std::optional<mock_interval_t> find_interval_value(register_t) const;
    std::optional<mock_interval_t> find_interval_at_loc(const reg_with_loc_t reg) const;
    std::optional<interval_cells_t> find_in_stack(uint64_t) const;
    void insert_in_registers(register_t, location_t, interval_t);
    void store_in_stack(uint64_t, mock_interval_t, int);
    void adjust_bb_for_types(location_t);
    void remove_overlap_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    void fill_values_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void reset_errors() { m_errors.clear(); }
    bool load_from_stack(register_t, uint64_t, location_t);
    void store_in_stack(const Mem&, uint64_t, int);
}; // end unsigned_interval_domain_t

} // namespace crab
