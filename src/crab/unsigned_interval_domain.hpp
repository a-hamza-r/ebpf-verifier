// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "array_domain.hpp"
#include "types.hpp"

namespace crab {

using check_require_func_t = std::function<bool(crab::domains::NumAbsDomain&, const crab::linear_constraint_t&, std::string)>;
using global_env_unsigned_registers_t = std::unordered_map<register_location_t, mock_interval_t>;

class unsigned_interval_registers_t {

    live_registers_t m_cur_register_def;
    std::shared_ptr<global_env_unsigned_registers_t> m_registers_env;
    bool m_is_bottom = false;

  public:
    unsigned_interval_registers_t() = default;
    unsigned_interval_registers_t(std::shared_ptr<global_env_unsigned_registers_t> registers_env)
        : m_registers_env(std::move(registers_env)) {}
    bool is_bottom() const;
    bool is_top() const;
    void set_to_bottom();
    void set_to_top();
    std::optional<mock_interval_t> find(register_location_t reg) const;
    std::optional<mock_interval_t> find(register_t key) const;
    void insert(register_t, const location_t&, interval_t);
    void operator-=(register_t);
    unsigned_interval_registers_t operator|(const unsigned_interval_registers_t& other) const;
    void adjust_bb_for_registers(location_t);
};

using unsigned_interval_stack_cell_t = std::pair<mock_interval_t, int>;    // intervals with width
using unsigned_interval_stack_cells_t = std::map<uint64_t, unsigned_interval_stack_cell_t>;

class unsigned_interval_stack_t {

    unsigned_interval_stack_cells_t m_cells;
    bool m_is_bottom = false;

  public:
    unsigned_interval_stack_t() = default;
    unsigned_interval_stack_t(unsigned_interval_stack_cells_t cells) : m_cells(std::move(cells)) {}
    bool is_bottom() const;
    bool is_top() const;
    void set_to_bottom();
    void set_to_top();
    static unsigned_interval_stack_t top();
    std::optional<unsigned_interval_stack_cell_t> find(uint64_t) const;
    void store(uint64_t, mock_interval_t, int);
    void operator-=(uint64_t);
    unsigned_interval_stack_t operator|(const unsigned_interval_stack_t& other) const;
    std::vector<uint64_t> get_keys() const;
    size_t size() const;
    void remove_overlap(const std::vector<uint64_t>&, uint64_t, int);
    void fill_values(const std::vector<uint64_t>&, uint64_t, int);
};

class unsigned_interval_domain_t final {
    unsigned_interval_registers_t m_registers;
    unsigned_interval_stack_t m_stack;
    std::vector<std::string> m_errors;
    bool m_is_bottom = false;

  public:

    unsigned_interval_domain_t() = default;
    unsigned_interval_domain_t(unsigned_interval_registers_t registers,
            unsigned_interval_stack_t stack) :
        m_registers(std::move(registers)), m_stack(std::move(stack)) {}
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
    void operator-=(register_t reg) { m_registers -= reg; }

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = location_t::top());
    void operator()(const Bin&, location_t loc = location_t::top());
    void operator()(const Un&, location_t loc = location_t::top());
    void operator()(const LoadMapFd&, location_t loc = location_t::top());
    void operator()(const Call&, location_t loc = location_t::top());
    void operator()(const Exit&, location_t loc = location_t::top());
    void operator()(const Jmp&, location_t loc = location_t::top());
    void operator()(const Mem&, location_t loc = location_t::top());
    void operator()(const Packet&, location_t loc = location_t::top());
    void operator()(const Assume&, location_t loc = location_t::top());
    void operator()(const Assert&, location_t loc = location_t::top());
    void operator()(const basic_block_t& bb);
    void write(std::ostream& os) const;
    crab::bound_t get_loop_count_upper_bound() const;
    void initialize_loop_counter(const label_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f);

    std::optional<mock_interval_t> find_interval_value(register_t) const;
    std::optional<mock_interval_t> find_interval_at_loc(const register_location_t reg) const;
    std::optional<unsigned_interval_stack_cell_t> find_in_stack(uint64_t) const;
    void insert_in_registers(register_t, location_t, interval_t);
    void store_in_stack(uint64_t, mock_interval_t, int);
    void adjust_bb_for_types(location_t);
    void remove_overlap_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    void fill_values_in_stack(const std::vector<uint64_t>&, uint64_t, int);
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void reset_errors() { m_errors.clear(); }
    bool load_from_stack(register_t, uint64_t, location_t);
    void store_in_stack(const Mem&, uint64_t, int);
    void print_ctx(std::ostream& o) const {};
    void print_stack(std::ostream& o) const {};
    void print_annotated_bb(std::ostream& o, const basic_block_t& bb) const {};
}; // end unsigned_interval_domain_t

} // namespace crab
