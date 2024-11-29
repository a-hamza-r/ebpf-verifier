// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "array_domain.hpp"
#include "refinement.hpp"

namespace crab {

using check_require_func_t = std::function<bool(crab::domains::NumAbsDomain&, const crab::linear_constraint_t&, std::string)>;
using global_env_offset_registers_t = std::unordered_map<register_location_t, refinement_t>;

class offset_registers_t {

    live_registers_t m_cur_register_def;
    std::shared_ptr<global_env_offset_registers_t> m_registers_env;
    std::shared_ptr<slacks_t> m_slacks;
    bool m_is_bottom = false;

    public:
        offset_registers_t(bool is_bottom = false) : m_registers_env(nullptr), m_slacks(nullptr),
        m_is_bottom(is_bottom) {}
        offset_registers_t(std::shared_ptr<global_env_offset_registers_t> registers_env,
                std::shared_ptr<slacks_t> slacks, bool is_bottom = false)
            : m_registers_env(registers_env), m_slacks(slacks), m_is_bottom(is_bottom) {}
        offset_registers_t(std::shared_ptr<global_env_offset_registers_t> registers_env,
                std::shared_ptr<slacks_t> slacks, const ebpf_context_descriptor_t* desc,
                bool is_bottom = false)
            : m_registers_env(registers_env), m_slacks(slacks), m_is_bottom(is_bottom) {

            location_t loc{label_t::entry, 0};
            if (desc->data >= 0) {
                insert(register_t{12}, loc, refinement_t::begin());
            }
        }

        explicit offset_registers_t(live_registers_t&& vars,
                std::shared_ptr<global_env_offset_registers_t> registers_env,
                std::shared_ptr<slacks_t> slacks, bool is_bottom = false)
            : m_cur_register_def(std::move(vars)), m_registers_env(registers_env), m_slacks(slacks),
            m_is_bottom(is_bottom) {}

        offset_registers_t operator|(const offset_registers_t&) const;
        void operator-=(register_t);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        void insert(register_t, const location_t&, refinement_t&&);
        void insert_slack_value(symbol_t, mock_interval_t);
        std::shared_ptr<slacks_t> get_slacks() const { return m_slacks; }
        std::optional<mock_interval_t> find_slack_value(symbol_t) const;
        std::optional<refinement_t> find(register_location_t reg) const;
        std::optional<refinement_t> find(register_t key) const;
        friend std::ostream& operator<<(std::ostream& o, const offset_registers_t& p);
        void adjust_bb_for_registers(location_t);
        void scratch_caller_saved_registers();
        void forget_packet_pointers(location_t);
};

using refinement_stack_cell_t = std::pair<refinement_t, int>;
using refinement_stack_cells_t = std::map<uint64_t, refinement_stack_cell_t>;

class offset_stack_t {
    refinement_stack_cells_t m_stack_cells;
    bool m_is_bottom = false;

    public:
        offset_stack_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        std::optional<refinement_stack_cell_t> find(uint64_t) const;
        void store(uint64_t, refinement_t, int);
        void operator-=(uint64_t);
        void operator-=(const std::vector<uint64_t>&);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        static offset_stack_t top();
        offset_stack_t operator|(const offset_stack_t&) const;
        explicit offset_stack_t(refinement_stack_cells_t&& cells, bool is_bottom = false)
            : m_stack_cells(std::move(cells)), m_is_bottom(is_bottom) {}
        std::vector<uint64_t> find_overlapping_cells(uint64_t, int) const;
        std::vector<uint64_t> get_keys() const;
};

class offset_ctx_t {
    using refinement_ctx_cells_t = std::unordered_map<uint64_t, refinement_t>;    // represents `cp[n] = rf;`
    refinement_ctx_cells_t m_ctx_cells;
    size_t m_size;

    public:
        offset_ctx_t(const ebpf_context_descriptor_t* desc);
        std::optional<refinement_t> find(uint64_t) const;
        size_t get_size() const { return m_size; }
        std::vector<uint64_t> get_keys() const;
};

class offset_domain_t final {

    bool m_is_bottom = false;
    offset_registers_t m_registers;
    offset_stack_t m_stack;
    std::shared_ptr<offset_ctx_t> m_ctx;
    std::vector<std::string> m_errors;

  public:
    offset_domain_t() = default;
    offset_domain_t(offset_domain_t&& o) = default;
    offset_domain_t(const offset_domain_t& o) = default;
    offset_domain_t& operator=(offset_domain_t&& o) = default;
    offset_domain_t& operator=(const offset_domain_t& o) = default;
    explicit offset_domain_t(offset_registers_t&& reg, offset_stack_t&& stack,
            std::shared_ptr<offset_ctx_t> ctx)
        : m_registers(std::move(reg)), m_stack(std::move(stack)), m_ctx(ctx) {}

    static offset_domain_t&& setup_entry();
    // bottom/top
    static offset_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    bool is_bottom() const;
    bool is_top() const;
    // inclusion
    bool operator<=(const offset_domain_t& other) const;
    // join
    void operator|=(const offset_domain_t& abs);
    void operator|=(offset_domain_t&& abs);
    offset_domain_t operator|(const offset_domain_t& other) const;
    offset_domain_t operator|(offset_domain_t&& abs) const;
    // meet
    offset_domain_t operator&(const offset_domain_t& other) const;
    // widening
    offset_domain_t widen(const offset_domain_t& other, bool);
    // narrowing
    offset_domain_t narrow(const offset_domain_t& other) const;
    //forget
    void operator-=(register_t reg) { m_registers -= reg; }

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = location_t::top());
    void operator()(const Bin&, location_t loc = location_t::top());
    void operator()(const Un&, location_t loc = location_t::top());
    void operator()(const LoadMapFd&, location_t loc = location_t::top());
    void operator()(const Atomic&, location_t loc = location_t::top()) {}
    void operator()(const Call&, location_t loc = location_t::top());
    void operator()(const Exit&, location_t loc = location_t::top());
    void operator()(const Jmp&, location_t loc = location_t::top());
    void operator()(const Mem&, location_t loc = location_t::top());
    void operator()(const Packet&, location_t loc = location_t::top());
    void operator()(const Assume&, location_t loc = location_t::top());
    void operator()(const Assert&, location_t loc = location_t::top());
    void operator()(const IncrementLoopCounter&, location_t loc = location_t::top()) {};
    void operator()(const basic_block_t& bb);
    void write(std::ostream& os) const;
    std::string domain_name() const;
    crab::bound_t get_loop_count_upper_bound() const;
    void initialize_loop_counter(const label_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f) {}

    void do_un(const Un&, interval_t, location_t);
    void do_load(const Mem&, const register_t&, std::optional<ptr_or_mapfd_t>, interval_t&&,
            location_t);
    void do_mem_store(const Mem&, std::optional<ptr_or_mapfd_t>&);
    void do_bin(const Bin&, const std::optional<interval_t>&,
            const std::optional<ptr_or_mapfd_t>&,
            const std::optional<interval_t>&,
            const std::optional<ptr_or_mapfd_t>&, mock_interval_t&&, location_t);
    void do_call(const Call&, const stack_cells_t&, location_t);
    bool check_packet_access(const Reg&, int, int, bool) const;
    void check_valid_access(const ValidAccess&, std::optional<ptr_or_mapfd_t>&, int);
    interval_t compute_packet_subtraction(register_t, register_t) const;

    std::vector<uint64_t> get_ctx_keys() const;
    std::optional<refinement_t> find_in_ctx(int) const;
    std::optional<refinement_stack_cell_t> find_in_stack(int) const;
    std::optional<refinement_t> find_refinement_at_loc(const register_location_t) const;
    std::optional<refinement_t> find_refinement_info(register_t reg) const;
    void insert_in_registers(register_t, location_t, refinement_t);
    void store_in_stack(uint64_t, refinement_t, int);
    void adjust_bb_for_types(location_t);
    void print_ctx(std::ostream& o) const {}
    void print_stack(std::ostream& o) const {}
    void print_annotated_bb(std::ostream& o, const basic_block_t& bb) const {}
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void reset_errors() { m_errors.clear(); }
}; // end offset_domain_t

} // end namespace crab
