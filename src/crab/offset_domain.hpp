// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "array_domain.hpp"
#include "refinement.hpp"

namespace crab {

using check_require_func_t = std::function<bool(crab::domains::NumAbsDomain&, const crab::linear_constraint_t&, std::string)>;
using global_env_offset_registers_t = std::unordered_map<register_location_t, refinement_t>;

constexpr uint8_t BEGIN_REG = 12;

class offset_registers_t {

    live_registers_t m_cur_register_def;
    std::shared_ptr<global_env_offset_registers_t> m_registers_env;
    bool m_is_bottom = false;

    public:
        offset_registers_t() :
            m_registers_env(std::make_shared<global_env_offset_registers_t>()) {
            insert(BEGIN_REG, location_t{label_t::entry, 0},
                   refinement_t::begin_with_constraints());
        }
        offset_registers_t(std::shared_ptr<global_env_offset_registers_t> env)
            : m_registers_env(env) {}
        offset_registers_t(const offset_registers_t& other)
            : m_cur_register_def(other.m_cur_register_def), m_is_bottom(other.m_is_bottom) {
            if (other.m_registers_env) {
                m_registers_env = std::make_shared<global_env_offset_registers_t>(*other.m_registers_env);
            }

        }
        offset_registers_t& operator=(const offset_registers_t& other) {
            if (this != &other) {
                m_cur_register_def = other.m_cur_register_def;
                m_is_bottom = other.m_is_bottom;
                if (other.m_registers_env) {
                    m_registers_env = std::make_shared<global_env_offset_registers_t>(*other.m_registers_env);
                }
            }
            return *this;
        }
        offset_registers_t(offset_registers_t&& other) = default;
        offset_registers_t& operator=(offset_registers_t&& other) = default;
        offset_registers_t operator|(const offset_registers_t&) const;
        bool operator<=(const offset_registers_t&) const;
        offset_registers_t widen(const offset_registers_t&) const;
        void operator-=(register_t);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        void insert(register_t, const location_t&, refinement_t);
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
        offset_stack_t() = default;
        explicit offset_stack_t(refinement_stack_cells_t cells)
            : m_stack_cells(std::move(cells)) {}
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
        bool operator<=(const offset_stack_t&) const;
        offset_stack_t widen(const offset_stack_t&) const;
        std::vector<uint64_t> find_overlapping_cells(uint64_t, int) const;
        std::vector<uint64_t> get_keys() const;
};

class offset_ctx_t {
    using refinement_ctx_cells_t = std::unordered_map<uint64_t, refinement_t>;    // represents `cp[n] = rf;`
    refinement_ctx_cells_t m_ctx_cells;
    size_t m_size = 0;

    public:
        offset_ctx_t(const ebpf_context_descriptor_t* desc, std::shared_ptr<slacks_t>);
        std::optional<refinement_t> find(uint64_t) const;
        size_t get_size() const { return m_size; }
        std::vector<uint64_t> get_keys() const;
};

class offset_domain_t final {

    bool m_is_bottom = false;
    std::shared_ptr<slacks_t> m_slacks;
    offset_registers_t m_registers;
    offset_stack_t m_stack;
    std::shared_ptr<offset_ctx_t> m_ctx;
    std::vector<std::string> m_errors;

  public:
    offset_domain_t() = default;
    offset_domain_t(offset_registers_t reg, offset_stack_t stack,
            std::shared_ptr<offset_ctx_t> ctx, std::shared_ptr<slacks_t> slacks) :
        m_slacks(std::move(slacks)), m_registers(std::move(reg)), m_stack(std::move(stack)),
        m_ctx(ctx) {}
    offset_domain_t(std::shared_ptr<slacks_t> slacks) : m_slacks(slacks) {}

    static offset_domain_t setup_entry(std::shared_ptr<slacks_t>);
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

    void do_load(const Mem&, const register_t&, std::optional<ptr_or_mapfd_t>, location_t);
    void do_mem_store(const Mem&, std::optional<ptr_or_mapfd_t>&);
    void do_bin(const Bin&, const std::optional<refinement_t>&, const std::optional<refinement_t>&,
                location_t);
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
