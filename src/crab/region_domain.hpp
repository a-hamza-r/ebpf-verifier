// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "array_domain.hpp"
#include "types.hpp"
#include "platform.hpp"

namespace crab {

using check_require_func_t = std::function<bool(crab::domains::NumAbsDomain&, const crab::linear_constraint_t&, std::string)>;
using shared_ptr_aliases_t = std::vector<std::set<int>>;
using ptr_or_mapfd_stack_cell_t = std::pair<ptr_or_mapfd_t, int>;

class region_ctx_t {
    std::vector<uint64_t> m_keys;
    std::size_t m_size = 0;
    bool m_is_bottom = false;

  public:
    region_ctx_t() = default;
    region_ctx_t(const ebpf_context_descriptor_t*);
    [[nodiscard]] const std::vector<uint64_t>& get_keys() const { return m_keys; }
    bool packet_ptr_at(uint64_t key) const;
    size_t get_size() const { return m_size; }
};

class region_stack_t {
    using ptr_or_mapfd_stack_cells_t = std::map<uint64_t, ptr_or_mapfd_stack_cell_t>;
    ptr_or_mapfd_stack_cells_t m_cells;
    bool m_is_bottom = false;

  public:
    region_stack_t() = default;
    region_stack_t operator|(const region_stack_t& other) const;
    bool operator<=(const region_stack_t& other) const;
    region_stack_t widen(const region_stack_t& other) const;
    void operator-=(uint64_t);
    void operator-=(const std::vector<uint64_t>&);
    void set_to_bottom();
    void set_to_top();
    static region_stack_t bottom();
    static region_stack_t top();
    bool is_bottom() const;
    bool is_top() const;
    const ptr_or_mapfd_stack_cells_t &get_cells() { return m_cells; }
    void store(uint64_t, ptr_or_mapfd_t, int);
    std::optional<ptr_or_mapfd_stack_cell_t> find(uint64_t) const;
    [[nodiscard]] std::vector<uint64_t> get_keys() const;
    std::vector<uint64_t> find_overlapping_cells(uint64_t, int) const;
    size_t size() const;
};

using global_env_region_registers_t = std::unordered_map<register_location_t, ptr_or_mapfd_t>;

class region_registers_t {

    live_registers_t m_cur_register_def;
    std::shared_ptr<global_env_region_registers_t> m_registers_env;
    bool m_is_bottom = false;

  public:
    region_registers_t() : m_registers_env(std::make_shared<global_env_region_registers_t>()) {}
    region_registers_t(const region_registers_t& other)
        : m_cur_register_def(other.m_cur_register_def), m_is_bottom(other.m_is_bottom) {
        if (other.m_registers_env) {
            m_registers_env = std::make_shared<global_env_region_registers_t>(*other.m_registers_env);
        }
    }
    region_registers_t(region_registers_t&& other) noexcept = default;
    region_registers_t& operator=(const region_registers_t& other) {
        if (this != &other) {
            m_cur_register_def = other.m_cur_register_def;
            if (other.m_registers_env) {
                m_registers_env = std::make_shared<global_env_region_registers_t>(*other.m_registers_env);
            }
            m_is_bottom = other.m_is_bottom;
        }
        return *this;
    }
    region_registers_t& operator=(region_registers_t&& other) noexcept = default;
    region_registers_t operator|(const region_registers_t& other) const;
    region_registers_t widen(const region_registers_t& other) const;
    bool operator<=(const region_registers_t& other) const;
    void operator-=(register_t var);
    void set_to_bottom();
    void set_to_top();
    bool is_bottom() const;
    bool is_top() const;
    void insert(register_t, const location_t&, const ptr_or_mapfd_t&);
    std::optional<ptr_or_mapfd_t> find(register_location_t reg) const;
    std::optional<ptr_or_mapfd_t> find(register_t key) const;
    [[nodiscard]] live_registers_t &get_vars() { return m_cur_register_def; }
    void forget_packet_ptrs();
    void scratch_caller_saved_registers();
    void adjust_bb_for_registers(location_t loc);
};

class region_domain_t final {

    bool m_is_bottom = false;
    region_stack_t m_stack;
    region_registers_t m_registers;
    std::shared_ptr<region_ctx_t> m_ctx;
    shared_ptr_aliases_t m_shared_ptr_aliases;
    std::vector<std::string> m_errors;

  public:

    region_domain_t() = default;
    region_domain_t(region_registers_t registers, region_stack_t stack,
            std::shared_ptr<region_ctx_t> ctx, shared_ptr_aliases_t shared_ptr_aliases = {})
            : m_stack(std::move(stack)), m_registers(std::move(registers)), m_ctx(ctx),
            m_shared_ptr_aliases(std::move(shared_ptr_aliases)) {}

    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static region_domain_t setup_entry(bool);
    // bottom/top
    static region_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    bool is_bottom() const;
    bool is_top() const;
    // inclusion
    bool operator<=(const region_domain_t& other) const;
    // join
    void operator|=(const region_domain_t& abs);
    void operator|=(region_domain_t&& abs);
    region_domain_t operator|(const region_domain_t& other) const;
    region_domain_t operator|(region_domain_t&& abs) const;
    // meet
    region_domain_t operator&(const region_domain_t& other) const;
    // widening
    region_domain_t widen(const region_domain_t& other, bool);
    // narrowing
    region_domain_t narrow(const region_domain_t& other) const;
    //forget
    void operator-=(register_t var) { m_registers -= var; }

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = location_t::top());
    void operator()(const Bin&, location_t loc = location_t::top());
    void operator()(const Un&, location_t loc = location_t::top());
    void operator()(const LoadMapFd&, location_t loc = location_t::top());
    void operator()(const Atomic&, location_t loc = location_t::top());
    void operator()(const Call&, location_t loc = location_t::top());
    void operator()(const Exit&, location_t loc = location_t::top());
    void operator()(const Jmp&, location_t loc = location_t::top());
    void operator()(const Mem&, location_t loc = location_t::top());
    void operator()(const Packet&, location_t loc = location_t::top());
    void operator()(const Assume&, location_t loc = location_t::top());
    void operator()(const Assert&, location_t loc = location_t::top());
    void operator()(const ValidAccess&, location_t loc = location_t::top());
    void operator()(const TypeConstraint&, location_t loc = location_t::top());
    void operator()(const ZeroCtxOffset&, location_t loc = location_t::top());
    void operator()(const IncrementLoopCounter&, location_t loc = location_t::top());
    void operator()(const basic_block_t& bb);
    void write(std::ostream& o) const {}
    crab::bound_t get_loop_count_upper_bound() const;
    void initialize_loop_counter(const label_t&);
    friend std::ostream& operator<<(std::ostream&, const region_domain_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f) {}

    interval_t get_map_value_size(const Reg&) const;
    bool get_map_fd_range(const Reg&, int32_t*, int32_t*) const;
    interval_t get_map_key_size(const Reg&) const;
    std::optional<uint32_t> get_map_type(const Reg&) const;
    std::optional<uint32_t> get_map_inner_map_fd(const Reg&) const;
    void check_type(const TypeConstraint&, bool);
    void do_load_mapfd(const register_t&, int, location_t);
    void do_load(const Mem&, const register_t&, bool, location_t);
    void do_mem_store(const Mem&);
    void do_bin(const Bin&, const std::optional<interval_t>&, const std::optional<interval_t>&,
                location_t);
    void do_call(const Call&, const stack_cells_t&, location_t);
    void check_valid_access(const ValidAccess &, int);
    void assume_cst(Condition::Op, ptr_with_off_t&&, int64_t, register_t, location_t);
    void update_ptr_or_mapfd(const ptr_or_mapfd_t&, const interval_t&, location_t, register_t);

    std::optional<crab::ptr_or_mapfd_t> find_ptr_or_mapfd_type(register_t) const;
    [[nodiscard]] size_t get_ctx_size() const;
    std::optional<crab::packet_ptr_t> find_in_ctx(uint64_t key) const;
    [[nodiscard]] const std::vector<uint64_t>& get_ctx_keys() const;
    std::optional<crab::ptr_or_mapfd_stack_cell_t> find_in_stack(uint64_t key) const;
    std::optional<crab::ptr_or_mapfd_t> find_ptr_or_mapfd_at_loc(const crab::register_location_t&) const;
    void insert_in_registers(register_t, location_t, const ptr_or_mapfd_t&);
    void store_in_stack(uint64_t, ptr_or_mapfd_t, int);
    void set_aliases(int, ptr_with_off_t&);
    [[nodiscard]] std::vector<uint64_t> get_stack_keys() const;
    void set_registers_to_top();
    void adjust_bb_for_types(location_t);
    void print_ctx(std::ostream& o) const {}
    void print_stack(std::ostream& o) const {}
    void print_annotated_bb(std::ostream& o, const basic_block_t& bb) const {};
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void reset_errors() { m_errors.clear(); }
}; // end region_domain_t

} // namespace crab
