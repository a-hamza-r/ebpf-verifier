// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "crab/refinement.hpp"

namespace crab {

using live_refinements_t = std::array<std::shared_ptr<reg_with_loc_t>, NUM_REGISTERS+3>;
using global_offset_env_t = std::unordered_map<reg_with_loc_t, refinement_t>;

class registers_state_t {

    live_refinements_t m_cur_def;
    std::shared_ptr<global_offset_env_t> m_offset_env;
    std::shared_ptr<slacks_t> m_slacks;
    bool m_is_bottom = false;

    public:
        registers_state_t(bool is_bottom = false) : m_offset_env(nullptr), m_slacks(nullptr),
        m_is_bottom(is_bottom) {}
        registers_state_t(std::shared_ptr<global_offset_env_t> offset_env,
                std::shared_ptr<slacks_t> slacks, bool is_bottom = false)
            : m_offset_env(offset_env), m_slacks(slacks), m_is_bottom(is_bottom) {}
        registers_state_t(std::shared_ptr<global_offset_env_t> offset_env,
                std::shared_ptr<slacks_t> slacks, const ebpf_context_descriptor_t* desc,
                bool is_bottom = false)
            : m_offset_env(offset_env), m_slacks(slacks), m_is_bottom(is_bottom) {

            auto loc = std::make_pair(label_t::entry, (unsigned int)0);
            if (desc->data >= 0) {
                insert(register_t{11}, loc, refinement_t::begin());
            }
            if (desc->end >= 0) {
                insert(register_t{12}, loc, refinement_t::end());
            }
            if (desc->meta >= 0) {
                insert(register_t{13}, loc, refinement_t::meta());
            }
        }

        explicit registers_state_t(live_refinements_t&& vars,
                std::shared_ptr<global_offset_env_t> offset_env,
                std::shared_ptr<slacks_t> slacks, bool is_bottom = false)
            : m_cur_def(std::move(vars)), m_offset_env(offset_env), m_slacks(slacks),
            m_is_bottom(is_bottom) {}

        registers_state_t operator|(const registers_state_t&) const;
        void operator-=(register_t);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        void insert(register_t, const location_t&, refinement_t&&);
        void insert_slack_value(symbol_t, mock_interval_t);
        std::shared_ptr<slacks_t> get_slacks() const { return m_slacks; }
        std::optional<mock_interval_t> find_slack_value(symbol_t) const;
        std::optional<refinement_t> find(reg_with_loc_t reg) const;
        std::optional<refinement_t> find(register_t key) const;
        friend std::ostream& operator<<(std::ostream& o, const registers_state_t& p);
        void adjust_bb_for_registers(location_t);
        void scratch_caller_saved_registers();
        void forget_packet_pointers(location_t);
};

using refinement_cells_t = std::pair<refinement_t, int>;
using stack_slot_refinements_t = std::map<uint64_t, refinement_cells_t>;
class stack_state_t {

    stack_slot_refinements_t m_slot_rfs;
    bool m_is_bottom = false;

    public:
        stack_state_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        std::optional<refinement_cells_t> find(uint64_t) const;
        void store(uint64_t, refinement_t, int);
        void operator-=(uint64_t);
        void operator-=(const std::vector<uint64_t>&);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        static stack_state_t top();
        stack_state_t operator|(const stack_state_t&) const;
        explicit stack_state_t(stack_slot_refinements_t&& stack_rfs, bool is_bottom = false)
            : m_slot_rfs(std::move(stack_rfs)), m_is_bottom(is_bottom) {}
        std::vector<uint64_t> find_overlapping_cells(uint64_t, int) const;
};

class ctx_offsets_t {
    using ctx_refinements_t = std::unordered_map<int, refinement_t>;    // represents `cp[n] = rf;`
    ctx_refinements_t m_rfs;
    int m_size;

    public:
        ctx_offsets_t(const ebpf_context_descriptor_t* desc);
        std::optional<refinement_t> find(int) const;
        int get_size() const;
};

class offset_domain_t final {

    bool m_is_bottom = false;
    registers_state_t m_reg_state;
    stack_state_t m_stack_state;
    std::shared_ptr<ctx_offsets_t> m_ctx_rfs;
    std::vector<std::string> m_errors;

  public:
    offset_domain_t() = default;
    offset_domain_t(offset_domain_t&& o) = default;
    offset_domain_t(const offset_domain_t& o) = default;
    offset_domain_t& operator=(offset_domain_t&& o) = default;
    offset_domain_t& operator=(const offset_domain_t& o) = default;
    explicit offset_domain_t(registers_state_t&& reg, stack_state_t&& stack,
            std::shared_ptr<ctx_offsets_t> ctx)
        : m_reg_state(std::move(reg)), m_stack_state(std::move(stack)), m_ctx_rfs(ctx) {}

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
    void operator-=(register_t reg) { m_reg_state -= reg; }

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = boost::none);
    void operator()(const Bin&, location_t loc = boost::none);
    void operator()(const Un&, location_t loc = boost::none);
    void operator()(const LoadMapFd&, location_t loc = boost::none);
    void operator()(const Atomic&, location_t loc = boost::none) {}
    void operator()(const Call&, location_t loc = boost::none);
    void operator()(const Callx&, location_t loc = boost::none);
    void operator()(const Exit&, location_t loc = boost::none);
    void operator()(const Jmp&, location_t loc = boost::none);
    void operator()(const Mem&, location_t loc = boost::none);
    void operator()(const Packet&, location_t loc = boost::none);
    void operator()(const Assume&, location_t loc = boost::none);
    void operator()(const Assert&, location_t loc = boost::none);
    void operator()(const IncrementLoopCounter&, location_t loc = boost::none) {};
    void operator()(const basic_block_t& bb, int print = 0);
    void write(std::ostream& os) const;
    std::string domain_name() const;
    crab::bound_t get_loop_count_upper_bound();
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

    std::optional<refinement_t> find_in_ctx(int) const;
    std::optional<refinement_cells_t> find_in_stack(int) const;
    std::optional<refinement_t> find_refinement_at_loc(const reg_with_loc_t) const;
    std::optional<refinement_t> find_refinement_info(register_t reg) const;
    void insert_in_registers(register_t, location_t, refinement_t);
    void store_in_stack(uint64_t, refinement_t, int);
    void adjust_bb_for_types(location_t);
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void reset_errors() { m_errors.clear(); }
}; // end offset_domain_t

} // end namespace crab
