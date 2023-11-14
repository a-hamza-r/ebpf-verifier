// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>

#include <boost/optional/optional_io.hpp>
#include "crab/abstract_domain.hpp"
#include "crab/common.hpp"

namespace crab {

constexpr int PACKET_END = -4100;
constexpr int PACKET_META = -1;

using weight_t = interval_t;
using slack_var_t = int;

enum class rop_t {
    R_GT,
    R_GE,
    R_LT,
    R_LE
};

struct dist_t {
    slack_var_t m_slack;
    weight_t m_dist;
    bool m_is_bottom = false;

    dist_t(weight_t d, slack_var_t s = -1, bool bottom = false)
        : m_slack(s), m_dist(d), m_is_bottom(bottom) {}
    dist_t() : m_slack(-1), m_dist(weight_t::top()), m_is_bottom(false) {}
    bool operator==(const dist_t& d) const;
    void write(std::ostream&) const;
    bool is_top() const;
    bool is_bottom() const;
    void set_to_top();
    void set_to_bottom();
    friend std::ostream& operator<<(std::ostream& o, const dist_t& d);
    bool is_meta_pointer() const;
    bool is_forward_pointer() const;
    bool is_backward_pointer() const;
    weight_t offset_from_reference() const;
};

struct inequality_t {
    slack_var_t m_slack;
    rop_t m_rel;
    weight_t m_value;
    bool m_is_bottom = false;

    inequality_t(slack_var_t slack, rop_t rel, weight_t val) : m_slack(slack), m_rel(rel)
                                                               , m_value(val) {}
    inequality_t() : m_slack(-1), m_value(weight_t::top()) {}
    bool is_top() const;
    bool is_bottom() const;
    void set_to_top();
    void set_to_bottom();
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const inequality_t&);
};    // represents `slack rel value;`, e.g., `s >= 0`

struct equality_t {
    dist_t m_lhs;
    dist_t m_rhs;
    bool m_is_bottom = false;

    equality_t(dist_t lhs, dist_t rhs) : m_lhs(lhs), m_rhs(rhs) {}
    equality_t() = default;
    bool is_top() const;
    bool is_bottom() const;
    void set_to_top();
    void set_to_bottom();
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const equality_t&);
};  // represents constraint `p[0] = p[1];`, e.g., `begin+8+s = end`

struct packet_constraint_t {
    equality_t m_eq;
    inequality_t m_ineq;
    bool m_is_meta_constraint;
    bool m_is_bottom = false;

    bool is_bottom() const;
    bool is_top() const;
    void set_to_bottom();
    void set_to_top();
    std::optional<bound_t> get_limit() const;
    packet_constraint_t() = default;
    packet_constraint_t(equality_t&& eq, inequality_t&& ineq, bool is_meta_constraint,
            bool is_bottom = false) : m_eq(eq), m_ineq(ineq),
    m_is_meta_constraint(is_meta_constraint), m_is_bottom(is_bottom) {}
    packet_constraint_t operator|(const packet_constraint_t&) const;
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const packet_constraint_t&);
};

using live_registers_t = std::array<std::shared_ptr<reg_with_loc_t>, NUM_REGISTERS>;
using global_offset_env_t = std::unordered_map<reg_with_loc_t, dist_t>;

class registers_state_t {

    live_registers_t m_cur_def;
    std::shared_ptr<global_offset_env_t> m_offset_env;
    bool m_is_bottom = false;

    public:
        registers_state_t(bool is_bottom = false) : m_offset_env(nullptr), m_is_bottom(is_bottom) {}
        registers_state_t(std::shared_ptr<global_offset_env_t> offset_env, bool is_bottom = false) : m_offset_env(offset_env), m_is_bottom(is_bottom) {}
        explicit registers_state_t(live_registers_t&& vars, std::shared_ptr<global_offset_env_t>
                offset_env, bool is_bottom = false)
            : m_cur_def(std::move(vars)), m_offset_env(std::move(offset_env)), m_is_bottom(is_bottom) {}

        registers_state_t operator|(const registers_state_t&) const;
        void operator-=(register_t);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        void insert(register_t, const location_t&, dist_t&&);
        std::optional<dist_t> find(reg_with_loc_t reg) const;
        std::optional<dist_t> find(register_t key) const;
        friend std::ostream& operator<<(std::ostream& o, const registers_state_t& p);
        void adjust_bb_for_registers(location_t);
        void scratch_caller_saved_registers();
        void forget_packet_pointers();
};

using dist_cells_t = std::pair<dist_t, int>;
using stack_slot_dists_t = std::map<uint64_t, dist_cells_t>;    // represents `sp[n] = dist;`, where n \belongs [0,511], e.g., `sp[508] = begin+16`

class stack_state_t {

    stack_slot_dists_t m_slot_dists;
    bool m_is_bottom = false;

    public:
        stack_state_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        std::optional<dist_cells_t> find(uint64_t) const;
        void store(uint64_t, dist_t, int);
        void operator-=(uint64_t);
        void operator-=(const std::vector<uint64_t>&);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        static stack_state_t top();
        stack_state_t operator|(const stack_state_t&) const;
        explicit stack_state_t(stack_slot_dists_t&& stack_dists, bool is_bottom = false)
            : m_slot_dists(std::move(stack_dists)), m_is_bottom(is_bottom) {}
        std::vector<uint64_t> find_overlapping_cells(uint64_t, int) const;
};

class extra_constraints_t {

    packet_constraint_t m_meta_and_begin;
    packet_constraint_t m_begin_and_end;
    bool m_is_bottom = false;

    public:
        extra_constraints_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        bool is_bottom() const;
        bool is_top() const;
        void set_to_top();
        void set_to_bottom();
        void normalize();
        std::optional<bound_t> get_end_limit() const;
        std::optional<bound_t> get_meta_limit() const;
        void add_meta_and_begin_constraint(equality_t&&, inequality_t&&);
        void add_begin_and_end_constraint(equality_t&&, inequality_t&&);
        extra_constraints_t operator|(const extra_constraints_t&) const;
        explicit extra_constraints_t(packet_constraint_t&& meta_and_begin,
                packet_constraint_t&& begin_and_end, bool is_bottom = false)
            : m_meta_and_begin(meta_and_begin), m_begin_and_end(begin_and_end),
            m_is_bottom(is_bottom) {}
};

class ctx_offsets_t {
    using ctx_dists_t = std::unordered_map<int, dist_t>;    // represents `cp[n] = dist;`
    ctx_dists_t m_dists;
    int m_size;

    public:
        ctx_offsets_t(const ebpf_context_descriptor_t* desc);
        std::optional<dist_t> find(int) const;
        int get_size() const;
};

class offset_domain_t final {

    bool m_is_bottom = false;
    registers_state_t m_reg_state;
    stack_state_t m_stack_state;
    extra_constraints_t m_extra_constraints;
    std::shared_ptr<ctx_offsets_t> m_ctx_dists;
    std::vector<std::string> m_errors;
    slack_var_t m_slack = 0;

  public:
    offset_domain_t() = default;
    offset_domain_t(offset_domain_t&& o) = default;
    offset_domain_t(const offset_domain_t& o) = default;
    offset_domain_t& operator=(offset_domain_t&& o) = default;
    offset_domain_t& operator=(const offset_domain_t& o) = default;
    explicit offset_domain_t(registers_state_t&& reg, stack_state_t&& stack,
            extra_constraints_t extra, std::shared_ptr<ctx_offsets_t> ctx, slack_var_t s = 0)
        : m_reg_state(std::move(reg)), m_stack_state(std::move(stack)),
        m_extra_constraints(std::move(extra)), m_ctx_dists(ctx), m_slack(s) {}

    explicit offset_domain_t(registers_state_t&& reg, stack_state_t&& stack,
            std::shared_ptr<ctx_offsets_t> ctx, slack_var_t s = 0) : m_reg_state(std::move(reg)),
    m_stack_state(std::move(stack)), m_ctx_dists(ctx), m_slack(s) {}

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
    void operator-=(variable_t var);
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
    void operator()(const ValidAccess&, location_t loc = boost::none);
    void operator()(const Comparable&, location_t loc = boost::none);
    void operator()(const Addable&, location_t loc = boost::none);
    void operator()(const ValidStore&, location_t loc = boost::none);
    void operator()(const TypeConstraint&, location_t loc = boost::none);
    void operator()(const ValidSize&, location_t loc = boost::none);
    void operator()(const ValidMapKeyValue&, location_t loc = boost::none);
    void operator()(const ZeroCtxOffset&, location_t loc = boost::none);
    void operator()(const ValidDivisor&, location_t loc = boost::none);
    void operator()(const FuncConstraint& s, location_t loc = boost::none) {};
    void operator()(const IncrementLoopCounter&, location_t loc = boost::none) {};
    void operator()(const basic_block_t& bb, int print = 0);
    void write(std::ostream& os) const;
    std::string domain_name() const;
    crab::bound_t get_loop_count_upper_bound();
    void initialize_loop_counter(const label_t&);
    string_invariant to_set();
    void set_require_check(check_require_func_t f) {}

    void do_load(const Mem&, const register_t&, std::optional<ptr_or_mapfd_t>, location_t loc);
    void do_mem_store(const Mem&, std::optional<ptr_or_mapfd_t>, std::optional<ptr_or_mapfd_t>&);
    interval_t do_bin(const Bin&, const std::optional<interval_t>&,
            const std::optional<interval_t>&,
            std::optional<ptr_or_mapfd_t>&,
            std::optional<ptr_or_mapfd_t>&, location_t);
    void do_call(const Call&, const stack_cells_t&, location_t);
    bool upper_bound_satisfied(const dist_t&, int, int, bool) const;
    bool lower_bound_satisfied(const dist_t&, int) const;
    bool check_packet_access(const Reg&, int, int, bool) const;
    void check_valid_access(const ValidAccess&, std::optional<ptr_or_mapfd_t>&, int);

    std::optional<dist_t> find_in_ctx(int) const;
    std::optional<dist_cells_t> find_in_stack(int) const;
    std::optional<dist_t> find_offset_at_loc(const reg_with_loc_t) const;
    std::optional<dist_t> find_offset_info(register_t reg) const;
    void update_offset_info(const dist_t&&, const interval_t&&, const location_t&,
            uint8_t, Bin::Op);
    dist_t update_offset(const dist_t&, const weight_t&, const interval_t&, Bin::Op);
    void insert_in_registers(register_t, location_t, dist_t);
    void store_in_stack(uint64_t, dist_t, int);
    void adjust_bb_for_types(location_t);
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void reset_errors() { m_errors.clear(); }
}; // end offset_domain_t

} // end namespace crab
