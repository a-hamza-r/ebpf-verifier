// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>

#include <boost/optional/optional_io.hpp>
#include "crab/abstract_domain.hpp"
#include "crab/region_domain.hpp"
#include "crab/cfg.hpp"
#include "linear_constraint.hpp"
#include "string_constraints.hpp"

constexpr int STACK_BEGIN = 0;
constexpr int CTX_BEGIN = 0;
constexpr int PACKET_BEGIN = 0;

using crab::ptr_or_mapfd_t;
using crab::ptr_t;
using crab::mapfd_t;
using crab::ptr_with_off_t;
using crab::ptr_no_off_t;
using crab::reg_with_loc_t;
using crab::interval_t;
using crab::bound_t;

//using symbol_t = register_t;    // a register with unknown value
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

    dist_t(weight_t d, slack_var_t s = -1) : m_slack(s), m_dist(d) {}
    dist_t() : m_slack(-1), m_dist(weight_t::top()) {}
    bool operator==(const dist_t& d) const;
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream& o, const dist_t& d);
};      // if dist is +ve, represents `begin+dist+slack;`, if dist is -ve, represents `end+dist+1`

struct inequality_t {
    slack_var_t m_slack;
    rop_t m_rel;
    weight_t m_value;

    inequality_t(slack_var_t slack, rop_t rel, weight_t val) : m_slack(slack), m_rel(rel)
                                                               , m_value(val) {}
    inequality_t() : m_slack(-1), m_value(weight_t::top()) {}
};    // represents `slack rel value;`, e.g., `s >= 0`

struct forward_and_backward_eq_t {
    dist_t m_forw;
    dist_t m_backw;

    forward_and_backward_eq_t(dist_t forw, dist_t backw) : m_forw(forw), m_backw(backw) {}
    forward_and_backward_eq_t() = default;
};  // represents constraint `p[0] = p[1];`, e.g., `begin+8+s = end`

using live_registers_t = std::array<std::shared_ptr<reg_with_loc_t>, 11>;
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
        void insert(register_t, const reg_with_loc_t&, const dist_t&);
        std::optional<dist_t> find(reg_with_loc_t reg) const;
        std::optional<dist_t> find(register_t key) const;
        friend std::ostream& operator<<(std::ostream& o, const registers_state_t& p);
        void adjust_bb_for_registers(location_t);
};

class stack_state_t {
    using stack_slot_dists_t = std::unordered_map<int, dist_t>;    // represents `sp[n] = dist;`, where n \belongs [0,511], e.g., `sp[508] = begin+16`

    stack_slot_dists_t m_slot_dists;
    bool m_is_bottom = false;

    public:
        stack_state_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        std::optional<dist_t> find(int) const;
        void store(int, dist_t);
        void operator-=(int);
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        static stack_state_t top();
        stack_state_t operator|(const stack_state_t&) const;
        explicit stack_state_t(stack_slot_dists_t&& stack_dists, bool is_bottom = false)
            : m_slot_dists(std::move(stack_dists)), m_is_bottom(is_bottom) {}
};

class extra_constraints_t {

    forward_and_backward_eq_t m_eq;
    inequality_t m_ineq;
    bool m_is_bottom = false;

    public:
        extra_constraints_t(bool is_bottom = false) : m_is_bottom(is_bottom) {}
        void set_to_top();
        void set_to_bottom();
        bool is_bottom() const;
        bool is_top() const;
        void add_equality(forward_and_backward_eq_t);
        void add_inequality(inequality_t);
        void normalize();
        bound_t get_limit() const;
        extra_constraints_t operator|(const extra_constraints_t&) const;
        explicit extra_constraints_t(forward_and_backward_eq_t&& fabeq, inequality_t ineq, bool is_bottom = false) : m_eq(fabeq), m_ineq(ineq), m_is_bottom(is_bottom) {}
};

class ctx_t {
    using ctx_dists_t = std::unordered_map<int, dist_t>;    // represents `cp[n] = dist;`
    ctx_dists_t m_dists;
    int m_size;

    public:
        ctx_t(const ebpf_context_descriptor_t* desc);
        std::optional<dist_t> find(int) const;
        int get_size() const;
};


class offset_domain_t final {

    bool m_is_bottom = false;
    registers_state_t m_reg_state;
    stack_state_t m_stack_state;
    extra_constraints_t m_extra_constraints;
    std::shared_ptr<ctx_t> m_ctx_dists;
    slack_var_t m_slack = 0;

  public:
    offset_domain_t() = default;
    offset_domain_t(offset_domain_t&& o) = default;
    offset_domain_t(const offset_domain_t& o) = default;
    offset_domain_t& operator=(offset_domain_t&& o) = default;
    offset_domain_t& operator=(const offset_domain_t& o) = default;
    explicit offset_domain_t(registers_state_t&& reg, stack_state_t&& stack,
            extra_constraints_t extra, std::shared_ptr<ctx_t> ctx, slack_var_t s = 0)
        : m_reg_state(std::move(reg)), m_stack_state(std::move(stack)),
        m_extra_constraints(std::move(extra)), m_ctx_dists(ctx), m_slack(s) {}

    explicit offset_domain_t(registers_state_t&& reg, stack_state_t&& stack,
            std::shared_ptr<ctx_t> ctx, slack_var_t s = 0) : m_reg_state(std::move(reg)),
    m_stack_state(std::move(stack)), m_ctx_dists(ctx), m_slack(s) {}

    static offset_domain_t setup_entry();
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
    offset_domain_t widen(const offset_domain_t& other) const;
    // narrowing
    offset_domain_t narrow(const offset_domain_t& other) const;
    //forget
    void operator-=(crab::variable_t var);

    //// abstract transformers
    void operator()(const Undefined &, location_t loc = boost::none, int print = 0);
    void operator()(const Bin &, location_t loc = boost::none, int print = 0);
    void operator()(const Un &, location_t loc = boost::none, int print = 0);
    void operator()(const LoadMapFd &, location_t loc = boost::none, int print = 0);
    void operator()(const Call &, location_t loc = boost::none, int print = 0);
    void operator()(const Exit &, location_t loc = boost::none, int print = 0);
    void operator()(const Jmp &, location_t loc = boost::none, int print = 0);
    void operator()(const Mem &, location_t loc = boost::none, int print = 0);
    void operator()(const Packet &, location_t loc = boost::none, int print = 0);
    void operator()(const LockAdd &, location_t loc = boost::none, int print = 0);
    void operator()(const Assume &, location_t loc = boost::none, int print = 0);
    void operator()(const Assert &, location_t loc = boost::none, int print = 0);
    void operator()(const ValidAccess&, location_t loc = boost::none, int print = 0) {}
    void operator()(const Comparable& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const Addable& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ValidStore& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const TypeConstraint& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ValidSize& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ValidMapKeyValue& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ZeroCtxOffset& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ValidDivisor& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const basic_block_t& bb, bool check_termination, int print = 0);
    void write(std::ostream& os) const;
    std::string domain_name() const;
    crab::bound_t get_instruction_count_upper_bound();
    string_invariant to_set();
    void set_require_check(check_require_func_t f) {}

    void do_load(const Mem&, const Reg&, std::optional<ptr_or_mapfd_t>&, location_t loc);
    void do_mem_store(const Mem&, const Reg&, std::optional<ptr_or_mapfd_t>&,
            std::optional<ptr_or_mapfd_t>&);
    void do_bin(const Bin&, std::optional<interval_t>, std::optional<ptr_or_mapfd_t>,
            std::optional<ptr_or_mapfd_t>, location_t);
    void check_valid_access(const ValidAccess&, std::optional<ptr_or_mapfd_t>&);

    std::optional<dist_t> find_in_ctx(int) const;
    std::optional<dist_t> find_in_stack(int) const;
    std::optional<dist_t> find_in_registers(const reg_with_loc_t) const;
    std::optional<dist_t> find_offset_info(register_t reg) const;
    void adjust_bb_for_types(location_t);
}; // end offset_domain_t
