// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>

#include "crab/abstract_domain.hpp"
#include "crab/region_domain.hpp"
#include "crab/cfg.hpp"
#include "linear_constraint.hpp"
#include "string_constraints.hpp"

using crab::ptr_t;
using crab::ptr_with_off_t;
using crab::ptr_no_off_t;

class type_domain_t final {
    region_domain_t m_region;
    bool m_is_bottom = false;

  public:

    type_domain_t() = default;
    type_domain_t(type_domain_t&& o) = default;
    type_domain_t(const type_domain_t& o) = default;
    explicit type_domain_t(region_domain_t&& reg, bool is_bottom = false) :
        m_region(reg), m_is_bottom(is_bottom) {}
    type_domain_t& operator=(type_domain_t&& o) = default;
    type_domain_t& operator=(const type_domain_t& o) = default;
    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static type_domain_t setup_entry();
    // bottom/top
    static type_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    bool is_bottom() const;
    bool is_top() const;
    // inclusion
    bool operator<=(const type_domain_t& other) const;
    // join
    void operator|=(const type_domain_t& abs);
    void operator|=(type_domain_t&& abs);
    type_domain_t operator|(const type_domain_t& other) const;
    type_domain_t operator|(type_domain_t&& abs) const;
    // meet
    type_domain_t operator&(const type_domain_t& other) const;
    // widening
    type_domain_t widen(const type_domain_t& other, bool);
    // narrowing
    type_domain_t narrow(const type_domain_t& other) const;
    //forget
    void operator-=(crab::variable_t var);

    //// abstract transformers
    void operator()(const Undefined &, location_t loc = boost::none, int print = 0);
    void operator()(const Bin &, location_t loc = boost::none, int print = 0);
    void operator()(const Un &, location_t loc = boost::none, int print = 0);
    void operator()(const LoadMapFd &, location_t loc = boost::none, int print = 0);
    void operator()(const Atomic&, location_t loc = boost::none, int print = 0);
    void operator()(const Call &, location_t loc = boost::none, int print = 0);
    void operator()(const Callx&, location_t loc = boost::none, int print = 0);
    void operator()(const Exit &, location_t loc = boost::none, int print = 0);
    void operator()(const Jmp &, location_t loc = boost::none, int print = 0);
    void operator()(const Mem &, location_t loc = boost::none, int print = 0);
    void operator()(const Packet &, location_t loc = boost::none, int print = 0);
    void operator()(const Assume &, location_t loc = boost::none, int print = 0);
    void operator()(const Assert &, location_t loc = boost::none, int print = 0);
    void operator()(const ValidAccess&, location_t loc = boost::none, int print = 0);
    void operator()(const Comparable& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const Addable& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ValidStore& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const TypeConstraint& s, location_t loc = boost::none, int print = 0);
    void operator()(const ValidSize& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ValidMapKeyValue& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ZeroCtxOffset& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const ValidDivisor& s, location_t loc = boost::none, int print = 0) {}
    void operator()(const FuncConstraint& s, location_t loc = boost::none, int print = 0) {};
    void operator()(const IncrementLoopCounter&, location_t loc = boost::none, int print = 0);
    void operator()(const basic_block_t& bb, int print = 0);
    void write(std::ostream& os) const;
    void initialize_loop_counter(label_t label);
    crab::bound_t get_loop_count_upper_bound();
    string_invariant to_set();
    void set_require_check(check_require_func_t f) {}

  private:

    void do_load(const Mem&, const Reg&, location_t, int print = 0);
    void do_mem_store(const Mem&, const Reg&, location_t, int print = 0);
    void print_initial_types() const;
    void report_type_error(std::string, location_t);
    void print_ctx() const;
    void print_stack() const;
    void print_initial_registers() const;

}; // end type_domain_t
