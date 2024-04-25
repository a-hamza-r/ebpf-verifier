// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "crab/region_domain.hpp"
#include "crab/interval_domain.hpp"
#include "crab/offset_domain.hpp"
#include "crab/common.hpp"
#include "crab/type_ostream.hpp"

namespace crab {

class type_domain_t final {
    region_domain_t m_region;
    offset_domain_t m_offset;
    interval_domain_t m_interval;
    bool m_is_bottom = false;
    std::vector<std::string> m_errors;

  public:

    type_domain_t() = default;
    type_domain_t(type_domain_t&& o) = default;
    type_domain_t(const type_domain_t& o) = default;
    explicit type_domain_t(region_domain_t&& reg, offset_domain_t&& off,
            interval_domain_t&& interval, bool is_bottom = false) :
        m_region(reg), m_offset(off), m_interval(interval), m_is_bottom(is_bottom) {}
    type_domain_t& operator=(type_domain_t&& o) = default;
    type_domain_t& operator=(const type_domain_t& o) = default;
    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static type_domain_t setup_entry(bool);
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

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = boost::none);
    void operator()(const Bin&, location_t loc = boost::none);
    void operator()(const Un&, location_t loc = boost::none);
    void operator()(const LoadMapFd&, location_t loc = boost::none);
    void operator()(const Atomic&, location_t loc = boost::none);
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
    void operator()(const FuncConstraint& s, location_t loc = boost::none);
    void operator()(const IncrementLoopCounter&, location_t loc = boost::none);
    void operator()(const basic_block_t& bb, int print = 0);
    void write(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& o, const type_domain_t& dom);
    void initialize_loop_counter(label_t label);
    crab::bound_t get_loop_count_upper_bound();
    string_invariant to_set() const;
    void set_require_check(check_require_func_t f) {}
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void print_ctx() const;
    void print_stack() const;
    std::optional<crab::ptr_or_mapfd_t> find_ptr_or_mapfd_at_loc(const crab::reg_with_loc_t&) const;
    std::optional<crab::refinement_t> find_refinement_at_loc(const crab::reg_with_loc_t&) const;
    std::optional<crab::mock_interval_t> find_signed_interval_at_loc(const crab::reg_with_loc_t&) const;
    std::optional<crab::mock_interval_t> find_unsigned_interval_at_loc(const crab::reg_with_loc_t&) const;
    static type_domain_t from_predefined_types(const std::set<std::string>&, bool);
    void insert_in_registers_in_region_domain(register_t, location_t, const ptr_or_mapfd_t&);
    void store_in_stack_in_region_domain(uint64_t, ptr_or_mapfd_t, int);
    void insert_in_registers_in_interval_domain(register_t, location_t, interval_t);
    void insert_in_registers_in_signed_interval_domain(register_t, location_t, interval_t);
    void insert_in_registers_in_unsigned_interval_domain(register_t, location_t, interval_t);
    void store_in_stack_in_interval_domain(uint64_t, mock_interval_t, int);
    void store_in_stack_in_signed_interval_domain(uint64_t, mock_interval_t, int);
    void store_in_stack_in_unsigned_interval_domain(uint64_t, mock_interval_t, int);
    void insert_in_registers_in_offset_domain(register_t, location_t, refinement_t);
    void store_in_stack_in_offset_domain(uint64_t, refinement_t, int);

  private:

    void do_load(const Mem&, const Reg&, bool, std::optional<ptr_or_mapfd_t>,
            location_t);
    void do_mem_store(const Mem&, std::optional<ptr_or_mapfd_t>&);
    void report_type_error(std::string, location_t);
    void print_registers() const;
    void adjust_bb_for_types(location_t);
    void operator+=(std::vector<std::string>& errs) {
        m_errors.insert(m_errors.end(), errs.begin(), errs.end());
    }
}; // end type_domain_t

} // namespace crab

void print_annotated(std::ostream&, const crab::type_domain_t&, const basic_block_t&, int);
