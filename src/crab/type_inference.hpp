// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "region_domain.hpp"
#include "interval_domain.hpp"
#include "offset_domain.hpp"
#include "type_ostream.hpp"

namespace crab {

class inference_domain_t final {
    std::shared_ptr<slacks_t> m_slacks;
    region_domain_t m_region;
    offset_domain_t m_offset;
    interval_domain_t m_interval;
    std::vector<std::string> m_errors;

  public:

    inference_domain_t() : m_slacks(std::make_shared<slacks_t>()),
        m_region(), m_offset(m_slacks), m_interval(m_slacks) {}
    inference_domain_t(region_domain_t region, offset_domain_t offset,
            interval_domain_t interval, std::shared_ptr<slacks_t> slacks) :
        m_slacks(std::move(slacks)), m_region(std::move(region)), m_offset(std::move(offset)),
        m_interval(std::move(interval)) {}
    inference_domain_t(const inference_domain_t& other) :
        m_slacks(other.m_slacks), m_region(other.m_region), m_offset(other.m_offset),
        m_interval(other.m_interval) {}
    inference_domain_t(inference_domain_t&& other) :
        m_slacks(std::move(other.m_slacks)), m_region(std::move(other.m_region)),
        m_offset(std::move(other.m_offset)), m_interval(std::move(other.m_interval)) {}
    inference_domain_t& operator=(const inference_domain_t& other) {
        if (this != &other) {
            m_slacks = other.m_slacks;
            m_region = other.m_region;
            m_offset = other.m_offset;
            m_interval = other.m_interval;
        }
        return *this;
    }
    inference_domain_t& operator=(inference_domain_t&& other) {
        if (this != &other) {
            m_slacks = std::move(other.m_slacks);
            m_region = std::move(other.m_region);
            m_offset = std::move(other.m_offset);
            m_interval = std::move(other.m_interval);
        }
        return *this;
    }

    // eBPF initialization: R1 points to ctx, R10 to stack, etc.
    static inference_domain_t setup_entry(bool);
    // bottom/top
    static inference_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    bool is_bottom() const;
    bool is_top() const;
    // inclusion
    bool operator<=(const inference_domain_t& other) const;
    // join
    void operator|=(const inference_domain_t& abs);
    void operator|=(inference_domain_t&& abs);
    inference_domain_t operator|(const inference_domain_t& other) const;
    inference_domain_t operator|(inference_domain_t&& abs) const;
    // meet
    inference_domain_t operator&(const inference_domain_t& other) const;
    // widening
    inference_domain_t widen(const inference_domain_t& other, bool);
    // narrowing
    inference_domain_t narrow(const inference_domain_t& other) const;

    //// abstract transformers
    void operator()(const Undefined&, location_t loc = location_t::top());
    void operator()(const Bin&, location_t loc = location_t::top());
    void operator()(const Un&, location_t loc = location_t::top());
    void operator()(const LoadMapFd&, location_t loc = location_t::top());
    void operator()(const Atomic&, location_t loc = location_t::top());
    void operator()(const Call&, location_t loc = location_t::top());
    void operator()(const CallLocal&, location_t loc = location_t::top()) {}
    void operator()(const Callx&, location_t loc = location_t::top());
    void operator()(const Exit&, location_t loc = location_t::top());
    void operator()(const Jmp&, location_t loc = location_t::top());
    void operator()(const Mem&, location_t loc = location_t::top());
    void operator()(const Packet&, location_t loc = location_t::top());
    void operator()(const Assume&, location_t loc = location_t::top());
    void operator()(const Assert&, location_t loc = location_t::top());
    void operator()(const ValidAccess&, location_t loc = location_t::top());
    void operator()(const Comparable&, location_t loc = location_t::top());
    void operator()(const Addable&, location_t loc = location_t::top());
    void operator()(const ValidStore&, location_t loc = location_t::top());
    void operator()(const TypeConstraint&, location_t loc = location_t::top());
    void operator()(const ValidSize&, location_t loc = location_t::top());
    void operator()(const ValidCall&, location_t loc = location_t::top()) {}
    void operator()(const ValidMapKeyValue&, location_t loc = location_t::top());
    void operator()(const ZeroCtxOffset&, location_t loc = location_t::top());
    void operator()(const ValidDivisor&, location_t loc = location_t::top());
    void operator()(const FuncConstraint& s, location_t loc = location_t::top());
    void operator()(const IncrementLoopCounter&, location_t loc = location_t::top());
    void operator()(const basic_block_t& bb);
    void write(std::ostream& os) const;
    friend std::ostream& operator<<(std::ostream& o, const inference_domain_t& dom);
    void initialize_loop_counter(label_t label);
    crab::bound_t get_loop_count_upper_bound() const;
    string_invariant to_set() const;
    void set_require_check(check_require_func_t f) {}
    [[nodiscard]] std::vector<std::string>& get_errors() { return m_errors; }
    void print_ctx(std::ostream&) const;
    void print_stack(std::ostream&) const;
    void print_annotated_bb(std::ostream&, const basic_block_t&) const;
    std::optional<crab::ptr_or_mapfd_t> find_ptr_or_mapfd_at_loc(const crab::register_location_t&) const;
    std::optional<crab::refinement_t> find_refinement_at_loc(const crab::register_location_t&) const;
    std::optional<crab::refinement_t> find_signed_interval_at_loc(const crab::register_location_t&) const;
    std::optional<crab::refinement_t> find_unsigned_interval_at_loc(const crab::register_location_t&) const;
    static inference_domain_t from_predefined_types(const std::set<std::string>&, bool);
    void insert_in_registers_in_region_domain(register_t, location_t, const ptr_or_mapfd_t&);
    void store_in_stack_in_region_domain(uint64_t, ptr_or_mapfd_t, int);
    void insert_in_registers_in_interval_domain(register_t, location_t, refinement_t);
    void insert_in_registers_in_signed_interval_domain(register_t, location_t, refinement_t);
    void insert_in_registers_in_unsigned_interval_domain(register_t, location_t, refinement_t);
    void store_in_stack_in_interval_domain(uint64_t, refinement_t, int);
    void store_in_stack_in_signed_interval_domain(uint64_t, refinement_t, int);
    void store_in_stack_in_unsigned_interval_domain(uint64_t, refinement_t, int);
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
}; // end inference_domain_t

} // namespace crab

