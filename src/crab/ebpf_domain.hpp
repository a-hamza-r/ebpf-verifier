// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

// This file is eBPF-specific, not derived from CRAB.

#include <functional>
#include <optional>

#include "crab/abstract_domain.hpp"
#include "crab/array_domain.hpp"
#include "crab/split_dbm.hpp"
#include "crab/type_domain.hpp"
#include "crab/variable.hpp"
#include "string_constraints.hpp"

namespace crab {

class ebpf_domain_t final {
  using location_t = boost::optional<std::pair<label_t, uint32_t>>;
  public:
    ebpf_domain_t();
    ebpf_domain_t(NumAbsDomain inv, domains::array_domain_t stack);

    // Generic abstract domain operations
    static ebpf_domain_t top();
    static ebpf_domain_t bottom();
    void set_to_top();
    void set_to_bottom();
    [[nodiscard]]
    bool is_bottom() const;
    [[nodiscard]]
    bool is_top() const;
    bool operator<=(const ebpf_domain_t& other);
    //bool operator==(const ebpf_domain_t& other) const;
    void operator|=(ebpf_domain_t&& other);
    void operator|=(const ebpf_domain_t& other);
    ebpf_domain_t operator|(ebpf_domain_t&& other) const;
    ebpf_domain_t operator|(const ebpf_domain_t& other) const;
    ebpf_domain_t operator&(const ebpf_domain_t& other) const;
    ebpf_domain_t widen(const ebpf_domain_t& other, bool to_constants) const;
    ebpf_domain_t widening_thresholds(const ebpf_domain_t& other, const thresholds_t& ts);
    ebpf_domain_t narrow(const ebpf_domain_t& other) const;

    typedef bool check_require_func_t(NumAbsDomain&, const linear_constraint_t&, std::string);
    void set_require_check(std::function<check_require_func_t> f);
    extended_number get_loop_count_upper_bound() const;
    static ebpf_domain_t setup_entry(bool init_r1);
    std::vector<std::string> get_errors() { return std::vector<std::string>(); }

    static ebpf_domain_t from_constraints(const std::set<std::string>& constraints, bool setup_constraints);
    string_invariant to_set() const;

    // abstract transformers
    void operator()(const basic_block_t& bb);

    void operator()(const Addable&, location_t loc = boost::none);
    void operator()(const Assert&, location_t loc = boost::none);
    void operator()(const Assume&, location_t loc = boost::none);
    void operator()(const Bin&, location_t loc = boost::none);
    void operator()(const Call&, location_t loc = boost::none);
    void operator()(const CallLocal&, location_t loc = boost::none);
    void operator()(const Callx&, location_t loc = boost::none);
    void operator()(const Comparable&, location_t loc = boost::none);
    void operator()(const Exit&, location_t loc = boost::none);
    void operator()(const FuncConstraint&, location_t loc = boost::none);
    void operator()(const Jmp&, location_t loc = boost::none);
    void operator()(const LoadMapFd&, location_t loc = boost::none);
    void operator()(const Atomic&, location_t loc = boost::none);
    void operator()(const Mem&, location_t loc = boost::none);
    void operator()(const ValidDivisor&, location_t loc = boost::none);
    void operator()(const Packet&, location_t loc = boost::none);
    void operator()(const TypeConstraint&, location_t loc = boost::none);
    void operator()(const Un&, location_t loc = boost::none);
    void operator()(const Undefined&, location_t loc = boost::none);
    void operator()(const ValidAccess&, location_t loc = boost::none);
    void operator()(const ValidCall&, location_t loc = boost::none);
    void operator()(const ValidMapKeyValue&, location_t loc = boost::none);
    void operator()(const ValidSize&, location_t loc = boost::none);
    void operator()(const ValidStore&, location_t loc = boost::none);
    void operator()(const ZeroCtxOffset&, location_t loc = boost::none);
    void operator()(const IncrementLoopCounter&, location_t loc = boost::none);

    // following operations are important to keep in ebpf_domain_t because of the parametric abstract domain
    void write(std::ostream& o) const;
    void print_ctx(const std::ostream& o) const {}
    void print_stack(const std::ostream& o) const {}
    void print_annotated_bb(std::ostream& o, const basic_block_t& bb) const {}

    void initialize_loop_counter(const label_t& label);
    static ebpf_domain_t calculate_constant_limits();

  private:
    // private generic domain functions
    void operator+=(const linear_constraint_t& cst);
    void operator-=(variable_t var);

    void assign(variable_t lhs, variable_t rhs);
    void assign(variable_t x, const linear_expression_t& e);
    void assign(variable_t x, int64_t e);

    void apply(arith_binop_t op, variable_t x, variable_t y, const number_t& z, int finite_width);
    void apply(arith_binop_t op, variable_t x, variable_t y, variable_t z, int finite_width);
    void apply(bitwise_binop_t op, variable_t x, variable_t y, variable_t z, int finite_width);
    void apply(bitwise_binop_t op, variable_t x, variable_t y, const number_t& k, int finite_width);
    void apply(binop_t op, variable_t x, variable_t y, const number_t& z, int finite_width);
    void apply(binop_t op, variable_t x, variable_t y, variable_t z, int finite_width);

    void add(const Reg& reg, int imm, int finite_width);
    void add(variable_t lhs, variable_t op2);
    void add(variable_t lhs, const number_t& op2);
    void sub(variable_t lhs, variable_t op2);
    void sub(variable_t lhs, const number_t& op2);
    void add_overflow(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void add_overflow(variable_t lhss, variable_t lhsu, const number_t& op2, int finite_width);
    void sub_overflow(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void sub_overflow(variable_t lhss, variable_t lhsu, const number_t& op2, int finite_width);
    void neg(variable_t lhss, variable_t lhsu, int finite_width);
    void mul(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void mul(variable_t lhss, variable_t lhsu, const number_t& op2, int finite_width);
    void sdiv(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void sdiv(variable_t lhss, variable_t lhsu, const number_t& op2, int finite_width);
    void udiv(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void udiv(variable_t lhss, variable_t lhsu, const number_t& op2, int finite_width);
    void srem(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void srem(variable_t lhss, variable_t lhsu, const number_t& op2, int finite_width);
    void urem(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void urem(variable_t lhss, variable_t lhsu, const number_t& op2, int finite_width);

    void bitwise_and(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void bitwise_and(variable_t lhss, variable_t lhsu, const number_t& op2);
    void bitwise_or(variable_t lhss, variable_t lhsu, variable_t op2, int finite_width);
    void bitwise_or(variable_t lhss, variable_t lhsu, const number_t& op2);
    void bitwise_xor(variable_t lhsss, variable_t lhsu, variable_t op2, int finite_width);
    void bitwise_xor(variable_t lhss, variable_t lhsu, const number_t& op2);
    void shl(const Reg& reg, int imm, int finite_width);
    void shl_overflow(variable_t lhss, variable_t lhsu, variable_t op2);
    void shl_overflow(variable_t lhss, variable_t lhsu, const number_t& op2);
    void lshr(const Reg& reg, int imm, int finite_width);
    void ashr(const Reg& dst_reg, const linear_expression_t& right_svalue, int finite_width);
    void sign_extend(const Reg& dst_reg, const linear_expression_t& right_svalue, int finite_width, Bin::Op op);

    void assume(const linear_constraint_t& cst);

    /// Forget everything we know about the value of a variable.
    void havoc(variable_t v);

    /// Forget everything about all offset variables for a given register.
    void havoc_offsets(const Reg& reg);

    static std::optional<variable_t> get_type_offset_variable(const Reg& reg, int type);
    [[nodiscard]]
    std::optional<variable_t> get_type_offset_variable(const Reg& reg, const NumAbsDomain& inv) const;
    [[nodiscard]]
    std::optional<variable_t> get_type_offset_variable(const Reg& reg) const;

    void scratch_caller_saved_registers();
    void save_callee_saved_registers(const std::string& prefix);
    void restore_callee_saved_registers(const std::string& prefix);
    [[nodiscard]]
    std::optional<uint32_t> get_map_type(const Reg& map_fd_reg) const;
    [[nodiscard]]
    std::optional<uint32_t> get_map_inner_map_fd(const Reg& map_fd_reg) const;
    [[nodiscard]]
    interval_t get_map_key_size(const Reg& map_fd_reg) const;
    [[nodiscard]]
    interval_t get_map_value_size(const Reg& map_fd_reg) const;
    [[nodiscard]]
    interval_t get_map_max_entries(const Reg& map_fd_reg) const;
    void forget_packet_pointers();
    void do_load_mapfd(const Reg& dst_reg, int mapfd, bool maybe_null);

    void assign_valid_ptr(const Reg& dst_reg, bool maybe_null);

    void require(domains::NumAbsDomain& inv, const linear_constraint_t& cst, const std::string& s) const;

    // memory check / load / store
    void check_access_stack(NumAbsDomain& inv, const linear_expression_t& lb, const linear_expression_t& ub) const;
    void check_access_context(NumAbsDomain& inv, const linear_expression_t& lb, const linear_expression_t& ub) const;
    void check_access_packet(NumAbsDomain& inv, const linear_expression_t& lb, const linear_expression_t& ub,
                             std::optional<variable_t> packet_size) const;
    void check_access_shared(NumAbsDomain& inv, const linear_expression_t& lb, const linear_expression_t& ub,
                             variable_t shared_region_size) const;

    void recompute_stack_numeric_size(NumAbsDomain& inv, const Reg& reg) const;
    void recompute_stack_numeric_size(NumAbsDomain& inv, variable_t type_variable) const;
    void do_load_stack(NumAbsDomain& inv, const Reg& target_reg, const linear_expression_t& addr, int width,
                       const Reg& src_reg);
    void do_load_ctx(NumAbsDomain& inv, const Reg& target_reg, const linear_expression_t& addr_vague, int width);
    void do_load_packet_or_shared(NumAbsDomain& inv, const Reg& target_reg, const linear_expression_t& addr, int width);
    void do_load(const Mem& b, const Reg& target_reg);

    template <typename X, typename Y, typename Z>
    void do_store_stack(domains::NumAbsDomain& inv, const number_t& width, const linear_expression_t& addr, X val_type,
                        Y val_svalue, Z val_uvalue, const std::optional<reg_pack_t>& opt_val_reg);

    template <typename Type, typename SValue, typename UValue>
    void do_mem_store(const Mem& b, Type val_type, SValue val_svalue, UValue val_uvalue,
                      const std::optional<reg_pack_t>& opt_val_reg);

    friend std::ostream& operator<<(std::ostream& o, const ebpf_domain_t& dom);

    static void initialize_packet(ebpf_domain_t& inv);

  private:
    /// Mapping from variables (including registers, types, offsets,
    /// memory locations, etc.) to numeric intervals or relationships
    /// to other variables.
    domains::NumAbsDomain m_inv;

    /// Represents the stack as a memory region, i.e., an array of bytes,
    /// allowing mapping to variable in the m_inv numeric domains
    /// while dealing with overlapping byte ranges.
    domains::array_domain_t stack;

    std::function<check_require_func_t> check_require{};
    bool get_map_fd_range(const Reg& map_fd_reg, int32_t* start_fd, int32_t* end_fd) const;

    TypeDomain type_inv;
    std::string current_assertion;
}; // end ebpf_domain_t

} // namespace crab
