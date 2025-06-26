// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include "constraint.hpp"

namespace crab {

enum class refinement_type_t {
    NUM,    // numeric refinement
    PACKET, // packet refinement
    ANY     // represents any refinement
};

class refinement_t {
    refinement_type_t _type;
    expression_t _value;
    constraint_t _meta_begin_constraint; // meta + sth <= begin
    constraint_t _begin_end_constraint; // begin + sth <= end
    bool _contains_pkt_constraints = false;

  public:
    refinement_t() : _type(refinement_type_t::ANY) {}
    refinement_t(refinement_type_t type, expression_t value,
                 constraint_t meta_begin_constraint, constraint_t begin_end_constraint,
                 bool contains_pkt_constraints) :
        _type(type), _value(value), _meta_begin_constraint(meta_begin_constraint),
        _begin_end_constraint(begin_end_constraint),
        _contains_pkt_constraints(contains_pkt_constraints) {}
    refinement_t(refinement_type_t type, expression_t value) :
        _type(type), _value(value), _meta_begin_constraint(constraint_t::true_constraint()),
        _begin_end_constraint(constraint_t::true_constraint()) {}
 
    [[nodiscard]] interval_t get_interval_value(std::shared_ptr<slacks_t>) const;
    [[nodiscard]] refinement_type_t get_type() const { return _type; }
    [[nodiscard]] expression_t get_value() const { return _value; }
    refinement_t operator+(int n) const;
    refinement_t operator+(interval_t) const;
    refinement_t operator-(interval_t) const;
    refinement_t add(const refinement_t &other, std::shared_ptr<slacks_t>) const;
    refinement_t subtract(const refinement_t &other, std::shared_ptr<slacks_t>) const;
    refinement_t join(const refinement_t &other, std::shared_ptr<slacks_t>) const;
    refinement_t widen(const refinement_t &other, std::shared_ptr<slacks_t>) const;
    bool inclusion(const refinement_t &other, std::shared_ptr<slacks_t>) const;
    constraint_t assume_le(const refinement_t &other) const;
    constraint_t assume_gt(const refinement_t &other) const;
    bool check_eq(const refinement_t &other, std::shared_ptr<slacks_t>) const;
    void write(std::ostream &o, std::shared_ptr<slacks_t>, bool) const;
    bool same_type(const refinement_t &other) const;
    std::pair<bool, bool> safe_access(const expression_t&, const expression_t&, bool,
                     std::shared_ptr<slacks_t>) const;
    bool check_consistent(const constraint_t&, std::shared_ptr<slacks_t>) const;
    constraint_t construct_meta_end_constraint() const;
    bool is_numeric_refinement() const;
    bool is_packet_refinement() const;
    interval_t simplify_for_subtraction(const symbol_t&, const symbol_t&,
                                        std::shared_ptr<slacks_t>) const;
    std::vector<symbol_t> get_all_slacks() const;
    void add_constraint(constraint_t, std::shared_ptr<slacks_t>);

    // construct a packet refinement type for begin of packet with initial constraints:
    // meta <= begin and begin <= end
    static refinement_t begin_with_constraints() {
        return refinement_t(refinement_type_t::PACKET,
                            expression_t::begin(),
                            constraint_t::meta_begin_init_constraint(),
                            constraint_t::begin_end_init_constraint(),
                            true);
    }

    static refinement_t begin() {
        return refinement_t(refinement_type_t::PACKET, expression_t::begin());
    }

    static refinement_t end() {
        return refinement_t(refinement_type_t::PACKET, expression_t::end());
    }

    static refinement_t meta() {
        return refinement_t(refinement_type_t::PACKET, expression_t::meta());
    }

    static refinement_t numeric_refinement(expression_t value) {
        return refinement_t(refinement_type_t::NUM, value);
    }

    static refinement_t numeric_refinement(interval_t i, std::shared_ptr<slacks_t> slacks) {
        symbol_t s = symbol_t::make();
        auto it = slacks->find(s);
        if (it != slacks->end()) {
            // Slack already exists, this should not happen.
            throw std::runtime_error("Newly-generated slack already exists.");
        }
        slacks->emplace(s, i);
        return numeric_refinement(expression_t(s));
    }

    static refinement_t numeric_refinement_top() {
        // expression_t::top() represents a top interval
        return numeric_refinement(expression_t::top());
    }
};


} // namespace crab
