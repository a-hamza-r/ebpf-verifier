// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include "constraint.hpp"

namespace crab {

enum class data_type_t {
    NUM,
    PACKET,
    ANY
};

class refinement_t {
    data_type_t _type;
    expression_t _value;
    std::vector<crab::constraint_t> _constraints;

  public:
    explicit refinement_t(data_type_t type, expression_t value,
            std::vector<crab::constraint_t> constraints = {})
        : _type(type), _value(value), _constraints(constraints) {}
    refinement_t() : _type(data_type_t::ANY) {}
 
    bool is_bottom(std::shared_ptr<slacks_t>);
    [[nodiscard]] data_type_t get_type() const { return _type; }
    [[nodiscard]] std::vector<crab::constraint_t> get_constraints() const { return _constraints; }
    [[nodiscard]] expression_t get_value() const { return _value; }
    refinement_t operator+(int n) const;
    refinement_t operator+(interval_t) const;
    refinement_t operator+(const refinement_t &other) const;
    refinement_t operator-(const refinement_t &other) const;
    refinement_t operator|(const refinement_t &other) const;
    constraint_t operator<=(const refinement_t &other) const;
    constraint_t operator>(const refinement_t &other) const;
    bool has_value(refinement_t&&) const;
    void write(std::ostream &o) const;
    bool same_type(const refinement_t &other) const;
    void add_constraint(const constraint_t&);
    void simplify();
    bool is_safe_with(refinement_t, std::shared_ptr<slacks_t>, bool) const;

    static refinement_t begin() {
        return refinement_t(data_type_t::PACKET, expression_t::begin());
    }

    static refinement_t end() {
        return refinement_t(data_type_t::PACKET, expression_t::end());
    }

    static refinement_t meta() {
        return refinement_t(data_type_t::PACKET, expression_t::meta());
    }

    static refinement_t make_slack() {
        return refinement_t(data_type_t::NUM, expression_t::make_slack());
    }
};

std::ostream &operator<<(std::ostream &, const refinement_t&);

} // namespace crab
