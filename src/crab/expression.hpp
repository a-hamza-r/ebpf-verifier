// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include <map>
#include <cassert>
#include "symbol.hpp"
#include "types.hpp"

namespace crab {
// An expression is of form: Ax + By + Cz + ... + I.
// x, y, z, ... are symbols, and A, B, C, ... are coefficients.
// I is an interval.
    
using slacks_t = std::map<symbol_t, interval_t>; // mapping from slack symbols to their values
using symbol_terms_t = std::map<symbol_t, int8_t>; // mapping from symbols to their coefficients
class expression_t {
  private:
    symbol_terms_t _symbol_terms;
    interval_t _constant_term;

  public:
    expression_t() : _constant_term(interval_t::top()) {};
    // We assume that the symbol is either a pkt symbol, or a slack variable whose interval is known.
    expression_t(symbol_t symbol) : _constant_term(interval_t{0}) {
        _symbol_terms[symbol] = 1;
    }
    expression_t(symbol_terms_t symbol_terms, interval_t interval)
        : _symbol_terms(std::move(symbol_terms)), _constant_term(interval) {}
    expression_t(interval_t interval) : _constant_term(interval) {}
    expression_t(int n) : _constant_term(interval_t{n}) {}
    expression_t get_equivalent_expression(std::shared_ptr<slacks_t>) const;
    bool check_eq(const expression_t &other, std::shared_ptr<slacks_t>) const;
    bool check_le(const expression_t &other, std::shared_ptr<slacks_t>) const;
    bool check_gt(const expression_t &other, std::shared_ptr<slacks_t>) const;
    void operator-=(const symbol_t& s) { _symbol_terms.erase(s); }
    bool is_constant() const { return _symbol_terms.empty(); }
    interval_t get_constant_term() const { return _constant_term; }
    expression_t operator+(const expression_t &other) const;
    expression_t operator+(interval_t) const;
    expression_t operator+(int n) const;
    expression_t operator-(const expression_t &other) const;
    expression_t join(const expression_t &other, std::shared_ptr<slacks_t>) const;
    expression_t widen(const expression_t &other, std::shared_ptr<slacks_t>) const;
    bool inclusion(const expression_t &other, std::shared_ptr<slacks_t>) const;
    void write(std::ostream &o) const;
    symbol_terms_t get_symbol_terms() const { return _symbol_terms; }
    int get_num_terms() const { return _symbol_terms.size(); }
    int8_t get_coefficient(const symbol_t &s) const;
    bool is_singleton() const;
    bool contains(const symbol_t &s) const;
    symbol_t get_singleton() const;
    std::vector<symbol_t> get_slacks() const;
    friend std::ostream& operator<<(std::ostream &, const expression_t&);

    static expression_t begin() {
        return expression_t(symbol_t::begin());
    }

    static expression_t end() {
        return expression_t(symbol_t::end());
    }

    static expression_t meta() {
        return expression_t(symbol_t::meta());
    }

    static expression_t top() {
        return expression_t(interval_t::top());
    }
};


}  // namespace crab
