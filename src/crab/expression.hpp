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
    
using slacks_t = std::map<symbol_t, mock_interval_t>;
using symbol_terms_t = std::map<symbol_t, int8_t>;
class expression_t {
  private:
    symbol_terms_t _symbol_terms;
    interval_t _constant_term;
    std::shared_ptr<slacks_t> _slacks = nullptr;

  public:
    expression_t() : _constant_term(interval_t::top()) {};
    expression_t(symbol_t symbol, std::shared_ptr<slacks_t> slacks = nullptr)
        : _constant_term(interval_t{0}), _slacks(slacks) {
        _symbol_terms[symbol] = 1;
    }
    expression_t(symbol_terms_t symbol_terms, interval_t interval,
            std::shared_ptr<slacks_t> slacks = nullptr)
        : _symbol_terms(symbol_terms), _constant_term(interval), _slacks(slacks) {}
    expression_t(interval_t interval, std::shared_ptr<slacks_t> slacks = nullptr)
        : _constant_term(interval), _slacks(slacks) {}
    expression_t(int n, std::shared_ptr<slacks_t> slacks = nullptr)
        : _constant_term(interval_t{n}), _slacks(slacks) {}
    expression_t get_equivalent_expression() const;
    bool operator==(const expression_t &other) const;
    bool operator<(const expression_t &other) const;
    bool operator<=(const expression_t &other) const;
    bool operator>(const expression_t &other) const;
    void operator-=(const symbol_t& s) { _symbol_terms.erase(s); }
    bool is_constant() const { return _symbol_terms.empty(); }
    expression_t operator+(const expression_t &other) const;
    expression_t operator+(interval_t) const;
    expression_t operator+(int n) const;
    expression_t operator-(const expression_t &other) const;
    expression_t operator|(const expression_t &other) const;
    void write(std::ostream &o) const;
    symbol_terms_t get_symbol_terms() const { return _symbol_terms; }
    int8_t get_coefficient(const symbol_t &s) const;
    bool is_singleton() const;
    bool contains(const symbol_t &s) const;
    bool only_has_interval() const;
    symbol_t get_singleton() const;
    interval_t get_constant_term() const { return _constant_term; }
    std::shared_ptr<slacks_t> get_slacks() const { return _slacks; }
    void set_slacks(std::shared_ptr<slacks_t> slacks) { _slacks = slacks; }
    std::map<symbol_t, mock_interval_t> get_slack_intervals() const;
    friend std::ostream& operator<<(std::ostream &, const expression_t&);

    static expression_t begin(std::shared_ptr<slacks_t> slacks = nullptr) {
        return expression_t(symbol_t::begin(), slacks);
    }

    static expression_t end(std::shared_ptr<slacks_t> slacks = nullptr) {
        return expression_t(symbol_t::end(), slacks);
    }

    static expression_t meta(std::shared_ptr<slacks_t> slacks = nullptr) {
        return expression_t(symbol_t::meta(), slacks);
    }

    static expression_t top(std::shared_ptr<slacks_t> slacks = nullptr) {
        return expression_t(interval_t::top(), slacks);
    }
};


}  // namespace crab
