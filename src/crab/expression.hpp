// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include <map>
#include <cassert>
#include "symbol.hpp"
#include "interval.hpp"
#include "common.hpp"

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
    expression_t(symbol_terms_t symbol_terms, interval_t interval = interval_t{0})
        : _symbol_terms(symbol_terms), _constant_term(interval) {}
    expression_t(symbol_t symbol, std::shared_ptr<slacks_t> slacks)
        : _constant_term(interval_t{0}), _slacks(slacks) {
        _symbol_terms[symbol] = 1;
    }
    expression_t(symbol_terms_t symbol_terms, interval_t interval,
            std::shared_ptr<slacks_t> slacks)
        : _symbol_terms(symbol_terms), _constant_term(interval), _slacks(slacks) {}
    expression_t(symbol_terms_t symbol_terms, std::shared_ptr<slacks_t> slacks)
        : _symbol_terms(symbol_terms), _constant_term{interval_t{0}}, _slacks(slacks) {}
    expression_t(interval_t interval)
        : _constant_term(interval) {}
    expression_t(int n)
        : _constant_term(interval_t{n}) {}
    expression_t get_equivalent_expression() const;
    bool is_equal(const expression_t &other) const;
    bool is_less_than(const expression_t &other) const;
    bool is_less_or_equal(const expression_t &other) const;
    bool is_greater_than(const expression_t &other) const;
    bool operator==(const expression_t &other) const;
    bool operator<(const expression_t &other) const;
    bool operator<=(const expression_t &other) const;
    bool operator>(const expression_t &other) const;
    void operator-=(const symbol_t& s) { _symbol_terms.erase(s); }
    bool is_constant() const { return _symbol_terms.empty(); }
    expression_t negate() const;
    expression_t operator+(const expression_t &other) const;
    expression_t operator+(interval_t) const;
    expression_t operator+(int n) const;
    expression_t operator|(const expression_t &other) const;
    void write(std::ostream &o) const;
    symbol_terms_t get_symbol_terms() const { return _symbol_terms; }
    bool is_singleton() const;
    expression_t substitute(const symbol_t &, const expression_t &) const;
    bool contains(const symbol_t &s) const;
    symbol_t get_singleton() const;
    interval_t get_constant_term() const { return _constant_term; }
    std::shared_ptr<slacks_t> get_slacks() const { return _slacks; }

    static expression_t begin() {
        return expression_t({std::make_pair(symbol_t::begin(), 1)});
    }

    static expression_t end() {
        return expression_t({std::make_pair(symbol_t::end(), 1)});
    }

    static expression_t meta() {
        return expression_t({std::make_pair(symbol_t::meta(), 1)});
    }

    static expression_t nu() {
        return expression_t({std::make_pair(symbol_t::nu(), 1)});
    }

    static expression_t pkt_symbol() {
        return expression_t({std::make_pair(symbol_t::pkt_symbol(), 1)});
    }

    static expression_t make_slack() {
        return expression_t({std::make_pair(symbol_t::make(), 1)});
    }
};

std::ostream& operator<<(std::ostream &, const expression_t&);

}  // namespace crab
