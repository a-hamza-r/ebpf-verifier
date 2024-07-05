// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "expression.hpp"

namespace crab {

std::ostream& operator<<(std::ostream& o, const expression_t& e) {
    e.write(o);
    return o;
}

// get values for all slack variables in the expression
std::map<symbol_t, mock_interval_t> expression_t::get_slack_intervals() const {
    std::map<symbol_t, mock_interval_t> slack_intervals;
    for (const auto &term : _symbol_terms) {
        if (term.first.is_slack()) {
            // assuming that there are no conflicting values
            slack_intervals[term.first] = (*_slacks)[term.first];
        }
    }
    return slack_intervals;
}

int8_t expression_t::get_coefficient(const symbol_t &s) const {
    auto it = _symbol_terms.find(s);
    if (it != _symbol_terms.end()) {
        return it->second;
    }
    return 0;
}

// get an equivalent expression with all slack variables replaced by their values
expression_t expression_t::get_equivalent_expression() const {
    interval_t value = _constant_term;
    symbol_terms_t symbol_terms;
    for (const auto &term : _symbol_terms) {
        if (term.first.is_slack()) {
            value = value + ((*_slacks)[term.first]).to_interval() * interval_t{static_cast<int>(term.second)};
        } else {
            symbol_terms[term.first] = term.second;
        }
    }
    return expression_t(symbol_terms, value, _slacks);
}

template <typename T, typename Op>
static bool check_op(const symbol_terms_t& s1, const symbol_terms_t& s2, const T& op1, const T& op2,
                     Op op) {
    return s1 == s2 && op(op1, op2);
}

// equality
bool expression_t::operator==(const expression_t &other) const {
    auto check = [](const interval_t &a, const interval_t &b) { return a == b; };
    if (check_op(_symbol_terms, other._symbol_terms, _constant_term, other._constant_term, check)) {
        // syntactic equality
        return true;
    }
    // check if two expressions are equal, by getting equivalent expressions and comparing them
    expression_t e1 = get_equivalent_expression();
    expression_t e2 = other.get_equivalent_expression();
    if (check_op(e1._symbol_terms, e2._symbol_terms, e1._constant_term, e2._constant_term, check)) {
        return true;
    }
    // expressions might still be equal, however, we can't prove it as it contains packet symbols
    return false;
}

// less than
bool expression_t::operator<(const expression_t &other) const {
    auto check = [](const bound_t &a, const bound_t &b) { return a < b; };
    if (check_op(_symbol_terms, other._symbol_terms, _constant_term.ub(), other._constant_term.lb(), check)) {
        return true;
    }
    expression_t e1 = get_equivalent_expression();
    expression_t e2 = other.get_equivalent_expression();
    if (check_op(e1._symbol_terms, e2._symbol_terms, e1._constant_term.ub(), e2._constant_term.lb(), check)) {
        return true;
    }
    // we cannot prove that one expression is less than the other, as they contain packet symbols
    return false;
}

// less than or equal
bool expression_t::operator<=(const expression_t &other) const {
    auto check = [](const bound_t &a, const bound_t &b) { return a <= b; };
    if (check_op(_symbol_terms, other._symbol_terms, _constant_term.ub(), other._constant_term.lb(), check)) {
        return true;
    }
    expression_t e1 = get_equivalent_expression();
    expression_t e2 = other.get_equivalent_expression();
    if (check_op(e1._symbol_terms, e2._symbol_terms, e1._constant_term.ub(), e2._constant_term.lb(), check)) {
        return true;
    }
    // we cannot prove that one expression is less than or equal to the other, as they contain packet symbols
    return false;
}

// greater than
bool expression_t::operator>(const expression_t &other) const {
    auto check = [](const bound_t &a, const bound_t &b) { return a > b; };
    if (check_op(_symbol_terms, other._symbol_terms, _constant_term.lb(), other._constant_term.ub(), check)) {
        return true;
    }
    expression_t e1 = get_equivalent_expression();
    expression_t e2 = other.get_equivalent_expression();
    if (check_op(e1._symbol_terms, e2._symbol_terms, e1._constant_term.lb(), e2._constant_term.ub(), check)) {
        return true;
    }
    // we cannot prove that one expression is greater than the other, as they contain packet symbols
    return false;
}

// check if an expression contains only an interval
bool expression_t::only_has_interval() const {
    return _symbol_terms.empty();
}

// check if an expression contains a single symbol, which can be slack or not
bool expression_t::is_singleton() const {
    return _symbol_terms.size() == 1 && _constant_term == interval_t{0};
}

bool expression_t::contains(const symbol_t &s) const {
    return _symbol_terms.find(s) != _symbol_terms.end();
}

symbol_t expression_t::get_singleton() const {
    if (!is_singleton()) {
        throw std::runtime_error("expression does not contain a single symbol");
    }
    return _symbol_terms.begin()->first;
}

// insert a symbol term into a list of symbol terms
static void insert(symbol_terms_t &terms, const std::pair<symbol_t, int8_t>& kv) {
    // add coefficient to term if it already exists
    if (terms.find(kv.first) != terms.end()) {
        terms[kv.first] += kv.second;
        // remove term if coefficient is zero
        if (terms[kv.first] == 0) {
            terms.erase(kv.first);
        }
    } else {
        terms[kv.first] = kv.second;
    }
}

// add two expressions, and return the result
expression_t expression_t::operator+(const expression_t &other) const {
    auto new_terms = _symbol_terms;
    for (const auto &term : other._symbol_terms) {
        insert(new_terms, term);
    }
    auto slacks = _slacks == nullptr ? other._slacks : _slacks;
    return expression_t(new_terms, _constant_term + other._constant_term, slacks);
}

// add a constant to an expression, and return the result
expression_t expression_t::operator+(interval_t constant) const {
    return expression_t(_symbol_terms, _constant_term + constant, _slacks);
}

// subtract an expression from another, and return the result
expression_t expression_t::operator-(const expression_t &other) const {
    auto new_terms = _symbol_terms;
    for (const auto &term : other._symbol_terms) {
        insert(new_terms, {term.first, -term.second});
    }
    auto slacks = _slacks == nullptr ? other._slacks : _slacks;
    return expression_t(new_terms, _constant_term - other._constant_term, slacks);
}

// add an integer to an expression, and return the result
expression_t expression_t::operator+(int n) const {
    return operator+(interval_t{n});
}

expression_t expression_t::operator|(const expression_t &other) const {
    auto slacks = _slacks == nullptr ? other._slacks : _slacks;
    if (*this == other) {
        // both expressions are equal
        return *this;
    }
    interval_t constant_term = _constant_term;
    interval_t other_constant_term = other._constant_term;
    symbol_terms_t new_terms;
    for (const auto &term : _symbol_terms) {
        if (other._symbol_terms.find(term.first) != other._symbol_terms.end()
            && term.second == other._symbol_terms.at(term.first)) {
            // for any symbol, if it exists in both expressions, and the coefficients are the same
            // then the term can be added to the new expression
            new_terms[term.first] = term.second;
        } else if (term.first.is_slack()) {
            // either the symbol was not found in the other expression, or the coefficients are different
            constant_term = constant_term +
                interval_t{static_cast<int>(term.second)} * (*slacks)[term.first].to_interval();
        }
        else {
            // in case where we have symbols other than slack variables, and either the symbol is not
            // found in the other expression, or the coefficients are different, we return top
            return expression_t();
        }
    }
    for (const auto &term : other._symbol_terms) {
        if (new_terms.find(term.first) != new_terms.end()) {
            // if the symbol is already handled in the previous loop, then we skip it
            continue;
        }
        if (term.first.is_slack()) {
            // the slack variable was not handled before, hence we resolve its value
            other_constant_term = other_constant_term +
                interval_t{static_cast<int>(term.second)} * (*slacks)[term.first].to_interval();
        }
        else {
            // in case where we have symbols other than slack variables, and either the symbol is not
            // found in the other expression, or the coefficients are different, we return top
            return expression_t();
        }
    }
    // in case there are different slack variables in the two expressions, we can't merge them,
    // so we construct a new slack variable and assign its value to the join of added intervals
    symbol_t new_slack = symbol_t::make();
    new_terms[new_slack] = 1;
    (*slacks)[new_slack] = constant_term | other_constant_term;
    return expression_t(new_terms, interval_t{0}, slacks);
}

expression_t expression_t::widen(const expression_t &other) const {
    // TODO: this needs more work
    if (_symbol_terms.size() != other._symbol_terms.size()) {
        return expression_t();
    }
    symbol_terms_t new_terms;
    for (auto &term : _symbol_terms) {
        auto it = other._symbol_terms.find(term.first);
        if (it == other._symbol_terms.end()) {
            // we need to know what is the other slack, but it's not accesible directly
                auto new_slack = symbol_t::make();
                new_terms[new_slack] = 1;
        }
        else {
            // TODO: some complex logic is needed here
            new_terms[term.first] = term.second;
        }
    }
    interval_t new_interval = _constant_term.widen(other._constant_term);
    return expression_t(new_terms, new_interval);
}

void expression_t::write(std::ostream &o) const {
    size_t i = 0;
    // print all symbol terms
    for (const auto &term : _symbol_terms) {
        if (term.second != 1) {
            o << static_cast<int>(term.second) << " * ";
        }
        o << term.first;
        if (i < _symbol_terms.size() - 1) {
            o << " + ";
        }
        i++;
    }
    // print constant term
    if (auto s = _constant_term.singleton()) {
        // print constant (singleton) term either if it is not zero,
        // of if it is zero but there are no other terms
        if ((cpp_int)*s != 0 || _symbol_terms.empty()) {
            if (!_symbol_terms.empty()) o << " + ";
            o <<  *s;
        }
    } else {
        if (!_symbol_terms.empty()) o << " + ";
        o << _constant_term;
    }
}

}  // namespace crab
