// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "expression.hpp"

namespace crab {

std::ostream& operator<<(std::ostream& o, const expression_t& e) {
    e.write(o);
    return o;
}

// get values for all slack variables in the expression
std::map<symbol_t, interval_t> expression_t::get_slack_intervals() const {
    std::map<symbol_t, interval_t> slack_intervals;
    for (const auto &term : _symbol_terms) {
        if (term.first.is_slack()) {
            // assuming that there are no conflicting values
            // TODO: check what to do in case value is not found
            auto it = _slacks->find(term.first);
            if (it != _slacks->end()) {
                slack_intervals.insert_or_assign(term.first, it->second);
            }
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
            auto it = _slacks->find(term.first);
            // TODO: check what to do in case value is not found
            if (it != _slacks->end()) {
                value = value + (it->second * interval_t{static_cast<int>(term.second)});
            }
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

// less than or equal, named check_le to avoid confusion with operator<= (inclusion operator)
bool expression_t::check_le(const expression_t &other) const {
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

using IntervalJoin = std::function<interval_t(const interval_t&, const interval_t&)>;

static inline expression_t join(const expression_t* e1, const expression_t &e2,
                                IntervalJoin join) {
    auto slacks = e1->get_slacks() == nullptr ? e2.get_slacks() : e1->get_slacks();
    if (*e1 == e2) {
        // both expressions are equal
        return *e1;
    }
    interval_t e1_constant_term = e1->get_constant_term();
    interval_t e2_constant_term = e2.get_constant_term();
    const symbol_terms_t& e1_terms = e1->get_symbol_terms();
    const symbol_terms_t& e2_terms = e2.get_symbol_terms();
    symbol_terms_t new_terms;
    for (const auto &term : e1_terms) {
        if (e2_terms.find(term.first) != e2_terms.end()
            && term.second == e2_terms.at(term.first)) {
            // for any symbol, if it exists in both expressions, and the coefficients are the same
            // then the term can be added to the new expression
            new_terms[term.first] = term.second;
        } else if (term.first.is_slack()) {
            // either the symbol was not found in the e2 expression, or the coefficients are different
            auto it = slacks->find(term.first);
            if (it != slacks->end()) {
                e1_constant_term = e1_constant_term +
                    (interval_t{static_cast<int>(term.second)} * it->second);
            }
        }
        else {
            // in case where we have symbols other than slack variables, and either the symbol is not
            // found in the e2 expression, or the coefficients are different, we return top
            return expression_t();
        }
    }
    for (const auto &term : e2_terms) {
        if (new_terms.find(term.first) != new_terms.end()) {
            // if the symbol is already handled in the previous loop, then we skip it
            continue;
        }
        if (term.first.is_slack()) {
            // the slack variable was not handled before, hence we resolve its value
            auto it = slacks->find(term.first);
            if (it != slacks->end()) {
                e2_constant_term = e2_constant_term +
                    (interval_t{static_cast<int>(term.second)} * it->second);
            }
        }
        else {
            // in case where we have symbols other than slack variables, and either the symbol is not
            // found in the e2 expression, or the coefficients are different, we return top
            return expression_t();
        }
    }
    // in case there are different slack variables in the two expressions, we can't merge them,
    // so we construct a new slack variable and assign its value to the join of added intervals
    symbol_t new_slack = symbol_t::make();
    new_terms[new_slack] = 1;
    // either join or widen
    slacks->insert_or_assign(new_slack, join(e1_constant_term, e2_constant_term));
    return expression_t(new_terms, interval_t{0}, slacks);
}

expression_t expression_t::operator|(const expression_t &other) const {
    return join(this, other, [](const interval_t& a, const interval_t& b) { return a | b; });
}

expression_t expression_t::widen(const expression_t &other) const {
    return join(this, other, [](const interval_t& a, const interval_t& b) { return a.widen(b); });
}


// inclusion operator
bool expression_t::operator<=(const expression_t &other) const {
    expression_t e1 = get_equivalent_expression();
    expression_t e2 = other.get_equivalent_expression();
    // if after removing slack variables, the expressions are equal, then the inclusion holds
    // otherwise, we can't prove the inclusion
    return e1._symbol_terms == e2._symbol_terms && e1._constant_term <= e2._constant_term;
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
