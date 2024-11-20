
#include "expression.hpp"

namespace crab {

std::ostream& operator<<(std::ostream& o, const expression_t& e) {
    e.write(o);
    return o;
}

// get values for all slack variables in the expression
std::vector<std::pair<symbol_t, interval_t>> expression_t::get_slack_intervals() const {
    std::vector<std::pair<symbol_t, interval_t>> slack_intervals;
    for (const auto &term : _symbol_terms) {
        if (term.first.is_slack()) {
            slack_intervals.push_back({term.first, (*_slacks)[term.first].to_interval()});
        }
    }
    return slack_intervals;
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

// check if two expressions are equal, by getting equivalent expressions and comparing them
bool expression_t::is_equal(const expression_t &other) const {
    return get_equivalent_expression() == other.get_equivalent_expression();
}

// check less than
bool expression_t::is_less_than(const expression_t &other) const {
    return get_equivalent_expression() < other.get_equivalent_expression();
}

// check less than or equal
bool expression_t::is_less_or_equal(const expression_t &other) const {
    return get_equivalent_expression() <= other.get_equivalent_expression();
}

// check greater than
bool expression_t::is_greater_than(const expression_t &other) const {
    return get_equivalent_expression() > other.get_equivalent_expression();
}

// syntactic equality
bool expression_t::operator==(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term == other._constant_term;
}

// syntactic less than
bool expression_t::operator<(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term.ub() < other._constant_term.lb();
}

// syntactic less than or equal
bool expression_t::operator<=(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term.ub() <= other._constant_term.lb();
}

// syntactic greater than
bool expression_t::operator>(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term.lb() > other._constant_term.ub();
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

// substitute a symbol in an expression with another expression
expression_t expression_t::substitute(const symbol_t &from, const expression_t &to) const {
    expression_t result = *this;
    auto it = result._symbol_terms.find(from);
    if (it != result._symbol_terms.end()) {
        result._symbol_terms.erase(it);
        result = result + to;
    }
    return result;
}

// negate an expression, and return it
expression_t expression_t::negate() const {
    symbol_terms_t new_terms;
    for (const auto &term : _symbol_terms) {
        new_terms[term.first] = -term.second;
    }
    return expression_t(new_terms, -_constant_term, _slacks);
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

// add an integer to an expression, and return the result
expression_t expression_t::operator+(int n) const {
    return operator+(interval_t{n});
}

expression_t expression_t::operator|(const expression_t &other) const {
    // TODO: support case where there are one slack variable in each expression,
    // then we can construct a new slack that is the interval of two.
    // we need to have intervals of slack variables to do this
    auto slacks = _slacks == nullptr ? other._slacks : _slacks;
    if (*this == other) {
        // both expressions are equal
        return expression_t(_symbol_terms, _constant_term, slacks);
    }
    else if (_symbol_terms == other._symbol_terms) {
        // both expressions have the same symbol terms, but different constant terms
        return expression_t(_symbol_terms, _constant_term | other._constant_term, slacks);
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
        if (term.first.is_slack()) {
            // if the symbol is a slack variable, it must not be in the first expression else it
            // would have been handled in the previous loop, hence we resolve its value
            other_constant_term = other_constant_term +
                interval_t{static_cast<int>(term.second)} * (*slacks)[term.first].to_interval();
        }
        else if (_symbol_terms.find(term.first) == _symbol_terms.end()
          || term.second != _symbol_terms.at(term.first)) {
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
