
#include "expression.hpp"

namespace crab {

std::ostream& operator<<(std::ostream& o, const expression_t& e) {
    e.write(o);
    return o;
}

expression_t expression_t::get_equivalent_expression() const {
    interval_t value = _constant_term;
    symbol_terms_t symbol_terms;
    for (const auto &term : _symbol_terms) {
        if (term.first.is_slack()) {
            value = value + ((*_slacks)[term.first]).to_interval() * interval_t{(int)term.second};
        } else {
            symbol_terms[term.first] = term.second;
        }
    }
    return expression_t(symbol_terms, value, _slacks);
}

bool expression_t::is_equal(const expression_t &other) const {
    return get_equivalent_expression() == other.get_equivalent_expression();
}

bool expression_t::is_less_than(const expression_t &other) const {
    return get_equivalent_expression() < other.get_equivalent_expression();
}

bool expression_t::is_less_or_equal(const expression_t &other) const {
    return get_equivalent_expression() <= other.get_equivalent_expression();
}

bool expression_t::is_greater_than(const expression_t &other) const {
    return get_equivalent_expression() > other.get_equivalent_expression();
}

bool expression_t::operator==(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term == other._constant_term;
}

bool expression_t::operator<(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term.ub() < other._constant_term.lb();
}

bool expression_t::operator<=(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term.ub() <= other._constant_term.lb();
}

bool expression_t::operator>(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _constant_term.lb() > other._constant_term.ub();
}

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

static void insert(symbol_terms_t &terms, const std::pair<symbol_t, int8_t>& kv) {
    if (terms.find(kv.first) != terms.end()) {
        terms[kv.first] += kv.second;
        if (terms[kv.first] == 0) {
            terms.erase(kv.first);
        }
    } else {
        terms[kv.first] = kv.second;
    }
}

expression_t expression_t::substitute(const symbol_t &from, const expression_t &to) const {
    expression_t result = *this;
    auto it = result._symbol_terms.find(from);
    if (it != result._symbol_terms.end()) {
        result._symbol_terms.erase(it);
        result = result + to;
    }
    return result;
}

expression_t expression_t::negate() const {
    symbol_terms_t new_terms;
    for (const auto &term : _symbol_terms) {
        new_terms[term.first] = -term.second;
    }
    return expression_t(new_terms, -_constant_term, _slacks);
}

expression_t expression_t::operator+(const expression_t &other) const {
    auto new_terms = _symbol_terms;
    for (const auto &term : other._symbol_terms) {
        insert(new_terms, term);
    }
    auto slacks = _slacks == nullptr ? other._slacks : _slacks;
    return expression_t(new_terms, _constant_term + other._constant_term, slacks);
}

expression_t expression_t::operator+(interval_t constant) const {
    return expression_t(_symbol_terms, _constant_term + constant, _slacks);
}

expression_t expression_t::operator+(int n) const {
    return operator+(interval_t{n});
}

expression_t expression_t::operator|(const expression_t &other) const {
    auto slacks = _slacks == nullptr ? other._slacks : _slacks;
    if (*this == other) {
        return expression_t(_symbol_terms, _constant_term, slacks);
    }
    interval_t constant_term = _constant_term;
    interval_t other_constant_term = other._constant_term;
    symbol_terms_t new_terms;
    for (const auto &term : _symbol_terms) {
        if (other._symbol_terms.find(term.first) != other._symbol_terms.end()) {
            new_terms[term.first] = term.second;
        } else {
            constant_term = constant_term +
                interval_t{(int)term.second} * (*slacks)[term.first].to_interval();
        }
    }
    for (const auto &term : other._symbol_terms) {
        if (_symbol_terms.find(term.first) == _symbol_terms.end()) {
            other_constant_term = other_constant_term +
                interval_t{(int)term.second} * (*slacks)[term.first].to_interval();
        }
    }
    symbol_t new_slack = symbol_t::make();
    new_terms[new_slack] = 1;
    (*slacks)[new_slack] = constant_term | other_constant_term;
    return expression_t(new_terms, interval_t{0}, slacks);
}

void expression_t::write(std::ostream &o) const {
    size_t i = 0;
    for (const auto &term : _symbol_terms) {
        if (term.second != 1) {
            o << (int)term.second << " * ";
        }
        o << term.first;
        if (i < _symbol_terms.size() - 1) {
            o << " + ";
        }
        i++;
    }
    if (auto s = _constant_term.singleton()) {
        if ((int)*s != 0 || _symbol_terms.empty()) {
            if (!_symbol_terms.empty()) o << " + ";
            o <<  *s;
        }
    } else {
        o << " + " << _constant_term;
    }
}

}  // namespace crab
