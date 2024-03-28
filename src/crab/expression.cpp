
#include "expression.hpp"

namespace crab {

std::ostream& operator<<(std::ostream& o, const expression_t& e) {
    e.write(o);
    return o;
}

bool expression_t::operator==(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _interval == other._interval;
}

bool expression_t::operator<(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _interval.ub() < other._interval.lb();
}

bool expression_t::operator<=(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _interval.ub() <= other._interval.lb();
}

bool expression_t::operator>(const expression_t &other) const {
    return _symbol_terms == other._symbol_terms && _interval.lb() > other._interval.ub();
}

bool expression_t::is_singleton() const {
    return _symbol_terms.size() == 1 && _interval == interval_t{0};
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
    return expression_t(new_terms, -_interval);
}

expression_t expression_t::operator+(const expression_t &other) const {
    auto new_terms = _symbol_terms;
    for (const auto &term : other._symbol_terms) {
        insert(new_terms, term);
    }
    return expression_t(new_terms, _interval + other._interval);
}

expression_t expression_t::operator+(interval_t interval) const {
    return expression_t(_symbol_terms, _interval + interval);
}

expression_t expression_t::operator+(int n) const {
    return operator+(interval_t{n});
}

expression_t expression_t::operator|(const expression_t &other) const {
    if (_symbol_terms.size() != other._symbol_terms.size()) {
        return expression_t();
    }
    symbol_terms_t new_terms;
    for (auto &term : _symbol_terms) {
        auto it = other._symbol_terms.find(term.first);
        if (it == other._symbol_terms.end()) {
            // we need to know what is the other slack, but it's not accesible directly
            //if (term.first.is_slack() && it->first.is_slack()) {
                auto new_slack = symbol_t::make();
                new_terms[new_slack] = 1;
            //}
            //else {
            //    return expression_t();
            //}
        }
        else {
            // TODO: some complex logic is needed here
            new_terms[term.first] = term.second;
        }
    }
    interval_t new_interval = _interval | other._interval;
    return expression_t(new_terms, new_interval);
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
    if (auto s = _interval.singleton()) {
        if ((int)*s != 0 || _symbol_terms.empty()) {
            if (!_symbol_terms.empty()) o << " + ";
            o <<  *s;
        }
    } else {
        o << " + " << _interval;
    }
}

}  // namespace crab
