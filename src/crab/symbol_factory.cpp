// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
/*
 * Factories for symbol names.
 */

#include "crab/symbol.hpp"

namespace crab {

int64_t symbol_t::count = 4;

void symbol_t::write(std::ostream& o) const {
    if (is_begin()) {
        o << "begin";
    } else if (is_end()) {
        o << "end";
    } else if (is_meta()) {
        o << "meta";
    } else if (is_nu()) {
        o << "v";
    } else if (is_pkt_symbol()) {
        o << "i";
    } else {
        o << "a_" << ((int)_id-4);
    }
}

std::ostream& operator<<(std::ostream& o, const symbol_t& s) {
    s.write(o);
    return o;
}

} // namespace crab
