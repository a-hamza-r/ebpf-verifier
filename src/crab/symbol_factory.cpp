// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
/*
 * Factories for symbol names.
 */

#include "symbol.hpp"

namespace crab {

uint64_t symbol_t::count = 3;

void symbol_t::write(std::ostream& o) const {
    if (is_begin()) {
        o << "begin";
    } else if (is_end()) {
        o << "end";
    } else if (is_meta()) {
        o << "meta";
    //} else if (is_nu()) {
    //    o << "v";
    } else {
        o << "s_" << ((uint64_t)_id-3);
    }
}

std::ostream& operator<<(std::ostream& o, const symbol_t& s) {
    s.write(o);
    return o;
}

} // namespace crab
