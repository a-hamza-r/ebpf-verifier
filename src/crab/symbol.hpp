// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>
#include <vector>

namespace crab {

using index_t = int64_t;

class symbol_t final {
    index_t _id;

    explicit symbol_t(index_t id) : _id(id) {}

  public:
    [[nodiscard]] std::size_t hash() const { return (size_t)_id; }
    bool operator==(symbol_t o) const { return _id == o._id; }
    bool operator!=(symbol_t o) const { return (!(operator==(o))); }
    bool operator<(symbol_t o) const { return _id < o._id; }
    bool operator>(symbol_t o) const { return _id > o._id; }

  private:
    static int64_t count;

  public:
    static symbol_t begin()     { return symbol_t(0); }
    static symbol_t end()       { return symbol_t(1); }
    static symbol_t meta()      { return symbol_t(2); }
    static symbol_t nu()        { return symbol_t(3); }
    static symbol_t pkt_symbol()     { return symbol_t(4); }
    static symbol_t make()      { count++;  return symbol_t(count); }
    bool is_nu() const { return *this == symbol_t::nu(); }
    bool is_pkt_symbol() const { return *this == symbol_t::pkt_symbol(); }
    bool is_meta() const { return *this == symbol_t::meta(); }
    bool is_end() const { return *this == symbol_t::end(); }
    bool is_begin() const { return *this == symbol_t::begin(); }
    bool is_slack() const { return _id >= 5; }
    void write(std::ostream& o) const;

    struct Hasher {
        std::size_t operator()(const symbol_t& s) const { return s.hash(); }
    };
};  // class symbol_t

std::ostream& operator<<(std::ostream& o, const symbol_t& s);

} // namespace crab
