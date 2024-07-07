// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include <boost/optional/optional_io.hpp>
#include <functional>
#include <optional>
#include <vector>
#include <unordered_map>

#include "string_constraints.hpp"
#include "asm_syntax.hpp"
#include "array_domain.hpp"

namespace crab {

using check_require_func_t = std::function<bool(crab::domains::NumAbsDomain&, const crab::linear_constraint_t&, std::string)>;

// 11 registers for the eBPF ISA, and 1 pseudo-register for the atomic operations
constexpr uint8_t NUM_REGISTERS = 12;

constexpr int STACK_BEGIN = 0;
constexpr int CTX_BEGIN = 0;
constexpr int PACKET_BEGIN = 0;
constexpr int SHARED_BEGIN = 0;
constexpr int MAX_PACKET_SIZE = 0xffff;

enum class region_t {
	T_CTX,
	T_STACK,
	T_PACKET,
	T_SHARED
};

inline std::string get_reg_ptr(const region_t& r) noexcept;

enum class nullness_t { MAYBE_NULL, NOT_NULL, _NULL };

class packet_ptr_t {
    region_t m_r = region_t::T_PACKET;

  public:
    packet_ptr_t() = default;
    packet_ptr_t(const packet_ptr_t &) = default;
    packet_ptr_t(packet_ptr_t &&) = default;
    packet_ptr_t &operator=(const packet_ptr_t &) = default;
    packet_ptr_t &operator=(packet_ptr_t &&) = default;

    friend std::ostream& operator<<(std::ostream& o, const packet_ptr_t& p);
    bool operator==(const packet_ptr_t&) const;
    bool operator!=(const packet_ptr_t&) const;
};

class mock_interval_t {
    bound_t _lb;
    bound_t _ub;

    public:
        static mock_interval_t top() {
            return mock_interval_t(bound_t::minus_infinity(), bound_t::plus_infinity());
        }
        [[nodiscard]] bound_t lb() const { return _lb; }
        [[nodiscard]] bound_t ub() const { return _ub; }
        mock_interval_t(bound_t lb, bound_t ub) : _lb(lb), _ub(ub) {};
        mock_interval_t() : _lb(bound_t::minus_infinity()), _ub(bound_t::plus_infinity()) {}
        mock_interval_t(const mock_interval_t& c) = default;
        mock_interval_t(const bound_t& b) : _lb(b), _ub(b) {}
        mock_interval_t& operator=(const mock_interval_t& o) = default;
        bool operator==(const mock_interval_t& o) const;
        mock_interval_t(const interval_t& i) : _lb(i.lb()), _ub(i.ub()) {}
        interval_t to_interval() const { return std::move(interval_t(_lb, _ub)); }
};

class ptr_with_off_t {
    region_t m_r;
    int m_id;
    mock_interval_t m_offset;
    nullness_t m_nullness = nullness_t::MAYBE_NULL;
    mock_interval_t m_region_size = mock_interval_t::top();

  public:
    ptr_with_off_t() = default;
    ptr_with_off_t(const ptr_with_off_t &) = default;
    ptr_with_off_t(ptr_with_off_t &&) = default;
    ptr_with_off_t &operator=(const ptr_with_off_t &) = default;
    ptr_with_off_t &operator=(ptr_with_off_t &&) = default;
    ptr_with_off_t(region_t _r, int _id, mock_interval_t _off,
            nullness_t _nullness = nullness_t::MAYBE_NULL,
            mock_interval_t _region_sz = mock_interval_t::top())
        : m_r(_r), m_id(_id), m_offset(_off), m_nullness(_nullness), m_region_size(_region_sz) {}
    ptr_with_off_t operator|(const ptr_with_off_t&) const;
    ptr_with_off_t widen(const ptr_with_off_t&) const;
    bool operator<=(const ptr_with_off_t&) const;
    [[nodiscard]] nullness_t get_nullness() const { return m_nullness; }
    void set_nullness(nullness_t);
    [[nodiscard]] int get_id() const { return m_id; }
    void set_id(int);
    [[nodiscard]] mock_interval_t get_region_size() const { return m_region_size; }
    void set_region_size(mock_interval_t);
    [[nodiscard]] mock_interval_t get_offset() const { return m_offset; }
    void set_offset(mock_interval_t);
    [[nodiscard]] region_t get_region() const { return m_r; }
    void set_region(region_t);
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream& o, const ptr_with_off_t& p);
    bool operator==(const ptr_with_off_t&) const;
    bool operator!=(const ptr_with_off_t&) const;
};

class mapfd_t {
    mock_interval_t m_mapfd;
    EbpfMapValueType m_value_type;

  public:
    mapfd_t(const mapfd_t&) = default;
    mapfd_t(mapfd_t&&) = default;
    mapfd_t &operator=(const mapfd_t&) = default;
    mapfd_t &operator=(mapfd_t&&) = default;
    mapfd_t operator|(const mapfd_t&) const;
    mapfd_t widen(const mapfd_t&) const;
    bool operator<=(const mapfd_t&) const;
    mapfd_t(mock_interval_t mapfd, EbpfMapValueType val_type)
        : m_mapfd(mapfd), m_value_type(val_type) {}
    friend std::ostream& operator<<(std::ostream&, const mapfd_t&);
    bool operator==(const mapfd_t&) const;
    bool operator!=(const mapfd_t&) const;
    void write(std::ostream&) const;

    bool has_type_map_programs() const;
    [[nodiscard]] EbpfMapValueType get_value_type() const { return m_value_type; }
    [[nodiscard]] mock_interval_t get_mapfd() const { return m_mapfd; }
};

using ptr_t = std::variant<packet_ptr_t, ptr_with_off_t>;
using register_t = uint8_t;
using location_t = boost::optional<std::pair<label_t, uint32_t>>;

class reg_with_loc_t {
    register_t m_reg;
    location_t m_loc;

  public:
    reg_with_loc_t(register_t _r, location_t _loc) : m_reg(_r), m_loc(_loc) {}
    bool operator==(const reg_with_loc_t& other) const;
    std::size_t hash() const;
    friend std::ostream& operator<<(std::ostream& o, const reg_with_loc_t& reg);
    void write(std::ostream& ) const;
};

using ptr_or_mapfd_t = std::variant<ptr_with_off_t, packet_ptr_t, mapfd_t>;

inline bool is_ptr_type(const std::optional<ptr_or_mapfd_t>& ptr_or_mapfd) {
    return (ptr_or_mapfd && (std::holds_alternative<ptr_with_off_t>(*ptr_or_mapfd)
                || std::holds_alternative<packet_ptr_t>(*ptr_or_mapfd)));
}

inline bool is_mapfd_type(const std::optional<ptr_or_mapfd_t>& ptr_or_mapfd) {
    return (ptr_or_mapfd && std::holds_alternative<mapfd_t>(*ptr_or_mapfd));
}

inline bool same_region(const ptr_or_mapfd_t& ptr1, const ptr_or_mapfd_t& ptr2) {
    if (std::holds_alternative<packet_ptr_t>(ptr1) && std::holds_alternative<packet_ptr_t>(ptr2))
        return true;
    return (std::holds_alternative<ptr_with_off_t>(ptr1)
            && std::holds_alternative<ptr_with_off_t>(ptr2)
            && std::get<ptr_with_off_t>(ptr1).get_region()
                == std::get<ptr_with_off_t>(ptr2).get_region());
}

inline bool is_stack_ptr(const std::optional<ptr_or_mapfd_t>& ptr) {
    return (ptr && std::holds_alternative<ptr_with_off_t>(*ptr)
            && std::get<ptr_with_off_t>(*ptr).get_region() == region_t::T_STACK);
}

inline bool is_ctx_ptr(const std::optional<ptr_or_mapfd_t>& ptr) {
    return (ptr && std::holds_alternative<ptr_with_off_t>(*ptr)
            && std::get<ptr_with_off_t>(*ptr).get_region() == region_t::T_CTX);
}

inline bool is_packet_ptr(const std::optional<ptr_or_mapfd_t>& ptr) {
    return (ptr && std::holds_alternative<packet_ptr_t>(*ptr));
}

inline bool is_shared_ptr(const std::optional<ptr_or_mapfd_t>& ptr) {
    return (ptr && std::holds_alternative<ptr_with_off_t>(*ptr)
            && std::get<ptr_with_off_t>(*ptr).get_region() == region_t::T_SHARED);
}

inline bool same_type(const std::optional<ptr_or_mapfd_t>& ptr_or_mapfd1,
        const std::optional<ptr_or_mapfd_t>& ptr_or_mapfd2,
        const std::optional<mock_interval_t>& interval1,
        const std::optional<mock_interval_t>& interval2) {
    if (is_mapfd_type(ptr_or_mapfd1) && is_mapfd_type(ptr_or_mapfd2)) return true;
    if (ptr_or_mapfd1 && ptr_or_mapfd2 && same_region(*ptr_or_mapfd1, *ptr_or_mapfd2))
        return true;
    if (interval1 && interval2) return true;
    return false;
}

inline std::ostream& operator<<(std::ostream& o, const region_t& t) {
    o << static_cast<std::underlying_type<region_t>::type>(t);
    return o;
}

using stack_cells_t = std::vector<std::pair<uint64_t, int>>;

} // namespace crab


namespace std {
    template <>
    struct hash<crab::reg_with_loc_t> {
        size_t operator()(const crab::reg_with_loc_t& reg) const { return reg.hash(); }
    };

    template <>
    struct equal_to<crab::ptr_t> {
        constexpr bool operator()(const crab::ptr_t& lhs, const crab::ptr_t& rhs) const {
            if (lhs.index() != rhs.index()) return false;
            return std::visit( overloaded
               {
                   []( const crab::ptr_with_off_t& x, const crab::ptr_with_off_t& y ){ return x == y;},
                   []( const crab::packet_ptr_t& x, const crab::packet_ptr_t& y ){ return x == y;},
                   []( auto& , auto& ) { return true;}
                }, lhs, rhs
            );
        }
    };

    template <>
    struct equal_to<crab::ptr_or_mapfd_t> {
        constexpr bool operator()(const crab::ptr_or_mapfd_t& lhs, const crab::ptr_or_mapfd_t& rhs) const {
            if (lhs.index() != rhs.index()) return false;
            return std::visit( overloaded
               {
                   []( const crab::ptr_with_off_t& x, const crab::ptr_with_off_t& y ){ return x == y;},
                   []( const crab::packet_ptr_t& x, const crab::packet_ptr_t& y ){ return x == y;},
                   []( const crab::mapfd_t& x, const crab::mapfd_t& y ){ return x == y;},
                   []( auto& , auto& ) { return true;}
                }, lhs, rhs
            );
        }
    };

    //crab::ptr_t get_ptr(const crab::ptr_or_mapfd_t& t);
}
