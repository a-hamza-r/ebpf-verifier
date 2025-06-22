// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include <optional>
#include <vector>

#include "register_utils.hpp"
#include "interval.hpp"

namespace crab {

constexpr int STACK_BEGIN = 0;
constexpr int STACK_END = 512;
constexpr int CTX_BEGIN = 0;
constexpr int PACKET_BEGIN = 0;
constexpr int SHARED_BEGIN = 0;
constexpr int MAX_PACKET_SIZE = 0xffff;
constexpr int MAX_META_SIZE = 4098;


enum class types_t {
    T_UNSPEC,
    T_NUM,
    T_STACK,
    T_CTX,
    T_PACKET,
    T_SHARED,
    T_MAP,
    T_MAP_PROGRAMS
};


enum class region_t {
	R_CTX,
	R_STACK,
	R_PACKET,
	R_SHARED
};

inline std::string region_to_string(const region_t& r) noexcept;

enum class nullness_t { MAYBE_NULL, NOT_NULL, _NULL };

class packet_ptr_t {
    region_t m_r = region_t::R_PACKET;

  public:
    friend std::ostream& operator<<(std::ostream& o, const packet_ptr_t& p);
    // because we only represent one packet pointer, we can always return true/false
    // TODO: clean this up in future
    bool operator==(const packet_ptr_t&) const { return true; }
    bool operator!=(const packet_ptr_t&) const { return false; }
    [[nodiscard]] region_t get_region() const { return m_r; }
};

class ptr_with_off_t {
    region_t m_r;
    interval_t m_offset;
    // following fields are used for shared pointers, default values used for other pointer types
    int m_id = -1;
    nullness_t m_nullness = nullness_t::MAYBE_NULL;
    interval_t m_region_size = interval_t::top();

  public:
    ptr_with_off_t(region_t _r, interval_t _off) : m_r(_r), m_offset(_off) {}
    ptr_with_off_t(region_t _r, interval_t _off, int _id, nullness_t _nullness,
                   interval_t _region_sz)
        : m_r(_r), m_offset(_off), m_id(_id), m_nullness(_nullness), m_region_size(_region_sz) {}
    static ptr_with_off_t shared_region_ptr(int _offset, int _id = -1,
                                            nullness_t _nullness = nullness_t::MAYBE_NULL,
                                            interval_t _region_sz = interval_t::top()) {
        interval_t offset_interval{number_t{_offset}};
        return ptr_with_off_t(region_t::R_SHARED, offset_interval, _id, _nullness, _region_sz);
    }
    static ptr_with_off_t ctx_region_ptr(int _offset) {
        interval_t offset_interval{number_t{_offset}};
        return ptr_with_off_t(region_t::R_CTX, offset_interval);
    }
    static ptr_with_off_t stack_region_ptr(int _offset) {
        interval_t offset_interval{number_t{_offset}};
        return ptr_with_off_t(region_t::R_STACK, offset_interval);
    }
    ptr_with_off_t operator|(const ptr_with_off_t&) const;
    ptr_with_off_t widen(const ptr_with_off_t&) const;
    bool operator<=(const ptr_with_off_t&) const;
    [[nodiscard]] nullness_t get_nullness() const { return m_nullness; }
    void set_nullness(nullness_t);
    [[nodiscard]] int get_id() const { return m_id; }
    void set_id(int);
    [[nodiscard]] interval_t get_region_size() const { return m_region_size; }
    void set_region_size(interval_t);
    [[nodiscard]] interval_t get_offset() const { return m_offset; }
    void set_offset(interval_t);
    [[nodiscard]] region_t get_region() const { return m_r; }
    void set_region(region_t);
    void write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream& o, const ptr_with_off_t& p);
    bool operator==(const ptr_with_off_t&) const;
    bool operator!=(const ptr_with_off_t&) const;
};

class mapfd_t {
    interval_t m_mapfd;
    EbpfMapValueType m_value_type;

  public:
    mapfd_t() = default;
    mapfd_t operator|(const mapfd_t&) const;
    mapfd_t widen(const mapfd_t&) const;
    bool operator<=(const mapfd_t&) const;
    mapfd_t(interval_t mapfd, EbpfMapValueType val_type)
        : m_mapfd(mapfd), m_value_type(val_type) {}
    friend std::ostream& operator<<(std::ostream&, const mapfd_t&);
    bool operator==(const mapfd_t&) const;
    bool operator!=(const mapfd_t&) const;
    void write(std::ostream&) const;

    bool has_type_map_programs() const;
    [[nodiscard]] EbpfMapValueType get_value_type() const { return m_value_type; }
    [[nodiscard]] interval_t get_mapfd() const { return m_mapfd; }
};

using ptr_t = std::variant<packet_ptr_t, ptr_with_off_t>;

using ptr_or_mapfd_t = std::variant<ptr_with_off_t, packet_ptr_t, mapfd_t>;

inline bool is_ptr_type(std::optional<ptr_or_mapfd_t> ptr_or_mapfd) {
    return (ptr_or_mapfd && !std::holds_alternative<mapfd_t>(*ptr_or_mapfd));
}

inline bool is_mapfd_type(std::optional<ptr_or_mapfd_t> ptr_or_mapfd) {
    return (ptr_or_mapfd && std::holds_alternative<mapfd_t>(*ptr_or_mapfd));
}

inline bool same_region(const ptr_or_mapfd_t& ptr1, const ptr_or_mapfd_t& ptr2) {
    if (std::holds_alternative<packet_ptr_t>(ptr1) && std::holds_alternative<packet_ptr_t>(ptr2))
        return true;
    if (std::holds_alternative<ptr_with_off_t>(ptr1) && std::holds_alternative<ptr_with_off_t>(ptr2)) {
        auto p1 = std::get<ptr_with_off_t>(ptr1);
        auto p2 = std::get<ptr_with_off_t>(ptr2);
        auto r1 = p1.get_region();
        auto r2 = p2.get_region();
        if (r1 == r2 && r1 != region_t::R_SHARED) return true;
    }
    return false;
}

inline bool is_stack_ptr(std::optional<ptr_or_mapfd_t> ptr) {
    return (ptr && std::holds_alternative<ptr_with_off_t>(*ptr)
            && std::get<ptr_with_off_t>(*ptr).get_region() == region_t::R_STACK);
}

inline bool is_ctx_ptr(std::optional<ptr_or_mapfd_t> ptr) {
    return (ptr && std::holds_alternative<ptr_with_off_t>(*ptr)
            && std::get<ptr_with_off_t>(*ptr).get_region() == region_t::R_CTX);
}

inline bool is_packet_ptr(std::optional<ptr_or_mapfd_t> ptr) {
    return (ptr && std::holds_alternative<packet_ptr_t>(*ptr));
}

inline bool is_shared_ptr(std::optional<ptr_or_mapfd_t> ptr) {
    return (ptr && std::holds_alternative<ptr_with_off_t>(*ptr)
            && std::get<ptr_with_off_t>(*ptr).get_region() == region_t::R_SHARED);
}

using stack_cells_t = std::vector<std::pair<uint64_t, int>>;

} // namespace crab


namespace std {
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
}
