// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/common.hpp"

/*
namespace std {
    crab::ptr_t get_ptr(const crab::ptr_or_mapfd_t& t) {
    return std::visit( overloaded
               {
                   []( const crab::ptr_with_off_t& x ){ return crab::ptr_t{x};},
                   []( const crab::packet_ptr_t& x ){ return crab::ptr_t{x};},
                   []( auto& ) { return crab::ptr_t{};}
                }, t
            );
    }
} // namespace std
*/

namespace crab {

inline std::string get_reg_ptr(const region_t& r) noexcept {
    switch (r) {
        case region_t::T_CTX:
            return "ctx_p";
        case region_t::T_STACK:
            return "stack_p";
        case region_t::T_SHARED:
            return "shared_p";
        default:
            return "packet_p";
    }
    __builtin_unreachable();
}

bool mock_interval_t::operator==(const mock_interval_t& other) const {
    return (to_interval() == other.to_interval());
}

bool ptr_with_off_t::operator==(const ptr_with_off_t& other) const {
    return (m_r == other.m_r && m_offset == other.m_offset
            && m_region_size == other.m_region_size);
}

bool ptr_with_off_t::operator!=(const ptr_with_off_t& other) const {
    return !(*this == other);
}

bool packet_ptr_t::operator==(const packet_ptr_t& other) const {
    return (m_r == other.m_r);
}

bool packet_ptr_t::operator!=(const packet_ptr_t& other) const {
    return !(*this == other);
}

bool mapfd_t::operator==(const mapfd_t& other) const {
    return (m_mapfd == other.m_mapfd);
}

bool mapfd_t::operator<=(const mapfd_t& other) const {
    bool map_type = m_value_type == other.m_value_type
        || other.m_value_type == EbpfMapValueType::ANY;
    return (m_mapfd.to_interval() <= other.m_mapfd.to_interval() && map_type);
}

mapfd_t mapfd_t::operator|(const mapfd_t& other) const {
    auto value_type = m_value_type == other.m_value_type ? m_value_type : EbpfMapValueType::ANY;
    const auto& mock_i = mock_interval_t(m_mapfd.to_interval() | other.m_mapfd.to_interval());
    return mapfd_t(std::move(mock_i), value_type);
}

mapfd_t mapfd_t::widen(const mapfd_t& other) const {
    auto value_type = m_value_type == other.m_value_type ? m_value_type : EbpfMapValueType::ANY;
    const auto& mock_i = mock_interval_t(m_mapfd.to_interval().widen(other.m_mapfd.to_interval()));
    return mapfd_t(std::move(mock_i), value_type);
}

void ptr_with_off_t::set_nullness(nullness_t n) { m_nullness = n; }

void ptr_with_off_t::set_id(int id) { m_id = id; }

void ptr_with_off_t::set_offset(mock_interval_t off) { m_offset = off; }

void ptr_with_off_t::set_region_size(mock_interval_t region_sz) { m_region_size = region_sz; }

void ptr_with_off_t::set_region(region_t r) { m_r = r; }

bool ptr_with_off_t::operator<=(const ptr_with_off_t& other) const {
    bool nullness = m_nullness == other.m_nullness || other.m_nullness == nullness_t::MAYBE_NULL;
    return (m_r == other.m_r && m_offset.to_interval() <= other.m_offset.to_interval()
            && m_region_size.to_interval() <= other.m_region_size.to_interval() && nullness);
}

ptr_with_off_t ptr_with_off_t::operator|(const ptr_with_off_t& other) const {
    auto&& mock_o = mock_interval_t(m_offset.to_interval() | other.m_offset.to_interval());
    auto&& mock_r_s = mock_interval_t(
            m_region_size.to_interval() | other.m_region_size.to_interval());
    auto&& nullness = m_nullness == other.m_nullness ? m_nullness : nullness_t::MAYBE_NULL;
    return ptr_with_off_t(m_r, -1, std::move(mock_o), std::move(nullness), std::move(mock_r_s));
}

ptr_with_off_t ptr_with_off_t::widen(const ptr_with_off_t& other) const {
    auto&& mock_o = mock_interval_t(m_offset.to_interval().widen(other.m_offset.to_interval()));
    auto&& mock_r_s = mock_interval_t(
            m_region_size.to_interval().widen(other.m_region_size.to_interval()));
    auto&& nullness = m_nullness == other.m_nullness ? m_nullness : nullness_t::MAYBE_NULL;
    return ptr_with_off_t(m_r, -1, std::move(mock_o), std::move(nullness), std::move(mock_r_s));
}

std::ostream& operator<<(std::ostream& o, const mapfd_t& m) {
    m.write(o);
    return o;
}

bool mapfd_t::has_type_map_programs() const {
    return (m_value_type == EbpfMapValueType::PROGRAM);
}

void mapfd_t::write(std::ostream& o) const {
    if (has_type_map_programs()) {
        o << "map_fd_programs ";
    }
    else {
        o << "map_fd ";
    }
    auto mapfd = m_mapfd.to_interval();
    if (auto mapfd_singleton = mapfd.singleton()) {
        o << *mapfd_singleton;
    }
    else {
        o << mapfd;
    }
}

void reg_with_loc_t::write(std::ostream& o) const {
    o << "r" << static_cast<unsigned int>(m_reg) << "@" << m_loc->second << " in " << m_loc->first << " ";
}

std::ostream& operator<<(std::ostream& o, const reg_with_loc_t& reg) {
    reg.write(o);
    return o;
}

bool reg_with_loc_t::operator==(const reg_with_loc_t& other) const {
    return (m_reg == other.m_reg && m_loc == other.m_loc);
}

std::size_t reg_with_loc_t::hash() const {
    // Similar to boost::hash_combine
    using std::hash;

    std::size_t seed = hash<register_t>()(m_reg);
    seed ^= hash<int>()(m_loc->first.from) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash<int>()(m_loc->first.to) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash<int>()(m_loc->second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

void ptr_with_off_t::write(std::ostream& o) const {
    o << get_reg_ptr(m_r);
    auto offset = m_offset.to_interval();
    auto region_size = m_region_size.to_interval();
    if (!offset.is_top()) {
        o << "<";
        if (auto off_singleton = offset.singleton()) {
            o << *off_singleton;
        }
        else {
            o << offset;
        }
        if (region_size.lb() >= number_t{0}) {
            o << ",";
            if (auto rs_singleton = region_size.singleton()) {
                o << *rs_singleton;
            }
            else {
                o << region_size;
            }
        }
        o << ">";
    }
}

std::ostream& operator<<(std::ostream& o, const ptr_with_off_t& p) {
    p.write(o);
    return o;
}

std::ostream& operator<<(std::ostream& o, const packet_ptr_t& p) {
    o << "packet_p";
    return o;
}

} // namespace crab
