// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "types.hpp"

namespace crab {

inline std::string region_to_string(const region_t& r) noexcept {
    switch (r) {
        case region_t::R_CTX:
            return "ctx_p";
        case region_t::R_STACK:
            return "stack_p";
        case region_t::R_SHARED:
            return "shared_p";
        case region_t::R_PACKET:
            return "packet_p";
        default:
            __builtin_unreachable();
    }
}

bool ptr_with_off_t::operator==(const ptr_with_off_t& other) const {
    return (m_r == other.m_r && m_offset == other.m_offset
            && m_region_size == other.m_region_size);
}

bool ptr_with_off_t::operator!=(const ptr_with_off_t& other) const {
    return !(*this == other);
}

bool mapfd_t::operator==(const mapfd_t& other) const {
    return (m_mapfd == other.m_mapfd);
}

bool mapfd_t::operator<=(const mapfd_t& other) const {
    bool map_type = m_value_type == other.m_value_type
        || other.m_value_type == EbpfMapValueType::ANY;
    return (m_mapfd <= other.m_mapfd && map_type);
}

mapfd_t mapfd_t::operator|(const mapfd_t& other) const {
    auto value_type = m_value_type == other.m_value_type ? m_value_type : EbpfMapValueType::ANY;
    return mapfd_t(m_mapfd | other.m_mapfd, value_type);
}

mapfd_t mapfd_t::widen(const mapfd_t& other) const {
    auto value_type = m_value_type == other.m_value_type ? m_value_type : EbpfMapValueType::ANY;
    return mapfd_t(m_mapfd.widen(other.m_mapfd), value_type);
}

void ptr_with_off_t::set_nullness(nullness_t n) { m_nullness = n; }

void ptr_with_off_t::set_id(int id) { m_id = id; }

void ptr_with_off_t::set_offset(interval_t off) { m_offset = off; }

void ptr_with_off_t::set_region_size(interval_t region_sz) { m_region_size = region_sz; }

void ptr_with_off_t::set_region(region_t r) { m_r = r; }

bool ptr_with_off_t::operator<=(const ptr_with_off_t& other) const {
    bool nullness = m_nullness == other.m_nullness || other.m_nullness == nullness_t::MAYBE_NULL;
    return (m_r == other.m_r && m_offset <= other.m_offset
            && m_region_size <= other.m_region_size && nullness);
}

ptr_with_off_t ptr_with_off_t::operator|(const ptr_with_off_t& other) const {
    auto nullness = m_nullness == other.m_nullness ? m_nullness : nullness_t::MAYBE_NULL;
    return ptr_with_off_t(m_r, m_offset | other.m_offset, -1, nullness,
                          m_region_size | other.m_region_size);
}

ptr_with_off_t ptr_with_off_t::widen(const ptr_with_off_t& other) const {
    auto nullness = m_nullness == other.m_nullness ? m_nullness : nullness_t::MAYBE_NULL;
    return ptr_with_off_t(m_r, m_offset.widen(other.m_offset), -1, nullness,
                          m_region_size.widen(other.m_region_size));
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
    if (auto mapfd_singleton = m_mapfd.singleton()) {
        o << *mapfd_singleton;
    }
    else {
        o << m_mapfd;
    }
}

void ptr_with_off_t::write(std::ostream& o) const {
    o << region_to_string(m_r);
    if (!m_offset.is_top()) {
        o << "<";
        if (auto off_singleton = m_offset.singleton()) {
            o << *off_singleton;
        }
        else {
            o << m_offset;
        }
        if (m_region_size.lb() >= number_t{0}) {
            o << ",";
            if (auto rs_singleton = m_region_size.singleton()) {
                o << *rs_singleton;
            }
            else {
                o << m_region_size;
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
