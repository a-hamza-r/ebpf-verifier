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

bool mapfd_t::operator==(const mapfd_t& other) const {
    return (m_mapfd == other.m_mapfd);
}

mapfd_t mapfd_t::operator|(const mapfd_t& other) const {
    auto value_type = m_value_type == other.m_value_type ? m_value_type : EbpfMapValueType::ANY;
    const auto& mock_i = mock_interval_t(m_mapfd.to_interval() | other.m_mapfd.to_interval());
    return mapfd_t(std::move(mock_i), value_type);
}

void ptr_with_off_t::set_nullness(nullness_t n) { m_nullness = n; }

void ptr_with_off_t::set_id(int id) { m_id = id; }

void ptr_with_off_t::set_offset(mock_interval_t off) { m_offset = off; }

void ptr_with_off_t::set_region_size(mock_interval_t region_sz) { m_region_size = region_sz; }

void ptr_with_off_t::set_region(region_t r) { m_r = r; }

ptr_with_off_t ptr_with_off_t::operator|(const ptr_with_off_t& other) const {
    auto&& mock_o = mock_interval_t(m_offset.to_interval() | other.m_offset.to_interval());
    auto&& mock_r_s = mock_interval_t(
            m_region_size.to_interval() | other.m_region_size.to_interval());
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

void ptr_with_off_t::write(std::ostream& o) const {
    o << region_to_string(m_r);
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
