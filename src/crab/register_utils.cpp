// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "register_utils.hpp"
#include <sstream>

namespace crab {

std::size_t register_t::hash() const {
    return std::hash<uint8_t>{}(m_val);
}

std::ostream& operator<<(std::ostream& o, const register_t& p) {
    if (p == register_t{R12_PKT_BEGIN}) {
        o << "r_offset";
    }
    else {
        o << "r" << (int)p;
    }
    return o;
}

bool location_t::operator==(const location_t& other) const {
    return m_bb_label == other.m_bb_label && m_line_num == other.m_line_num;
}

void location_t::write(std::ostream& o) const {
    o << "instr#" << m_line_num << " in " << m_bb_label << " ";
}

std::string location_t::to_string() const {
    std::ostringstream oss;
    write(oss);
    return oss.str();
}

std::ostream& operator<<(std::ostream& o, const location_t& loc) {
    loc.write(o);
    return o;
}

std::size_t register_location_t::hash() const {
    // Similar to boost::hash_combine
    using std::hash;

    std::size_t seed = m_reg.hash();
    seed ^= hash<uint32_t>()(m_loc.m_line_num) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash<int>()(m_loc.m_bb_label.from) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash<int>()(m_loc.m_bb_label.to) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

bool register_location_t::operator==(const register_location_t& other) const {
    return m_reg == other.m_reg && m_loc == other.m_loc;
}

void register_location_t::write(std::ostream& o) const {
    o << m_reg << "@" << m_loc;
}

std::ostream& operator<<(std::ostream& o, const register_location_t& reg) {
    reg.write(o);
    return o;
}

} // namespace crab
