// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "asm_syntax.hpp"

namespace crab {

// 11 registers for the eBPF ISA, 1 register used in offset domain to represent r_begin,
// and 1 pseudo-register for the atomic operations
constexpr uint8_t NUM_REGISTERS = 13;

// Represents a register, e.g., r0, r1, ..., r_begin, etc.
class register_t {
    uint8_t m_val;

  public:
    register_t(uint8_t val) : m_val(val) {}
    register_t() : m_val(NUM_REGISTERS+1) {}
    std::size_t hash() const;
    uint8_t operator()() const { return m_val; }
    operator int() const { return static_cast<int>(m_val); }
    friend std::ostream& operator<<(std::ostream& o, const register_t& p);
};


// Represents a location in the eBPF program, i.e., a basic block label and a line number
// it is shown as 'line_num in bb_label'
class location_t {
  public:
    label_t m_bb_label;
    uint32_t m_line_num; // uint32_t is used to represent the line number in the eBPF program
                       // which should be sufficient for ebpf
    location_t(label_t _label, uint32_t _line_num) : m_bb_label(_label), m_line_num(_line_num) {}
    location_t() : m_bb_label(-2, -2), m_line_num(0) {} // not a valid location
    static location_t top() { return location_t(); }
    bool operator==(const location_t& other) const;
    friend std::ostream& operator<<(std::ostream& o, const location_t& loc);
    void write(std::ostream& o) const;
};

// Represents a register and a location in the eBPF program
// it is shown as 'r0@line_num in bb_label'
class register_location_t {
    register_t m_reg;
    location_t m_loc;

  public:
    register_location_t(register_t _r, location_t _loc) : m_reg(_r), m_loc(_loc) {}
    std::size_t hash() const;
    bool operator==(const register_location_t& other) const;
    friend std::ostream& operator<<(std::ostream& o, const register_location_t& reg);
    void write(std::ostream& o) const;
};

using live_registers_t = std::array<std::shared_ptr<register_location_t>, NUM_REGISTERS>;

} // namespace crab


namespace std {
    template <>
    struct hash<crab::register_location_t> {
        std::size_t operator()(const crab::register_location_t& reg) const {
            return reg.hash();
        }
    };
} // namespace std
