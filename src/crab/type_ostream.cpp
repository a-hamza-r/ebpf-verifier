// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/type_ostream.hpp"

void print_non_numeric_memory_cell(std::ostream& o, int start, int end,
        const crab::ptr_or_mapfd_t& ptr, std::optional<crab::dist_t> d) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << "[" << start << "-" << end << "] : " << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            o << "[" << start << "-" << end << "] : " <<
                std::get<crab::packet_ptr_t>(ptr) << "<" << *d << ">";
        }
        else {
            o << "[" << start << "-" << end << "] : " << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << "[" << start << "-" << end << "] : " << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_memory_cell(std::ostream& o, int start, int end, crab::interval_t n) {
    if (n.is_top()) {
        o << "[" << start << "-" << end << "] : number";
    }
    else {
        if (auto n_singleton = n.singleton()) {
            o << "[" << start << "-" << end << "] : number<" << *n_singleton << ">";
        }
        else {
            o << "[" << start << "-" << end << "] : number<" << n << ">";
        }
    }
}

void print_memory_cell(std::ostream& o, int start, int end,
        const std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::dist_t> d,
        std::optional<crab::mock_interval_t> n) {
    if (n) print_numeric_memory_cell(o, start, end, n->to_interval());
    else if (p) print_non_numeric_memory_cell(o, start, end, *p, d);
}

void print_non_numeric_register(std::ostream& o, Reg r, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::dist_t> d) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << r << " : " << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            o << r << " : " << std::get<crab::packet_ptr_t>(ptr) << "<" << *d << ">";
        }
        else {
            o << r << " : " << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << r << " : " << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_register(std::ostream& o, Reg r, crab::interval_t n) {
    if (n.is_top()) {
        o << r << " : number";
    }
    else {
        if (auto n_singleton = n.singleton()) {
            o << r << " : number<" << *n_singleton << ">";
        }
        else {
            o << r << " : number<" << n << ">";
        }
    }
}

void print_register(std::ostream& o, Reg r, const std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::dist_t> d, std::optional<crab::mock_interval_t> n) {
    if (n) print_numeric_register(o, r, n->to_interval());
    else if (p) print_non_numeric_register(o, r, *p, d);
}

inline std::string size_(int w) { return std::string("u") + std::to_string(w * 8); }

void print_annotated(std::ostream& o, const Call& call, std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::mock_interval_t>& n) {
    o << "  ";
    print_register(o, Reg{(uint8_t)R0_RETURN_VALUE}, p, std::nullopt, n);
    o << " = " << call.name << ":" << call.func << "(...)\n";
}

void print_annotated(std::ostream& o, const Bin& b, std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::dist_t>& d, std::optional<crab::mock_interval_t>& n) {
    o << "  ";
    print_register(o, b.dst, p, d, n);
    o << " " << b.op << "= " << b.v << "\n";
}

void print_annotated(std::ostream& o, const LoadMapFd& u, std::optional<crab::ptr_or_mapfd_t>& p) {
    o << "  ";
    print_register(o, u.dst, p, std::nullopt, std::nullopt);
    o << " = map_fd " << u.mapfd << "\n";
}

void print_annotated(std::ostream& o, const Mem& b, std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::dist_t>& d, std::optional<crab::mock_interval_t>& n) {
    o << "  ";
    print_register(o, std::get<Reg>(b.value), p, d, n);
    o << " = ";
    std::string sign = b.access.offset < 0 ? " - " : " + ";
    int offset = std::abs(b.access.offset);
    o << "*(" << size_(b.access.width) << " *)";
    o << "(" << b.access.basereg << sign << offset << ")\n";
}

std::string op(Un::Op op) {
    switch (op) {
        case Un::Op::NEG:
            return "-";
        case Un::Op::BE16:
            return "be16";
        case Un::Op::BE32:
            return "be32";
        case Un::Op::BE64:
            return "be64";
        case Un::Op::LE16:
            return "le16";
        case Un::Op::LE32:
            return "le32";
        case Un::Op::LE64:
            return "le64";
        default:
            return "unknown";
    }
}

void print_annotated(std::ostream& o, const Un& b, std::optional<crab::mock_interval_t>& n) {
    o << "  ";
    print_register(o, b.dst, std::nullopt, std::nullopt, n);
    o << " = " << op(b.op) << " " << b.dst << "\n";
}
