// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/type_ostream.hpp"

void print_ptr_type(std::ostream& o, const crab::ptr_or_mapfd_t& ptr, std::optional<crab::dist_t> d) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        o << std::get<crab::packet_ptr_t>(ptr) << "<" << (d ? *d : crab::dist_t{}) << ">";
    }
}

void print_number(std::ostream& o, crab::interval_t n) {
    o << "number";
    if (!n.is_top()) {
        if (auto n_singleton = n.singleton()) {
            o << "<" << *n_singleton << ">";
        }
        else {
            o << "<" << n << ">";
        }
    }
}

void print_ptr_or_mapfd_type(std::ostream& o, const crab::ptr_or_mapfd_t& ptr_or_mapfd, std::optional<crab::dist_t> d) {
    if (std::holds_alternative<crab::mapfd_t>(ptr_or_mapfd)) {
        o << std::get<crab::mapfd_t>(ptr_or_mapfd);
    }
    else {
        print_ptr_type(o, ptr_or_mapfd, d);
    }
}

void print_register(std::ostream& o, const Reg& r, const std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::dist_t> d, std::optional<crab::mock_interval_t> n) {
    o << r << " : ";
    if (p) {
        print_ptr_or_mapfd_type(o, *p, d);
    }
    else if (n) {
        print_number(o, n->to_interval());
    }
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