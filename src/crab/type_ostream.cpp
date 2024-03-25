// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/type_ostream.hpp"

void print_non_numeric_memory_cell(std::ostream& o, int start, int end,
        const crab::ptr_or_mapfd_t& ptr, std::optional<crab::refinement_t> d) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << "[" << start << "-" << end << "] : " << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            o << "[" << start << "-" << end << "] : " << *d;
        }
        else {
            o << "[" << start << "-" << end << "] : " << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << "[" << start << "-" << end << "] : " << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_memory_cell(std::ostream& o, int start, int end, crab::interval_t n,
        bool is_signed) {
    if (n.is_top()) {
        if (is_signed) {
            o << "[" << start << "-" << end << "] : snumber";
        }
        else {
            o << "[" << start << "-" << end << "] : unumber";
        }
    }
    else {
        if (auto n_singleton = n.singleton()) {
            if (is_signed) {
                o << "[" << start << "-" << end << "] : snumber<" << *n_singleton << ">";
            }
            else {
                o << "[" << start << "-" << end << "] : unumber<" << *n_singleton << ">";
            }
        }
        else {
            if (is_signed) {
                o << "[" << start << "-" << end << "] : snumber<" << n << ">";
            }
            else {
                o << "[" << start << "-" << end << "] : unumber<" << n << ">";
            }
        }
    }
}

void print_memory_cell(std::ostream& o, int start, int end,
        const std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::refinement_t> d
        , std::optional<crab::mock_interval_t> signed_interval,
        std::optional<crab::mock_interval_t> unsigned_interval) {
    if (signed_interval) {
        print_numeric_memory_cell(o, start, end, signed_interval->to_interval(), true);
    }
    if (unsigned_interval) {
        print_numeric_memory_cell(o, start, end, unsigned_interval->to_interval(), false);
    }
    else if (p) print_non_numeric_memory_cell(o, start, end, *p, d);
}

void print_non_numeric_register(std::ostream& o, Reg r, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::refinement_t> d) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << r << " : " << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            o << r << " : " << *d;
        }
        else {
            o << r << " : " << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << r << " : " << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_register(std::ostream& o, Reg r, crab::interval_t n, bool is_signed) {
    if (n.is_top()) {
        if (is_signed) {
            o << r << " : snumber";
        }
        else {
            o << r << " : unumber";
        }
    }
    else {
        if (auto n_singleton = n.singleton()) {
            if (is_signed) {
                o << r << " : snumber<" << *n_singleton << ">";
            }
            else {
                o << r << " : unumber<" << *n_singleton << ">";
            }
        }
        else {
            if (is_signed) {
                o << r << " : snumber<" << n << ">";
            }
            else {
                o << r << " : unumber<" << n << ">";
            }
        }
    }
}

void print_register(std::ostream& o, Reg r, const std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::refinement_t> d, std::optional<crab::mock_interval_t> interval,
        bool is_signed) {
    if (interval) print_numeric_register(o, r, interval->to_interval(), is_signed);
    else if (p) print_non_numeric_register(o, r, *p, d);
}

inline std::string size_(int w) { return std::string("u") + std::to_string(w * 8); }

void print_annotated(std::ostream& o, const Call& call, std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::refinement_t>& d,
        std::optional<crab::mock_interval_t>& n, bool is_signed) {
    o << "  ";
    print_register(o, Reg{(uint8_t)R0_RETURN_VALUE}, p, d, n, is_signed);
    o << " = " << call.name << ":" << call.func << "(...)\n";
    if (d) {
        o << "  " << *d << "\n";
    }
}

void print_annotated(std::ostream& o, const Bin& b, std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::refinement_t>& d, std::optional<crab::mock_interval_t>& n,
        bool is_signed) {
    o << "  ";
    print_register(o, b.dst, p, d, n, is_signed);
    o << " " << b.op << "= " << b.v << "\n";
    if (d) {
        o << "  " << *d << "\n";
    }
}

void print_annotated(std::ostream& o, const LoadMapFd& u, std::optional<crab::ptr_or_mapfd_t>& p) {
    o << "  ";
    print_register(o, u.dst, p, std::nullopt, std::nullopt, false);
    o << " = map_fd " << u.mapfd << "\n";
}

void print_annotated(std::ostream& o, const Mem& b, std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::refinement_t>& d, std::optional<crab::mock_interval_t>& n, bool is_signed) {
    o << "  ";
    print_register(o, std::get<Reg>(b.value), p, d, n, is_signed);
    o << " = ";
    std::string sign = b.access.offset < 0 ? " - " : " + ";
    int offset = std::abs(b.access.offset);
    o << "*(" << size_(b.access.width) << " *)";
    o << "(" << b.access.basereg << sign << offset << ")\n";
    if (d) {
        o << "  " << *d << "\n";
    }
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

void print_annotated(std::ostream& o, const Un& b, std::optional<crab::mock_interval_t>& n,
        bool is_signed) {
    o << "  ";
    print_register(o, b.dst, std::nullopt, std::nullopt, n, is_signed);
    o << " = " << op(b.op) << " " << b.dst << "\n";
}
