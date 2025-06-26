// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/type_ostream.hpp"

void print_non_numeric_memory_cell(std::ostream& o, int start, int end,
        const crab::ptr_or_mapfd_t& ptr, std::optional<crab::refinement_t> d,
                                   std::shared_ptr<crab::slacks_t> slacks) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << "[" << start << "-" << end << "] : " << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            o << "[" << start << "-" << end << "] : ";
            d->write(o, slacks);
        }
        else {
            o << "[" << start << "-" << end << "] : " << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << "[" << start << "-" << end << "] : " << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_memory_cell(std::ostream& o, int start, int end, crab::refinement_t n,
        bool is_signed, std::shared_ptr<crab::slacks_t> slacks) {
    crab::interval_t i = n.get_interval_value(slacks);
    if (i.is_bottom()) {
        o << "[" << start << "-" << end << "] : bottom";
        return;
    }
    if (i.is_top()) {
        if (is_signed) {
            o << "[" << start << "-" << end << "] : snumber";
        }
        else {
            o << "[" << start << "-" << end << "] : unumber";
        }
    }
    else {
        // TODO: differentiate between signed and unsigned
        o << "[" << start << "-" << end << "] : ";
        n.write(o, slacks);
    }
}

void print_memory_cell(std::ostream& o, int start, int end,
        const std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::refinement_t> d
        , std::optional<crab::refinement_t> signed_numeric,
        std::optional<crab::refinement_t> unsigned_numeric,
                       std::shared_ptr<crab::slacks_t> slacks) {
    if (signed_numeric) {
        print_numeric_memory_cell(o, start, end, *signed_numeric, true, slacks);
    }
    if (unsigned_numeric) {
        print_numeric_memory_cell(o, start, end, *unsigned_numeric, false, slacks);
    }
    else if (p) {
        print_non_numeric_memory_cell(o, start, end, *p, d, slacks);
    }
}

void print_non_numeric_register(std::ostream& o, Reg r, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::refinement_t> d, std::shared_ptr<crab::slacks_t> slacks) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << r << " : " << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            o << r << " : ";
            d->write(o, slacks);
        }
        else {
            o << r << " : " << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << r << " : " << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_register(std::ostream& o, Reg r, crab::refinement_t n, bool is_signed,
        std::shared_ptr<crab::slacks_t> slacks) {
    crab::interval_t i = n.get_interval_value(slacks);
    if (i.is_bottom()) {
        o << r << " : bottom";
        return;
    }
    if (i.is_top()) {
        if (is_signed) {
            o << r << " : snumber";
        }
        else {
            o << r << " : unumber";
        }
    }
    else {
        o << r << " : ";
        n.write(o, slacks);
    }
}

void print_register(std::ostream& o, Reg r, const std::optional<crab::ptr_or_mapfd_t>& p,
        const std::optional<crab::refinement_t>& d, const std::optional<crab::refinement_t>& numeric,
        bool is_signed, std::shared_ptr<crab::slacks_t> slacks) {
    if (numeric) print_numeric_register(o, r, *numeric, is_signed, slacks);
    else if (p) print_non_numeric_register(o, r, *p, d, slacks);
    else o << r << " : unknown";
}

inline std::string size_(int w) { return std::string("u") + std::to_string(w * 8); }

void print_annotated(std::ostream& o, const Call& call, std::optional<crab::ptr_or_mapfd_t>& p,
                     std::optional<crab::refinement_t>& d,
                     std::optional<crab::refinement_t>& n, bool is_signed,
                     std::shared_ptr<crab::slacks_t> slacks) {
    o << "  ";
    print_register(o, Reg{(uint8_t)R0_RETURN_VALUE}, p, d, n, is_signed, slacks);
    o << " = " << call.name << ":" << call.func << "(...)\n";
}

void print_annotated(std::ostream& o, const Bin& b, std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::refinement_t>& d, std::optional<crab::refinement_t>& n,
        bool is_signed, std::shared_ptr<crab::slacks_t> slacks) {
    o << "  ";
    print_register(o, b.dst, p, d, n, is_signed, slacks);
    o << " " << b.op << "= " << b.v << "\n";
}

void print_annotated(std::ostream& o, const LoadMapFd& u, std::optional<crab::ptr_or_mapfd_t>& p) {
    o << "  ";
    print_register(o, u.dst, p, std::nullopt, std::nullopt, false, nullptr);
    o << " = map_fd " << u.mapfd << "\n";
}

void print_annotated(std::ostream& o, const Mem& b, std::optional<crab::ptr_or_mapfd_t>& p,
        std::optional<crab::refinement_t>& d, std::optional<crab::refinement_t>& n, bool is_signed,
                     std::shared_ptr<crab::slacks_t> slacks) {
    o << "  ";
    print_register(o, std::get<Reg>(b.value), p, d, n, is_signed, slacks);
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

void print_annotated(std::ostream& o, const Un& b, std::optional<crab::refinement_t>& n,
        bool is_signed, std::shared_ptr<crab::slacks_t> slacks) {
    o << "  ";
    print_register(o, b.dst, std::nullopt, std::nullopt, n, is_signed, slacks);
    o << " = " << op(b.op) << " " << b.dst << "\n";
}

void print_bb(std::ostream& o, const basic_block_t& bb) {
    o << bb << "\n";
}

void print_instr(std::ostream& o, const Instruction& i) {
    std::visit([&](auto&& i) { o << "  " << i << "\n"; }, i);
}
