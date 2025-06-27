// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/type_ostream.hpp"

void print_region(std::ostream& o, crab::region_t region) {
    if (region == crab::region_t::R_STACK) {
        o << "stack";
    }
    else if (region == crab::region_t::R_CTX) {
        o << "ctx";
    }
    else if (region == crab::region_t::R_PACKET) {
        o << "packet";
    }
    else if (region == crab::region_t::R_SHARED) {
        o << "shared";
    }
}

void print_non_numeric_memory_cell(std::ostream& o, int start, int end,
                                   const crab::ptr_or_mapfd_t& ptr,
                                   std::optional<crab::refinement_t> d,
                                   std::shared_ptr<crab::slacks_t> slacks,
                                   crab::region_t region) {
    print_region(o, region);
    o << "[" << start << "-" << end << "] : ";
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            d->write(o, slacks, false); // pkt pointers are not un/signed
        }
        else {
            o << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_memory_cell(std::ostream& o, int start, int end, crab::refinement_t n,
        bool is_signed, std::shared_ptr<crab::slacks_t> slacks, crab::region_t region) {
    print_region(o, region);
    crab::interval_t i = n.get_interval_value(slacks);
    o << "[" << start << "-" << end << "] : ";
    if (i.is_bottom()) {
        o << "bottom";
        return;
    }
    if (i.is_top()) {
        if (is_signed) {
            o << "snumber<-oo, +oo>";
        }
        else {
            o << "unumber<-oo, +oo>";
        }
    }
    else {
        n.write(o, slacks, is_signed);
    }
}

void print_memory_cell(std::ostream& o, int start, int end,
                       std::optional<crab::ptr_or_mapfd_t> p,
                       std::optional<crab::refinement_t> d,
                       std::optional<crab::refinement_t> signed_numeric,
                       std::optional<crab::refinement_t> unsigned_numeric,
                       std::shared_ptr<crab::slacks_t> slacks, crab::region_t region) {
    if (signed_numeric) {
        print_numeric_memory_cell(o, start, end, *signed_numeric, true, slacks, region);
    }
    if (unsigned_numeric) {
        print_numeric_memory_cell(o, start, end, *unsigned_numeric, false, slacks, region);
    }
    else if (p) {
        print_non_numeric_memory_cell(o, start, end, *p, d, slacks, region);
    }
}

void print_non_numeric_register(std::ostream& o, crab::register_t r, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::refinement_t> d, std::shared_ptr<crab::slacks_t> slacks) {
    if (std::holds_alternative<crab::ptr_with_off_t>(ptr)) {
        o << r << " : " << std::get<crab::ptr_with_off_t>(ptr);
    }
    else if (std::holds_alternative<crab::packet_ptr_t>(ptr)) {
        if (d) {
            o << r << " : ";
            d->write(o, slacks, false); // pkt pointers are not un/signed
        }
        else {
            o << r << " : " << std::get<crab::packet_ptr_t>(ptr);
        }
    }
    else {
        o << r << " : " << std::get<crab::mapfd_t>(ptr);
    }
}

void print_numeric_register(std::ostream& o, crab::register_t r, crab::refinement_t n, bool is_signed,
        std::shared_ptr<crab::slacks_t> slacks) {
    crab::interval_t i = n.get_interval_value(slacks);
    if (i.is_bottom()) {
        o << r << " : bottom";
        return;
    }
    if (i.is_top()) {
        if (is_signed) {
            o << r << " : snumber<-oo, +oo>";
        }
        else {
            o << r << " : unumber<-oo, +oo>";
        }
    }
    else {
        o << r << " : ";
        n.write(o, slacks, is_signed);
    }
}

void print_register(std::ostream& o, Reg r,
                    std::optional<crab::ptr_or_mapfd_t> ptr_or_mapfd,
                    std::optional<crab::refinement_t> offset,
                    std::optional<crab::refinement_t> numeric,
                    bool is_signed,
                    std::shared_ptr<crab::slacks_t> slacks) {
    crab::register_t reg(r.v);
    if (ptr_or_mapfd) {
        print_non_numeric_register(o, reg, *ptr_or_mapfd, offset, slacks);
    }
    else if (numeric) {
        print_numeric_register(o, reg, *numeric, is_signed, slacks);
    }
    else {
        o << reg << " : unknown";
    }
}

inline std::string size_(int w) { return std::string("u") + std::to_string(w * 8); }

void print_annotated(std::ostream& o, const Call& call,
                     std::optional<crab::ptr_or_mapfd_t> ptr_or_mapfd,
                     std::optional<crab::refinement_t> offset,
                     std::optional<crab::refinement_t> signed_numeric,
                     std::optional<crab::refinement_t> unsigned_numeric,
                     std::shared_ptr<crab::slacks_t> slacks) {
    o << "  ";
    Reg r0 = Reg{R0_RETURN_VALUE};
    print_register(o, r0, ptr_or_mapfd, offset, signed_numeric, true, slacks);
    o << " = " << call.name << ":" << call.func << "(...)";
    if (unsigned_numeric) {
        o << "\n\t\t\t";
        print_register(o, r0, {}, {}, unsigned_numeric, false, slacks);
    }
    o << "\n";
}

void print_annotated(std::ostream& o, const Bin& b,
                     std::optional<crab::ptr_or_mapfd_t> ptr_or_mapfd,
                     std::optional<crab::refinement_t> offset,
                     std::optional<crab::refinement_t> signed_numeric,
                     std::optional<crab::refinement_t> unsigned_numeric,
                     std::shared_ptr<crab::slacks_t> slacks) {
    o << "  ";
    print_register(o, b.dst, ptr_or_mapfd, offset, signed_numeric, true, slacks);
    o << " " << b.op << "= " << b.v;
    if (unsigned_numeric) {
        o << "\n\t\t\t";
        print_register(o, b.dst, {}, {}, unsigned_numeric, false, slacks);
    }
    o << "\n";
}

void print_annotated(std::ostream& o, const LoadMapFd& u,
                     std::optional<crab::ptr_or_mapfd_t> ptr_or_mapfd) {
    o << "  ";
    print_register(o, u.dst, ptr_or_mapfd, {}, {}, false, nullptr);
    o << " = map_fd " << u.mapfd << "\n";
}

void print_annotated(std::ostream& o, const Mem& b,
                     std::optional<crab::ptr_or_mapfd_t> ptr_or_mapfd,
                     std::optional<crab::refinement_t> offset,
                     std::optional<crab::refinement_t> signed_numeric,
                     std::optional<crab::refinement_t> unsigned_numeric,
                     std::shared_ptr<crab::slacks_t> slacks) {
    // This is load operation
    o << "  ";
    print_register(o, std::get<Reg>(b.value), ptr_or_mapfd, offset, signed_numeric, true, slacks);
    o << " = ";
    std::string sign = b.access.offset < 0 ? " - " : " + ";
    int offset_int = std::abs(b.access.offset);
    o << "*(" << size_(b.access.width) << " *)";
    o << "(" << b.access.basereg << sign << offset_int << ")";
    if (unsigned_numeric) {
        o << "\n\t\t\t";
        print_register(o, std::get<Reg>(b.value), {}, {}, unsigned_numeric, false, slacks);
    }
    o << "\n";
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

void print_annotated(std::ostream& o, const Un& b,
                     std::optional<crab::refinement_t> signed_numeric,
                     std::optional<crab::refinement_t> unsigned_numeric,
                     std::shared_ptr<crab::slacks_t> slacks) {
    o << "  ";
    print_register(o, b.dst, {}, {}, signed_numeric, true, slacks);
    o << " = " << op(b.op) << " " << b.dst;
    if (unsigned_numeric) {
        o << "\n\t\t\t";
        print_register(o, b.dst, {}, {}, unsigned_numeric, false, slacks);
    }
    o << "\n";
}

void print_bb(std::ostream& o, const basic_block_t& bb) {
    o << bb << "\n";
}

void print_instr(std::ostream& o, const Instruction& i) {
    std::visit([&](auto&& i) { o << "  " << i << "\n"; }, i);
}
