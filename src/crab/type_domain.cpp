// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include <unordered_map>

#include "crab/type_domain.hpp"

namespace std {
    static ptr_t get_ptr(const ptr_or_mapfd_t& t) {
    return std::visit( overloaded
               {
                   []( const ptr_with_off_t& x ){ return ptr_t{x};},
                   []( const ptr_no_off_t& x ){ return ptr_t{x};},
                   []( auto& ) { return ptr_t{};}
                }, t
            );
    }
}

static std::string size(int w) { return std::string("u") + std::to_string(w * 8); }

static void print_ptr_no_off_type(const ptr_no_off_t& ptr, std::optional<dist_t>& dist) {
    std::cout << ptr;
    if (dist) {
        std::cout << "<" << dist.value() << ">";
    }
}

static void print_ptr_type(const ptr_t& ptr, std::optional<dist_t>& dist) {
    if (std::holds_alternative<ptr_with_off_t>(ptr)) {
        ptr_with_off_t ptr_with_off = std::get<ptr_with_off_t>(ptr);
        std::cout << ptr_with_off;
    }
    else {
        ptr_no_off_t ptr_no_off = std::get<ptr_no_off_t>(ptr);
        print_ptr_no_off_type(ptr_no_off, dist);
    }
}

static void print_ptr_or_mapfd_type(const ptr_or_mapfd_t& ptr_or_mapfd, std::optional<dist_t>& d) {
    if (std::holds_alternative<mapfd_t>(ptr_or_mapfd)) {
        std::cout << std::get<mapfd_t>(ptr_or_mapfd);
    }
    else {
        auto ptr = get_ptr(ptr_or_mapfd);
        print_ptr_type(ptr, d);
    }
}

static void print_number(interval_t num) {
    std::cout << "number";
    if (!num.is_top()) {
        std::cout << "<" << num << ">";
    }
}

static void print_register(Reg r, std::optional<ptr_or_mapfd_t>& p, std::optional<dist_t>& d,
        std::optional<interval_t> n) {
    std::cout << r << " : ";
    if (p) {
        print_ptr_or_mapfd_type(p.value(), d);
    }
    if (n) {
        print_number(n.value());
    }
}

static void print_annotated(std::ostream& o, const Call& call, std::optional<ptr_or_mapfd_t>& p,
        std::optional<dist_t>& d, std::optional<interval_t>& n) {
    o << "  ";
    print_register(Reg{(uint8_t)R0_RETURN_VALUE}, p, d, n);
    o << " = " << call.name << ":" << call.func << "(...)\n";
}

static void print_annotated(std::ostream& o, const Bin& b, std::optional<ptr_or_mapfd_t>& p,
        std::optional<dist_t>& d, std::optional<interval_t>& n) {
    o << "  ";
    print_register(b.dst, p, d, n);
    o << " " << b.op << "= " << b.v << "\n";
}

static void print_annotated(std::ostream& o, const LoadMapFd& u, std::optional<ptr_or_mapfd_t>& p) {
    o << "  ";
    std::optional<dist_t> d;
    std::optional<interval_t> n;
    print_register(u.dst, p, d, n);
    o << " = map_fd " << u.mapfd << "\n";
}

static void print_annotated(std::ostream& o, const Mem& b, std::optional<ptr_or_mapfd_t>& p,
        std::optional<dist_t>& d, std::optional<interval_t>& n) {
    o << "  ";
    print_register(std::get<Reg>(b.value), p, d, n);
    o << " = ";
    std::string sign = b.access.offset < 0 ? " - " : " + ";
    int offset = std::abs(b.access.offset);
    o << "*(" << size(b.access.width) << " *)";
    o << "(" << b.access.basereg << sign << offset << ")\n";
}

bool type_domain_t::is_bottom() const {
    return m_is_bottom;
}

bool type_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_region.is_top() && m_offset.is_top() && m_interval.is_top());
}

type_domain_t type_domain_t::bottom() {
    type_domain_t typ;
    typ.set_to_bottom();
    return typ;
}

void type_domain_t::set_to_bottom() {
    m_is_bottom = true;
}

void type_domain_t::set_to_top() {
    m_region.set_to_top();
    m_offset.set_to_top();
    m_interval.set_to_top();
}

bool type_domain_t::operator<=(const type_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return true;
}

void type_domain_t::operator|=(const type_domain_t& abs) {
    type_domain_t tmp{abs};
    operator|=(std::move(tmp));
}

void type_domain_t::operator|=(type_domain_t&& abs) {
    if (is_bottom()) {
        *this = abs;
        return;
    }
    *this = *this | std::move(abs);
}

type_domain_t type_domain_t::operator|(const type_domain_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return type_domain_t(m_region | other.m_region, m_offset | other.m_offset,
            m_interval | other.m_interval);
}

type_domain_t type_domain_t::operator|(type_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return type_domain_t(m_region | std::move(other.m_region),
            m_offset | std::move(other.m_offset), m_interval | std::move(other.m_interval));
}

type_domain_t type_domain_t::operator&(const type_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

type_domain_t type_domain_t::widen(const type_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

type_domain_t type_domain_t::narrow(const type_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

std::string type_domain_t::domain_name() const {
    return "type_domain";
}

crab::bound_t type_domain_t::get_instruction_count_upper_bound() {
    /* WARNING: The operation is not implemented yet.*/
    return crab::bound_t{crab::number_t{0}};
}

string_invariant type_domain_t::to_set() {
    return string_invariant{};
}

void type_domain_t::operator()(const Undefined & u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const Un &u, location_t loc, int print) {
    m_interval(u, loc);
}

void type_domain_t::operator()(const LoadMapFd &u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const Call &u, location_t loc, int print) {

    using interval_values_stack_t = std::map<uint64_t, interval_cells_t>;
    interval_values_stack_t stack_values;
    for (ArgPair param : u.pairs) {
        if (param.kind == ArgPair::Kind::PTR_TO_WRITABLE_MEM) {
            auto maybe_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(param.mem.v);
            auto maybe_width_interval = m_interval.find_interval_value(param.size.v);
            if (!maybe_ptr_or_mapfd || !maybe_width_interval) continue;
            auto ptr_or_mapfd = maybe_ptr_or_mapfd.value();
            auto width_interval = maybe_width_interval.value();
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd)) {
                auto ptr_with_off = std::get<ptr_with_off_t>(ptr_or_mapfd);
                if (ptr_with_off.get_region() == region_t::T_STACK) {
                    auto offset_singleton = ptr_with_off.get_offset().singleton();
                    if (!offset_singleton) {
                        //std::cout << "type error: storing at an unknown offset in stack\n";
                        m_errors.push_back("storing at an unknown offset in stack");
                        continue;
                    }
                    auto offset = (uint64_t)offset_singleton.value();
                    if (auto single_width = width_interval.singleton(); single_width) {
                        int width = (int)single_width.value();
                        stack_values[offset] = std::make_pair(interval_t::top(), width);
                    }
                }
            }
        }
    }
    m_region(u, loc);
    m_offset(u, loc);
    m_interval.do_call(u, stack_values, loc);
}

void type_domain_t::operator()(const Exit &u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const Jmp &u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const Packet & u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const LockAdd &u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const Assume &u, location_t loc, int print) {
    Condition cond = u.cond;
    auto ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(cond.left.v);
    if (ptr_or_mapfd) {
        auto ptr_or_mapfd_type = ptr_or_mapfd.value();
        if (std::holds_alternative<ptr_no_off_t>(ptr_or_mapfd_type))
            m_offset(u, loc);
    }
    m_interval(u, loc);
}

void type_domain_t::operator()(const ValidAccess& s, location_t loc, int print) {
    auto reg_type = m_region.find_ptr_or_mapfd_type(s.reg.v);
    auto interval_type = m_interval.find_interval_value(s.reg.v);
    std::optional<interval_t> width_interval = {};
    if (std::holds_alternative<Reg>(s.width)) {
        width_interval = m_interval.find_interval_value(std::get<Reg>(s.width).v);
    }
    m_offset.check_valid_access(s, reg_type, interval_type, width_interval);
}

void type_domain_t::operator()(const TypeConstraint& s, location_t loc, int print) {
    m_region(s, loc);
}

void type_domain_t::operator()(const Assert &u, location_t loc, int print) {
    std::visit([this, loc, print](const auto& v) { std::apply(*this, std::make_tuple(v, loc, print)); }, u.cst);
}

static bool is_mapfd_type(const ptr_or_mapfd_t& ptr_or_mapfd) {
    return (std::holds_alternative<mapfd_t>(ptr_or_mapfd));
}

static region_t get_region(const ptr_t& ptr) {
    if (std::holds_alternative<ptr_with_off_t>(ptr)) {
        return std::get<ptr_with_off_t>(ptr).get_region();
    }
    else {
        return std::get<ptr_no_off_t>(ptr).get_region();
    }
}

void type_domain_t::operator()(const Comparable& u, location_t loc, int print) {

    auto maybe_ptr_or_mapfd1 = m_region.find_ptr_or_mapfd_type(u.r1.v);
    auto maybe_ptr_or_mapfd2 = m_region.find_ptr_or_mapfd_type(u.r2.v);
    auto maybe_interval1 = m_interval.find_interval_value(u.r1.v);
    auto maybe_interval2 = m_interval.find_interval_value(u.r2.v);
    if (maybe_ptr_or_mapfd1 && maybe_ptr_or_mapfd2) {
        if (!maybe_interval1 && !maybe_interval2) {
            // an extra check just to make sure registers are not labelled both ptrs and numbers
            auto ptr_or_mapfd1 = maybe_ptr_or_mapfd1.value();
            auto ptr_or_mapfd2 = maybe_ptr_or_mapfd1.value();
            if (is_mapfd_type(ptr_or_mapfd1) && is_mapfd_type(ptr_or_mapfd2)) {
                return;
            }
            else if (!is_mapfd_type(ptr_or_mapfd1) && !is_mapfd_type(ptr_or_mapfd2)) {
                auto ptr1 = get_ptr(ptr_or_mapfd1);
                auto ptr2 = get_ptr(ptr_or_mapfd2);
                if (get_region(ptr1) == get_region(ptr2)) {
                    return;
                }
            }
        }
    }
    else if (!maybe_ptr_or_mapfd1 && !maybe_ptr_or_mapfd2) {
        // all other cases when we do not have a ptr or mapfd, the type is a number
        return;
    }
    //std::cout << "type error: Non-comparable types\n";
    m_errors.push_back("Non-comparable types");
}

void type_domain_t::operator()(const Addable& u, location_t loc, int print) {
    m_region(u, loc);
}

void type_domain_t::operator()(const ValidStore& u, location_t loc, int print) {
    m_region(u, loc);
}


void type_domain_t::operator()(const ValidSize& u, location_t loc, int print) {
    m_interval(u, loc);
}

void type_domain_t::operator()(const ValidMapKeyValue& u, location_t loc, int print) {

    int width;
    auto maybe_ptr_or_mapfd_basereg = m_region.find_ptr_or_mapfd_type(u.access_reg.v);
    auto maybe_mapfd = m_region.find_ptr_or_mapfd_type(u.map_fd_reg.v);
    if (maybe_ptr_or_mapfd_basereg && maybe_mapfd) {
        auto ptr_or_mapfd_basereg = maybe_ptr_or_mapfd_basereg.value();
        auto mapfd = maybe_mapfd.value();
        if (is_mapfd_type(mapfd)) {
            auto mapfd_type = std::get<mapfd_t>(mapfd);
            if (u.key) {
                width = (int)mapfd_type.get_key_size();
            }
            else {
                width = (int)mapfd_type.get_value_size();
            }
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd_basereg)) {
                auto ptr_with_off = std::get<ptr_with_off_t>(ptr_or_mapfd_basereg);
                if (ptr_with_off.get_region() == region_t::T_STACK) {
                    auto offset_singleton = ptr_with_off.get_offset().singleton();
                    if (!offset_singleton) {
                        //std::cout << "type error: reading the stack at an unknown offset\n";
                        m_errors.push_back("reading the stack at an unknown offset");
                        return;
                    }
                    auto offset_to_check = (uint64_t)offset_singleton.value();
                    auto it = m_interval.all_numeric_in_stack(offset_to_check, width);
                    if (it) return;
                    auto it2 = m_region.find_in_stack(offset_to_check);
                    if (it2) {
                        //std::cout << "type error: map update with a non-numerical value\n";
                        m_errors.push_back("map update with a non-numerical value");
                    }
                }
            }
            else if (std::holds_alternative<ptr_no_off_t>(ptr_or_mapfd_basereg)) {
                auto ptr_no_off = std::get<ptr_no_off_t>(ptr_or_mapfd_basereg);
                if (ptr_no_off.get_region() == region_t::T_PACKET) {
                    if (m_offset.check_packet_access(u.access_reg, width, 0, true)) return;
                }
            }
        }
    }
    //std::cout << "type error: valid map key value assertion failed\n";
    m_errors.push_back("valid map key value assertion failed");
}

void type_domain_t::operator()(const ZeroCtxOffset& u, location_t loc, int print) {

    auto maybe_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(u.reg.v);
    if (maybe_ptr_or_mapfd) {
        if (std::holds_alternative<ptr_with_off_t>(maybe_ptr_or_mapfd.value())) {
            auto ptr_type_with_off = std::get<ptr_with_off_t>(maybe_ptr_or_mapfd.value());
            if (ptr_type_with_off.get_offset() == interval_t{crab::number_t{0}}) return;
        }
        auto maybe_dist = m_offset.find_offset_info(u.reg.v);
        if (maybe_dist) {
            auto dist_val = maybe_dist.value().m_dist;
            auto single_val = dist_val.singleton();
            if (single_val) {
                auto dist_value = single_val.value();
                if (dist_value == crab::number_t{0}) return;
            }
        }
    }
    //std::cout << "type error: Zero Offset assertion fail\n";
    m_errors.push_back("Zero Offset assertion fail");
}

type_domain_t type_domain_t::setup_entry() {
    region_domain_t reg = region_domain_t::setup_entry();
    offset_domain_t off = offset_domain_t::setup_entry();
    interval_prop_domain_t cp = interval_prop_domain_t::setup_entry();
    type_domain_t typ(std::move(reg), std::move(off), std::move(cp));
    return typ;
}

void type_domain_t::operator()(const Bin& bin, location_t loc, int print) {

    auto dst_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(bin.dst.v);
    auto dst_interval = m_interval.find_interval_value(bin.dst.v);
    // both region domain and interval domain should not have a definition for the dst register
    assert(!dst_ptr_or_mapfd || !dst_interval);

    std::optional<ptr_or_mapfd_t> src_ptr_or_mapfd;
    std::optional<interval_t> src_interval;
    if (std::holds_alternative<Reg>(bin.v)) {
        Reg r = std::get<Reg>(bin.v);
        src_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(r.v);
        src_interval = m_interval.find_interval_value(r.v);
        // both region domain and interval domain should not have a definition for the src register
        assert(!src_ptr_or_mapfd || !src_interval);
    }
    else {
        auto imm = std::get<Imm>(bin.v);
        src_interval = interval_t{crab::number_t{static_cast<int>(imm.v)}};
    }

    using Op = Bin::Op;
    // for all operations except mov, add, sub, the src and dst should be numbers
    if ((src_ptr_or_mapfd || dst_ptr_or_mapfd)
            && (bin.op != Op::MOV && bin.op != Op::ADD && bin.op != Op::SUB)) {
        //std::cout << "type error: operation on pointers not allowed\n";
        m_errors.push_back("operation on pointers not allowed");
        m_region -= bin.dst.v;
        m_offset -= bin.dst.v;
        m_interval -= bin.dst.v;
        return;
    }

    interval_t subtracted_reg =
        m_region.do_bin(bin, src_interval, dst_interval, src_ptr_or_mapfd, dst_ptr_or_mapfd, loc);
    interval_t subtracted_off =
        m_offset.do_bin(bin, src_interval, dst_interval, src_ptr_or_mapfd, dst_ptr_or_mapfd, loc);
    auto subtracted = subtracted_reg.is_bottom() ? subtracted_off : subtracted_reg;
    m_interval.do_bin(bin, src_interval, dst_interval, src_ptr_or_mapfd, dst_ptr_or_mapfd,
            subtracted, loc);
}

void type_domain_t::do_load(const Mem& b, const Reg& target_reg, location_t loc, int print) {

    Reg basereg = b.access.basereg;
    auto basereg_type = m_region.find_ptr_or_mapfd_type(basereg.v);

    m_region.do_load(b, target_reg, loc);

    auto updated_targetreg_type = m_region.find_ptr_or_mapfd_type(target_reg.v);

    m_interval.do_load(b, target_reg, basereg_type, updated_targetreg_type, loc);
    m_offset.do_load(b, target_reg, basereg_type, loc);
}

void type_domain_t::do_mem_store(const Mem& b, const Reg& target_reg, location_t loc, int print) {

    Reg basereg = b.access.basereg;
    auto basereg_type = m_region.find_ptr_or_mapfd_type(basereg.v);
    auto targetreg_type = m_region.find_ptr_or_mapfd_type(target_reg.v);

    m_region.do_mem_store(b, target_reg, loc);
    m_interval.do_mem_store(b, target_reg, basereg_type);
    m_offset.do_mem_store(b, target_reg, basereg_type, targetreg_type);
}

void type_domain_t::operator()(const Mem& b, location_t loc, int print) {
    if (std::holds_alternative<Reg>(b.value)) {
        if (b.is_load) {
            do_load(b, std::get<Reg>(b.value), loc, print);
        } else {
            do_mem_store(b, std::get<Reg>(b.value), loc, print);
        }
    } else {
        std::string s = std::to_string(static_cast<unsigned int>(std::get<Imm>(b.value).v));
        std::string desc = std::string("\tEither loading to a number (not allowed) or storing a number (not allowed yet) - ") + s + "\n";
        //std::cout << desc;
        m_errors.push_back(desc);
        return;
    }
}

// the method does not work well as it requires info about the label of basic block we are in
// this info is not available when we are only printing any state
// but it is available when we are processing a basic block for all its instructions:w
//
void type_domain_t::print_registers() const {
    std::cout << "  register types: {\n";
    for (size_t i = 0; i < NUM_REGISTERS; i++) {
        register_t reg = (register_t)i;
        auto maybe_ptr_or_mapfd_type = m_region.find_ptr_or_mapfd_type(reg);
        auto maybe_offset = m_offset.find_offset_info(reg);
        auto maybe_interval = m_interval.find_interval_value(reg);
        if (maybe_ptr_or_mapfd_type || maybe_interval) {
            std::cout << "    ";
            print_register(Reg{(uint8_t)reg}, maybe_ptr_or_mapfd_type, maybe_offset,
                    maybe_interval);
            std::cout << "\n";
        }
    }
    std::cout << "  }\n";
}

void type_domain_t::print_ctx() const {
    std::vector<uint64_t> ctx_keys = m_region.get_ctx_keys();
    std::cout << "  ctx: {\n";
    for (auto const& k : ctx_keys) {
        auto ptr = m_region.find_in_ctx(k);
        auto dist = m_offset.find_in_ctx(k);
        if (ptr) {
            std::cout << "    " << k << ": ";
            print_ptr_type(ptr.value(), dist);
            std::cout << ",\n";
        }
    }
    std::cout << "  }\n";
}

void type_domain_t::print_stack() const {
    std::vector<uint64_t> stack_keys_region = m_region.get_stack_keys();
    std::vector<uint64_t> stack_keys_interval = m_interval.get_stack_keys();
    std::cout << "  stack: {\n";
    for (auto const& k : stack_keys_region) {
        auto maybe_ptr_or_mapfd_cells = m_region.find_in_stack(k);
        auto dist = m_offset.find_in_stack(k);
        if (maybe_ptr_or_mapfd_cells) {
            auto ptr_or_mapfd_cells = maybe_ptr_or_mapfd_cells.value();
            int width = ptr_or_mapfd_cells.second;
            auto ptr_or_mapfd = ptr_or_mapfd_cells.first;
            std::cout << "    [" << k << "-" << k+width-1 << "] : ";
            print_ptr_or_mapfd_type(ptr_or_mapfd, dist);
            std::cout << ",\n";
        }
    }
    for (auto const& k : stack_keys_interval) {
        auto maybe_interval_cells = m_interval.find_in_stack(k);
        if (maybe_interval_cells) {
            auto interval_cells = maybe_interval_cells.value();
            int width = interval_cells.second;
            std::cout << "    [" << k << "-" << k+width-1 << "] : ";
            print_number(interval_cells.first);
            std::cout << ",\n";
        }
    }
    std::cout << "  }\n";
}

void type_domain_t::adjust_bb_for_types(location_t loc) {
    m_region.adjust_bb_for_types(loc);
    m_offset.adjust_bb_for_types(loc);
    m_interval.adjust_bb_for_types(loc);
}

void type_domain_t::operator()(const basic_block_t& bb, bool check_termination, int print) {

    if (print != 0) {
        write(std::cout, bb, print);
        return;
    }

    auto label = bb.label();
    uint32_t curr_pos = 0;
    location_t loc = location_t(std::make_pair(label, curr_pos));
    if (print == 0)
        adjust_bb_for_types(loc);

    for (const Instruction& statement : bb) {
        loc = location_t(std::make_pair(label, ++curr_pos));
        std::visit([this, loc, print](const auto& v) { std::apply(*this, std::make_tuple(v, loc, print)); }, statement);
    }

    operator+=(m_region.get_errors());
    operator+=(m_offset.get_errors());
    operator+=(m_interval.get_errors());
}

void type_domain_t::write(std::ostream& o, const basic_block_t& bb, int print) const {
    if (is_bottom()) {
        o << bb << "\n";
        return;
    }
    if (print < 0) {
        o << "state of stack and ctx in program:\n";
        print_ctx();
        print_stack();
        o << "\n";
        return;
    }

    o << bb.label() << ":\n";
    uint32_t curr_pos = 0;
    for (const Instruction& statement : bb) {
        ++curr_pos;
        location_t loc = location_t(std::make_pair(bb.label(), curr_pos));
        o << "   " << curr_pos << ".";
        if (std::holds_alternative<Call>(statement)) {
            auto r0_reg = reg_with_loc_t(register_t{R0_RETURN_VALUE}, loc);
            auto region = m_region.find_ptr_or_mapfd_at_loc(r0_reg);
            auto offset = m_offset.find_offset_at_loc(r0_reg);
            auto interval = m_interval.find_interval_at_loc(r0_reg);
            print_annotated(o, std::get<Call>(statement), region, offset, interval);
        }
        else if (std::holds_alternative<Bin>(statement)) {
            auto b = std::get<Bin>(statement);
            auto reg_with_loc = reg_with_loc_t(b.dst.v, loc);
            auto region = m_region.find_ptr_or_mapfd_at_loc(reg_with_loc);
            auto offset = m_offset.find_offset_at_loc(reg_with_loc);
            auto interval = m_interval.find_interval_at_loc(reg_with_loc);
            print_annotated(o, b, region, offset, interval);
        }
        else if (std::holds_alternative<Mem>(statement)) {
            auto u = std::get<Mem>(statement);
            if (u.is_load) {
                auto target_reg = std::get<Reg>(u.value);
                auto target_reg_loc = reg_with_loc_t(target_reg.v, loc);
                auto region = m_region.find_ptr_or_mapfd_at_loc(target_reg_loc);
                auto offset = m_offset.find_offset_at_loc(target_reg_loc);
                auto interval = m_interval.find_interval_at_loc(target_reg_loc);
                print_annotated(o, u, region, offset, interval);
            }
            else o << "  " << u << "\n";
        }
        else if (std::holds_alternative<LoadMapFd>(statement)) {
            auto u = std::get<LoadMapFd>(statement);
            auto reg = reg_with_loc_t(u.dst.v, loc);
            auto region = m_region.find_ptr_or_mapfd_at_loc(reg);
            print_annotated(o, u, region);
        }
        else o << "  " << statement << "\n";
    }

    auto [it, et] = bb.next_blocks();
    if (it != et) {
        o << "  " << "goto ";
        for (; it != et;) {
            o << *it;
            ++it;
            if (it == et) {
                o << ";";
            } else {
                o << ",";
            }
        }
    }
    o << "\n\n";
}

std::ostream& operator<<(std::ostream& o, const type_domain_t& typ) {
    typ.write(o);
    return o;
}
