// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/type_domain.hpp"

namespace crab {

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

void type_domain_t::operator()(const Undefined& u, location_t loc, int print) {
    // nothing to do here
}

void type_domain_t::operator()(const Un& u, location_t loc, int print) {
    m_interval(u, loc);
}

void type_domain_t::operator()(const LoadMapFd& u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const Call& u, location_t loc, int print) {

    using interval_values_stack_t = std::map<uint64_t, interval_cells_t>;
    interval_values_stack_t stack_values;
    for (ArgPair param : u.pairs) {
        if (param.kind == ArgPair::Kind::PTR_TO_WRITABLE_MEM) {
            auto maybe_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(param.mem.v);
            auto maybe_width_interval = m_interval.find_interval_value(param.size.v);
            if (!maybe_ptr_or_mapfd || !maybe_width_interval) continue;
            if (is_stack_ptr(maybe_ptr_or_mapfd)) {
                auto ptr_with_off = std::get<ptr_with_off_t>(*maybe_ptr_or_mapfd);
                auto width_interval = maybe_width_interval->to_interval();

                auto offset_singleton = ptr_with_off.get_offset().to_interval().singleton();
                if (!offset_singleton) {
                    //std::cout << "type error: storing at an unknown offset in stack\n";
                    m_errors.push_back("storing at an unknown offset in stack");
                    continue;
                }
                auto offset = (uint64_t)offset_singleton.value();
                if (auto single_width = width_interval.singleton()) {
                    int width = (int)single_width.value();

                    // Removing overlapping cells
                    // TODO: do it for all domains
                    auto overlapping_cells
                        = m_interval.find_overlapping_cells_in_stack(offset, width);
                    m_interval.remove_overlap_in_stack(overlapping_cells, offset, width);

                    stack_values[offset] = std::make_pair(interval_t::top(), width);
                }
            }
        }
    }
    m_region(u, loc);
    m_offset(u, loc);
    m_interval.do_call(u, stack_values, loc);
}

void type_domain_t::operator()(const Exit& u, location_t loc, int print) {
    // nothing to do here
}

void type_domain_t::operator()(const Jmp& u, location_t loc, int print) {}

void type_domain_t::operator()(const LockAdd& u, location_t loc, int print) {
    // nothing to do here
}

void type_domain_t::operator()(const Packet& u, location_t loc, int print) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void type_domain_t::operator()(const Assume& s, location_t loc, int print) {
    Condition cond = s.cond;
    const auto& maybe_left_type = m_region.find_ptr_or_mapfd_type(cond.left.v);
    if (std::holds_alternative<Reg>(cond.right)) {
        const auto& right_reg = std::get<Reg>(cond.right);
        const auto& maybe_right_type = m_region.find_ptr_or_mapfd_type(right_reg.v);
        if (maybe_left_type && maybe_right_type) {
            // both pointers
            if (is_packet_ptr(maybe_left_type) && is_packet_ptr(maybe_right_type)) {
                m_offset(s, loc, print);
            }
        }
        else if (!maybe_left_type && !maybe_right_type) {
            // both numbers
            auto left_interval = m_interval.find_interval_value(cond.left.v);
            m_interval.assume_cst(cond.op, cond.is64, cond.left.v, cond.right, loc);
        }
        else {
            // We should only reach here if `--assume-assert` is off
            assert(!thread_local_options.assume_assertions || is_bottom());
            // be sound in any case, it happens to flush out bugs:
            m_region.set_registers_to_top();
        }
    }
    else {
        m_interval.assume_cst(cond.op, cond.is64, cond.left.v, cond.right, loc);
    }
}

void type_domain_t::operator()(const ValidDivisor& s, location_t loc, int print) {
    /* WARNING: The operation is not implemented yet.*/
}

void type_domain_t::operator()(const ValidAccess& s, location_t loc, int print) {
    auto reg_type = m_region.find_ptr_or_mapfd_type(s.reg.v);
    auto mock_interval_type = m_interval.find_interval_value(s.reg.v);
    auto interval_type = 
        mock_interval_type ? mock_interval_type->to_interval() : std::optional<interval_t>{};
    if (reg_type) {
        m_region(s, loc);
        std::optional<mock_interval_t> width_mock_interval;
        if (std::holds_alternative<Reg>(s.width)) {
            width_mock_interval = m_interval.find_interval_value(std::get<Reg>(s.width).v);
        }
        else {
            auto imm = std::get<Imm>(s.width); 
            width_mock_interval = mock_interval_t{interval_t(number_t{imm.v}, number_t{imm.v})};
        }
        auto width_interval = 
            width_mock_interval ? width_mock_interval->to_interval() : std::optional<interval_t>{};
        if (is_packet_ptr(reg_type)) {
            m_offset.check_valid_access(s, reg_type, interval_type, width_interval);
        }
    }
    else {
        m_interval(s, loc);
    }
}

void type_domain_t::operator()(const TypeConstraint& s, location_t loc, int print) {
    auto reg_type = m_region.find_ptr_or_mapfd_type(s.reg.v);
    auto mock_interval_type = m_interval.find_interval_value(s.reg.v);
    assert(!reg_type || !mock_interval_type);
    m_region(s, loc);
}

void type_domain_t::operator()(const Assert& u, location_t loc, int print) {
    std::visit([this, loc, print](const auto& v) { std::apply(*this, std::make_tuple(v, loc, print)); }, u.cst);
}

void type_domain_t::operator()(const Comparable& u, location_t loc, int print) {

    auto maybe_ptr_or_mapfd1 = m_region.find_ptr_or_mapfd_type(u.r1.v);
    auto maybe_ptr_or_mapfd2 = m_region.find_ptr_or_mapfd_type(u.r2.v);
    auto maybe_num_type1 = m_interval.find_interval_value(u.r1.v);
    auto maybe_num_type2 = m_interval.find_interval_value(u.r2.v);
    assert(!maybe_ptr_or_mapfd1 || !maybe_num_type1);
    assert(!maybe_ptr_or_mapfd2 || !maybe_num_type2);
    if (maybe_ptr_or_mapfd1 && maybe_ptr_or_mapfd2) {
        if (is_mapfd_type(maybe_ptr_or_mapfd1) && is_mapfd_type(maybe_ptr_or_mapfd2)) return;
        if (!is_shared_ptr(maybe_ptr_or_mapfd1)
                && same_region(*maybe_ptr_or_mapfd1, *maybe_ptr_or_mapfd2)) return;
    }
    else if (!maybe_ptr_or_mapfd2) {
        // TODO: interval check here
        // two numbers can be compared
        // if r1 is a pointer, r2 must be a number
        return;
    }
    //std::cout << "type error: Non-comparable types\n";
    m_errors.push_back("Non-comparable types");
}

void type_domain_t::operator()(const Addable& u, location_t loc, int print) {
    auto maybe_ptr_or_mapfd1 = m_region.find_ptr_or_mapfd_type(u.ptr.v);
    auto maybe_ptr_or_mapfd2 = m_region.find_ptr_or_mapfd_type(u.num.v);
    auto maybe_num_type1 = m_interval.find_interval_value(u.ptr.v);
    auto maybe_num_type2 = m_interval.find_interval_value(u.num.v);
    assert(!maybe_ptr_or_mapfd1 || !maybe_num_type1);
    assert(!maybe_ptr_or_mapfd2 || !maybe_num_type2);
    m_region(u, loc);
    // TODO: move the definition from the region domain to the type domain
}

void type_domain_t::operator()(const ValidStore& u, location_t loc, int print) {
    auto maybe_ptr_or_mapfd1 = m_region.find_ptr_or_mapfd_type(u.mem.v);
    auto maybe_ptr_or_mapfd2 = m_region.find_ptr_or_mapfd_type(u.val.v);
    auto maybe_num_type1 = m_interval.find_interval_value(u.mem.v);
    auto maybe_num_type2 = m_interval.find_interval_value(u.val.v);
    assert(!maybe_ptr_or_mapfd1 || !maybe_num_type1);
    assert(!maybe_ptr_or_mapfd2 || !maybe_num_type2);
    m_region(u, loc);
    // TODO: move the definition from the region domain to the type domain
}

void type_domain_t::operator()(const ValidSize& u, location_t loc, int print) {
    auto maybe_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(u.reg.v);
    auto maybe_num_type = m_interval.find_interval_value(u.reg.v);
    assert(!maybe_ptr_or_mapfd || !maybe_num_type);
    m_interval(u, loc);
}

void type_domain_t::operator()(const ValidMapKeyValue& u, location_t loc, int print) {

    // TODO: move map-related function to common
    //auto fd_type = m_region.get_map_type(u.map_fd_reg);

    int width;
    if (u.key) {
        auto key_size = m_region.get_map_key_size(u.map_fd_reg).singleton();
        if (!key_size.has_value()) {
            m_errors.push_back("Map key size is not singleton");
            return;
        }
        width = (int)key_size.value();
    } else {
        auto value_size = m_region.get_map_value_size(u.map_fd_reg).singleton();
        if (!value_size.has_value()) {
            m_errors.push_back("Map value size is not singleton");
            return;
        }
        width = (int)value_size.value();
    }
    auto maybe_ptr_or_mapfd_basereg = m_region.find_ptr_or_mapfd_type(u.access_reg.v);
    auto maybe_mapfd = m_region.find_ptr_or_mapfd_type(u.map_fd_reg.v);
    if (maybe_ptr_or_mapfd_basereg && maybe_mapfd) {
        auto mapfd = maybe_mapfd.value();
        if (is_mapfd_type(maybe_mapfd)) {
            auto ptr_or_mapfd_basereg = maybe_ptr_or_mapfd_basereg.value();
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd_basereg)) {
                auto ptr_with_off = std::get<ptr_with_off_t>(ptr_or_mapfd_basereg);
                if (ptr_with_off.get_region() == region_t::T_STACK) {
                    auto offset_singleton = ptr_with_off.get_offset().to_interval().singleton();
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
                    return;
                }
            }
            else if (std::holds_alternative<ptr_no_off_t>(ptr_or_mapfd_basereg)) {
                if (m_offset.check_packet_access(u.access_reg, width, 0, true)) return;
            }
            else {
                m_errors.push_back("Only stack or packet can be used as a parameter");
            }
        }
    }
    //std::cout << "type error: valid map key value assertion failed\n";
    m_errors.push_back("valid map key value assertion failed");
}

void type_domain_t::operator()(const ZeroCtxOffset& u, location_t loc, int print) {
    m_region(u, loc);
}

type_domain_t type_domain_t::setup_entry() {
    auto&& reg = crab::region_domain_t::setup_entry();
    auto&& off = offset_domain_t::setup_entry();
    auto&& interval = interval_prop_domain_t::setup_entry();
    type_domain_t typ(std::move(reg), std::move(off), std::move(interval));
    return typ;
}

void type_domain_t::operator()(const Bin& bin, location_t loc, int print) {
    if (is_bottom()) return;
    std::optional<ptr_or_mapfd_t> src_ptr_or_mapfd;
    std::optional<interval_t> src_interval;

    if (std::holds_alternative<Reg>(bin.v)) {
        Reg r = std::get<Reg>(bin.v);
        src_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(r.v);
        auto src_mock_interval = m_interval.find_interval_value(r.v);
        if (src_mock_interval) {
            src_interval = src_mock_interval->to_interval();
        }
    }
    else {
        int64_t imm;
        if (bin.is64) {
            imm = static_cast<int64_t>(std::get<Imm>(bin.v).v);
        }
        else {
            imm = static_cast<int>(std::get<Imm>(bin.v).v);
        }
        src_interval = interval_t{number_t{imm}};
    }
    assert(!src_ptr_or_mapfd || !src_interval);

    auto dst_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(bin.dst.v);
    auto dst_mock_interval = m_interval.find_interval_value(bin.dst.v);
    auto dst_interval = dst_mock_interval ? dst_mock_interval->to_interval() :
        std::optional<interval_t>();
    assert(!dst_ptr_or_mapfd.has_value() || !dst_interval.has_value());

    using Op = Bin::Op;
    // for all operations except mov, add, sub, the src and dst should be numbers
    if ((src_ptr_or_mapfd || dst_ptr_or_mapfd)
            && (bin.op != Op::MOV && bin.op != Op::ADD && bin.op != Op::SUB)) {
        m_errors.push_back("operation on pointers not allowed");
        //std::cout << "type error: operation on pointers not allowed\n";
        m_region -= bin.dst.v;
        m_offset -= bin.dst.v;
        m_interval -= bin.dst.v;
        return;
    }

    interval_t subtracted_reg =
        m_region.do_bin(bin, src_interval, src_ptr_or_mapfd, dst_ptr_or_mapfd, loc);
    interval_t subtracted_off =
        m_offset.do_bin(bin, src_interval, dst_interval, src_ptr_or_mapfd, dst_ptr_or_mapfd, loc);
    auto subtracted = subtracted_reg.is_bottom() ? subtracted_off : subtracted_reg;
    m_interval.do_bin(bin, src_interval, dst_interval, src_ptr_or_mapfd, dst_ptr_or_mapfd,
            subtracted, loc);
}

void type_domain_t::do_load(const Mem& b, const Reg& target_reg, bool unknown_ptr,
        std::optional<ptr_or_mapfd_t> basereg_opt, location_t loc) {
    m_region.do_load(b, target_reg, unknown_ptr, loc);
    m_offset.do_load(b, target_reg, basereg_opt, loc);
    // TODO: replace with a bool value returned from region do_load
    auto load_in_region = m_region.find_ptr_or_mapfd_type(target_reg.v).has_value();
    m_interval.do_load(b, target_reg, basereg_opt, load_in_region, loc);
}

void type_domain_t::do_mem_store(const Mem& b, std::optional<ptr_or_mapfd_t> target_opt,
        std::optional<ptr_or_mapfd_t>& basereg_opt, location_t loc) {
    m_region.do_mem_store(b, loc);
    m_interval.do_mem_store(b, basereg_opt);
    // TODO: replace target_opt with a bool value representing whether we have a packet pointer,
    // because that is the case target_opt is needed for
    m_offset.do_mem_store(b, target_opt, basereg_opt);
}

void type_domain_t::operator()(const Mem& b, location_t loc, int print) {
    auto basereg = b.access.basereg;
    auto base_ptr_or_mapfd_opt = m_region.find_ptr_or_mapfd_type(basereg.v);
    bool unknown_ptr = !base_ptr_or_mapfd_opt.has_value();
    if (unknown_ptr) {
        std::string s = std::to_string(static_cast<unsigned int>(basereg.v));
        m_errors.push_back(
                std::string("load/store using an unknown pointer, or number - r") + s);
    }
    if (std::holds_alternative<Reg>(b.value)) {
        auto targetreg = std::get<Reg>(b.value);
        auto targetreg_type = m_region.find_ptr_or_mapfd_type(targetreg.v);
        if (b.is_load) do_load(b, targetreg, unknown_ptr, base_ptr_or_mapfd_opt, loc);
        else if (!unknown_ptr) do_mem_store(b, targetreg_type, base_ptr_or_mapfd_opt, loc);
    }
    else if (!unknown_ptr && !b.is_load) {
        do_mem_store(b, std::nullopt, base_ptr_or_mapfd_opt, loc);
    }
}

void type_domain_t::print_ctx() const {
    std::vector<uint64_t> ctx_keys = m_region.get_ctx_keys();
    std::cout << "\tctx: {\n";
    for (auto const& k : ctx_keys) {
        auto ptr = m_region.find_in_ctx(k);
        auto dist = m_offset.find_in_ctx(k);
        if (ptr) {
            std::cout << "\t\t" << k << ": ";
            print_ptr_type(std::cout, *ptr, dist);
            std::cout << ",\n";
        }
    }
    std::cout << "\t}\n";
}

void type_domain_t::print_stack() const {
    std::vector<uint64_t> stack_keys_region = m_region.get_stack_keys();
    std::vector<uint64_t> stack_keys_interval = m_interval.get_stack_keys();
    std::cout << "\tstack: {\n";
    for (auto const& k : stack_keys_region) {
        auto maybe_ptr_or_mapfd_cells = m_region.find_in_stack(k);
        auto dist = m_offset.find_in_stack(k);
        if (maybe_ptr_or_mapfd_cells) {
            auto ptr_or_mapfd_cells = maybe_ptr_or_mapfd_cells.value();
            int width = ptr_or_mapfd_cells.second;
            auto ptr_or_mapfd = ptr_or_mapfd_cells.first;
            std::cout << "\t\t[" << k << "-" << k+width-1 << "] : ";
            if (dist) 
                print_ptr_or_mapfd_type(std::cout, ptr_or_mapfd, dist->first);
            else 
                print_ptr_or_mapfd_type(std::cout, ptr_or_mapfd, std::nullopt);
            std::cout << ",\n";
        }
    }
    for (auto const& k : stack_keys_interval) {
        auto maybe_interval_cells = m_interval.find_in_stack(k);
        if (maybe_interval_cells) {
            auto interval_cells = maybe_interval_cells.value();
            std::cout << "\t\t" << "[" << k << "-" << k+interval_cells.second-1 << "] : ";
            print_number(std::cout, interval_cells.first.to_interval());
            std::cout << ",\n";
        }
    }
    std::cout << "\t}\n";
}

void type_domain_t::adjust_bb_for_types(location_t loc) {
    m_region.adjust_bb_for_types(loc);
    m_offset.adjust_bb_for_types(loc);
    m_interval.adjust_bb_for_types(loc);
}

void type_domain_t::operator()(const basic_block_t& bb, bool check_termination, int print) {

    if (print != 0) {
        print_annotated(std::cout, *this, bb, print);
        return;
    }

    // A temporary fix to avoid printing errors for multiple basic blocks
    m_errors.clear();
    m_region.reset_errors();
    m_offset.reset_errors();
    m_interval.reset_errors();

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

std::optional<crab::ptr_or_mapfd_t>
type_domain_t::find_ptr_or_mapfd_at_loc(const crab::reg_with_loc_t& loc) const {
    return m_region.find_ptr_or_mapfd_at_loc(loc);
}

std::optional<crab::dist_t>
type_domain_t::find_offset_at_loc(const crab::reg_with_loc_t& loc) const {
    return m_offset.find_offset_at_loc(loc);
}

std::optional<crab::mock_interval_t>
type_domain_t::find_interval_at_loc(const crab::reg_with_loc_t& loc) const {
    return m_interval.find_interval_at_loc(loc);
}

std::ostream& operator<<(std::ostream& o, const type_domain_t& typ) {
    typ.write(o);
    return o;
}

} // namespace crab

void print_annotated(std::ostream& o, const crab::type_domain_t& typ,
        const basic_block_t& bb, int print) {
    if (typ.is_bottom()) {
        o << bb << "\n";
        return;
    }
    if (print < 0) {
        o << "state of stack and ctx in program:\n";
        typ.print_ctx();
        typ.print_stack();
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
            auto r0_reg = crab::reg_with_loc_t(register_t{R0_RETURN_VALUE}, loc);
            auto region = typ.find_ptr_or_mapfd_at_loc(r0_reg);
            auto interval = typ.find_interval_at_loc(r0_reg);
            print_annotated(o, std::get<Call>(statement), region, interval);
        }
        else if (std::holds_alternative<Bin>(statement)) {
            auto b = std::get<Bin>(statement);
            auto reg_with_loc = crab::reg_with_loc_t(b.dst.v, loc);
            auto region = typ.find_ptr_or_mapfd_at_loc(reg_with_loc);
            auto offset = typ.find_offset_at_loc(reg_with_loc);
            auto interval = typ.find_interval_at_loc(reg_with_loc);
            print_annotated(o, b, region, offset, interval);
        }
        else if (std::holds_alternative<Mem>(statement)) {
            auto u = std::get<Mem>(statement);
            if (u.is_load) {
                auto target_reg = std::get<Reg>(u.value);
                auto target_reg_loc = crab::reg_with_loc_t(target_reg.v, loc);
                auto region = typ.find_ptr_or_mapfd_at_loc(target_reg_loc);
                auto offset = typ.find_offset_at_loc(target_reg_loc);
                auto interval = typ.find_interval_at_loc(target_reg_loc);
                print_annotated(o, u, region, offset, interval);
            }
            else o << "  " << u << "\n";
        }
        else if (std::holds_alternative<LoadMapFd>(statement)) {
            auto u = std::get<LoadMapFd>(statement);
            auto reg = crab::reg_with_loc_t(u.dst.v, loc);
            auto region = typ.find_ptr_or_mapfd_at_loc(reg);
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

