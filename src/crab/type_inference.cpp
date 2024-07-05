// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/type_inference.hpp"
#include <regex>

namespace crab {

bool inference_domain_t::is_bottom() const {
    return (m_region.is_bottom() || m_offset.is_bottom() || m_interval.is_bottom());
}

bool inference_domain_t::is_top() const {
    return (m_region.is_top() && m_offset.is_top() && m_interval.is_top());
}

inference_domain_t inference_domain_t::bottom() {
    inference_domain_t typ;
    typ.set_to_bottom();
    return typ;
}

void inference_domain_t::set_to_bottom() {
    m_region.set_to_bottom();
    m_offset.set_to_bottom();
    m_interval.set_to_bottom();
}

void inference_domain_t::set_to_top() {
    m_region.set_to_top();
    m_offset.set_to_top();
    m_interval.set_to_top();
}

bool inference_domain_t::operator<=(const inference_domain_t& abs) const {
    if (abs.is_top() || is_bottom()) return true;
    if (is_top() || abs.is_bottom()) return false;
    return (m_region <= abs.m_region && m_offset <= abs.m_offset && m_interval <= abs.m_interval);
}

inference_domain_t inference_domain_t::widen(const inference_domain_t& other, bool to_constants) {
    return inference_domain_t(m_region.widen(other.m_region, to_constants),
            m_offset.widen(other.m_offset, to_constants),
            m_interval.widen(other.m_interval, to_constants));
}

void inference_domain_t::operator|=(const inference_domain_t& abs) {
    inference_domain_t tmp{abs};
    operator|=(std::move(tmp));
}

void inference_domain_t::operator|=(inference_domain_t&& other) {
    if (is_bottom()) {
        *this = std::move(other);
        return;
    }
    if (other.is_bottom()) return;
    *this = *this | std::move(other);
}

inference_domain_t inference_domain_t::operator|(const inference_domain_t& other) const {
    return inference_domain_t(m_region | other.m_region, m_offset | other.m_offset,
            m_interval | other.m_interval, m_slacks);
}

inference_domain_t inference_domain_t::operator|(inference_domain_t&& other) const {
    return inference_domain_t(m_region | std::move(other.m_region),
            m_offset | std::move(other.m_offset),
            m_interval | std::move(other.m_interval),
            std::move(other.m_slacks));
}

inference_domain_t inference_domain_t::operator&(const inference_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

inference_domain_t inference_domain_t::narrow(const inference_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

void inference_domain_t::initialize_loop_counter(label_t label) {
    // WARNING: Not implemented yet
}

crab::bound_t inference_domain_t::get_loop_count_upper_bound() const {
    // WARNING: Not implemented yet
    return crab::bound_t{crab::number_t{0}};
}

string_invariant inference_domain_t::to_set() const {
    if (is_top()) return string_invariant::top();
    std::set<std::string> result;
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        auto maybe_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(register_t{i});
        auto maybe_rf = m_offset.find_refinement_info(register_t{i});
        if (maybe_ptr_or_mapfd.has_value()) {
            std::stringstream elem;
            print_register(elem, Reg{i}, maybe_ptr_or_mapfd, maybe_rf, {}, false);
            result.insert(elem.str());
        }
        auto maybe_signed_interval = m_interval.find_signed_interval_value(register_t{i});
        if (maybe_signed_interval.has_value()) {
            std::stringstream elem;
            print_register(elem, Reg{i}, maybe_ptr_or_mapfd, maybe_rf, maybe_signed_interval, true);
            result.insert(elem.str());
        }
        auto maybe_unsigned_interval = m_interval.find_unsigned_interval_value(register_t{i});
        if (maybe_unsigned_interval.has_value()) {
            std::stringstream elem;
            print_register(elem, Reg{i}, maybe_ptr_or_mapfd, maybe_rf,
                    maybe_unsigned_interval, false);
            result.insert(elem.str());
        }
    }
   const std::vector<uint64_t>& stack_keys_region = m_region.get_stack_keys();
    for (auto const& k : stack_keys_region) {
        std::stringstream elem;
        auto maybe_ptr_or_mapfd_cells = m_region.find_in_stack(k);
        auto rf = m_offset.find_in_stack(k);
        if (maybe_ptr_or_mapfd_cells.has_value()) {
            auto ptr_or_mapfd_cells = maybe_ptr_or_mapfd_cells.value();
            int width = ptr_or_mapfd_cells.second;
            auto ptr_or_mapfd = ptr_or_mapfd_cells.first;
            elem << "stack";
            if (rf) {
                print_non_numeric_memory_cell(elem, k, k+width-1, ptr_or_mapfd, rf->first);
            }
            else {
                print_non_numeric_memory_cell(elem, k, k+width-1, ptr_or_mapfd);
            }
        }
        result.insert(elem.str());
    }

    const std::vector<uint64_t>& stack_keys_interval = m_interval.get_stack_keys();
    for (auto const& k : stack_keys_interval) {
        auto maybe_interval_cells_signed = m_interval.find_in_stack_signed(k);
        if (maybe_interval_cells_signed.has_value()) {
            std::stringstream elem;
            auto signed_interval_cells = maybe_interval_cells_signed.value();
            elem << "stack";
            print_numeric_memory_cell(elem, k, k+signed_interval_cells.second,
                    signed_interval_cells.first, true);
            result.insert(elem.str());
        }
        auto maybe_interval_cells_unsigned = m_interval.find_in_stack_unsigned(k);
        if (maybe_interval_cells_unsigned.has_value()) {
            std::stringstream elem;
            auto unsigned_interval_cells = maybe_interval_cells_unsigned.value();
            elem << "stack";
            print_numeric_memory_cell(elem, k, k+unsigned_interval_cells.second,
                    unsigned_interval_cells.first, false);
            result.insert(elem.str());
        }
    }
    return string_invariant{result};
}

void inference_domain_t::operator()(const Undefined& u, location_t loc) {
    // nothing to do here
}

void inference_domain_t::operator()(const Un& u, location_t loc) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

void inference_domain_t::operator()(const LoadMapFd& u, location_t loc) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

// Construct a Bin operation that does the main operation that a given Atomic operation does atomically.
static Bin atomic_to_bin(const Atomic& a) {
    Bin bin{
        .dst = Reg{R11_ATOMIC_SCRATCH}, .v = a.valreg, .is64 = (a.access.width == sizeof(uint64_t)), .lddw = false};
    switch (a.op) {
    case Atomic::Op::ADD: bin.op = Bin::Op::ADD; break;
    case Atomic::Op::OR: bin.op = Bin::Op::OR; break;
    case Atomic::Op::AND: bin.op = Bin::Op::AND; break;
    case Atomic::Op::XOR: bin.op = Bin::Op::XOR; break;
    case Atomic::Op::XCHG:
    case Atomic::Op::CMPXCHG: bin.op = Bin::Op::MOV; break;
    default: throw std::exception();
    }
    return bin;
}

void inference_domain_t::operator()(const Atomic &u, location_t loc) {
    // WARNING: Not implemented yet
    if (is_bottom()) return;
    std::optional<ptr_or_mapfd_t> base_reg_opt
        = m_region.find_ptr_or_mapfd_type(u.access.basereg.v);
    std::optional<refinement_t> value_reg_opt = m_interval.find_interval_value(u.valreg.v);
    if (!base_reg_opt || !value_reg_opt) return;
    if (is_stack_ptr(base_reg_opt)) {
        if (u.op == Atomic::Op::CMPXCHG) {
            m_region -= register_t{R0_RETURN_VALUE};
            m_offset -= register_t{R0_RETURN_VALUE};
            insert_in_registers_in_interval_domain(register_t{R0_RETURN_VALUE},
                    loc, refinement_t::numeric_refinement_top());
        }
        else if (u.fetch) {
            insert_in_registers_in_interval_domain(u.valreg.v, loc,
                                    refinement_t::numeric_refinement_top());
        }
        return;
    }
    // Fetch the current value into the R11 pseudo-register.
    const Reg r11{11};
    (*this)(Mem{.access = u.access, .value = r11, .is_load = true}, loc);

    // Compute the new value in R11.
    (*this)(atomic_to_bin(u), loc);

    if (u.op == Atomic::Op::CMPXCHG) {
        // For CMPXCHG, store the original value in r0.
        (*this)(Mem{.access = u.access, .value = Reg{R0_RETURN_VALUE}, .is_load = true}, loc);

        // For the destination, there are 3 possibilities:
        // 1) dst.value == r0.value : set R11 to valreg
        // 2) dst.value != r0.value : don't modify R11
        // 3) dst.value may or may not == r0.value : set R11 to the union of R11 and valreg
        // For now we just havoc the value of R11.
        m_region -= register_t{11};
        m_offset -= register_t{11};
        insert_in_registers_in_interval_domain(register_t{11}, loc,
                refinement_t::numeric_refinement_top());
    } else if (u.fetch) {
        // For other FETCH operations, store the original value in the src register.
        (*this)(Mem{.access = u.access, .value = u.valreg, .is_load = true}, loc);
    }

    // Store the new value back in the original shared memory location.
    // Note that do_mem_store() currently doesn't track shared memory values,
    // but stack memory values are tracked and are legal here.
    (*this)(Mem{.access = u.access, .value = r11, .is_load = false}, loc);

    // Clear the R11 pseudo-register.
    m_region -= register_t{11};
    m_offset -= register_t{11};
    m_interval -= register_t{11};
}

void inference_domain_t::operator()(const IncrementLoopCounter &u, location_t loc) {
    // WARNING: Not implemented yet
}

void inference_domain_t::operator()(const Call& u, location_t loc) {

    stack_cells_t stack_values;
    for (ArgPair param : u.pairs) {
        if (param.kind == ArgPair::Kind::PTR_TO_WRITABLE_MEM) {
            auto maybe_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(param.mem.v);
            auto maybe_width_rf = m_interval.find_signed_interval_value(param.size.v);
            if (!maybe_ptr_or_mapfd || !maybe_width_rf) continue;
            if (is_stack_ptr(maybe_ptr_or_mapfd)) {
                auto ptr_with_off = std::get<ptr_with_off_t>(*maybe_ptr_or_mapfd);
                auto width_interval = maybe_width_rf->get_interval_value();

                auto offset_singleton = ptr_with_off.get_offset().to_interval().singleton();
                if (!offset_singleton) {
                    //std::cout << "type error: storing at an unknown offset in stack\n";
                    m_errors.push_back("storing at an unknown offset in stack");
                    continue;
                }
                auto offset = offset_singleton.value().cast_to<int>();
                if (auto single_width = width_interval.singleton()) {
                    int width = single_width.value().cast_to<int>();
                    stack_values.push_back(std::make_pair(offset, width));
                }
            }
        }
    }
    m_region.do_call(u, stack_values, loc);
    m_offset.do_call(u, stack_values, loc);
    m_interval.do_call(u, stack_values, loc);
}

void inference_domain_t::operator()(const Callx &u, location_t loc) {
    // WARNING: Not implemented yet
    // This operation is not needed for current benchmarks,
    // TODO: implement this if needed
}

void inference_domain_t::operator()(const Exit& u, location_t loc) {
    // nothing to do here
}

void inference_domain_t::operator()(const Jmp& u, location_t loc) {
    // nothing to do here
}

void inference_domain_t::operator()(const Packet& u, location_t loc) {
    m_region(u, loc);
    m_offset(u, loc);
    m_interval(u, loc);
}

static inline bool same_type(const std::optional<ptr_or_mapfd_t>& ptr_or_mapfd1,
        const std::optional<ptr_or_mapfd_t>& ptr_or_mapfd2,
        const std::optional<refinement_t>& interval1,
        const std::optional<refinement_t>& interval2) {
    if (is_mapfd_type(ptr_or_mapfd1) && is_mapfd_type(ptr_or_mapfd2)) return false;
    if (ptr_or_mapfd1 && ptr_or_mapfd2 && same_region(*ptr_or_mapfd1, *ptr_or_mapfd2))
        return true;
    if (interval1 && interval2) return true;
    return false;
}

void inference_domain_t::operator()(const Assume& s, location_t loc) {
    Condition cond = s.cond;
    const auto& maybe_left_type = m_region.find_ptr_or_mapfd_type(cond.left.v);
    const auto& maybe_left_rf = m_interval.find_interval_value(cond.left.v);
    assert(!maybe_left_type.has_value() || !maybe_left_rf.has_value());
    if (std::holds_alternative<Reg>(cond.right)) {
        const auto& right_reg = std::get<Reg>(cond.right);
        const auto& maybe_right_type = m_region.find_ptr_or_mapfd_type(right_reg.v);
        const auto& maybe_right_rf = m_interval.find_interval_value(right_reg.v);
        assert(!maybe_right_type.has_value() || !maybe_right_rf.has_value());
        // TODO: it does not handle for mapfd yet
        if (same_type(maybe_left_type, maybe_right_type, maybe_left_rf, maybe_right_rf)) {
            if (maybe_left_rf) {
                // both numbers
                m_interval.assume_cst(cond.op, cond.is64, register_t{cond.left.v},
                        cond.right, loc);
            }
            else if (maybe_left_type) {
                if (is_packet_ptr(maybe_left_type)) {
                    // both packet pointers
                    m_offset(s, loc);
                }
                else {
                    // other cases, not implemented yet
                }
            }
        }
        else {
            // We should only reach here if `--assume-assert` is off
            assert(!thread_local_options.assume_assertions || is_bottom());
            // be sound in any case, it happens to flush out bugs:
            m_region.set_registers_to_top();
        }
    }
    else {
        if (is_shared_ptr(maybe_left_type)) {
            // left is a shared pointer
            int64_t imm = static_cast<int64_t>(std::get<Imm>(cond.right).v);
            auto shared_ptr = std::get<ptr_with_off_t>(*maybe_left_type);
            m_region.assume_cst(cond.op, std::move(shared_ptr), imm, cond.left.v, loc);
        }
        if (is_mapfd_type(maybe_left_type)) {
            // left is  a mapfd
            // TODO: need to work with values
        }
        else if (maybe_left_rf) {
            m_interval.assume_cst(cond.op, cond.is64, register_t{cond.left.v}, cond.right, loc);
        }
    }
}

void inference_domain_t::operator()(const FuncConstraint& s, location_t loc) {
    // WARNING: Not implemented yet
    // This operation is not needed for current benchmarks,
    // TODO: implement this if needed
}

void inference_domain_t::operator()(const ValidDivisor& u, location_t loc) {
    auto maybe_ptr_or_mapfd_reg = m_region.find_ptr_or_mapfd_type(u.reg.v);
    auto maybe_num_type_reg = m_interval.find_unsigned_interval_value(u.reg.v);
    assert(!maybe_ptr_or_mapfd_reg.has_value() || !maybe_num_type_reg.has_value());

    if (is_ptr_type(maybe_ptr_or_mapfd_reg)) {
        m_errors.push_back("Only numbers can be used as divisors");
    }
    else if (maybe_num_type_reg.has_value() && !thread_local_options.allow_division_by_zero) {
        auto num_type_reg = maybe_num_type_reg->get_interval_value();
        if (interval_t{number_t{0}} <= num_type_reg) {
            m_errors.push_back("Possible division by zero");
        }
    }
}

void inference_domain_t::operator()(const ValidAccess& s, location_t loc) {
    auto reg_type = m_region.find_ptr_or_mapfd_type(s.reg.v);
    if (reg_type) {
        interval_t width_interval = interval_t::bottom();
        if (std::holds_alternative<Reg>(s.width)) {
            auto width_rf = m_interval.find_interval_value(std::get<Reg>(s.width).v);
            if (!width_rf) {
                m_errors.push_back("width is unknown for valid access");
                return;
            }
            width_interval = width_rf->get_interval_value();
        }
        else {
            auto imm = std::get<Imm>(s.width); 
            width_interval = interval_t{number_t{imm.v}};
        }
        if (auto width_number = width_interval.ub().number()) {
            int width = width_number->cast_to<int>();
            m_region.check_valid_access(s, width);
            if (is_packet_ptr(reg_type)) {
                m_offset.check_valid_access(s, reg_type, width);
            }
            if (s.access_type == AccessType::read && is_stack_ptr(reg_type)) {
                auto stack_ptr = std::get<ptr_with_off_t>(*reg_type);
                auto offset_ptr = stack_ptr.get_offset().to_interval();
                m_interval.check_valid_access(s, std::move(offset_ptr), width, true);
            }
        }
        else {
            m_errors.push_back("width is unknown for valid access");
        }
    }
    else {
        auto rf_type = m_interval.find_interval_value(s.reg.v);
        if (rf_type) {
            m_interval.check_valid_access(s, rf_type->get_interval_value());
        }
        else {
            m_errors.push_back("valid access on unknown register");
        }
    }
}

void inference_domain_t::operator()(const TypeConstraint& s, location_t loc) {
    auto reg_type = m_region.find_ptr_or_mapfd_type(s.reg.v);
    auto rf_type = m_interval.find_interval_value(s.reg.v);
    assert(!reg_type.has_value() || !rf_type.has_value());
    m_region.check_type(s, rf_type.has_value());
}

void inference_domain_t::operator()(const Assert& u, location_t loc) {
    std::visit([this, loc](const auto& v) { std::apply(*this, std::make_tuple(v, loc)); }, u.cst);
}

void inference_domain_t::operator()(const Comparable& u, location_t loc) {
    auto maybe_ptr_or_mapfd1 = m_region.find_ptr_or_mapfd_type(u.r1.v);
    auto maybe_ptr_or_mapfd2 = m_region.find_ptr_or_mapfd_type(u.r2.v);
    auto maybe_num_type1 = m_interval.find_interval_value(u.r1.v);
    auto maybe_num_type2 = m_interval.find_interval_value(u.r2.v);
    assert(!maybe_ptr_or_mapfd1.has_value() || !maybe_num_type1.has_value());
    assert(!maybe_ptr_or_mapfd2.has_value() || !maybe_num_type2.has_value());
    if (maybe_ptr_or_mapfd1 && maybe_ptr_or_mapfd2) {
        if (is_mapfd_type(maybe_ptr_or_mapfd1) && is_mapfd_type(maybe_ptr_or_mapfd2)) return;
        if (same_region(*maybe_ptr_or_mapfd1, *maybe_ptr_or_mapfd2)) return;
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

void inference_domain_t::operator()(const Addable& u, location_t loc) {
    auto maybe_ptr_or_mapfd_ptr = m_region.find_ptr_or_mapfd_type(u.ptr.v);
    auto maybe_ptr_or_mapfd_num = m_region.find_ptr_or_mapfd_type(u.num.v);
    auto maybe_num_type_ptr = m_interval.find_interval_value(u.ptr.v);
    auto maybe_num_type_num = m_interval.find_interval_value(u.num.v);
    assert(!maybe_ptr_or_mapfd_ptr.has_value() || !maybe_num_type_ptr.has_value());
    assert(!maybe_ptr_or_mapfd_num.has_value() || !maybe_num_type_num.has_value());

    // a -> b <-> !a || b
    // is_ptr(ptr) -> is_num(num) <-> !is_ptr(ptr) || is_num(num)
    if (!is_ptr_type(maybe_ptr_or_mapfd_ptr) ||
      (!maybe_ptr_or_mapfd_num.has_value() || maybe_num_type_num.has_value())) {
        return;
    }
    m_errors.push_back("Addable assertion fail");
}

void inference_domain_t::operator()(const ValidStore& u, location_t loc) {
    auto maybe_ptr_or_mapfd_mem = m_region.find_ptr_or_mapfd_type(u.mem.v);
    auto maybe_ptr_or_mapfd_val = m_region.find_ptr_or_mapfd_type(u.val.v);
    auto maybe_num_type_mem = m_interval.find_interval_value(u.mem.v);
    auto maybe_num_type_val = m_interval.find_interval_value(u.val.v);
    assert(!maybe_ptr_or_mapfd_mem.has_value() || !maybe_num_type_mem.has_value());
    assert(!maybe_ptr_or_mapfd_val.has_value() || !maybe_num_type_val.has_value());

    // a -> b <-> !a || b
    // !is_stack_ptr(mem) -> is_num(val) <-> is_stack_ptr(mem) || is_num(val)
    if (is_stack_ptr(maybe_ptr_or_mapfd_mem) ||
            (!maybe_ptr_or_mapfd_val.has_value() || maybe_num_type_val.has_value())) {
        return;
    }
    m_errors.push_back("Valid store assertion fail");
}

void inference_domain_t::operator()(const ValidSize& u, location_t loc) {
    auto maybe_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(u.reg.v);
    auto maybe_num_type = m_interval.find_interval_value(u.reg.v);
    assert(!maybe_ptr_or_mapfd || !maybe_num_type);

    if (maybe_num_type) {
        auto reg_value = maybe_num_type->get_interval_value();
        if ((u.can_be_zero && reg_value.lb() >= bound_t{number_t{0}})
                || (!u.can_be_zero && reg_value.lb() > bound_t{number_t{0}})) {
            return;
        }
    }
    m_errors.push_back("Valid Size assertion fail");
}

void inference_domain_t::operator()(const ValidMapKeyValue& u, location_t loc) {

    // TODO: move map-related function to common
    //auto fd_type = m_region.get_map_type(u.map_fd_reg);

    int width;
    if (u.key) {
        auto key_size = m_region.get_map_key_size(u.map_fd_reg).singleton();
        if (!key_size.has_value()) {
            m_errors.push_back("Map key size is not singleton");
            return;
        }
        width = key_size.value().cast_to<int>();
    } else {
        auto value_size = m_region.get_map_value_size(u.map_fd_reg).singleton();
        if (!value_size.has_value()) {
            m_errors.push_back("Map value size is not singleton");
            return;
        }
        width = value_size.value().cast_to<int>();
    }
    auto maybe_ptr_or_mapfd_basereg = m_region.find_ptr_or_mapfd_type(u.access_reg.v);
    auto maybe_mapfd = m_region.find_ptr_or_mapfd_type(u.map_fd_reg.v);
    if (maybe_ptr_or_mapfd_basereg && maybe_mapfd) {
        auto mapfd = maybe_mapfd.value();
        if (is_mapfd_type(maybe_mapfd)) {
            if (is_stack_ptr(maybe_ptr_or_mapfd_basereg)) {
                auto ptr_with_off = std::get<ptr_with_off_t>(*maybe_ptr_or_mapfd_basereg);
                auto offset_singleton = ptr_with_off.get_offset().to_interval().singleton();
                if (!offset_singleton) {
                    //std::cout << "type error: reading the stack at an unknown offset\n";
                    m_errors.push_back("reading the stack at an unknown offset");
                    return;
                }
                auto offset_to_check = offset_singleton.value().cast_to<uint64_t>();
                auto it = m_interval.all_numeric_in_stack(offset_to_check, width);
                if (it) return;
            }
            else if (is_packet_ptr(maybe_ptr_or_mapfd_basereg)) {
                if (m_offset.check_packet_access(u.access_reg, width, 0, true)) return;
            }
            else {
                m_errors.push_back("Only stack or packet can be used as a parameter");
                return;
            }
        }
    }
    m_errors.push_back("map update with a non-numerical value");
}

void inference_domain_t::operator()(const ZeroCtxOffset& u, location_t loc) {
    m_region(u, loc);
}

inference_domain_t inference_domain_t::setup_entry(bool init_r1) {
    std::shared_ptr<slacks_t> slacks = std::make_shared<slacks_t>();
    return inference_domain_t{
        region_domain_t::setup_entry(init_r1),
        offset_domain_t::setup_entry(slacks),
        interval_domain_t::setup_entry(slacks),
        slacks
    };
}

void inference_domain_t::operator()(const Bin& bin, location_t loc) {
    std::optional<ptr_or_mapfd_t> src_ptr_or_mapfd;
    std::optional<refinement_t> src_signed_rf;
    std::optional<interval_t> dst_signed_interval, src_signed_interval;

    auto dst_register = register_t{bin.dst.v};
    if (std::holds_alternative<Reg>(bin.v)) {
        Reg r = std::get<Reg>(bin.v);
        src_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(r.v);
        src_signed_rf = m_interval.find_signed_interval_value(r.v);
        if (src_signed_rf) {
            src_signed_interval = src_signed_rf->get_interval_value();
        }
    }
    auto dst_ptr_or_mapfd = m_region.find_ptr_or_mapfd_type(dst_register);
    auto dst_signed_rf = m_interval.find_signed_interval_value(dst_register);
    if (dst_signed_rf) {
        dst_signed_interval = dst_signed_rf->get_interval_value();
    }

    std::optional<interval_t> subtracted;
    using Op = Bin::Op;
    // ptr -= ptr
    if (std::holds_alternative<Reg>(bin.v) && bin.op == Op::SUB) {
        if (dst_ptr_or_mapfd && src_ptr_or_mapfd) {
            auto dst_ptr = *dst_ptr_or_mapfd;
            auto src_ptr = *src_ptr_or_mapfd;
            if (std::holds_alternative<mapfd_t>(dst_ptr)
                    && std::holds_alternative<mapfd_t>(src_ptr)) {
                m_errors.push_back("mapfd registers subtraction not defined");
            }
            else if (same_region(dst_ptr, src_ptr)) {
                if (std::holds_alternative<ptr_with_off_t>(dst_ptr)) {
                    auto dst_ptr_with_off = std::get<ptr_with_off_t>(dst_ptr);
                    auto src_ptr_with_off = std::get<ptr_with_off_t>(src_ptr);
                    subtracted = dst_ptr_with_off.get_offset().to_interval() -
                                src_ptr_with_off.get_offset().to_interval();
                }
                else if (std::holds_alternative<packet_ptr_t>(dst_ptr)) {
                    register_t src_reg = std::get<Reg>(bin.v).v;
                    register_t dst_reg = bin.dst.v;
                    subtracted = m_offset.compute_packet_subtraction(dst_reg, src_reg);
                }
                else {
                    // This should not happen as same_region only allows non-shared pointers
                    m_errors.push_back("subtraction between pointers of different region");
                }
            }
            m_region -= dst_register;
            m_offset -= dst_register;
        }
    }

    m_interval.do_bin(bin, subtracted, loc);
    m_region.do_bin(bin, dst_signed_interval, src_signed_interval, loc);
    m_offset.do_bin(bin, dst_signed_rf, src_signed_rf, loc);
}

void inference_domain_t::do_load(const Mem& b, const Reg& target_reg, bool unknown_ptr,
        std::optional<ptr_or_mapfd_t> basereg_opt, location_t loc) {
    m_region.do_load(b, register_t{target_reg.v}, unknown_ptr, loc);
    // TODO: replace with a bool value returned from region do_load
    auto load_in_region = m_region.find_ptr_or_mapfd_type(target_reg.v).has_value();
    m_interval.do_load(b, register_t{target_reg.v}, basereg_opt, load_in_region, loc);
    m_offset.do_load(b, register_t{target_reg.v}, basereg_opt, loc);
}

void inference_domain_t::do_mem_store(const Mem& b, std::optional<ptr_or_mapfd_t>& basereg_opt) {
    m_region.do_mem_store(b);
    // TODO: Before storing into interval domain, check if there is no overlap with packet pointers
    m_interval.do_mem_store(b, basereg_opt);
    m_offset.do_mem_store(b, basereg_opt);
}

void inference_domain_t::operator()(const Mem& b, location_t loc) {
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
        if (b.is_load) do_load(b, targetreg, unknown_ptr, base_ptr_or_mapfd_opt, loc);
        else if (!unknown_ptr) do_mem_store(b, base_ptr_or_mapfd_opt);
    }
    else if (!unknown_ptr && !b.is_load) {
        do_mem_store(b, base_ptr_or_mapfd_opt);
    }
}

void inference_domain_t::print_ctx(std::ostream& o) const {
    const std::vector<uint64_t>& ctx_keys = m_region.get_ctx_keys();
    o << "\tctx: {";
    for (auto const& k : ctx_keys) {
        auto dist = m_offset.find_in_ctx(k);
        if (dist) {
            o << "\t\t";
            print_non_numeric_memory_cell(o, k, k+3, packet_ptr_t{}, dist);
            o << ",\n";
        }
    }
    o << "\t}\n";
}

void inference_domain_t::print_stack(std::ostream& o) const {
    const std::vector<uint64_t>& stack_keys_region = m_region.get_stack_keys();
    const std::vector<uint64_t>& stack_keys_interval = m_interval.get_stack_keys();
    o << "\tstack: {\n";
    for (auto const& k : stack_keys_region) {
        auto maybe_ptr_or_mapfd_cells = m_region.find_in_stack(k);
        auto dist = m_offset.find_in_stack(k);
        if (maybe_ptr_or_mapfd_cells) {
            auto ptr_or_mapfd_cells = maybe_ptr_or_mapfd_cells.value();
            int width = ptr_or_mapfd_cells.second;
            auto ptr_or_mapfd = ptr_or_mapfd_cells.first;
            o << "\t\t";
            if (dist) {
                print_non_numeric_memory_cell(o, k, k+width-1, std::move(ptr_or_mapfd),
                        std::optional<refinement_t>(dist->first));
            }
            else {
                print_non_numeric_memory_cell(o, k, k+width-1, std::move(ptr_or_mapfd));
            }
            o << ",\n";
        }
    }
    for (auto const& k : stack_keys_interval) {
        auto maybe_signed_interval_cells = m_interval.find_in_stack_signed(k);
        if (maybe_signed_interval_cells) {
            auto interval_cells = maybe_signed_interval_cells.value();
            o << "\t\t";
            print_numeric_memory_cell(o, k, k+interval_cells.second-1,
                    interval_cells.first, true);
            o << ",\n";
        }
        auto maybe_unsigned_interval_cells = m_interval.find_in_stack_unsigned(k);
        if (maybe_unsigned_interval_cells) {
            auto interval_cells = maybe_unsigned_interval_cells.value();
            o << "\t\t";
            print_numeric_memory_cell(o, k, k+interval_cells.second-1,
                    interval_cells.first, false);
            o << ",\n";
        }
    }
    o << "\t}\n";
}

void inference_domain_t::adjust_bb_for_types(location_t loc) {
    m_region.adjust_bb_for_types(loc);
    m_offset.adjust_bb_for_types(loc);
    m_interval.adjust_bb_for_types(loc);
}

void inference_domain_t::operator()(const basic_block_t& bb) {

    // A temporary fix to avoid printing errors for multiple basic blocks
    m_errors.clear();
    m_region.reset_errors();
    m_offset.reset_errors();
    m_interval.reset_errors();

    auto label = bb.label();
    uint32_t curr_pos = 0;
    location_t loc{label, curr_pos};
    adjust_bb_for_types(loc);

    for (const Instruction& statement : bb) {
        loc = location_t(label, ++curr_pos);
        std::visit([this, loc](const auto& v) { std::apply(*this, std::make_tuple(v, loc)); }, statement);
    }

    operator+=(m_region.get_errors());
    operator+=(m_offset.get_errors());
    operator+=(m_interval.get_errors());
}

std::optional<crab::ptr_or_mapfd_t>
inference_domain_t::find_ptr_or_mapfd_at_loc(const crab::register_location_t& loc) const {
    return m_region.find_ptr_or_mapfd_at_loc(loc);
}

std::optional<crab::refinement_t>
inference_domain_t::find_refinement_at_loc(const crab::register_location_t& loc) const {
    return m_offset.find_refinement_at_loc(loc);
}

std::optional<crab::refinement_t>
inference_domain_t::find_signed_interval_at_loc(const crab::register_location_t& loc) const {
    return m_interval.find_signed_interval_at_loc(loc);
}

std::optional<crab::refinement_t>
inference_domain_t::find_unsigned_interval_at_loc(const crab::register_location_t& loc) const {
    return m_interval.find_unsigned_interval_at_loc(loc);
}

static inline region_t string_to_region(const std::string& s) {
    static std::map<std::string, region_t> string_to_region{
        {std::string("ctx"), region_t::R_CTX},
        {std::string("stack"), region_t::R_STACK},
        {std::string("packet"), region_t::R_PACKET},
        {std::string("shared"), region_t::R_SHARED},
    };
    if (string_to_region.count(s)) {
        return string_to_region[s];
    }
    throw std::runtime_error(std::string("Unsupported region name: ") + s);
}

void inference_domain_t::insert_in_registers_in_interval_domain(register_t r, location_t loc,
                                                            refinement_t rf) {
    m_interval.insert_in_registers(r, loc, rf);
}

void inference_domain_t::insert_in_registers_in_signed_interval_domain(register_t r,
                                                            location_t loc, refinement_t rf) {
    m_interval.insert_in_registers_signed(r, loc, rf);
}

void inference_domain_t::insert_in_registers_in_unsigned_interval_domain(register_t r,
                                                            location_t loc, refinement_t rf) {
    m_interval.insert_in_registers_unsigned(r, loc, rf);
}

void inference_domain_t::store_in_stack_in_interval_domain(uint64_t key, refinement_t p, int width) {
    m_interval.store_in_stack(key, p, width);
}

void inference_domain_t::store_in_stack_in_signed_interval_domain(uint64_t key, refinement_t p,
                                                            int width) {
    m_interval.store_in_stack_signed(key, p, width);
}

void inference_domain_t::store_in_stack_in_unsigned_interval_domain(uint64_t key, refinement_t p,
                                                            int width) {
    m_interval.store_in_stack_unsigned(key, p, width);
}

void inference_domain_t::insert_in_registers_in_offset_domain(register_t r, location_t loc,
                                                             refinement_t d) {
    m_offset.insert_in_registers(r, loc, d);
}

void inference_domain_t::store_in_stack_in_offset_domain(uint64_t key, refinement_t d, int width) {
    m_offset.store_in_stack(key, d, width);
}

void inference_domain_t::insert_in_registers_in_region_domain(register_t r, location_t loc,
                                                             const ptr_or_mapfd_t& p) {
    m_region.insert_in_registers(r, loc, p);
}

void inference_domain_t::store_in_stack_in_region_domain(uint64_t key, ptr_or_mapfd_t p, int width) {
    m_region.store_in_stack(key, p, width);
}

inference_domain_t inference_domain_t::from_predefined_types(const std::set<std::string>& types,
                                                             bool setup_constraints) {
    // TODO: redo the method according to the new offset domain
    // also, need to store the intervals in the offset domain
    using std::regex;
    using std::regex_match;

    #define NUMERIC R"_(\s*\[?([-+]?(?:\d+|oo))(?:,\s*([-+]?(?:\d+|oo))\])?\s*)_"
    #define NUMERIC_ENCLOSED "<" NUMERIC ">"
    #define NUMERIC_NUMERIC_ENCLOSED "<" NUMERIC ",\\s*" NUMERIC ">"
    #define REG R"_(\s*r(\d\d?)\s*)_"
    #define STACK_CELL R"_(\s*stack\[(\d+)-(\d+)\]\s*)_"
    #define SHARED_PTR "\\s*shared_p(?:" NUMERIC_NUMERIC_ENCLOSED ")?\\s*"
    #define CTX_OR_STACK_PTR "\\s*(ctx|stack)_p(?:" NUMERIC_ENCLOSED ")?\\s*"
    #define PACKET_PTR "\\s*packet_p(?:<(begin|end|meta)\\+" NUMERIC ">)?\\s*"
    #define SNUMBER "\\s*snumber(?:" NUMERIC_ENCLOSED ")?\\s*"
    #define UNUMBER "\\s*unumber(?:" NUMERIC_ENCLOSED ")?\\s*"
    #define MAPFD "\\s*(map_fd|map_fd_programs)" NUMERIC "\\s*"

    auto create_interval = [](std::string lb, std::string ub) {
        if (lb == "" && ub == "") {
            return crab::mock_interval_t::top();
        }
        bound_t lb_num = bound_t::minus_infinity();
        if (lb != "-oo") {
            try {
                lb_num = bound_t{number_t{static_cast<int64_t>(std::stoll(lb))}};
            } catch (std::out_of_range& e) {
                // TODO: Separate handling for such cases
                lb_num = bound_t{number_t{static_cast<uint64_t>(std::stoull(lb))}};
            }
        }
        auto ub_num = lb_num;
        if (ub != "" && ub != "+oo") {
            try {
                ub_num = bound_t{number_t{static_cast<int64_t>(std::stoll(ub))}};
            }
            catch (std::out_of_range& e) {
                // TODO: Separate handling for such cases
                ub_num = bound_t{number_t{static_cast<uint64_t>(std::stoull(ub))}};
            }
        }
        else if (ub == "+oo") ub_num = bound_t::plus_infinity();
        return crab::mock_interval_t{lb_num, ub_num};
    };

    auto create_ptr = [create_interval](std::string region, std::string off_lb,
            std::string off_ub, std::string region_sz_lb = "", std::string region_sz_ub = "") {
        auto region_type = string_to_region(region);
        auto mock_offset = create_interval(off_lb, off_ub);
        auto mock_region_size = create_interval(region_sz_lb, region_sz_ub);
        return crab::ptr_with_off_t{region_type, -1, mock_offset, nullness_t::MAYBE_NULL,
            mock_region_size};
    };

    auto create_mapfd = [create_interval](std::string mapfd_type, std::string lb_mapfd,
            std::string ub_mapfd) {
        auto interval = create_interval(lb_mapfd, ub_mapfd);
        if (mapfd_type == "map_fd_programs") {
            return crab::mapfd_t(interval, EbpfMapValueType::PROGRAM);
        }
        else {
            return crab::mapfd_t(interval, EbpfMapValueType::MAP);
        }
    };

    /*
    auto create_pkt_offset = [create_interval](std::string offset_type, std::string offset_lb,
            std::string offset_ub) {
        auto offset = create_interval(offset_lb, offset_ub).to_interval();
        if (offset_type == "begin") {
            return dist_t{offset};
        }
        else if (offset_type == "end") {
            auto packet_end = crab::interval_t{number_t{PACKET_END}};
            return dist_t{packet_end - offset};
        }
        else {
            auto packet_meta = crab::interval_t{number_t{PACKET_META}};
            return dist_t{packet_meta - offset};
        }
    };
    */

    // TODO: Incomplete implementation, needs work
    inference_domain_t typ;
    if (setup_constraints) {
        typ = inference_domain_t::setup_entry(false);
    }
    else {
        typ.set_to_top();
    }
    location_t loc{label_t::entry, 0};
    for (const auto& t : types) {
        std::smatch m;
        if (regex_match(t, m, regex(REG ":" CTX_OR_STACK_PTR))) {
            auto reg = register_t{static_cast<uint8_t>(std::stoul(m[1]))};
            auto ptr = create_ptr(m[2], m[3], m[4]);
            typ.insert_in_registers_in_region_domain(reg, loc, ptr);
        }
        else if (regex_match(t, m, regex(REG ":" PACKET_PTR))) {
            /*
            auto reg = register_t{static_cast<uint8_t>(std::stoul(m[1]))};
            auto ptr = packet_ptr_t{};
            auto offset = create_pkt_offset(m[2], m[3], m[4]);
            typ.insert_in_registers_in_region_domain(reg, loc, ptr);
            typ.insert_in_registers_in_offset_domain(reg, loc, offset);
            */
        }
        else if (regex_match(t, m, regex(REG ":" SHARED_PTR))) {
            auto reg = register_t{static_cast<uint8_t>(std::stoul(m[1]))};
            auto ptr = create_ptr("shared", m[2], m[3], m[4], m[5]);
            typ.insert_in_registers_in_region_domain(reg, loc, ptr);
        }
        else if (regex_match(t, m, regex(REG ":" SNUMBER))) {
            auto reg = register_t{static_cast<uint8_t>(std::stoul(m[1]))};
            auto num = create_interval(m[2], m[3]).to_interval();
            //typ.insert_in_registers_in_signed_interval_domain(reg, loc, num);
        }
        else if (regex_match(t, m, regex(REG ":" UNUMBER))) {
            auto reg = register_t{static_cast<uint8_t>(std::stoul(m[1]))};
            auto num = create_interval(m[2], m[3]).to_interval();
            //typ.insert_in_registers_in_unsigned_interval_domain(reg, loc, num);
        }
        else if (regex_match(t, m, regex(REG ":" MAPFD))) {
            auto reg = register_t{static_cast<uint8_t>(std::stoul(m[1]))};
            auto mapfd = create_mapfd(m[2], m[3], m[4]);
            typ.insert_in_registers_in_region_domain(reg, loc, mapfd);
        }
        else if (regex_match(t, m, regex(STACK_CELL ":" CTX_OR_STACK_PTR))) {
            auto stack_cell_start = static_cast<uint64_t>(std::stoul(m[1]));
            auto stack_cell_end = static_cast<uint64_t>(std::stoul(m[2]));
            auto ptr = create_ptr(m[3], m[4], m[5]);
            typ.store_in_stack_in_region_domain(stack_cell_start, ptr,
                    stack_cell_end-stack_cell_start);
        }
        else if (regex_match(t, m, regex(STACK_CELL ":" PACKET_PTR))) {
            /*
            auto stack_cell_start = static_cast<uint64_t>(std::stoul(m[1]));
            auto stack_cell_end = std::stoi(m[2]);
            auto ptr = packet_ptr_t{};
            auto pkt_offset = create_pkt_offset(m[3], m[4], m[5]);
            int width = stack_cell_end - stack_cell_start;
            typ.store_in_stack_in_region_domain(stack_cell_start, ptr, width);
            typ.store_in_stack_in_offset_domain(stack_cell_start, pkt_offset, width);
            */
        }
        else if (regex_match(t, m, regex(STACK_CELL ":" SHARED_PTR))) {
            auto stack_cell_start = static_cast<uint64_t>(std::stoul(m[1]));
            auto stack_cell_end = std::stoi(m[2]);
            auto ptr = create_ptr("shared", m[3], m[4], m[5], m[6]);
            typ.store_in_stack_in_region_domain(stack_cell_start, ptr,
                    stack_cell_end-stack_cell_start);
        }
        else if (regex_match(t, m, regex(STACK_CELL ":" SNUMBER))) {
            auto stack_cell_start = static_cast<uint64_t>(std::stoul(m[1]));
            auto stack_cell_end = std::stoi(m[2]);
            auto num = create_interval(m[3], m[4]);
            //typ.store_in_stack_in_signed_interval_domain(stack_cell_start, num,
            //        stack_cell_end-stack_cell_start);
        }
        else if (regex_match(t, m, regex(STACK_CELL ":" UNUMBER))) {
            auto stack_cell_start = static_cast<uint64_t>(std::stoul(m[1]));
            auto stack_cell_end = std::stoi(m[2]);
            auto num = create_interval(m[3], m[4]);
            //typ.store_in_stack_in_unsigned_interval_domain(stack_cell_start, num,
            //        stack_cell_end-stack_cell_start);
        }
        else if (regex_match(t, m, regex(STACK_CELL ":" MAPFD))) {
            auto stack_cell_start = static_cast<uint64_t>(std::stoul(m[1]));
            auto stack_cell_end = std::stoi(m[2]);
            auto mapfd = create_mapfd(m[3], m[4], m[5]);
            typ.store_in_stack_in_region_domain(stack_cell_start, mapfd,
                    stack_cell_end-stack_cell_start);
        }
        else {
            std::cout << "type not recognized: " << t << "\n";
        }
    }
    return typ;
}

void inference_domain_t::write(std::ostream& os) const {
    os << to_set();
}

std::ostream& operator<<(std::ostream& o, const inference_domain_t& typ) {
    typ.write(o);
    return o;
}

void inference_domain_t::print_annotated_bb(std::ostream& o, const basic_block_t& bb) const {
    if (is_bottom()) {
        print_bb(o, bb);
        return;
    }

    o << bb.label() << ":\n";
    uint32_t curr_pos = 0;
    // TODO: add support for printing unsigned intervals as well
    for (const Instruction& statement : bb) {
        ++curr_pos;
        crab::location_t loc{bb.label(), curr_pos};
        o << "   " << curr_pos << ".";
        // TODO: print unsigned intervals in a proper way
        if (std::holds_alternative<Call>(statement)) {
            auto r0_reg = crab::register_location_t(register_t{R0_RETURN_VALUE}, loc);
            auto region = find_ptr_or_mapfd_at_loc(r0_reg);
            auto rf = find_refinement_at_loc(r0_reg);
            auto signed_interval = find_signed_interval_at_loc(r0_reg);
            print_annotated(o, std::get<Call>(statement), region, rf, signed_interval, true);
            auto unsigned_interval = find_unsigned_interval_at_loc(r0_reg);
            //print_annotated(o, std::get<Call>(statement), region, unsigned_interval, false);
        }
        else if (std::holds_alternative<Bin>(statement)) {
            auto b = std::get<Bin>(statement);
            auto register_location = crab::register_location_t(b.dst.v, loc);
            auto region = find_ptr_or_mapfd_at_loc(register_location);
            auto rf = find_refinement_at_loc(register_location);
            auto signed_interval = find_signed_interval_at_loc(register_location);
            print_annotated(o, b, region, rf, signed_interval, true);
            auto unsigned_interval = find_unsigned_interval_at_loc(register_location);
            //print_annotated(o, b, region, rf, unsigned_interval, false);
        }
        else if (std::holds_alternative<Mem>(statement)) {
            auto u = std::get<Mem>(statement);
            if (u.is_load) {
                auto target_reg = std::get<Reg>(u.value);
                auto target_reg_loc = crab::register_location_t(target_reg.v, loc);
                auto region = find_ptr_or_mapfd_at_loc(target_reg_loc);
                auto rf = find_refinement_at_loc(target_reg_loc);
                auto signed_interval = find_signed_interval_at_loc(target_reg_loc);
                print_annotated(o, u, region, rf, signed_interval, true);
                auto unsigned_interval = find_unsigned_interval_at_loc(target_reg_loc);
                //print_annotated(o, u, region, rf, unsigned_interval, false);
            }
            else print_instr(o, u);
        }
        else if (std::holds_alternative<LoadMapFd>(statement)) {
            auto u = std::get<LoadMapFd>(statement);
            auto reg = crab::register_location_t(u.dst.v, loc);
            auto region = find_ptr_or_mapfd_at_loc(reg);
            print_annotated(o, u, region);
        }
        else if (std::holds_alternative<Un>(statement)) {
            auto u = std::get<Un>(statement);
            auto reg = crab::register_location_t(u.dst.v, loc);
            auto signed_interval = find_signed_interval_at_loc(reg);
            print_annotated(o, u, signed_interval, true);
            auto unsigned_interval = find_unsigned_interval_at_loc(reg);
            print_annotated(o, u, unsigned_interval, false);
        }
        else print_instr(o, statement);
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

} // namespace crab
