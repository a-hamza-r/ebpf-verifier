// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "offset_domain.hpp"

namespace crab {

void offset_registers_t::insert(register_t reg, const location_t& loc, refinement_t rf) {
    register_location_t register_location{reg, loc};
    m_registers_env->insert_or_assign(register_location, rf);
    m_cur_register_def[reg] = std::make_shared<register_location_t>(register_location);
}

std::optional<refinement_t> offset_registers_t::find(register_location_t reg) const {
    auto it = m_registers_env->find(reg);
    if (it == m_registers_env->end()) return {};
    return it->second;
}

std::optional<refinement_t> offset_registers_t::find(register_t key) const {
    if (m_cur_register_def[key] == nullptr) return {};
    return find(*(m_cur_register_def[key]));
}

std::vector<uint64_t> offset_stack_t::find_overlapping_cells(uint64_t start, int width) const {
    std::vector<uint64_t> overlapping_cells;
    auto it = m_stack_cells.begin();
    while (it != m_stack_cells.end() && it->first < start) {
        it++;
    }
    if (it != m_stack_cells.begin()) {
        it--;
        auto key = it->first;
        auto width_key = it->second.second;
        if (key < start && key+width_key > start) overlapping_cells.push_back(key);
    }

    for (; it != m_stack_cells.end(); it++) {
        auto key = it->first;
        if (key >= start && key < start+width) overlapping_cells.push_back(key);
        if (key >= start+width) break;
    }
    return overlapping_cells;
}

void offset_registers_t::set_to_top() {
    m_registers_env = std::make_shared<global_env_offset_registers_t>();
    m_cur_register_def = live_registers_t{nullptr};
    m_is_bottom = false;
}

void offset_registers_t::set_to_bottom() {
    m_cur_register_def = live_registers_t{nullptr};
    m_is_bottom = true;
}

bool offset_registers_t::is_top() const {
    if (m_is_bottom) return false;
    if (m_registers_env == nullptr) return true;
    for (auto &it : m_cur_register_def) {
        if (it != nullptr) return false;
    }
    return true;
}

bool offset_registers_t::is_bottom() const {
    return m_is_bottom;
}

void offset_registers_t::operator-=(register_t to_forget) {
    if (is_bottom()) {
        return;
    }
    m_cur_register_def[to_forget] = nullptr;
}

bool offset_registers_t::inclusion(const offset_registers_t& other,
                                   std::shared_ptr<slacks_t> slacks) const {
    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (other.m_cur_register_def[i] == nullptr) continue;
        if (m_cur_register_def[i] == nullptr) return false;
        auto it1 = find(*(m_cur_register_def[i]));
        auto it2 = other.find(*(other.m_cur_register_def[i]));
        if (it1 && it2) {
            if (!(it1->inclusion(*it2, slacks))) return false;
        }
    }
    return true;
}

offset_registers_t offset_registers_t::join(const offset_registers_t& other,
                                            std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }

    offset_registers_t joined_state;
    location_t loc = location_t::top();

    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (m_cur_register_def[i] == nullptr || other.m_cur_register_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_register_def[i]));
        auto it2 = other.find(*(other.m_cur_register_def[i]));
        if (it1 && it2) {
            auto rf1 = *it1, rf2 = *it2;
            if (rf1.same_type(rf2)) {
                joined_state.insert(register_t{i}, loc, rf1.join(rf2, slacks));
            }
        }
    }
    return joined_state;
}

offset_registers_t offset_registers_t::widen(const offset_registers_t& other,
                                             std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }

    offset_registers_t joined_state;
    location_t loc = location_t::top();

    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (other.m_cur_register_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_register_def[i]));
        auto it2 = other.find(*(other.m_cur_register_def[i]));
        if (it1 && it2) {
            auto rf1 = *it1, rf2 = *it2;
            if (rf1.same_type(rf2)) {
                joined_state.insert(register_t{i}, loc, rf1.widen(rf2, slacks));
            }
        }
    }
    return joined_state;
}

void offset_registers_t::adjust_bb_for_registers(location_t loc) {
    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (auto it = find(register_t{i})) {
            insert(register_t{i}, loc, *it);
        }
    }
}

void offset_registers_t::scratch_caller_saved_registers() {
    for (uint8_t r = R1_ARG; r <= R5_ARG; r++) {
        operator-=(register_t{r});
    }
}

void offset_registers_t::forget_packet_pointers(location_t loc) {
    for (uint8_t r = R0_RETURN_VALUE; r < NUM_REGISTERS-2; r++) {
        if (auto it = find(register_t{r})) {
            if (it->is_packet_refinement()) {
                operator-=(register_t{r});
            }
        }
    }
    // Reset the constraints on the packet
    insert(R12_PKT_BEGIN, loc, refinement_t::begin_with_constraints());
}

void offset_stack_t::set_to_top() {
    m_stack_cells.clear();
    m_is_bottom = false;
}

void offset_stack_t::set_to_bottom() {
    m_stack_cells.clear();
    m_is_bottom = true;
}

bool offset_stack_t::is_top() const {
    if (m_is_bottom) return false;
    return m_stack_cells.empty();
}

bool offset_stack_t::is_bottom() const {
    return m_is_bottom;
}

offset_stack_t offset_stack_t::top() {
    return offset_stack_t();
}

std::optional<refinement_stack_cell_t> offset_stack_t::find(uint64_t key) const {
    auto it = m_stack_cells.find(key);
    if (it == m_stack_cells.end()) return {};
    return it->second;
}

void offset_stack_t::store(uint64_t key, refinement_t d, int width) {
    m_stack_cells[key] = std::make_pair(d, width);
}

std::vector<uint64_t> offset_stack_t::get_keys() const {
    std::vector<uint64_t> keys;
    keys.reserve(m_stack_cells.size());

    for (auto const& kv : m_stack_cells) {
        keys.push_back(kv.first);
    }
    return keys;
}

void offset_stack_t::operator-=(uint64_t to_erase) {
    if (is_bottom()) {
        return;
    }
    m_stack_cells.erase(to_erase);
}

void offset_stack_t::operator-=(const std::vector<uint64_t>& keys) {
    for (auto &key : keys) {
       *this -= key;
    }
}

bool offset_stack_t::inclusion(const offset_stack_t& other, std::shared_ptr<slacks_t> slacks) const {
    size_t size1 = m_stack_cells.size();
    size_t size2 = other.m_stack_cells.size();
    if (size2 > size1) return false;
    for (auto const &kv : other.m_stack_cells) {
        auto it = m_stack_cells.find(kv.first);
        if (it == m_stack_cells.end()) return false;
        auto rf1 = it->second.first;
        auto rf2 = kv.second.first;
        auto width1 = it->second.second;
        auto width2 = kv.second.second;
        if (!(rf1.inclusion(rf2, slacks) && width1 == width2)) return false;
    }
    return true;
}

offset_stack_t offset_stack_t::join(const offset_stack_t& other,
                                    std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }

    refinement_stack_cells_t out_stack_rfs;
    // We do not join rf cells because different rf values different types of offsets
    for (auto const&kv: m_stack_cells) {
        auto maybe_rf_cells = other.find(kv.first);
        if (maybe_rf_cells) {
            auto rf_cells1 = kv.second;
            auto rf_cells2 = *maybe_rf_cells;
            auto rf1 = rf_cells1.first;
            auto rf2 = rf_cells2.first;
            int width1 = rf_cells1.second;
            int width2 = rf_cells2.second;
            if (rf1.same_type(rf2) && width1 == width2) {
                out_stack_rfs.insert({kv.first, std::make_pair(rf1.join(rf2, slacks), width1)});
            }
        }
    }
    return offset_stack_t(std::move(out_stack_rfs));
}

offset_stack_t offset_stack_t::widen(const offset_stack_t& other,
                                     std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }

    refinement_stack_cells_t out_stack_rfs;
    // We do not join rf cells because different rf values different types of offsets
    for (auto const&kv: m_stack_cells) {
        auto maybe_rf_cells = other.find(kv.first);
        if (maybe_rf_cells) {
            auto rf_cells1 = kv.second;
            auto rf_cells2 = *maybe_rf_cells;
            auto rf1 = rf_cells1.first;
            auto rf2 = rf_cells2.first;
            int width1 = rf_cells1.second;
            int width2 = rf_cells2.second;
            if (rf1.same_type(rf2) && width1 == width2) {
                out_stack_rfs.insert({kv.first, std::make_pair(rf1.widen(rf2, slacks), width1)});
            }
        }
    }
    return offset_stack_t(std::move(out_stack_rfs));
}

offset_ctx_t::offset_ctx_t(const ebpf_context_descriptor_t* desc) {
    if (desc == nullptr) return;
    if (desc->data >= 0) {
        m_ctx_cells.insert_or_assign(desc->data, refinement_t::begin());
    }
    if (desc->end >= 0) {
        m_ctx_cells.insert_or_assign(desc->end, refinement_t::end());
    }
    if (desc->meta >= 0) {
        m_ctx_cells.insert_or_assign(desc->meta, refinement_t::meta());
    }
    m_size = std::max(0, desc->size);
}

std::vector<uint64_t> offset_ctx_t::get_keys() const {
    std::vector<uint64_t> keys;
    keys.reserve(m_ctx_cells.size());

    for (auto const& kv : m_ctx_cells) {
        keys.push_back(kv.first);
    }
    return keys;
}

std::optional<refinement_t> offset_ctx_t::find(uint64_t key) const {
    auto it = m_ctx_cells.find(key);
    if (it == m_ctx_cells.end()) return {};
    return it->second;
}

offset_domain_t offset_domain_t::setup_entry(std::shared_ptr<slacks_t> slacks) {
    location_t loc{label_t::entry, 0};
    offset_registers_t offset_regs;
    const auto r12_pkt_begin = register_t{R12_PKT_BEGIN};
    const auto type_begin_rf = refinement_t::begin_with_constraints();
    offset_regs.insert(r12_pkt_begin, loc, type_begin_rf);

    return offset_domain_t{std::move(offset_regs), offset_stack_t::top(),
            std::make_shared<offset_ctx_t>(global_program_info->type.context_descriptor),
            slacks};
}

offset_domain_t offset_domain_t::bottom() {
    offset_domain_t off;
    off.set_to_bottom();
    return off;
}

void offset_domain_t::set_to_top() {
    m_is_bottom = false;
    m_registers.set_to_top();
    m_stack.set_to_top();
}

void offset_domain_t::set_to_bottom() {
    m_is_bottom = true;
    m_registers.set_to_bottom();
    m_stack.set_to_bottom();
}

bool offset_domain_t::is_bottom() const {
    if (m_is_bottom) return true;
    return (m_registers.is_bottom() || m_stack.is_bottom());
}

bool offset_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_registers.is_top() && m_stack.is_top());
}

// inclusion
bool offset_domain_t::operator<=(const offset_domain_t& other) const {
    return (m_registers.inclusion(other.m_registers, m_slacks)
    && m_stack.inclusion(other.m_stack, m_slacks));
}

// join
void offset_domain_t::operator|=(const offset_domain_t& abs) {
    offset_domain_t tmp{abs};
    operator|=(std::move(tmp));
}

void offset_domain_t::operator|=(offset_domain_t&& abs) {
    if (is_bottom()) {
        *this = abs;
        return;
    }
    *this = *this | std::move(abs);
}

offset_domain_t offset_domain_t::operator|(const offset_domain_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return offset_domain_t(m_registers.join(other.m_registers, m_slacks),
                           m_stack.join(other.m_stack, m_slacks), m_ctx, m_slacks);
}

offset_domain_t offset_domain_t::operator|(offset_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    // TODO: check if slacks need to be moved
    return offset_domain_t(
        m_registers.join(std::move(other.m_registers), m_slacks),
        m_stack.join(std::move(other.m_stack), m_slacks),
        std::move(other.m_ctx),
        std::move(other.m_slacks)
    );
}

// meet
offset_domain_t offset_domain_t::operator&(const offset_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

// widening
offset_domain_t offset_domain_t::widen(const offset_domain_t& other, bool to_constants) {
    if (is_bottom() || other.is_top()) {
        return other;
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return offset_domain_t(m_registers.widen(other.m_registers, m_slacks),
            m_stack.widen(other.m_stack, m_slacks), m_ctx, m_slacks);
}

// narrowing
offset_domain_t offset_domain_t::narrow(const offset_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

void offset_domain_t::write(std::ostream& os) const {}

std::string offset_domain_t::domain_name() const {
    return "offset_domain";
}

crab::bound_t offset_domain_t::get_loop_count_upper_bound() const {
    /* WARNING: The operation is not implemented yet.*/
    return crab::bound_t{crab::number_t{0}};
}

void offset_domain_t::initialize_loop_counter(const label_t& label) {
    /* WARNING: The operation is not implemented yet.*/
}

string_invariant offset_domain_t::to_set() { return string_invariant{}; }

void offset_domain_t::operator()(const Assume &b, location_t loc) {
    Condition cond = b.cond;
    if (std::holds_alternative<Reg>(cond.right)) {
        auto right_reg = std::get<Reg>(cond.right).v;
        auto rf_left = m_registers.find(cond.left.v);
        auto rf_right = m_registers.find(right_reg);
        if (!rf_left || !rf_right) {
            // this should not happen, comparison between a packet pointer and either
            // other region's pointers or numbers; possibly raise type error
            m_errors.push_back(loc.to_string() + ": Both registers are not packet pointers.");
            return;
        }
        register_t pkt_begin{R12_PKT_BEGIN};
        auto begin = m_registers.find(pkt_begin);
        if (cond.op == Condition::Op::LE) {
            begin->add_constraint(rf_left->assume_le(*rf_right), m_slacks);
            m_registers.insert(pkt_begin, loc, std::move(*begin));
        }
        else if (cond.op == Condition::Op::GT) {
            begin->add_constraint(rf_left->assume_gt(*rf_right), m_slacks);
            m_registers.insert(pkt_begin, loc, std::move(*begin));
        }
        else if (cond.op == Condition::Op::GE) {
            begin->add_constraint(rf_right->assume_le(*rf_left), m_slacks);
            m_registers.insert(pkt_begin, loc, std::move(*begin));
        }
        else if (cond.op == Condition::Op::LT) {
            begin->add_constraint(rf_right->assume_gt(*rf_left), m_slacks);
            m_registers.insert(pkt_begin, loc, std::move(*begin));
        }
        // other comparisons not supported

        // To keep the info about the packet pointer, we store the types here again
        m_registers.insert(register_t{cond.left.v}, loc, *rf_left);
        m_registers.insert(register_t{right_reg}, loc, *rf_right);
    }
}

interval_t offset_domain_t::compute_packet_subtraction(register_t dst, register_t src) const {
    auto dst_rf = m_registers.find(dst);
    auto src_rf = m_registers.find(src);
    if (!dst_rf || !src_rf) return interval_t::bottom();
    expression_t dst_expr = dst_rf->get_value();
    expression_t src_expr = src_rf->get_value();
    refinement_t result_rf = dst_rf->subtract(*src_rf, m_slacks);
    expression_t result_expr = result_rf.get_value().get_equivalent_expression(m_slacks);
    if (result_expr.is_constant()) {
        return result_expr.get_constant_term();
    }
    // with non-singleton expressions, we might be able to compute subtraction,
    // but it might be complicated
    if (!dst_expr.is_singleton() || !src_expr.is_singleton()) return interval_t::top();
    auto dst_symbol = dst_expr.get_singleton();
    auto src_symbol = src_expr.get_singleton();
    std::optional<refinement_t> begin_rf = m_registers.find(register_t{R12_PKT_BEGIN});
    return begin_rf->simplify_for_subtraction(dst_symbol, src_symbol, m_slacks);
}

void offset_domain_t::do_bin(const Bin& bin, std::optional<refinement_t> dst_rf_numeric_opt,
                             std::optional<refinement_t> src_rf_numeric_opt, location_t loc) {

    using Op = Bin::Op;

    auto dst_register = register_t{bin.dst.v};

    if (auto pimm = std::get_if<Imm>(&bin.v)) {
        int64_t imm;
        if (bin.is64) {
            // Use the full signed value.
            imm = to_signed(pimm->v);
        } else {
            // Use only the low 32 bits of the value.
            imm = gsl::narrow_cast<int32_t>(imm);
        }
        auto imm_interval = interval_t{imm};
        switch (bin.op) {
            case Op::MOV: {
                // ra = imm
                m_registers -= dst_register;
                break;
            }
            case Op::ADD: {
                // ra += imm
                if (imm == 0) break;
                if (auto dst_rf_ptr_opt = m_registers.find(dst_register)) {
                    m_registers.insert(dst_register, loc, *dst_rf_ptr_opt + imm_interval);
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::SUB: {
                // ra -= imm
                if (imm == 0) break;
                if (auto dst_rf_ptr_opt = m_registers.find(dst_register)) {
                    m_registers.insert(dst_register, loc, *dst_rf_ptr_opt + (-imm_interval));
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            default: {
                // no other operations supported for packet pointers in the offset domain
                m_registers -= dst_register;
                break;
            }
        }
    }
    else {
        auto src = std::get<Reg>(bin.v);
        switch (bin.op) {
            case Op::MOV: {
                // ra = rb
                if (auto src_rf_ptr_opt = m_registers.find(src.v)) {
                    m_registers.insert(dst_register, loc, *src_rf_ptr_opt);
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::ADD: {
                // ra += rb
                auto dst_rf_ptr_opt = m_registers.find(dst_register);
                auto src_rf_ptr_opt = m_registers.find(src.v);
                if (dst_rf_ptr_opt.has_value() && src_rf_ptr_opt.has_value()) {
                    // possibly adding two pointers
                    set_to_bottom();
                }
                else if (dst_rf_ptr_opt.has_value() && src_rf_numeric_opt.has_value()) {
                    m_registers.insert(dst_register, loc, dst_rf_ptr_opt->add(*src_rf_numeric_opt, m_slacks));
                }
                else if (dst_rf_numeric_opt.has_value() && src_rf_ptr_opt.has_value()) {
                    m_registers.insert(dst_register, loc, dst_rf_numeric_opt->add(*src_rf_ptr_opt, m_slacks));
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::SUB: {
                // ra -= rb
                auto dst_rf_ptr_opt = m_registers.find(dst_register);
                auto src_rf_ptr_opt = m_registers.find(src.v);
                if (dst_rf_ptr_opt.has_value() && src_rf_ptr_opt.has_value()) {
                    // possibly subtracting two pointers
                    m_registers -= dst_register;
                }
                else if (dst_rf_ptr_opt.has_value() && src_rf_numeric_opt.has_value()) {
                    m_registers.insert(dst_register, loc, dst_rf_ptr_opt->subtract(*src_rf_numeric_opt, m_slacks));
                }
                else { // when dst is numeric, we keep no information about subtraction
                    m_registers -= dst_register;
                }
                break;
            }
            default: {
                // no other operations supported for packet pointers in the offset domain
                m_registers -= dst_register;
                break;
            }
        }
    }
}

void offset_domain_t::operator()(const Bin& bin, location_t loc) {
    // nothing to do here
}

void offset_domain_t::operator()(const Undefined& u, location_t loc) {
    // nothing to do here
}

void offset_domain_t::operator()(const Un& u, location_t loc) {
    m_registers -= register_t{u.dst.v};
}

void offset_domain_t::operator()(const LoadMapFd& u, location_t loc) {
    m_registers -= register_t{u.dst.v};
}

void offset_domain_t::operator()(const LoadVariable& u, location_t loc) {
    m_registers -= register_t{u.dst.v};
}

void offset_domain_t::do_call(const Call& u, const stack_cells_t& cells, location_t loc) {
    for (const auto& [offset, width] : cells) {
        m_stack -= m_stack.find_overlapping_cells(offset, width);
    }
    m_registers -= register_t{R0_RETURN_VALUE};
    m_registers.scratch_caller_saved_registers();
    if (u.reallocate_packet) {
        m_registers.forget_packet_pointers(loc);
    }
}

void offset_domain_t::operator()(const Call& u, location_t loc) {
    // nothing to do here
}
void offset_domain_t::operator()(const Exit& u, location_t loc) {}

void offset_domain_t::operator()(const Jmp& u, location_t loc) {
    // nothing to do here
}

void offset_domain_t::operator()(const Packet& u, location_t loc) {
    m_registers -= register_t{R0_RETURN_VALUE};
    m_registers.scratch_caller_saved_registers();
}

bool offset_domain_t::check_packet_access(const Reg& r, int width, int offset,
        bool is_comparison_check, location_t loc) {
    auto begin = m_registers.find(register_t{R12_PKT_BEGIN});
    std::string loc_str = loc.to_string();
    auto reg = m_registers.find(r.v);
    if (!reg) {
        m_errors.push_back(loc_str + ": Checking valid access on unknown register");
        return false;
    }
    auto check_lb = *reg + offset;
    auto check_ub = check_lb + width;
    auto [check_lb_satisfied, check_ub_satisfied] = begin->safe_access(check_lb.get_value(),
                                            check_ub.get_value(), is_comparison_check, m_slacks);
    if (!check_lb_satisfied) {
        m_errors.push_back(loc_str + ": Lower bound must be at least meta_offset");
        return false;
    }
    if (!check_ub_satisfied) {
        if (is_comparison_check) {
            m_errors.push_back(loc_str + ": Upper bound must be at most " +
                               std::to_string(MAX_PACKET_SIZE));
        } else {
            m_errors.push_back(loc_str + ": Upper bound must be at most packet size");
        }
        return false;
    }
    return true;
}

void offset_domain_t::check_valid_access(const ValidAccess& s, int w, location_t loc) {
    bool is_comparison_check = s.width == (Value)Imm{0};
    check_packet_access(s.reg, w, s.offset, is_comparison_check, loc);
}

void offset_domain_t::operator()(const Assert &u, location_t loc) {
    // nothing to do here
}

void offset_domain_t::operator()(const basic_block_t& bb) {
    // nothing to do here
}

void offset_domain_t::do_mem_store(const Mem& b,
        std::optional<ptr_or_mapfd_t>& maybe_basereg_type) {
    std::optional<refinement_t> rf_info = std::nullopt;
    if (!is_stack_ptr(maybe_basereg_type)) return;

    int offset = b.access.offset;
    int width = b.access.width;
    auto basereg_with_off = std::get<ptr_with_off_t>(*maybe_basereg_type);
    auto offset_reg = basereg_with_off.get_offset();
    if (auto finite = offset_reg.finite_size()) {
        int finite_size = finite->cast_to<int>();
        const number_t lb = offset_reg.lb().number().value();
        uint64_t lb_n = lb.cast_to<uint64_t>();
        uint64_t store_at = lb_n + offset;
        m_stack -= m_stack.find_overlapping_cells(store_at, width + finite_size);

        if (auto offset_singleton = offset_reg.singleton()) {
            if (auto target_reg = std::get_if<Reg>(&b.value)) {
                rf_info = m_registers.find(target_reg->v);
                if (rf_info) m_stack.store(store_at, *rf_info, width);
            }
        }
    }
}

bool offset_domain_t::do_load(const Mem& b, const register_t& target_register,
        std::optional<ptr_or_mapfd_t> basereg_type, location_t loc) {

    bool is_stack_p = is_stack_ptr(basereg_type);
    bool is_ctx_p = is_ctx_ptr(basereg_type);

    if (!is_stack_p && !is_ctx_p) {
        m_registers -= target_register;
        return false;
    }

    std::string loc_str = loc.to_string();
    int width = b.access.width;
    int offset = b.access.offset;
    auto type_with_off = std::get<ptr_with_off_t>(*basereg_type);
    auto p_offset = type_with_off.get_offset();
    auto offset_singleton = p_offset.singleton();
    if (is_stack_p) {
        if (!offset_singleton) {
            m_registers -= target_register;
        }
        else {
            if (width != 1 && width != 2 && width != 4 && width != 8) {
                m_registers -= target_register;
                return false;
            }
            auto ptr_offset = offset_singleton.value();
            auto load_at = (ptr_offset + offset).cast_to<uint64_t>();

            auto loaded = m_stack.find(load_at);
            if (!loaded) {
                // no field at loaded offset in stack
                m_registers -= target_register;
                return false;
            }
            m_registers.insert(target_register, loc, loaded->first);
        }
    }
    else {
        if (offset_singleton) {
            auto ptr_offset = offset_singleton.value();
            auto load_at = (ptr_offset + offset).cast_to<uint64_t>();

            auto loaded = m_ctx->find(load_at);
            if (loaded && width == 4) {
                m_registers.insert(target_register, loc, *loaded);
                return false;
            }
            m_registers -= target_register;
        }
        // These checks are important to ensure that either we read a complete ptr (loaded before),
        // or a numeric value (when cells contains nothing).
        for (auto const& k : m_ctx->get_keys()) {
            // The value 4 should be dynamic, however, pkt pointers are always stored as 4-byte
            auto start = p_offset.lb() + bound_t{offset};
            auto end = start + bound_t{width-1};
            if (end < bound_t{k} || start > bound_t{k + 4 - 1}) {
                // no overlap with stored range
                continue;
            }
            m_errors.push_back(loc_str + ": Load range in ctx contains pointers");
            return true;
        }
    }
    return false;
}

void offset_domain_t::operator()(const Mem& b, location_t loc) {
    // nothing to do here
}

std::vector<uint64_t> offset_domain_t::get_ctx_keys() const {
    return m_ctx->get_keys();
}

std::optional<refinement_t> offset_domain_t::find_refinement_at_loc(const register_location_t reg) const {
    return m_registers.find(reg);
}

std::optional<refinement_t> offset_domain_t::find_in_ctx(int key) const {
    return m_ctx->find(key);
}

std::optional<refinement_stack_cell_t> offset_domain_t::find_in_stack(int key) const {
    return m_stack.find(key);
}

std::optional<refinement_t> offset_domain_t::find_refinement_info(register_t reg) const {
    return m_registers.find(reg);
}

void offset_domain_t::insert_in_registers(register_t reg, location_t loc, refinement_t rf) {
    m_registers.insert(reg, loc, rf);
}

void offset_domain_t::store_in_stack(uint64_t key, refinement_t d, int width) {
    m_stack.store(key, d, width);
}

void offset_domain_t::adjust_bb_for_types(location_t loc) {
    m_registers.adjust_bb_for_registers(loc);
}

} // namespace crab
