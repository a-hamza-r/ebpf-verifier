// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "offset_domain.hpp"

namespace crab {

void offset_registers_t::insert(register_t reg, const location_t& loc, refinement_t&& rf) {
    register_location_t register_location{reg, loc};
    (*m_registers_env)[register_location] = std::move(rf);
    m_cur_register_def[reg] = std::make_shared<register_location_t>(register_location);
}

void offset_registers_t::insert_slack_value(symbol_t sym, mock_interval_t in) {
    (*m_slacks)[sym] = std::move(in);
}

std::optional<refinement_t> offset_registers_t::find(register_location_t reg) const {
    auto it = m_registers_env->find(reg);
    if (it == m_registers_env->end()) return {};
    return it->second;
}

std::optional<mock_interval_t> offset_registers_t::find_slack_value(symbol_t sym) const {
    auto it = m_slacks->find(sym);
    if (it == m_slacks->end()) return {};
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

offset_registers_t offset_registers_t::operator|(const offset_registers_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }

    offset_registers_t joined_state(m_registers_env, m_slacks);
    location_t loc = location_t::top();

    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (m_cur_register_def[i] == nullptr || other.m_cur_register_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_register_def[i]));
        auto it2 = other.find(*(other.m_cur_register_def[i]));
        if (it1 && it2) {
            auto rf1 = *it1, rf2 = *it2;
            if (rf1.same_type(rf2)) {
                joined_state.insert(register_t{i}, loc, rf1 | rf2);
            }
        }
    }
    return joined_state;
}

void offset_registers_t::adjust_bb_for_registers(location_t loc) {
    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (auto it = find(register_t{i})) {
            insert(register_t{i}, loc, std::move(*it));
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
            if (it->get_type() == refinement_type_t::PACKET) {
                operator-=(register_t{r});
            }
        }
    }
    insert(BEGIN_REG, loc, refinement_t::begin(true));
    // TODO: verify if this is all needed
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

offset_stack_t offset_stack_t::operator|(const offset_stack_t& other) const {
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
            // TODO: for numerical values, the width does not have to be the same
            // hence, handle accordingly
            if (rf1.same_type(rf2) && width1 == width2) {
                out_stack_rfs.insert({kv.first, std::make_pair(rf1 | rf2, width1)});
            }
        }
    }
    return offset_stack_t(std::move(out_stack_rfs));
}

offset_ctx_t::offset_ctx_t(const ebpf_context_descriptor_t* desc) {
    if (desc->data >= 0) {
        m_ctx_cells[desc->data] = refinement_t::begin();
    }
    if (desc->end >= 0) {
        m_ctx_cells[desc->end] = refinement_t::end();
    }
    if (desc->meta >= 0) {
        m_ctx_cells[desc->meta] = refinement_t::meta();
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

offset_domain_t offset_domain_t::setup_entry() {
    offset_registers_t regs(std::make_shared<global_env_offset_registers_t>(),
                            std::make_shared<slacks_t>(),
                            global_program_info->type.context_descriptor->data);

    return offset_domain_t{std::move(regs), offset_stack_t::top(),
                    std::make_shared<offset_ctx_t>(global_program_info->type.context_descriptor)};
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
bool offset_domain_t::operator<=(const offset_domain_t& other) const { return true; }

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
    return offset_domain_t(m_registers | other.m_registers, m_stack | other.m_stack, m_ctx);
}

offset_domain_t offset_domain_t::operator|(offset_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return offset_domain_t(
        m_registers | std::move(other.m_registers),
        m_stack | std::move(other.m_stack),
        std::move(m_ctx)
    );
}

// meet
offset_domain_t offset_domain_t::operator&(const offset_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

// widening
offset_domain_t offset_domain_t::widen(const offset_domain_t& other, bool to_constants) {
    /* WARNING: The operation is not implemented yet.*/
    return other;
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
            m_errors.push_back("one of the pointers being compared isn't packet pointer");
            return;
        }
        auto begin = m_registers.find(register_t{BEGIN_REG});
        if (cond.op == Condition::Op::LE) {
            begin->add_constraint(*rf_left <= *rf_right);
            m_registers.insert(register_t{BEGIN_REG}, loc, std::move(*begin));
        }
        else if (cond.op == Condition::Op::GT) {
            begin->add_constraint(*rf_left > *rf_right);
            m_registers.insert(register_t{BEGIN_REG}, loc, std::move(*begin));
        }
        // other comparisons not supported
    }
}

interval_t offset_domain_t::compute_packet_subtraction(register_t dst, register_t src) const {
    auto dst_rf = m_registers.find(dst);
    auto src_rf = m_registers.find(src);
    if (!dst_rf || !src_rf) return interval_t::bottom();
    expression_t dst_expr = dst_rf->get_value();
    expression_t src_expr = src_rf->get_value();
    refinement_t result_rf = *dst_rf - *src_rf;
    expression_t result_expr = result_rf.get_value().get_equivalent_expression();
    if (result_expr.is_constant()) {
        return result_expr.get_constant_term();
    }
    // with non-singleton expressions, we might be able to compute subtraction,
    // but it might be complicated
    if (!dst_expr.is_singleton() || !src_expr.is_singleton()) return interval_t::top();
    auto dst_symbol = dst_expr.get_singleton();
    auto src_symbol = src_expr.get_singleton();
    std::optional<refinement_t> begin_rf = m_registers.find(register_t{BEGIN_REG});
    return begin_rf->simplify_for_subtraction(dst_symbol, src_symbol);
}

static void create_numeric_refinement(offset_registers_t& reg_state, mock_interval_t&& interval,
        location_t loc, register_t reg) {
    symbol_t s = symbol_t::make();
    reg_state.insert_slack_value(s, std::move(interval));
    expression_t value = expression_t(s, reg_state.get_slacks());
    reg_state.insert(reg, loc, refinement_t::numeric_refinement(std::move(value)));
}

void offset_domain_t::do_bin(const Bin& bin,
        const std::optional<interval_t>& src_signed_interval_opt,
        const std::optional<ptr_or_mapfd_t>& src_ptr_or_mapfd_opt,
        const std::optional<interval_t>& dst_signed_interval_opt,
        const std::optional<ptr_or_mapfd_t>& dst_ptr_or_mapfd_opt,
        mock_interval_t &&interval_result, location_t loc) {

    using Op = Bin::Op;

    auto dst_register = register_t{bin.dst.v};

    if (std::holds_alternative<Imm>(bin.v)) {
        int64_t imm;
        if (bin.is64) {
            // Use the full signed value.
            imm = static_cast<int64_t>(std::get<Imm>(bin.v).v);
        } else {
            // Use only the low 32 bits of the value.
            imm = static_cast<int>(std::get<Imm>(bin.v).v);
        }
        auto imm_interval = interval_t{imm};
        switch (bin.op) {
            case Op::MOV: {
                // ra = imm
                // we just get the value of the immediate as an interval from interval domain
                create_numeric_refinement(m_registers, std::move(interval_result),
                        loc, dst_register);
                break;
            }
            case Op::ADD: {
                // ra += imm
                if (imm == 0) break;
                if (auto dst_rf_opt = m_registers.find(dst_register)) {
                    auto rf = *dst_rf_opt + imm_interval;
                    m_registers.insert(dst_register, loc, std::move(rf));
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::SUB: {
                // ra -= imm
                if (imm == 0) break;
                if (auto dst_rf_opt = m_registers.find(dst_register)) {
                    auto rf = *dst_rf_opt + (-imm_interval);
                    m_registers.insert(dst_register, loc, std::move(rf));
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            default: {
                if (dst_signed_interval_opt) {
                    create_numeric_refinement(m_registers, std::move(interval_result),
                            loc, dst_register);
                }
                else {
                    // no other operations supported for packet pointers in the offset domain
                    m_registers -= dst_register;
                }
                break;
            }
        }
    }
    else {
        auto src = std::get<Reg>(bin.v);
        switch (bin.op) {
            case Op::MOV: {
                // ra = rb
                if (auto src_rf_opt = m_registers.find(src.v)) {
                    m_registers.insert(dst_register, loc, std::move(*src_rf_opt));
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::ADD: {
                // ra += rb
                if (is_packet_ptr(src_ptr_or_mapfd_opt) && is_packet_ptr(dst_ptr_or_mapfd_opt)) {
                    // possibly adding two pointers
                    set_to_bottom();
                }
                else {
                    if (auto src_rf_opt = m_registers.find(src.v)) {
                        if (auto dst_rf_opt = m_registers.find(dst_register)) {
                            auto rf = *dst_rf_opt + *src_rf_opt;
                            m_registers.insert(dst_register, loc, std::move(rf));
                        }
                        else {
                            m_registers -= dst_register;
                        }
                    }
                    else {
                        m_registers -= dst_register;
                    }
                }
                break;
            }
            case Op::SUB: {
                // ra -= rb
                if (is_packet_ptr(src_ptr_or_mapfd_opt) && is_packet_ptr(dst_ptr_or_mapfd_opt)) {
                    create_numeric_refinement(m_registers, std::move(interval_result),
                            loc, dst_register);
                    return;
                }
                else {
                    if (auto src_rf_opt = m_registers.find(src.v)) {
                        if (auto dst_rf_opt = m_registers.find(dst_register)) {
                            auto rf = *dst_rf_opt - *src_rf_opt;
                            m_registers.insert(dst_register, loc, std::move(rf));
                        }
                        else {
                            m_registers -= dst_register;
                        }
                    }
                    else {
                        m_registers -= dst_register;
                    }
                }
                break;
            }
            default: {
                if (dst_ptr_or_mapfd_opt || src_ptr_or_mapfd_opt) {
                    // no other operations supported for packet pointers in the offset domain
                    m_registers -= dst_register;
                }
                else {
                    create_numeric_refinement(m_registers, std::move(interval_result),
                            loc, dst_register);
                }
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
    // nothing to do here
}

void offset_domain_t::do_un(const Un& u, interval_t interval, location_t loc) {
    if (interval == interval_t::bottom()) {
        m_registers -= u.dst.v;
    }
    else {
        create_numeric_refinement(m_registers, std::move(interval), loc, register_t{u.dst.v});
    }
}

void offset_domain_t::operator()(const LoadMapFd& u, location_t loc) {
    m_registers -= u.dst.v;
}

void offset_domain_t::do_call(const Call& u, const stack_cells_t& cells, location_t loc) {
    for (const auto& kv : cells) {
        auto rf = kv.first;
        auto width = kv.second;
        auto overlapping_cells = m_stack.find_overlapping_cells(rf, width);
        m_stack -= overlapping_cells;
    }
    m_registers.scratch_caller_saved_registers();
    register_t r0{R0_RETURN_VALUE};
    if (u.reallocate_packet) {
        m_registers -= r0;
        m_registers.forget_packet_pointers(loc);
    }
    else if (u.is_map_lookup) {
        m_registers -= r0;
    }
    else {
        // slack needs to be fixed, as it can have any value
        create_numeric_refinement(m_registers, mock_interval_t::top(), loc, r0);
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
    create_numeric_refinement(m_registers, mock_interval_t::top(), loc,
            register_t{R0_RETURN_VALUE});
    m_registers.scratch_caller_saved_registers();
}

bool offset_domain_t::check_packet_access(const Reg& r, int width, int offset,
        bool is_comparison_check) const {
    auto begin = m_registers.find(register_t{BEGIN_REG});
    auto reg = m_registers.find(r.v);
    if (!reg) return false;
    auto toCheck_lb = (*reg + offset).get_value();
    auto toCheck_ub = (*reg + offset + width).get_value();
    return begin->safe_access(toCheck_lb, toCheck_ub, is_comparison_check);
}

void offset_domain_t::check_valid_access(const ValidAccess& s,
        std::optional<ptr_or_mapfd_t>& reg_type, int w) {
    if (w == 0 || !reg_type) return;

    bool is_comparison_check = s.width == (Value)Imm{0};
    if (check_packet_access(s.reg, w, s.offset, is_comparison_check)) return;
    m_errors.push_back("valid access check failed");
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

    if (std::holds_alternative<Reg>(b.value)) {
        auto target_reg = std::get<Reg>(b.value);
        rf_info = m_registers.find(target_reg.v);
    }
    else {
        symbol_t s = symbol_t::make();
        expression_t value = expression_t(s, m_registers.get_slacks());
        rf_info = refinement_t::numeric_refinement(std::move(value));
        interval_t interval = interval_t{number_t{static_cast<uint64_t>(std::get<Imm>(b.value).v)}};
        m_registers.insert_slack_value(s, std::move(interval));
    }
    if (!rf_info) return;

    int offset = b.access.offset;
    int width = b.access.width;
    auto basereg_with_off = std::get<ptr_with_off_t>(*maybe_basereg_type);
    auto basereg_off_singleton = basereg_with_off.get_offset().to_interval().singleton();
    if (!basereg_off_singleton) return;
    auto store_at = (*basereg_off_singleton + offset).cast_to<uint64_t>();
    auto overlapping_cells = m_stack.find_overlapping_cells(store_at, width);
    m_stack -= overlapping_cells;
    m_stack.store(store_at, *rf_info, width);
}

void offset_domain_t::do_load(const Mem& b, const register_t& target_register,
        std::optional<ptr_or_mapfd_t> basereg_type, interval_t &&interval_result, location_t loc) {

    bool is_stack_p = is_stack_ptr(basereg_type);
    bool is_ctx_p = is_ctx_ptr(basereg_type);
    bool is_packet_p = is_packet_ptr(basereg_type);
    bool is_shared_p = is_shared_ptr(basereg_type);

    if (interval_result != interval_t::bottom()) {
        if (is_ctx_p || is_shared_p || is_packet_p) {
            create_numeric_refinement(m_registers, std::move(interval_result), loc,
                    target_register);
            return;
        }
    }

    if (!is_stack_p && !is_ctx_p) {
        m_registers -= target_register;
        return;
    }

    int width = b.access.width;
    int offset = b.access.offset;
    auto type_with_off = std::get<ptr_with_off_t>(*basereg_type);
    auto p_offset = type_with_off.get_offset();
    auto offset_singleton = p_offset.to_interval().singleton();
    if (is_stack_p) {
        if (!offset_singleton) {
            for (auto const& k : m_stack.get_keys()) {
                auto start = p_offset.lb();
                auto end = p_offset.ub()+number_t{offset+width-1};
                interval_t range{start, end};
                // TODO: fix this
                /*
                if (range[number_t{(int)k}]) {
                    //std::cout << "stack load at unknown offset, and offset range contains pointers\n";
                    m_errors.push_back("stack load at unknown offset, and offset range contains pointers");
                    break;
                }
                */
            }
            m_registers -= target_register;
        }
        else {
            if (width != 1 && width != 2 && width != 4 && width != 8) {
                m_registers -= target_register;
                return;
            }
            auto ptr_offset = offset_singleton.value();
            auto load_at = (ptr_offset + offset).cast_to<uint64_t>();

            auto loaded = m_stack.find(load_at);
            if (!loaded) {
                // no field at loaded offset in stack
                m_registers -= target_register;
                return;
            }
            m_registers.insert(target_register, loc, std::move(loaded->first));
        }
    }
    else {
        if (!offset_singleton) {
            for (auto const& k : m_ctx->get_keys()) {
                auto start = p_offset.lb();
                auto end = p_offset.ub()+crab::bound_t{offset+width-1};
                interval_t range{start, end};
                // TODO: fix this
                /*
                if (range[number_t{(int)k}]) {
                    //std::cout << "ctx load at unknown offset, and offset range contains pointers\n";
                    m_errors.push_back("ctx load at unknown offset, and offset range contains pointers");
                    break;
                }
                */
            }
            m_registers -= target_register;
        }
        else {
            auto ptr_offset = offset_singleton.value();
            auto load_at = (ptr_offset + offset).cast_to<uint64_t>();

            auto loaded = m_ctx->find(load_at);
            if (!loaded) {
                // no field at loaded offset in ctx
                m_registers -= target_register;
                return;
            }
            m_registers.insert(target_register, loc, std::move(*loaded));
        }
    }
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
    m_registers.insert(reg, loc, std::move(rf));
}

void offset_domain_t::store_in_stack(uint64_t key, refinement_t d, int width) {
    m_stack.store(key, d, width);
}

void offset_domain_t::adjust_bb_for_types(location_t loc) {
    m_registers.adjust_bb_for_registers(loc);
}

} // namespace crab
