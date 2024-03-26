// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/offset_domain.hpp"

namespace crab {

void registers_state_t::insert(register_t reg, const location_t& loc, refinement_t&& rf) {
    reg_with_loc_t reg_with_loc{reg, loc};
    (*m_offset_env)[reg_with_loc] = std::move(rf);
    m_cur_def[reg] = std::make_shared<reg_with_loc_t>(reg_with_loc);
}

void registers_state_t::insert_slack_value(symbol_t sym, mock_interval_t in) {
    (*m_slacks)[sym] = std::move(in);
}

std::optional<refinement_t> registers_state_t::find(reg_with_loc_t reg) const {
    auto it = m_offset_env->find(reg);
    if (it == m_offset_env->end()) return {};
    return it->second;
}

std::optional<mock_interval_t> registers_state_t::find_slack_value(symbol_t sym) const {
    auto it = m_slacks->find(sym);
    if (it == m_slacks->end()) return {};
    return it->second;
}

std::optional<refinement_t> registers_state_t::find(register_t key) const {
    if (m_cur_def[key] == nullptr) return {};
    return find(*(m_cur_def[key]));
}

std::vector<uint64_t> stack_state_t::find_overlapping_cells(uint64_t start, int width) const {
    std::vector<uint64_t> overlapping_cells;
    auto it = m_slot_rfs.begin();
    while (it != m_slot_rfs.end() && it->first < start) {
        it++;
    }
    if (it != m_slot_rfs.begin()) {
        it--;
        auto key = it->first;
        auto width_key = it->second.second;
        if (key < start && key+width_key > start) overlapping_cells.push_back(key);
    }

    for (; it != m_slot_rfs.end(); it++) {
        auto key = it->first;
        if (key >= start && key < start+width) overlapping_cells.push_back(key);
        if (key >= start+width) break;
    }
    return overlapping_cells;
}

void registers_state_t::set_to_top() {
    m_offset_env = std::make_shared<global_offset_env_t>();
    m_cur_def = live_refinements_t{nullptr};
    m_is_bottom = false;
}

void registers_state_t::set_to_bottom() {
    m_cur_def = live_refinements_t{nullptr};
    m_is_bottom = true;
}

bool registers_state_t::is_top() const {
    if (m_is_bottom) return false;
    if (m_offset_env == nullptr) return true;
    for (auto &it : m_cur_def) {
        if (it != nullptr) return false;
    }
    return true;
}

bool registers_state_t::is_bottom() const {
    return m_is_bottom;
}

void registers_state_t::operator-=(register_t to_forget) {
    if (is_bottom()) {
        return;
    }
    m_cur_def[to_forget] = nullptr;
}

registers_state_t registers_state_t::operator|(const registers_state_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }

    registers_state_t joined_state(m_offset_env, m_slacks);
    location_t loc = location_t(std::make_pair(label_t(-2, -2), 0));

    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (m_cur_def[i] == nullptr || other.m_cur_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_def[i]));
        auto it2 = other.find(*(other.m_cur_def[i]));
        if (it1 && it2) {
            auto rf1 = *it1, rf2 = *it2;
            if (rf1.same_type(rf2)) {
                expression_t rf1_value = rf1.get_value();
                expression_t rf2_value = rf2.get_value();
                if (rf1_value == rf2_value) {
                    joined_state.insert(register_t{i}, loc, std::move(rf1));
                }
                else {
                    //joined_state.insert(register_t{i}, loc, std::move(rf1 | rf2));
                    symbol_t s = symbol_t::make();
                    symbol_terms_t terms;
                    bool added = false;
                    auto symbol_terms1 = rf1_value.get_symbol_terms();
                    auto symbol_terms2 = rf2_value.get_symbol_terms();
                    auto interval_rf1_value = rf1_value.get_interval();
                    auto interval_rf2_value = rf2_value.get_interval();
                    for (auto it = symbol_terms1.begin(); it != symbol_terms1.end(); it++) {
                        // assuming slack variables are at same position in both states
                        if (!added && it->first.is_slack()) {
                            auto slack1 = it->first;
                            auto dist = std::distance(symbol_terms1.begin(), it);
                            auto it2 = std::next(symbol_terms2.begin(), dist);
                            auto slack2 = it2->first;
                            auto slack1_value = find_slack_value(slack1);
                            auto slack2_value = other.find_slack_value(slack2);
                            if (slack1_value && slack2_value) {
                                auto interval1 = slack1_value->to_interval() + interval_rf1_value;
                                auto interval2 = slack2_value->to_interval() + interval_rf2_value;
                                auto interval = interval1 | interval2;
                                terms[s] = 1;
                                (*m_slacks)[s] = mock_interval_t{interval};
                                added = true;
                                continue;
                            }
                        }
                        terms[it->first] = it->second;
                    }
                    auto rf = refinement_t(rf1.get_type(), expression_t(terms));
                    joined_state.insert(register_t{i}, loc, std::move(rf));
                }
            }
        }
    }
    // need special handling for the registers v_begin, v_end, and v_meta
    for (uint8_t i = NUM_REGISTERS; i < NUM_REGISTERS+3; i++) {
        if (m_cur_def[i] == nullptr || other.m_cur_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_def[i]));
        auto it2 = other.find(*(other.m_cur_def[i]));
        if (it1 && it2) {
            auto rf1 = *it1, rf2 = *it2;
            auto rf_joined = rf1 | rf2;
            joined_state.insert(register_t{i}, loc, std::move(rf_joined));
        }
    }
    return joined_state;
}

void registers_state_t::adjust_bb_for_registers(location_t loc) {
    for (uint8_t i = 0; i < NUM_REGISTERS+3; i++) {
        if (auto it = find(register_t{i})) {
            insert(register_t{i}, loc, std::move(*it));
        }
    }
}

void registers_state_t::scratch_caller_saved_registers() {
    for (uint8_t r = R1_ARG; r <= R5_ARG; r++) {
        operator-=(register_t{r});
    }
}

void registers_state_t::forget_packet_pointers(location_t loc) {
    for (uint8_t r = R0_RETURN_VALUE; r < NUM_REGISTERS; r++) {
        if (auto it = find(register_t{r})) {
            if (it->get_type() == data_type_t::PACKET) {
                operator-=(register_t{r});
            }
        }
    }
    insert(register_t{11}, loc, refinement_t::begin());
    // TODO: verify if this is all needed
}

void stack_state_t::set_to_top() {
    m_slot_rfs.clear();
    m_is_bottom = false;
}

void stack_state_t::set_to_bottom() {
    m_slot_rfs.clear();
    m_is_bottom = true;
}

bool stack_state_t::is_top() const {
    if (m_is_bottom) return false;
    return m_slot_rfs.empty();
}

bool stack_state_t::is_bottom() const {
    return m_is_bottom;
}

stack_state_t stack_state_t::top() {
    return stack_state_t(false);
}

std::optional<refinement_cells_t> stack_state_t::find(uint64_t key) const {
    auto it = m_slot_rfs.find(key);
    if (it == m_slot_rfs.end()) return {};
    return it->second;
}

void stack_state_t::store(uint64_t key, refinement_t d, int width) {
    m_slot_rfs[key] = std::make_pair(d, width);
}

void stack_state_t::operator-=(uint64_t to_erase) {
    if (is_bottom()) {
        return;
    }
    m_slot_rfs.erase(to_erase);
}

void stack_state_t::operator-=(const std::vector<uint64_t>& keys) {
    for (auto &key : keys) {
       *this -= key;
    }
}

stack_state_t stack_state_t::operator|(const stack_state_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }

    stack_slot_refinements_t out_stack_rfs;
    // We do not join rf cells because different rf values different types of offsets
    for (auto const&kv: m_slot_rfs) {
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
    return stack_state_t(std::move(out_stack_rfs));
}

ctx_offsets_t::ctx_offsets_t(const ebpf_context_descriptor_t* desc) {
    if (desc->data >= 0) {
        m_rfs[desc->data] = refinement_t::begin();
    }
    if (desc->end >= 0) {
        m_rfs[desc->end] = refinement_t::end();
    }
    if (desc->meta >= 0) {
        m_rfs[desc->meta] = refinement_t::meta();
    }
    if (desc->size >= 0) {
        m_size = desc->size;
    }
}

int ctx_offsets_t::get_size() const {
    return m_size;
}

std::optional<refinement_t> ctx_offsets_t::find(int key) const {
    auto it = m_rfs.find(key);
    if (it == m_rfs.end()) return {};
    return it->second;
}

offset_domain_t&& offset_domain_t::setup_entry() {
    std::shared_ptr<ctx_offsets_t> ctx
        = std::make_shared<ctx_offsets_t>(global_program_info->type.context_descriptor);
    registers_state_t regs(std::make_shared<global_offset_env_t>(),
            std::make_shared<slacks_t>(), global_program_info->type.context_descriptor);

    static offset_domain_t off_d(std::move(regs), stack_state_t::top(), ctx);
    return std::move(off_d);
}

offset_domain_t offset_domain_t::bottom() {
    offset_domain_t off;
    off.set_to_bottom();
    return off;
}

void offset_domain_t::set_to_top() {
    m_reg_state.set_to_top();
    m_stack_state.set_to_top();
    m_is_bottom = false;
}

void offset_domain_t::set_to_bottom() {
    m_is_bottom = true;
    m_reg_state.set_to_bottom();
    m_stack_state.set_to_bottom();
}

bool offset_domain_t::is_bottom() const {
    if (m_is_bottom) return true;
    return (m_reg_state.is_bottom() || m_stack_state.is_bottom());
}

bool offset_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_reg_state.is_top() && m_stack_state.is_top());
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
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return offset_domain_t(
            m_reg_state | other.m_reg_state,
            m_stack_state | other.m_stack_state,
            m_ctx_rfs
    );
}

offset_domain_t offset_domain_t::operator|(offset_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return offset_domain_t(
            m_reg_state | std::move(other.m_reg_state),
            m_stack_state | std::move(other.m_stack_state),
            m_ctx_rfs
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

crab::bound_t offset_domain_t::get_loop_count_upper_bound() {
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
        auto rf_left = m_reg_state.find(cond.left.v);
        auto rf_right = m_reg_state.find(right_reg);
        if (!rf_left || !rf_right) {
            // this should not happen, comparison between a packet pointer and either
            // other region's pointers or numbers; possibly raise type error
            m_errors.push_back("one of the pointers being compared isn't packet pointer");
            return;
        }
        if (cond.op == Condition::Op::LE) {
            auto b = m_reg_state.find(register_t{11});
            auto le_rf = *rf_left <= *rf_right;
            if (b) {
                b->add_constraint(le_rf);
                if (rf_right->has_value(refinement_t::end())) {
                    b->add_constraint(refinement_t::meta() <= refinement_t::begin());
                }
                else if (rf_right->has_value(refinement_t::begin())) {
                    b->add_constraint(refinement_t::begin() <= refinement_t::end());
                }
                b->simplify();
                auto b_copy = *b;
                if (b_copy.is_bottom(m_reg_state.get_slacks())) {
                    set_to_bottom();
                }
                else {
                    m_reg_state.insert(register_t{11}, loc, std::move(*b));
                }
            }
        }
        else if (cond.op == Condition::Op::GT) {
            auto b = m_reg_state.find(register_t{11});
            auto gt_rf = *rf_left > *rf_right;
            if (b) {
                b->add_constraint(gt_rf);
                b->simplify();
                auto b_copy = *b;
                if (b_copy.is_bottom(m_reg_state.get_slacks())) {
                    set_to_bottom();
                }
                else {
                    m_reg_state.insert(register_t{11}, loc, std::move(*b));
                }
            }

        }
        // other comparisons not supported
    }
}

static void create_numeric_refinement(registers_state_t& reg_state, mock_interval_t&& interval,
        location_t loc, register_t reg) {
    symbol_t s = symbol_t::make();
    refinement_t rf = refinement_t(data_type_t::NUM, expression_t(s));
    reg_state.insert(reg, loc, std::move(rf));
    reg_state.insert_slack_value(s, std::move(interval));
}

void offset_domain_t::do_bin(const Bin& bin,
        const std::optional<interval_t>& src_signed_interval_opt,
        const std::optional<ptr_or_mapfd_t>& src_ptr_or_mapfd_opt,
        const std::optional<interval_t>& dst_signed_interval_opt,
        const std::optional<ptr_or_mapfd_t>& dst_ptr_or_mapfd_opt,
        mock_interval_t &&interval_result, location_t loc) {

    using Op = Bin::Op;

    auto dst_register = register_t{bin.dst.v};

    //if (!is_packet_ptr(src_ptr_or_mapfd_opt) && !is_packet_ptr(dst_ptr_or_mapfd_opt) &&
    //        (!src_signed_interval_opt && bin.op != Op::MOV) && !dst_signed_interval_opt) {
    //    return interval_t::bottom();
    //}

    if (std::holds_alternative<Imm>(bin.v)) {
        int64_t imm;
        if (bin.is64) {
            // Use the full signed value.
            imm = static_cast<int64_t>(std::get<Imm>(bin.v).v);
        } else {
            // Use only the low 32 bits of the value.
            imm = static_cast<int>(std::get<Imm>(bin.v).v);
        }
        auto imm_interval = interval_t{number_t{imm}};
        switch (bin.op) {
            case Op::MOV: {
                // ra = imm, we forget the type in the offset domain
                create_numeric_refinement(m_reg_state, std::move(interval_result),
                        loc, dst_register);
                break;
            }
            case Op::ADD: {
                // ra += imm
                if (imm == 0) break;
                if (auto dst_rf_opt = m_reg_state.find(dst_register)) {
                    auto rf = *dst_rf_opt + imm;
                    m_reg_state.insert(dst_register, loc, std::move(rf));
                }
                else {
                    m_reg_state -= dst_register;
                }
                break;
            }
            case Op::SUB: {
                // ra -= imm
                if (imm == 0) break;
                if (auto dst_rf_opt = m_reg_state.find(dst_register)) {
                    auto rf = *dst_rf_opt + (-imm);
                    m_reg_state.insert(dst_register, loc, std::move(rf));
                }
                else {
                    m_reg_state -= dst_register;
                }
                break;
            }
            default: {
                if (dst_signed_interval_opt) {
                    symbol_t s = symbol_t::make();
                    refinement_t rf = refinement_t(data_type_t::NUM, expression_t(s));
                    m_reg_state.insert(dst_register, loc, std::move(rf));
                    m_reg_state.insert_slack_value(s, interval_result);
                }
                else {
                    // no other operations supported for packet pointers in the offset domain
                    m_reg_state -= dst_register;
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
                if (auto src_rf_opt = m_reg_state.find(src.v)) {
                    auto rf = *src_rf_opt;
                    m_reg_state.insert(dst_register, loc, std::move(rf));
                }
                else {
                    m_reg_state -= dst_register;
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
                    if (auto src_rf_opt = m_reg_state.find(src.v)) {
                        if (auto dst_rf_opt = m_reg_state.find(dst_register)) {
                            auto rf = *dst_rf_opt + *src_rf_opt;
                            m_reg_state.insert(dst_register, loc, std::move(rf));
                        }
                        else {
                            m_reg_state -= dst_register;
                        }
                    }
                    else {
                        m_reg_state -= dst_register;
                    }
                }
                break;
            }
            case Op::SUB: {
                // ra -= rb
                // TODO: be precise with ptr -= ptr, and possibly assign a slack to the result
                if (is_packet_ptr(src_ptr_or_mapfd_opt) && is_packet_ptr(dst_ptr_or_mapfd_opt)) {
                    create_numeric_refinement(m_reg_state, std::move(interval_result),
                            loc, dst_register);
                    return;
                }
                else {
                    if (auto src_rf_opt = m_reg_state.find(src.v)) {
                        if (auto dst_rf_opt = m_reg_state.find(dst_register)) {
                            auto rf = *dst_rf_opt - *src_rf_opt;
                            m_reg_state.insert(dst_register, loc, std::move(rf));
                        }
                        else {
                            m_reg_state -= dst_register;
                        }
                    }
                    else {
                        m_reg_state -= dst_register;
                    }
                }
                break;
            }
            default: {
                if (dst_ptr_or_mapfd_opt || src_ptr_or_mapfd_opt) {
                    // no other operations supported for packet pointers in the offset domain
                    m_reg_state -= dst_register;
                }
                else {
                    create_numeric_refinement(m_reg_state, std::move(interval_result),
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
        m_reg_state -= u.dst.v;
    }
    else {
        create_numeric_refinement(m_reg_state, std::move(interval), loc, register_t{u.dst.v});
    }
}

void offset_domain_t::operator()(const LoadMapFd& u, location_t loc) {
    m_reg_state -= u.dst.v;
}

void offset_domain_t::do_call(const Call& u, const stack_cells_t& cells, location_t loc) {
    for (const auto& kv : cells) {
        auto rf = kv.first;
        auto width = kv.second;
        auto overlapping_cells = m_stack_state.find_overlapping_cells(rf, width);
        m_stack_state -= overlapping_cells;
    }
    m_reg_state.scratch_caller_saved_registers();
    register_t r0{R0_RETURN_VALUE};
    if (u.reallocate_packet) {
        m_reg_state -= r0;
        m_reg_state.forget_packet_pointers(loc);
    }
    else if (u.is_map_lookup) {
        m_reg_state -= r0;
    }
    else {
        // slack needs to be fixed, as it can have any value
        create_numeric_refinement(m_reg_state, mock_interval_t::top(), loc, r0);
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
    create_numeric_refinement(m_reg_state, mock_interval_t::top(), loc,
            register_t{R0_RETURN_VALUE});
    m_reg_state.scratch_caller_saved_registers();
}

bool offset_domain_t::lower_bound_satisfied(const dist_t& dist, int offset) const {
    auto meta_limit = m_extra_constraints.get_meta_limit();
    auto end_limit = m_extra_constraints.get_end_limit();

    dist_t dist1 = dist;
    if (dist.is_meta_pointer()) {
        dist1 = dist_t(dist.offset_from_reference() + (meta_limit ?
                    weight_t{*meta_limit-number_t{PACKET_META}} : weight_t{number_t{0}}));
    }
    if (dist.is_backward_pointer()) {
        dist1 = dist_t(dist.offset_from_reference()
                + (end_limit ? weight_t{*end_limit} : weight_t{number_t{0}}));
    }

    bound_t lb = meta_limit ? *meta_limit-number_t{PACKET_META} : bound_t{number_t{0}};
    return (dist1.m_dist.lb()+number_t{offset} >= lb);
}

bool offset_domain_t::upper_bound_satisfied(const dist_t& dist, int offset, int width,
        bool is_comparison_check) const {
    auto meta_limit = m_extra_constraints.get_meta_limit();
    auto end_limit = m_extra_constraints.get_end_limit();

    dist_t dist1 = dist;
    if (dist.is_meta_pointer()) {
        dist1 = dist_t(dist.offset_from_reference() + (meta_limit ?
                    weight_t{*meta_limit-number_t{PACKET_META}} : weight_t{number_t{0}}));
    }
    if (dist.is_backward_pointer()) {
        dist1 = dist_t(dist.offset_from_reference()
                + (end_limit ? weight_t{*end_limit} :
                    weight_t{number_t{is_comparison_check ? MAX_PACKET_SIZE : 0}}));
    }

    bound_t ub = is_comparison_check ? bound_t{MAX_PACKET_SIZE}
        : (end_limit ? *end_limit : number_t{0});
    return (dist1.m_dist.ub()+number_t{offset+width} <= ub);
}

bool offset_domain_t::check_packet_access(const Reg& r, int width, int offset,
        bool is_comparison_check) const {
    auto begin = m_reg_state.find(register_t{11});
    if (!begin) return false;
    auto reg = m_reg_state.find(r.v);
    if (!reg) return false;
    auto toCheck = *reg + (offset+width);
    return toCheck.is_safe_with(*begin, m_reg_state.get_slacks(), is_comparison_check);
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

void offset_domain_t::operator()(const basic_block_t& bb, int print) {
    // nothing to do here
}

void offset_domain_t::do_mem_store(const Mem& b,
        std::optional<ptr_or_mapfd_t>& maybe_basereg_type) {
    auto target_reg = std::get<Reg>(b.value);
    auto rf_info = m_reg_state.find(target_reg.v);
    if (!rf_info) return;

    int offset = b.access.offset;
    int width = b.access.width;
    if (is_stack_ptr(maybe_basereg_type)) {
        auto basereg_with_off = std::get<ptr_with_off_t>(*maybe_basereg_type);
        auto basereg_off_singleton = basereg_with_off.get_offset().to_interval().singleton();
        if (!basereg_off_singleton) return;
        auto store_at = (uint64_t)(*basereg_off_singleton + offset);
        auto overlapping_cells = m_stack_state.find_overlapping_cells(store_at, width);
        m_stack_state -= overlapping_cells;

        m_stack_state.store(store_at, *rf_info, width);
    }
}

void offset_domain_t::do_load(const Mem& b, const register_t& target_register,
        std::optional<ptr_or_mapfd_t> basereg_type, interval_t &&interval_result, location_t loc) {

    bool is_stack_p = is_stack_ptr(basereg_type);
    bool is_ctx_p = is_ctx_ptr(basereg_type);
    bool is_packet_p = is_packet_ptr(basereg_type);
    bool is_shared_p = is_shared_ptr(basereg_type);

    if (interval_result != interval_t::bottom()) {
        if (is_ctx_p || is_shared_p || is_packet_p) {
            create_numeric_refinement(m_reg_state, std::move(interval_result), loc,
                    target_register);
            return;
        }
    }

    if (!is_stack_p && !is_ctx_p) {
        m_reg_state -= target_register;
        return;
    }

    int offset = b.access.offset;
    auto type_with_off = std::get<ptr_with_off_t>(*basereg_type);
    auto p_offset = type_with_off.get_offset();
    auto offset_singleton = p_offset.to_interval().singleton();
    if (!offset_singleton) {
        m_reg_state -= target_register;
        return;
    }
    auto ptr_offset = *offset_singleton;
    auto load_at = (uint64_t)(ptr_offset + offset);

    if (is_stack_p) {
        auto it = m_stack_state.find(load_at);
        if (!it) {
            m_reg_state -= target_register;
            return;
        }
        m_reg_state.insert(target_register, loc, std::move(it->first));
    }
    else if (is_ctx_p) {
        auto it = m_ctx_rfs->find(load_at);
        if (!it) {
            m_reg_state -= target_register;
            return;
        }
        m_reg_state.insert(target_register, loc, std::move(*it));
    }
}

void offset_domain_t::operator()(const Mem& b, location_t loc) {
    // nothing to do here
}

std::optional<refinement_t> offset_domain_t::find_refinement_at_loc(const reg_with_loc_t reg) const {
    return m_reg_state.find(reg);
}

std::optional<refinement_t> offset_domain_t::find_in_ctx(int key) const {
    return m_ctx_rfs->find(key);
}

std::optional<refinement_cells_t> offset_domain_t::find_in_stack(int key) const {
    return m_stack_state.find(key);
}

std::optional<refinement_t> offset_domain_t::find_refinement_info(register_t reg) const {
    return m_reg_state.find(reg);
}

void offset_domain_t::insert_in_registers(register_t reg, location_t loc, refinement_t rf) {
    m_reg_state.insert(reg, loc, std::move(rf));
}

void offset_domain_t::store_in_stack(uint64_t key, refinement_t d, int width) {
    m_stack_state.store(key, d, width);
}

void offset_domain_t::adjust_bb_for_types(location_t loc) {
    m_reg_state.adjust_bb_for_registers(loc);
}

} // namespace crab
