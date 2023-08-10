// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/offset_domain.hpp"

namespace crab {

bool dist_t::operator==(const dist_t& d) const {
    return (m_dist == d.m_dist && m_slack == d.m_slack);
}

weight_t dist_t::offset_from_reference() const {
    if (is_meta_pointer()) {
        return (-m_dist+interval_t{number_t{PACKET_META}});
    }
    if (is_backward_pointer()) {
        return (m_dist-interval_t{number_t{PACKET_END}});
    }
    return m_dist;
}

void dist_t::write(std::ostream& o) const {
    if (m_slack != -1)
        o << "s" << m_slack << "+";
    if (is_forward_pointer()) o << "begin+";
    else if (is_meta_pointer()) o << "meta+";
    else if (is_backward_pointer()) o << "end+";
    auto offset = offset_from_reference();
    auto singleton_val = offset.singleton();
    if (singleton_val) o << singleton_val.value();
    else o << offset;
}

bool dist_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_slack == -1 && m_dist.is_top());
}

bool dist_t::is_bottom() const {
    return m_is_bottom;
}

void dist_t::set_to_top() {
    m_slack = -1;
    m_dist = interval_t::top();
    m_is_bottom = false;
}

void dist_t::set_to_bottom() {
    m_is_bottom = true;
}

bool dist_t::is_meta_pointer() const {
    return (m_dist.lb() > number_t{PACKET_END} && m_dist.ub() <= number_t{PACKET_META});
}
bool dist_t::is_forward_pointer() const {
    return (m_dist.lb() >= number_t{PACKET_BEGIN});
}
bool dist_t::is_backward_pointer() const {
    return (m_dist.ub() <= number_t{PACKET_END});
}

std::ostream& operator<<(std::ostream& o, const dist_t& d) {
    d.write(o);
    return o;
}

bool inequality_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_slack == -1 && m_value.is_top());
}

bool inequality_t::is_bottom() const {
    return m_is_bottom;
}

void inequality_t::set_to_top() {
    m_value = interval_t::top();
    m_slack = -1;
    m_is_bottom = false;
}

void inequality_t::set_to_bottom() {
    m_is_bottom = true;
}


std::ostream& operator<<(std::ostream& o, const inequality_t& ineq) {
    ineq.write(o);
    return o;
}

void inequality_t::write(std::ostream& o) const {
    o << m_slack << (m_rel == rop_t::R_GT ? ">" :
            m_rel == rop_t::R_GE ? ">=" :
            m_rel == rop_t::R_LT ? "<" : "<=")
        << m_value;
}

bool equality_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_lhs.is_top() && m_rhs.is_top());
}

bool equality_t::is_bottom() const {
    return m_is_bottom;
}

void equality_t::set_to_top() {
    m_lhs.set_to_top();
    m_rhs.set_to_top();
    m_is_bottom = false;
}

void equality_t::set_to_bottom() {
    m_is_bottom = true;
}

std::ostream& operator<<(std::ostream& o, const equality_t& eq) {
    eq.write(o);
    return o;
}

void equality_t::write(std::ostream& o) const {
    o << m_lhs << " = " << m_rhs;
}

void registers_state_t::insert(register_t reg, const reg_with_loc_t& reg_with_loc,
        const dist_t& dist) {
    (*m_offset_env)[reg_with_loc] = dist;
    m_cur_def[reg] = std::make_shared<reg_with_loc_t>(reg_with_loc);
}

std::optional<dist_t> registers_state_t::find(reg_with_loc_t reg) const {
    auto it = m_offset_env->find(reg);
    if (it == m_offset_env->end()) return {};
    return it->second;
}

std::optional<dist_t> registers_state_t::find(register_t key) const {
    if (m_cur_def[key] == nullptr) return {};
    const reg_with_loc_t& reg = *(m_cur_def[key]);
    return find(reg);
}

std::vector<uint64_t> stack_state_t::find_overlapping_cells(uint64_t start, int width) const {
    std::vector<uint64_t> overlapping_cells;
    auto it = m_slot_dists.begin();
    while (it != m_slot_dists.end() && it->first < start) {
        it++;
    }
    if (it != m_slot_dists.begin()) {
        it--;
        auto key = it->first;
        auto width_key = it->second.second;
        if (key < start && key+width_key > start) overlapping_cells.push_back(key);
    }

    for (; it != m_slot_dists.end(); it++) {
        auto key = it->first;
        if (key >= start && key < start+width) overlapping_cells.push_back(key);
        if (key >= start+width) break;
    }
    return overlapping_cells;
}

void registers_state_t::set_to_top() {
    m_cur_def = live_registers_t{nullptr};
    m_is_bottom = false;
}

void registers_state_t::set_to_bottom() {
    m_cur_def = live_registers_t{nullptr};
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
    live_registers_t out_vars;
    location_t loc = location_t(std::make_pair(label_t(-2, -2), 0));

    for (size_t i = 0; i < m_cur_def.size(); i++) {
        if (m_cur_def[i] == nullptr || other.m_cur_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_def[i]));
        auto it2 = other.find(*(other.m_cur_def[i]));
        if (it1 && it2) {
            dist_t dist1 = it1.value(), dist2 = it2.value();
            auto reg = reg_with_loc_t((register_t)i, loc);
            if (dist1.m_slack != dist2.m_slack) continue;
            auto dist_joined = dist_t(std::move(dist1.m_dist | dist2.m_dist), dist1.m_slack);
            out_vars[i] = std::make_shared<reg_with_loc_t>(reg);
            (*m_offset_env)[reg] = dist_joined;
        }
    }
    return registers_state_t(std::move(out_vars), m_offset_env, false);
}

void registers_state_t::adjust_bb_for_registers(location_t loc) {
    location_t old_loc = location_t(std::make_pair(label_t(-2, -2), 0));
    for (size_t i = 0; i < m_cur_def.size(); i++) {
        auto new_reg = reg_with_loc_t((register_t)i, loc);
        auto old_reg = reg_with_loc_t((register_t)i, old_loc);

        auto it = find((register_t)i);
        if (!it) continue;

        if (*m_cur_def[i] == old_reg)
            m_offset_env->erase(old_reg);

        m_cur_def[i] = std::make_shared<reg_with_loc_t>(new_reg);
        (*m_offset_env)[new_reg] = it.value();
    }
}

void stack_state_t::set_to_top() {
    m_slot_dists.clear();
    m_is_bottom = false;
}

void stack_state_t::set_to_bottom() {
    m_slot_dists.clear();
    m_is_bottom = true;
}

bool stack_state_t::is_top() const {
    if (m_is_bottom) return false;
    return m_slot_dists.empty();
}

bool stack_state_t::is_bottom() const {
    return m_is_bottom;
}

stack_state_t stack_state_t::top() {
    return stack_state_t(false);
}

std::optional<dist_cells_t> stack_state_t::find(uint64_t key) const {
    auto it = m_slot_dists.find(key);
    if (it == m_slot_dists.end()) return {};
    return it->second;
}

void stack_state_t::store(uint64_t key, dist_t d, int width) {
    m_slot_dists[key] = std::make_pair(d, width);
}

void stack_state_t::operator-=(uint64_t to_erase) {
    if (is_bottom()) {
        return;
    }
    m_slot_dists.erase(to_erase);
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
    stack_slot_dists_t out_stack_dists;
    // We do not join dist cells because different dist values different types of offsets
    for (auto const&kv: m_slot_dists) {
        auto maybe_dist_cells = other.find(kv.first);
        if (maybe_dist_cells) {
            auto dist_cells1 = kv.second;
            auto dist_cells2 = *maybe_dist_cells;
            auto dist1 = dist_cells1.first;
            auto dist2 = dist_cells2.first;
            int width1 = dist_cells1.second;
            int width2 = dist_cells2.second;
            if (dist1 == dist2 && width1 == width2) {
                out_stack_dists.insert(kv);
            }
        }
    }
    return stack_state_t(std::move(out_stack_dists), false);
}

bool extra_constraints_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_meta_and_begin.is_top() && m_begin_and_end.is_top());
}

bool extra_constraints_t::is_bottom() const {
    return m_is_bottom;
}

void extra_constraints_t::set_to_top() {
    m_meta_and_begin.set_to_top();
    m_begin_and_end.set_to_top();
    m_is_bottom = false;
}

void extra_constraints_t::set_to_bottom() {
    m_is_bottom = true;
}

void extra_constraints_t::add_meta_and_begin_constraint(equality_t&& eq,
        inequality_t&& ineq) {
    m_meta_and_begin = packet_constraint_t(std::move(eq), std::move(ineq), true);
}

void extra_constraints_t::add_begin_and_end_constraint(equality_t&& eq,
        inequality_t&& ineq) {
    m_begin_and_end = packet_constraint_t(std::move(eq), std::move(ineq), false);
}
/*
void extra_constraints_t::normalize() {
    weight_t dist_lhs = m_eq.m_lhs.m_dist - m_eq.m_rhs.m_dist - 4099;
    weight_t dist_rhs = -4099;
    slack_var_t s = m_eq.m_lhs.m_slack;
    dist_lhs += m_ineq.m_value;
    weight_t ineq_val = 0;
    rop_t ineq_rel = m_ineq.m_rel;

    m_eq = equality_t(dist_t(dist_lhs, s), dist_t(dist_rhs));
    m_ineq = inequality_t(s, ineq_rel, ineq_val);
}
*/

packet_constraint_t packet_constraint_t::operator|(const packet_constraint_t& other) const {
    //normalize();
    //other.normalize();

    weight_t dist1 = m_eq.m_lhs.m_dist;
    weight_t dist2 = other.m_eq.m_lhs.m_dist;
    slack_var_t s = m_eq.m_lhs.m_slack;

    dist_t lhs = dist_t(dist1 | dist2, s);
    dist_t rhs;
    if (m_is_meta_constraint) rhs = dist_t(weight_t{number_t{PACKET_BEGIN}});
    else rhs = dist_t(weight_t{number_t{PACKET_END}});

    equality_t out_eq(lhs, rhs);
    inequality_t out_ineq(s, m_ineq.m_rel, weight_t{number_t{0}});
    return packet_constraint_t(std::move(out_eq), std::move(out_ineq), m_is_meta_constraint);
        // have to handle case for different slack vars
}

std::ostream& operator<<(std::ostream& o, const packet_constraint_t& p) {
    p.write(o);
    return o;
}

void packet_constraint_t::write(std::ostream& o) const {
    o << m_eq << "\n";
    o << m_ineq << "\n";
}

void packet_constraint_t::set_to_top() {
    m_eq.set_to_top();
    m_ineq.set_to_top();
    m_is_bottom = false;
}

void packet_constraint_t::set_to_bottom() {
    m_is_bottom = true;
}

bool packet_constraint_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_eq.is_top() && m_ineq.is_top());
}

bool packet_constraint_t::is_bottom() const {
    return m_is_bottom;
}

std::optional<bound_t> packet_constraint_t::get_limit() const {
    // TODO: normalize constraint, if required
    auto dist = m_eq.m_lhs.m_dist;
    if (dist.is_top()) return {};
    return dist.ub();
}

extra_constraints_t extra_constraints_t::operator|(const extra_constraints_t& other) const {
    auto meta_and_begin = m_meta_and_begin | other.m_meta_and_begin;
    auto begin_and_end = m_begin_and_end | other.m_begin_and_end;
    return extra_constraints_t(std::move(meta_and_begin), std::move(begin_and_end), false);
}

std::optional<bound_t> extra_constraints_t::get_end_limit() const {
    return m_begin_and_end.get_limit();
}

std::optional<bound_t> extra_constraints_t::get_meta_limit() const {
    return m_meta_and_begin.get_limit();
}

ctx_offsets_t::ctx_offsets_t(const ebpf_context_descriptor_t* desc) {
    if (desc->data >= 0) {
        m_dists[desc->data] = dist_t(weight_t{number_t{PACKET_BEGIN}});
    }
    if (desc->end >= 0) {
        m_dists[desc->end] = dist_t(weight_t{number_t{PACKET_END}});
    }
    if (desc->meta >= 0) {
        m_dists[desc->meta] = dist_t(weight_t{number_t{PACKET_META}});
    }
    if (desc->size >= 0) {
        m_size = desc->size;
    }
}

int ctx_offsets_t::get_size() const {
    return m_size;
}

std::optional<dist_t> ctx_offsets_t::find(int key) const {
    auto it = m_dists.find(key);
    if (it == m_dists.end()) return {};
    return it->second;
}

offset_domain_t offset_domain_t::setup_entry() {
    std::shared_ptr<ctx_offsets_t> ctx = std::make_shared<ctx_offsets_t>(global_program_info->type.context_descriptor);
    std::shared_ptr<global_offset_env_t> all_types = std::make_shared<global_offset_env_t>();
    registers_state_t regs(all_types);

    offset_domain_t off_d(std::move(regs), stack_state_t::top(), ctx);
    return off_d;
}

offset_domain_t offset_domain_t::bottom() {
    offset_domain_t off;
    off.set_to_bottom();
    return off;
}

void offset_domain_t::set_to_top() {
    m_reg_state.set_to_top();
    m_stack_state.set_to_top();
    m_extra_constraints.set_to_top();
}

void offset_domain_t::set_to_bottom() {
    m_is_bottom = true;
}

bool offset_domain_t::is_bottom() const {
    return m_is_bottom;
}

bool offset_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_reg_state.is_top() && m_stack_state.is_top() && m_extra_constraints.is_top());
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
    return offset_domain_t(m_reg_state | other.m_reg_state, m_stack_state | other.m_stack_state, m_extra_constraints | other.m_extra_constraints, m_ctx_dists, std::max(m_slack, other.m_slack));
}

offset_domain_t offset_domain_t::operator|(offset_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return offset_domain_t(m_reg_state | std::move(other.m_reg_state), m_stack_state | std::move(other.m_stack_state), m_extra_constraints | std::move(other.m_extra_constraints), m_ctx_dists, std::max(m_slack, other.m_slack));
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

//forget
void offset_domain_t::operator-=(variable_t var) {}

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

void offset_domain_t::operator()(const Assume &b, location_t loc, int print) {
    Condition cond = b.cond;
    if (cond.op == Condition::Op::LE) {
        if (std::holds_alternative<Reg>(cond.right)) {
            auto right_reg = std::get<Reg>(cond.right).v;
            auto dist_left = m_reg_state.find(cond.left.v);
            auto dist_right = m_reg_state.find(right_reg);
            if (!dist_left || !dist_right) {
                // this should not happen, comparison between a packet pointer and either
                // other region's pointers or numbers; possibly raise type error
                m_errors.push_back("one of the pointers being compared isn't packet pointer");
                //std::cout << "type_error: one of the pointers being compared isn't packet pointer\n";
                return;
            }
            dist_t left_reg_dist = dist_left.value();
            dist_t right_reg_dist = dist_right.value();
            slack_var_t s = m_slack++;
            dist_t f = dist_t(left_reg_dist.m_dist, s);
            dist_t b = dist_t(right_reg_dist.m_dist);
            auto eq = equality_t(f, b);
            auto ineq = inequality_t(s, rop_t::R_GE, weight_t{number_t{0}});
            if (f.is_meta_pointer() && b.is_forward_pointer()) {
                m_extra_constraints.add_meta_and_begin_constraint(std::move(eq), std::move(ineq));
            }
            else if (f.is_forward_pointer() && b.is_backward_pointer()) {
                m_extra_constraints.add_begin_and_end_constraint(std::move(eq), std::move(ineq));
            }
        }
    }
    else {}     //we do not need to deal with other cases
}

void offset_domain_t::update_offset_info(const dist_t&& dist, const interval_t&& change,
        const reg_with_loc_t& reg_with_loc, uint8_t reg, Bin::Op op) {
    auto offset = dist.m_dist;
    if (op == Bin::Op::ADD) {
        if (dist.is_forward_pointer()) offset += change;
        else if (dist.is_backward_pointer()) offset -= change;
        else offset -= change;
    }
    else if (op == Bin::Op::SUB) {
        // TODO: needs precise handling of subtraction
        offset = interval_t::top();
    }
    m_reg_state.insert(reg, reg_with_loc, dist_t(offset));
}

interval_t offset_domain_t::do_bin(const Bin &bin,
        const std::optional<interval_t>& src_interval_opt,
        const std::optional<interval_t>& dst_interval_opt,
        std::optional<ptr_or_mapfd_t>& src_ptr_or_mapfd_opt,
        std::optional<ptr_or_mapfd_t>& dst_ptr_or_mapfd_opt, location_t loc) {

    using Op = Bin::Op;
    // if both src and dst are numbers, nothing to do in offset domain
    // if we are doing a move, where src is a number and dst is not set, nothing to do
    if ((dst_interval_opt && src_interval_opt)
            || (src_interval_opt && !dst_ptr_or_mapfd_opt && bin.op == Op::MOV))
        return interval_t::bottom();
    // offset domain only handles packet pointers
    if (!is_packet_ptr(src_ptr_or_mapfd_opt) && !is_packet_ptr(dst_ptr_or_mapfd_opt))
        return interval_t::bottom();

    interval_t src_interval = interval_t::bottom(), dst_interval = interval_t::bottom();
    if (src_interval_opt) src_interval = std::move(src_interval_opt.value());
    if (dst_interval_opt) dst_interval = std::move(dst_interval_opt.value());

    Reg src;
    if (std::holds_alternative<Reg>(bin.v)) src = std::get<Reg>(bin.v);

    auto reg_with_loc = reg_with_loc_t(bin.dst.v, loc);
    switch (bin.op)
    {
        // ra = rb;
        case Op::MOV: {
            if (!is_packet_ptr(src_ptr_or_mapfd_opt)) {
                m_reg_state -= bin.dst.v;
                return interval_t::bottom();
            }
            auto src_offset_opt = m_reg_state.find(src.v);
            if (!src_offset_opt) {
                m_errors.push_back("src is a packet_pointer and no offset info found");
                //std::cout << "type_error: src is a packet_pointer and no offset info found\n";
                return interval_t::bottom();
            }
            m_reg_state.insert(bin.dst.v, reg_with_loc, src_offset_opt.value());
            break;
        }
        // ra += rb
        case Op::ADD: {
            dist_t dist_to_update;
            interval_t interval_to_add = interval_t::bottom();
            if (is_packet_ptr(dst_ptr_or_mapfd_opt)
                    && is_packet_ptr(src_ptr_or_mapfd_opt)) {
                m_reg_state -= bin.dst.v;
                return interval_t::bottom();
            }
            else if (is_packet_ptr(dst_ptr_or_mapfd_opt) && src_interval_opt) {
                auto dst_offset_opt = m_reg_state.find(bin.dst.v);
                if (!dst_offset_opt) {
                    m_errors.push_back("dst is a packet_pointer and no offset info found");
                    //std::cout << "type_error: dst is a packet_pointer and no offset info found\n";
                    m_reg_state -= bin.dst.v;
                    return interval_t::bottom();
                }
                dist_to_update = std::move(dst_offset_opt.value());
                interval_to_add = std::move(src_interval_opt.value());
            }
            else {
                auto src_offset_opt = m_reg_state.find(src.v);
                if (!src_offset_opt) {
                    m_errors.push_back("src is a packet_pointer and no offset info found");
                    //std::cout << "type_error: src is a packet_pointer and no offset info found\n";
                    m_reg_state -= bin.dst.v;
                    return interval_t::bottom();
                }
                dist_to_update = std::move(src_offset_opt.value());
                interval_to_add = std::move(dst_interval_opt.value());
            }
            update_offset_info(std::move(dist_to_update), std::move(interval_to_add),
                    reg_with_loc, bin.dst.v, bin.op);
            break;
        }
        // ra -= rb
        case Op::SUB: {
            dist_t dist_to_update;
            interval_t interval_to_sub = interval_t::bottom();
            if (is_packet_ptr(dst_ptr_or_mapfd_opt)
                    && is_packet_ptr(src_ptr_or_mapfd_opt)) {
                m_reg_state -= bin.dst.v;
                return interval_t::top();
            }
            else if (is_packet_ptr(dst_ptr_or_mapfd_opt) && src_interval_opt) {
                auto dst_offset_opt = m_reg_state.find(bin.dst.v);
                if (!dst_offset_opt) {
                    m_errors.push_back("dst is a packet_pointer and no offset info found");
                    //std::cout << "type_error: dst is a packet_pointer and no offset info found\n";
                    m_reg_state -= bin.dst.v;
                    return interval_t::bottom();
                }
                dist_to_update = std::move(dst_offset_opt.value());
                interval_to_sub = std::move(src_interval_opt.value());
            }
            else {
                auto src_offset_opt = m_reg_state.find(src.v);
                if (!src_offset_opt) {
                    m_errors.push_back("src is a packet_pointer and no offset info found");
                    //std::cout << "type_error: src is a packet_pointer and no offset info found\n";
                    m_reg_state -= bin.dst.v;
                    return interval_t::bottom();
                }
                dist_to_update = std::move(src_offset_opt.value());
                interval_to_sub = std::move(dst_interval_opt.value());
            }
            update_offset_info(std::move(dist_to_update), std::move(interval_to_sub),
                    reg_with_loc, bin.dst.v, bin.op);
            break;
        }
        default: {
            m_reg_state -= bin.dst.v;
            break;
        }
    }
    return interval_t::bottom();
}

void offset_domain_t::operator()(const Bin& bin, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const Undefined& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const Un& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const LoadMapFd& u, location_t loc, int print) {
    m_reg_state -= u.dst.v;
}

void offset_domain_t::operator()(const Call& u, location_t loc, int print) {
    m_reg_state -= register_t{R0_RETURN_VALUE};
}

void offset_domain_t::operator()(const Exit& u, location_t loc, int print) {}

void offset_domain_t::operator()(const Jmp& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const Packet& u, location_t loc, int print) {
    m_reg_state -= register_t{R0_RETURN_VALUE};
}

void offset_domain_t::operator()(const ValidDivisor& u, location_t loc, int print) {
    /* WARNING: This operation is not implemented yet. */
}

void offset_domain_t::operator()(const ValidAccess& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const Comparable& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const Addable& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const ValidStore& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const TypeConstraint& u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const ValidSize& u, location_t loc, int print) {
    /* WARNING: This operation is not implemented yet. */
}

void offset_domain_t::operator()(const ValidMapKeyValue& u, location_t loc, int print) {
    /* WARNING: This operation is not implemented yet. */
}

void offset_domain_t::operator()(const ZeroCtxOffset&, location_t loc, int print) {
    // nothing to do here
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
    auto it = m_reg_state.find(r.v);
    if (!it) return false;
    dist_t dist = it.value();

    return (lower_bound_satisfied(dist, offset)
            && upper_bound_satisfied(dist, offset, width, is_comparison_check));
}

void offset_domain_t::check_valid_access(const ValidAccess& s,
        std::optional<ptr_or_mapfd_t>& reg_type,
        std::optional<interval_t>, std::optional<interval_t>) {
    if (std::holds_alternative<Reg>(s.width)) return;
    int w = std::get<Imm>(s.width).v;
    if (w == 0 || !reg_type) return;

    bool is_comparison_check = s.width == (Value)Imm{0};
    if (check_packet_access(s.reg, w, s.offset, is_comparison_check)) return;
    m_errors.push_back("valid access check failed");
    //std::cout << "type_error: valid access assert fail\n";
    //exit(1);
}

void offset_domain_t::operator()(const Assert &u, location_t loc, int print) {
    // nothing to do here
}

void offset_domain_t::operator()(const basic_block_t& bb, int print) {
    // nothing to do here
}

void offset_domain_t::do_mem_store(const Mem& b, std::optional<ptr_or_mapfd_t> maybe_targetreg_type, std::optional<ptr_or_mapfd_t>& maybe_basereg_type) {
    bool target_is_reg = std::holds_alternative<Reg>(b.value);
    if (target_is_reg) {
    }

    int offset = b.access.offset;
    int width = b.access.width;

    if (is_stack_ptr(maybe_basereg_type)) {
        auto basereg_with_off = std::get<ptr_with_off_t>(*maybe_basereg_type);
        auto basereg_off_singleton = basereg_with_off.get_offset().to_interval().singleton();
        if (!basereg_off_singleton) return;
        auto store_at = (uint64_t)(*basereg_off_singleton + offset);
        auto overlapping_cells = m_stack_state.find_overlapping_cells(store_at, width);
        m_stack_state -= overlapping_cells;

        if (!is_packet_ptr(maybe_targetreg_type)) return;
        auto target_reg = std::get<Reg>(b.value);
        auto offset_info = m_reg_state.find(target_reg.v);
        if (!offset_info) {
            m_errors.push_back("register is a packet_pointer and no offset info found");
            //std::cout << "type_error: register is a packet_pointer and no offset info found\n";
                return;
        }
        m_stack_state.store(store_at, *offset_info, width);
    }
}

void offset_domain_t::do_load(const Mem& b, const Reg& target_reg,
        std::optional<ptr_or_mapfd_t> basereg_type, location_t loc) {
    if (!basereg_type) {
        m_reg_state -= target_reg.v;
        return;
    }
    auto basereg_ptr_type = basereg_type.value();
    int offset = b.access.offset;
    
    if (std::holds_alternative<ptr_with_off_t>(basereg_ptr_type)) {
        auto p_with_off = std::get<ptr_with_off_t>(basereg_ptr_type);
        auto offset_singleton = p_with_off.get_offset().to_interval().singleton();
        if (!offset_singleton) {
            m_errors.push_back("basereg offset is not a singleton");
            //std::cout << "type_error: basereg offset is not a singleton\n";
            m_reg_state -= target_reg.v;
            return;
        }

        if (p_with_off.get_region() == crab::region_t::T_CTX) {
            if (!offset_singleton) {
                m_reg_state -= target_reg.v;
                return;
            }
            auto load_at = (uint64_t)offset_singleton.value() + (uint64_t)offset;
            auto it = m_ctx_dists->find(load_at);
            if (!it) {
                m_reg_state -= target_reg.v;
                return;
            }
            dist_t d = it.value();
            auto reg = reg_with_loc_t(target_reg.v, loc);
            m_reg_state.insert(target_reg.v, reg, dist_t(d));
        }
        else if (p_with_off.get_region() == crab::region_t::T_STACK) {
            if (!offset_singleton) {
                m_reg_state -= target_reg.v;
                return;
            }
            auto ptr_offset = offset_singleton.value();
            auto load_at = (uint64_t)(ptr_offset + offset);
            auto it = m_stack_state.find(load_at);

            if (!it) {
                m_reg_state -= target_reg.v;
                return;
            }
            dist_t d = it->first;
            auto reg = reg_with_loc_t(target_reg.v, loc);
            m_reg_state.insert(target_reg.v, reg, d);
        }
        else {  // shared
            m_reg_state -= target_reg.v;
        }
    }
    else {  // we are loading from packet, or we have mapfd
        m_reg_state -= target_reg.v;
    }
}

void offset_domain_t::operator()(const Mem& b, location_t loc, int print) {
    // nothing to do here
}

std::optional<dist_t> offset_domain_t::find_offset_at_loc(const reg_with_loc_t reg) const {
    return m_reg_state.find(reg);
}

std::optional<dist_t> offset_domain_t::find_in_ctx(int key) const {
    return m_ctx_dists->find(key);
}

std::optional<dist_cells_t> offset_domain_t::find_in_stack(int key) const {
    return m_stack_state.find(key);
}

std::optional<dist_t> offset_domain_t::find_offset_info(register_t reg) const {
    return m_reg_state.find(reg);
}

void offset_domain_t::adjust_bb_for_types(location_t loc) {
    m_reg_state.adjust_bb_for_registers(loc);
}

} // namespace crab
