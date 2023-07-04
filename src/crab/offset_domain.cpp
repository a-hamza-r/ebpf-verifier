// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/offset_domain.hpp"

#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)

#define STACK_BEGIN 0
#define STACK_END 512
#define CTX_BEGIN 0
#define PACKET_BEGIN 0

namespace std {
    template <>
    struct hash<crab::reg_with_loc_t> {
        size_t operator()(const crab::reg_with_loc_t& reg) const { return reg.hash(); }
    };
}

bool dist_t::operator==(const dist_t& d) const {
    return (m_dist == d.m_dist && m_slack == d.m_slack);
}

void dist_t::write(std::ostream& o) const {
    if (m_slack != -1)
        o << "s" << m_slack << "+";
    if (m_dist >= 0)
        o << "begin+" << m_dist;
    else if (m_dist >= -4098)
        o << "meta";
    else
        o << "end-" << (-1)*m_dist-4099;
}

std::ostream& operator<<(std::ostream& o, const dist_t& d) {
    d.write(o);
    return o;
}

void registers_state_t::insert(register_t reg, const reg_with_loc_t& reg_with_loc, const dist_t& dist) {
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
    for (size_t i = 0; i < m_cur_def.size(); i++) {
        if (m_cur_def[i] == nullptr || other.m_cur_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_def[i]));
        auto it2 = other.find(*(other.m_cur_def[i]));
        if (it1 && it2) {
            dist_t dist1 = it1.value(), dist2 = it2.value();
            if (dist1 == dist2) {
                out_vars[i] = m_cur_def[i];
            }
        }
    }

    return registers_state_t(std::move(out_vars), m_offset_env, false);
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

std::optional<dist_t> stack_state_t::find(int key) const {
    auto it = m_slot_dists.find(key);
    if (it == m_slot_dists.end()) return {};
    return it->second;
}

void stack_state_t::store(int key, dist_t d) {
    m_slot_dists[key] = d;
}

void stack_state_t::operator-=(int to_erase) {
    if (is_bottom()) {
        return;
    }
    m_slot_dists.erase(to_erase);
}

stack_state_t stack_state_t::operator|(const stack_state_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    stack_slot_dists_t out_stack_dists;
    for (auto const&kv: m_slot_dists) {
        auto it = other.m_slot_dists.find(kv.first);
        if (it != m_slot_dists.end() && kv.second == it->second)
            out_stack_dists.insert(kv);
    }
    return stack_state_t(std::move(out_stack_dists), false);
}

void extra_constraints_t::set_to_top() {
    add_equality(forward_and_backward_eq_t());
    add_inequality(inequality_t());
}

void extra_constraints_t::set_to_bottom() {
    m_is_bottom = true;
}

bool extra_constraints_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_eq.m_forw.m_slack == -1 && m_ineq.m_slack == -1);
}

bool extra_constraints_t::is_bottom() const {
    return m_is_bottom;
}

void extra_constraints_t::add_equality(forward_and_backward_eq_t fabeq) {
    m_eq = std::move(fabeq);
}

void extra_constraints_t::add_inequality(inequality_t ineq) {
    m_ineq = std::move(ineq);
}

weight_t extra_constraints_t::get_limit() const {
    return m_eq.m_forw.m_dist;
}

void extra_constraints_t::normalize() {
    weight_t dist_forw = m_eq.m_forw.m_dist - m_eq.m_backw.m_dist - 4099;
    weight_t dist_backw = -4099;
    slack_var_t s = m_eq.m_forw.m_slack;
    dist_forw += m_ineq.m_value;
    weight_t ineq_val = 0;
    rop_t ineq_rel = m_ineq.m_rel;

    m_eq = forward_and_backward_eq_t(dist_t(dist_forw, s), dist_t(dist_backw));
    m_ineq = inequality_t(s, ineq_rel, ineq_val);
}

extra_constraints_t extra_constraints_t::operator|(const extra_constraints_t& other) const {
    //normalize();
    //other.normalize();

    weight_t dist1 = m_eq.m_forw.m_dist;
    weight_t dist2 = other.m_eq.m_forw.m_dist;
    slack_var_t s = m_eq.m_forw.m_slack;

    dist_t f = dist_t(min(dist1, dist2), s);
    dist_t b = dist_t(-4099);

    forward_and_backward_eq_t out_eq(f, b);
    inequality_t out_ineq(s, m_ineq.m_rel, 0);

    return extra_constraints_t(std::move(out_eq), std::move(out_ineq), false);
        // have to handle case for different slack vars
}

ctx_t::ctx_t(const ebpf_context_descriptor_t* desc) {
    if (desc->data != -1) {
        m_dists[desc->data] = dist_t(0);
    }
    if (desc->end != -1) {
        m_dists[desc->end] = dist_t(-4099);
    }
    if (desc->meta != -1) {
        m_dists[desc->meta] = dist_t(-1);
    }
    m_size = desc->size;
}

int ctx_t::get_size() const {
    return m_size;
}

std::optional<dist_t> ctx_t::find(int key) const {
    auto it = m_dists.find(key);
    if (it == m_dists.end()) return {};
    return it->second;
}

offset_domain_t offset_domain_t::setup_entry() {
    std::shared_ptr<ctx_t> ctx = std::make_shared<ctx_t>(global_program_info.get().type.context_descriptor);
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
    return offset_domain_t(m_reg_state | other.m_reg_state, m_stack_state | other.m_stack_state, m_extra_constraints | other.m_extra_constraints, m_ctx_dists, max(m_slack, other.m_slack));
}

offset_domain_t offset_domain_t::operator|(offset_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return offset_domain_t(m_reg_state | std::move(other.m_reg_state), m_stack_state | std::move(other.m_stack_state), m_extra_constraints | std::move(other.m_extra_constraints), m_ctx_dists, max(m_slack, other.m_slack));
}

// meet
offset_domain_t offset_domain_t::operator&(const offset_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

// widening
offset_domain_t offset_domain_t::widen(const offset_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

// narrowing
offset_domain_t offset_domain_t::narrow(const offset_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

//forget
void offset_domain_t::operator-=(crab::variable_t var) {}

void offset_domain_t::write(std::ostream& os) const {}

std::string offset_domain_t::domain_name() const {
    return "offset_domain";
}

crab::bound_t offset_domain_t::get_instruction_count_upper_bound() {
    /* WARNING: The operation is not implemented yet.*/
    return crab::bound_t(crab::number_t(0));
}

string_invariant offset_domain_t::to_set() { return string_invariant{}; }

void offset_domain_t::operator()(const Assume &b, location_t loc, int print) {
    Condition cond = b.cond;
    if (cond.op == Condition::Op::LE) {
        if (std::holds_alternative<Reg>(cond.right)) {
            auto right_reg = std::get<Reg>(cond.right).v;
            auto dist_left = m_reg_state.find(cond.left.v);
            auto dist_right = m_reg_state.find(right_reg);
            if (!dist_left && !dist_right) {
                return;
            }
            else if (!dist_left || !dist_right) {
                // this should not happen, comparison between a packet pointer and either
                // other region's pointers or numbers; possibly raise type error
                //exit(1);
                std::cout << "type_error: one of the pointers being compared isn't packet pointer\n";
                return;
            }
            dist_t left_reg_dist = dist_left.value();
            dist_t right_reg_dist = dist_right.value();
            slack_var_t s = m_slack++;
            dist_t f = dist_t(left_reg_dist.m_dist, s);
            dist_t b = dist_t(right_reg_dist.m_dist, slack_var_t{-1});
            m_extra_constraints.add_equality(forward_and_backward_eq_t(f, b));
            m_extra_constraints.add_inequality(inequality_t(s, rop_t::R_GE, 0));
        }
    }
    else {}     //we do not need to deal with other cases
}

bool is_packet_pointer(std::optional<ptr_t>& type) {
    if (!type) {    // not a pointer
        return false;
    }
    ptr_t ptr_type = type.value();
    if (std::holds_alternative<ptr_no_off_t>(ptr_type)
        && std::get<ptr_no_off_t>(ptr_type).get_region() == crab::region::T_PACKET) {
        return true;
    }
    return false;
}

bool is_stack_pointer(std::optional<ptr_t>& type) {
    if (!type) {    // not a pointer
        return false;
    }
    ptr_t ptr_type = type.value();
    if (std::holds_alternative<ptr_with_off_t>(ptr_type)
        && std::get<ptr_with_off_t>(ptr_type).get_region() == crab::region::T_STACK) {
        return true;
    }
    return false;
}

void offset_domain_t::do_bin(const Bin &bin, std::shared_ptr<int> src_const_value,
        std::optional<ptr_t> src_type, std::optional<ptr_t> dst_type, location_t loc) {
    if (is_bottom()) return;

    auto reg_with_loc = reg_with_loc_t(bin.dst.v, loc);
    if (std::holds_alternative<Reg>(bin.v)) {
        Reg src = std::get<Reg>(bin.v);
        switch (bin.op)
        {
            // ra = rb;
            case Bin::Op::MOV: {
                if (!is_packet_pointer(src_type)) {
                    m_reg_state -= bin.dst.v;
                    return;
                }
                auto it = m_reg_state.find(src.v);
                if (!it) {
                    std::cout << "type_error: src is a packet_pointer and no offset info found\n";
                    //exit(1);
                    return;
                }
                m_reg_state.insert(bin.dst.v, reg_with_loc, it.value());
                //std::cout << "offset: " << (*m_reg_state.get(bin.dst.v)).m_dist << "\n";
                break;
            }
            // ra += rb
            case Bin::Op::ADD: {
                if (!is_packet_pointer(dst_type)) {
                    m_reg_state -= bin.dst.v;
                    return;
                }
                auto it = m_reg_state.find(bin.dst.v);
                if (!it) {
                    std::cout << "type_error: dst is a packet_pointer and no offset info found\n";
                    //exit(1);
                    return;
                }
                auto dst_dist = it.value();
                if (src_const_value) {
                    weight_t updated_dist;
                    if (dst_dist.m_dist >= 0) {
                        updated_dist = dst_dist.m_dist + (*src_const_value);
                    }
                    else if (dst_dist.m_dist >= -4098) {
                        // TODO: special handling of meta pointer required
                        updated_dist = dst_dist.m_dist - (*src_const_value);
                    }
                    else {
                        updated_dist = dst_dist.m_dist - (*src_const_value);
                    }
                    m_reg_state.insert(bin.dst.v, reg_with_loc, dist_t(updated_dist));
                    //std::cout << "offset: " << (*m_reg_state.get(bin.dst.v)).m_dist << "\n";
                }
                else {
                    m_reg_state -= bin.dst.v;
                }
                break;
            }

            default: {
                m_reg_state -= bin.dst.v;
                break;
            }
        }
    }
    else {
        int imm = static_cast<int>(std::get<Imm>(bin.v).v);
        auto it = m_reg_state.find(bin.dst.v);
        switch (bin.op)
        {
            case Bin::Op::ADD: {
                if (!is_packet_pointer(dst_type)) {
                    m_reg_state -= bin.dst.v;
                    return;
                }
                if (!it) {
                    std::cout << "type_error: dst is a packet_pointer and no offset info found\n";
                    //exit(1);
                    return;
                }
                auto dst_dist = it.value();

                weight_t updated_dist;
                if (dst_dist.m_dist >= 0) {
                    updated_dist = dst_dist.m_dist + imm;
                }
                else if (dst_dist.m_dist >= -4098) {
                    // TODO: special handling of meta pointer required
                    updated_dist = dst_dist.m_dist - imm;
                }
                else {
                    updated_dist = dst_dist.m_dist - imm;
                }
                m_reg_state.insert(bin.dst.v, reg_with_loc, dist_t(updated_dist));
                //std::cout << "offset: " << (*m_reg_state.get(bin.dst.v)).m_dist << "\n";
                break;
            }

            default: {
                m_reg_state -= bin.dst.v;
                break;
            }
        }
    }
}

void offset_domain_t::operator()(const Bin &bin, location_t loc, int print) {
    do_bin(bin, nullptr, {}, {}, loc);
}

void offset_domain_t::operator()(const Undefined &, location_t loc, int print) {}

void offset_domain_t::operator()(const Un &, location_t loc, int print) {}

void offset_domain_t::operator()(const LoadMapFd &, location_t loc, int print) {}

void offset_domain_t::operator()(const Call &, location_t loc, int print) {}

void offset_domain_t::operator()(const Exit &, location_t loc, int print) {}

void offset_domain_t::operator()(const Jmp &, location_t loc, int print) {}

void offset_domain_t::operator()(const Packet &, location_t loc, int print) {}

void offset_domain_t::operator()(const LockAdd &, location_t loc, int print) {}

void offset_domain_t::check_valid_access(const ValidAccess& s, std::optional<ptr_t>& reg_type) {
    if (std::holds_alternative<Imm>(s.width)) {
        int w = std::get<Imm>(s.width).v;
        if (w == 0 || !reg_type) return;

        m_extra_constraints.normalize();
        ptr_t reg_ptr_type = reg_type.value();
        if (std::holds_alternative<ptr_with_off_t>(reg_ptr_type)) {
            auto reg_with_off_ptr_type = std::get<ptr_with_off_t>(reg_ptr_type);
            int offset = reg_with_off_ptr_type.get_offset();
            int offset_to_check = offset+s.offset;
            if (reg_with_off_ptr_type.get_region() == crab::region::T_STACK) {
                if (offset_to_check >= STACK_BEGIN && offset_to_check+w <= STACK_END) return;
            }
            else {
                if (offset_to_check >= CTX_BEGIN && offset_to_check+w <= m_ctx_dists->get_size())
                    return;
            }
        }
        else {
            auto reg_no_off_ptr_type = std::get<ptr_no_off_t>(reg_ptr_type);
            if (reg_no_off_ptr_type.get_region() == crab::region::T_PACKET) {
                auto it = m_reg_state.find(s.reg.v);
                int limit = m_extra_constraints.get_limit();
                if (it) {
                    dist_t dist = it.value();
                    // TODO: handle meta and end pointers separately
                    if (dist.m_dist >= PACKET_BEGIN && dist.m_dist+w <= limit) return;
                }
            }
            else {
                return;
            }
        }
    }
    else {
        return;
    }
    std::cout << "valid access assert fail\n";
    //exit(1);
}

void offset_domain_t::operator()(const Assert &u, location_t loc, int print) {
    std::visit(*this, u.cst);
}

void offset_domain_t::operator()(const basic_block_t& bb, bool check_termination, int print) {
    for (const Instruction& statement : bb) {
        location_t loc = boost::none;
        std::visit([this, loc, print](const auto& v) { std::apply(*this, std::make_tuple(v, loc, print)); }, statement);
    }
}

void offset_domain_t::do_mem_store(const Mem& b, const Reg& target_reg, std::optional<ptr_t>& basereg_type, std::optional<ptr_t>& targetreg_type) {
    int offset = b.access.offset;

    if (is_stack_pointer(basereg_type)) {
        auto basereg_with_off = std::get<ptr_with_off_t>(basereg_type.value());
        int store_at = basereg_with_off.get_offset() + offset;
        if (is_packet_pointer(targetreg_type)) {
            ptr_t targetreg_ptr_type = targetreg_type.value();
            auto it = m_reg_state.find(target_reg.v);
            if (!it) {
                std::cout << "type_error: register is a packet_pointer and no offset info found\n";
                return;
            }
            m_stack_state.store(store_at, it.value());
        }
        else {
            m_stack_state -= store_at;
        }
    }
    else {}  // in the rest cases, we do not store
}

void offset_domain_t::do_load(const Mem& b, const Reg& target_reg, std::optional<ptr_t>& basereg_type, location_t loc) {
    if (!basereg_type) {
        m_reg_state -= target_reg.v;
        return;
    }
    ptr_t basereg_ptr_type = basereg_type.value();
    int offset = b.access.offset;
    
    if (std::holds_alternative<ptr_with_off_t>(basereg_ptr_type)) {
        auto p_with_off = std::get<ptr_with_off_t>(basereg_ptr_type);
        int to_load = p_with_off.get_offset() + offset;

        if (p_with_off.get_region() == crab::region::T_CTX) {
            auto it = m_ctx_dists->find(to_load);
            if (!it) {
                m_reg_state -= target_reg.v;
                return;
            }
            dist_t d = it.value();
            auto reg = reg_with_loc_t(target_reg.v, loc);
            m_reg_state.insert(target_reg.v, reg, dist_t(d));
            //std::cout << "offset: " << (*m_reg_state.get(target_reg.v)).m_dist << "\n";
        }
        else if (p_with_off.get_region() == crab::region::T_STACK) {
            auto it = m_stack_state.find(to_load);

            if (!it) {
                m_reg_state -= target_reg.v;
                return;
            }
            dist_t d = it.value();
            auto reg = reg_with_loc_t(target_reg.v, loc);
            m_reg_state.insert(target_reg.v, reg, dist_t(d));
            //std::cout << "offset: " << (*m_reg_state.get(target_reg.v)).m_dist << "\n";
        }
    }
    else {  // we are loading from packet or shared
        m_reg_state -= target_reg.v;
    }
}

void offset_domain_t::operator()(const Mem &b, location_t loc, int print) {
}

std::optional<dist_t> offset_domain_t::find_in_registers(const reg_with_loc_t reg) const {
    return m_reg_state.find(reg);
}

std::optional<dist_t> offset_domain_t::find_in_ctx(int key) const {
    return m_ctx_dists->find(key);
}

std::optional<dist_t> offset_domain_t::find_in_stack(int key) const {
    return m_stack_state.find(key);
}
