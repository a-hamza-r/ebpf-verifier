// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include <unordered_map>

#include "crab/type_domain.hpp"

using crab::___print___;
using crab::ptr_t;
using crab::ptr_with_off_t;
using crab::ptr_no_off_t;
using crab::ctx_t;
using crab::global_type_env_t;
using crab::reg_with_loc_t;
using crab::live_registers_t;
using crab::register_types_t;

static std::string size(int w) { return std::string("u") + std::to_string(w * 8); }


namespace std {
    template <>
    struct hash<crab::reg_with_loc_t> {
        size_t operator()(const crab::reg_with_loc_t& reg) const { return reg.hash(); }
    };

    // does not seem to work for me
    /*
    template <>
    struct equal_to<crab::ptr_t> {
        constexpr bool operator()(const crab::ptr_t& p1, const crab::ptr_t& p2) const {
            if (p1.index() != p2.index()) return false;
            if (std::holds_alternative<crab::ptr_no_off_t>(p1)) {
                auto ptr_no_off1 = std::get<crab::ptr_no_off_t>(p1);
                auto ptr_no_off2 = std::get<crab::ptr_no_off_t>(p2);
                return (ptr_no_off1.get_region() == ptr_no_off2.get_region());
            }
            else {
                auto ptr_with_off1 = std::get<crab::ptr_with_off_t>(p1);
                auto ptr_with_off2 = std::get<crab::ptr_with_off_t>(p2);
                return (ptr_with_off1.get_region() == ptr_with_off2.get_region() && ptr_with_off1.get_offset() == ptr_with_off2.get_offset());
            }
        }
    };

    template <>
    struct equal_to<crab::ptr_with_off_t> {
        constexpr bool operator()(const crab::ptr_with_off_t& p1, const crab::ptr_with_off_t& p2) const {
            return (p1.get_region() == p2.get_region() && p1.get_offset() == p2.get_offset());
        }
    };

    template <>
    struct equal_to<crab::ptr_no_off_t> {
        constexpr bool operator()(const crab::ptr_no_off_t& p1, const crab::ptr_no_off_t& p2) const {
            return (p1.get_region() == p2.get_region());
        }
    };
    */
}


static void print_ptr_type(const ptr_t& p) {
    if (std::holds_alternative<ptr_with_off_t>(p)) {
        auto t = std::get<ptr_with_off_t>(p);
        std::cout << t;
    }
    else {
        auto t = std::get<ptr_no_off_t>(p);
        std::cout << t;
    }
}

static void print_type(register_t r, const ptr_t& p) {
    std::cout << "r" << static_cast<unsigned int>(r) << " : ";
    print_ptr_type(p);
}

static void print_annotated(Mem const& b, const ptr_t& p, std::ostream& os_) {
    if (b.is_load) {
        os_ << "  ";
        print_type(std::get<Reg>(b.value).v, p);
        os_ << " = ";
    }
    std::string sign = b.access.offset < 0 ? " - " : " + ";
    int offset = std::abs(b.access.offset);
    os_ << "*(" << size(b.access.width) << " *)";
    os_ << "(" << b.access.basereg << sign << offset << ")\n";
}

static void print_annotated(Call const& call, const ptr_t& p, std::ostream& os_) {
    os_ << "  ";
    print_type(0, p);
    os_ << " = " << call.name << ":" << call.func << "(...)\n";
}

static void print_annotated(Bin const& b, const ptr_t& p, std::ostream& os_) {
    os_ << "  ";
    print_type(b.dst.v, p);
    // add better checks as we add more support
    if (std::holds_alternative<Reg>(b.v)) {
        if (b.op == Bin::Op::MOV)
            os_ << " = r" << static_cast<unsigned int>(std::get<Reg>(b.v).v) << ";";
        else if (b.op == Bin::Op::ADD)
            os_ << " += r" << static_cast<unsigned int>(std::get<Reg>(b.v).v) << ";";
    }
    else {
        if (b.op == Bin::Op::ADD)
            os_ << " += " << static_cast<int>(std::get<Imm>(b.v).v) << ";";
    }
}

namespace crab {

inline std::string get_reg_ptr(const region& r) {
    switch (r) {
        case region::T_CTX:
            return "ctx_p";
        case region::T_STACK:
            return "stack_p";
        case region::T_PACKET:
            return "packet_p";
        default:
            return "shared_p";
    }
}

inline std::ostream& operator<<(std::ostream& o, const region& t) {
    o << static_cast<std::underlying_type<region>::type>(t);
    return o;
}

bool operator==(const ptr_with_off_t& p1, const ptr_with_off_t& p2) {
    return (p1.get_region() == p2.get_region() && p1.get_offset() == p2.get_offset());
}

bool ptr_with_off_t::operator!=(const ptr_with_off_t& p2) {
    return !(*this == p2);
}

void ptr_with_off_t::write(std::ostream& o) const {
    o << get_reg_ptr(m_r) << "<" << m_offset << ">";
}

std::ostream& operator<<(std::ostream& o, const ptr_with_off_t& p) {
    p.write(o);
    return o;
}

void ptr_with_off_t::set_offset(int off) { m_offset = off; }

constexpr int ptr_with_off_t::get_offset() const { return m_offset; }

void ptr_with_off_t::set_region(region r) { m_r = r; }

constexpr region ptr_with_off_t::get_region() const { return m_r; }

bool operator==(const ptr_no_off_t& p1, const ptr_no_off_t& p2) {
    return (p1.get_region() == p2.get_region());
}

bool ptr_no_off_t::operator!=(const ptr_no_off_t& p2) {
    return !(*this == p2);
}

void ptr_no_off_t::write(std::ostream& o) const {
    o << get_reg_ptr(get_region());
}

std::ostream& operator<<(std::ostream& o, const ptr_no_off_t& p) {
    p.write(o);
    return o;
}

void ptr_no_off_t::set_region(region r) { m_r = r; }

constexpr region ptr_no_off_t::get_region() const { return m_r; }

void reg_with_loc_t::write(std::ostream& o) const {
    o << "r" << static_cast<unsigned int>(m_reg) << "@" << m_loc->second << " in " << m_loc->first << " ";
}

std::ostream& operator<<(std::ostream& o, const reg_with_loc_t& reg) {
    reg.write(o);
    return o;
}

bool reg_with_loc_t::operator==(const reg_with_loc_t& other) const {
    return (m_reg == other.m_reg && m_loc == other.m_loc);
}

std::size_t reg_with_loc_t::hash() const {
    // Similar to boost::hash_combine
    using std::hash;

    std::size_t seed = hash<register_t>()(m_reg);
    seed ^= hash<int>()(m_loc->first.from) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash<int>()(m_loc->first.to) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash<int>()(m_loc->second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

std::ostream& operator<<(std::ostream& o, const stack_t& st) {
    o << "Stack: ";
    if (st.is_bottom())
        o << "_|_\n";
    else {
        o << "{";
        for (auto s : st.m_ptrs) {
            o << s.first << ": ";
            print_ptr_type(s.second);
            o << ", ";
        }
        o << "}";
    }
    return o;
}

std::ostream& operator<<(std::ostream& o, const ctx_t& _ctx) {

    o << "type of context: " << (_ctx.m_packet_ptrs.empty() ? "top" : "") << "\n";
    for (const auto& it : _ctx.m_packet_ptrs) {
        o << "  stores at " << it.first << ": " << it.second << "\n";
    }
    return o;
}

ctx_t::ctx_t(const ebpf_context_descriptor_t* desc)
{
    if (desc->data != -1)
        m_packet_ptrs[desc->data] = crab::ptr_no_off_t(crab::region::T_PACKET);
    if (desc->end != -1)
        m_packet_ptrs[desc->end] = crab::ptr_no_off_t(crab::region::T_PACKET);
    if (desc->meta != -1)
        m_packet_ptrs[desc->meta] = crab::ptr_no_off_t(crab::region::T_PACKET);
}

std::optional<ptr_no_off_t> ctx_t::find(int key) const {
    auto it = m_packet_ptrs.find(key);
    if (it == m_packet_ptrs.end()) return {};
    return it->second;
}


std::ostream& operator<<(std::ostream& o, const register_types_t& typ) {
    if (typ.is_bottom())
        o << "_|_\n";
    else {
        for (const auto& v : *(typ.m_reg_type_env)) {
            o << v.first << ": ";
            print_ptr_type(v.second);
            o << "\n";
        }
    }
    return o;
}

register_types_t register_types_t::operator|(const register_types_t& other) const {
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
        if (it1 && it2 && it1.value() == it2.value()) {
            out_vars[i] = m_cur_def[i];
        }
    }

    return register_types_t(std::move(out_vars), m_reg_type_env, false);
}

void register_types_t::operator-=(register_t var) {
    if (is_bottom()) {
        return;
    }
    m_cur_def[var] = nullptr;
}

void register_types_t::set_to_bottom() {
    m_cur_def = live_registers_t{nullptr};
    m_is_bottom = true;
}

void register_types_t::set_to_top() {
    m_cur_def = live_registers_t{nullptr};
    m_is_bottom = false;
}

bool register_types_t::is_bottom() const { return m_is_bottom; }

bool register_types_t::is_top() const {
    if (m_is_bottom) { return false; }
    if (m_reg_type_env == nullptr) return true;
    for (auto it : m_cur_def) {
        if (it != nullptr) return false;
    }
    return true;
}

void register_types_t::insert(register_t reg, const reg_with_loc_t& reg_with_loc, const ptr_t& type) {
    (*m_reg_type_env)[reg_with_loc] = type;
    m_cur_def[reg] = std::make_shared<reg_with_loc_t>(reg_with_loc);
}

std::optional<ptr_t> register_types_t::find(reg_with_loc_t reg) const {
    auto it = m_reg_type_env->find(reg);
    if (it == m_reg_type_env->end()) return {};
    return it->second;
}

std::optional<ptr_t> register_types_t::find(register_t key) const {
    if (m_cur_def[key] == nullptr) return {};
    const reg_with_loc_t& reg = *(m_cur_def[key]);
    return find(reg);
}

stack_t stack_t::operator|(const stack_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    offset_to_ptr_t out_ptrs;
    for (auto const&kv: m_ptrs) {
        auto it = other.find(kv.first);
        if (it && kv.second == it.value())
            out_ptrs.insert(kv);
    }
    return stack_t(std::move(out_ptrs), false);
}

void stack_t::operator-=(int key) {
    auto it = find(key);
    if (it)
        m_ptrs.erase(key);
}

void stack_t::set_to_bottom() {
    m_ptrs.clear();
    m_is_bottom = true;
}

void stack_t::set_to_top() {
    m_ptrs.clear();
    m_is_bottom = false;
}

stack_t stack_t::bottom() { return stack_t(true); }

stack_t stack_t::top() { return stack_t(false); }

bool stack_t::is_bottom() const { return m_is_bottom; }

bool stack_t::is_top() const {
    if (m_is_bottom)
        return false;
    return m_ptrs.empty();
}

void stack_t::insert(int key, ptr_t value) {
    m_ptrs[key] = value;
}

std::optional<ptr_t> stack_t::find(int key) const {
    auto it = m_ptrs.find(key);
    if (it == m_ptrs.end()) return {};
    return it->second;
}

}

bool type_domain_t::is_bottom() const {
    if (m_is_bottom) return true;
    return (m_stack.is_bottom() || m_registers.is_bottom());
}

bool type_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_stack.is_top() && m_registers.is_top());
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
    m_stack.set_to_top();
    m_registers.set_to_top();
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
    return type_domain_t(m_registers | other.m_registers, m_stack | other.m_stack, other.m_ctx);
}

type_domain_t type_domain_t::operator|(type_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return type_domain_t(m_registers | std::move(other.m_registers), m_stack | std::move(other.m_stack), other.m_ctx);
}

type_domain_t type_domain_t::operator&(const type_domain_t& abs) const {
    return abs;
}

type_domain_t type_domain_t::widen(const type_domain_t& abs) const {
    return abs;
}

type_domain_t type_domain_t::narrow(const type_domain_t& other) const {
    return other;
}

void type_domain_t::write(std::ostream& os) const { 
    os << m_registers;
    os << m_stack << "\n";
}

std::string type_domain_t::domain_name() const {
    return "type_domain";
}

int type_domain_t::get_instruction_count_upper_bound() {
    return 0;
}

string_invariant type_domain_t::to_set() {
    return string_invariant{};
}

void type_domain_t::operator()(const Undefined & u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0)
        std::cout << "  " << u << ";\n";
}
void type_domain_t::operator()(const Un &u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0)
        std::cout << "  " << u << ";\n";
}
void type_domain_t::operator()(const LoadMapFd &u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0) {
        std::cout << "  " << u << ";\n";
        return;
    }
    m_registers -= u.dst.v;
}
void type_domain_t::operator()(const Call &u, location_t loc, int print) {
    if (is_bottom()) return;
    register_t r0_reg{R0_RETURN_VALUE};
    auto r0 = reg_with_loc_t(r0_reg, loc);
    if (print > 0) {
        if (u.is_map_lookup) {
            auto it = m_registers.find(r0);
            if (it) {
                print_annotated(u, it.value(), std::cout);
            }
        }
        else
            std::cout << "  " << u << ";\n";
        return;
    }
    if (u.is_map_lookup) {
        auto type = ptr_no_off_t(crab::region::T_SHARED);
        m_registers.insert(r0_reg, r0, type);
    }
    else {
        m_registers -= r0_reg;
    }
}
void type_domain_t::operator()(const Exit &u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0)
        std::cout << "  " << u << ";\n";
}
void type_domain_t::operator()(const Jmp &u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0)
        std::cout << "  " << u << ";\n";
}
void type_domain_t::operator()(const Packet & u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0) {
        std::cout << "  " << u << ";\n";
        return;
    }
    m_registers -= register_t{0};
}
void type_domain_t::operator()(const LockAdd &u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0)
        std::cout << "  " << u << ";\n";
}
void type_domain_t::operator()(const Assume &u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0)
        std::cout << "  " << u << ";\n";
}
void type_domain_t::operator()(const Assert &u, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0)
        std::cout << "  " << u << ";\n";
}

type_domain_t type_domain_t::setup_entry() {

    std::shared_ptr<ctx_t> ctx = std::make_shared<ctx_t>(global_program_info.type.context_descriptor);
    std::shared_ptr<global_type_env_t> all_types = std::make_shared<global_type_env_t>();

    live_registers_t vars;
    register_types_t typ(std::move(vars), all_types);

    auto r1 = reg_with_loc_t(R1_ARG, std::make_pair(label_t::entry, static_cast<unsigned int>(0)));
    auto r10 = reg_with_loc_t(R10_STACK_POINTER, std::make_pair(label_t::entry, static_cast<unsigned int>(0)));

    typ.insert(R1_ARG, r1, ptr_with_off_t(crab::region::T_CTX, 0));
    typ.insert(R10_STACK_POINTER, r10, ptr_with_off_t(crab::region::T_STACK, 512));

    type_domain_t inv(std::move(typ), crab::stack_t::top(), ctx);
    return inv;
}

void type_domain_t::report_type_error(std::string s, location_t loc) {
    std::cout << "type_error at line " << loc->second << " in bb " << loc->first << "\n";
    std::cout << s;
    error_location = loc;
    set_to_bottom();
}

void type_domain_t::operator()(const Bin& bin, location_t loc, int print) {
    if (is_bottom()) return;
    if (print > 0) {
        if (print == 2) {
            if ((std::holds_alternative<Reg>(bin.v) && (bin.op == Bin::Op::MOV || bin.op == Bin::Op::ADD)) ||
            (std::holds_alternative<Imm>(bin.v) && bin.op == Bin::Op::ADD)) {
                auto reg_with_loc = reg_with_loc_t(bin.dst.v, loc);
                auto it = m_registers.find(reg_with_loc);
                if (it) {
                    print_annotated(bin, it.value(), std::cout);
                    std::cout << "\n";
                    return;
                }
            }
        }
        std::cout << "  " << bin << ";\n";
        return;
    }

    ptr_t dst_reg;
    if (bin.op == Bin::Op::ADD) {
        auto it = m_registers.find(bin.dst.v);
        if (!it) {
            m_registers -= bin.dst.v;
            return;
        }
        dst_reg = it.value();
    }

    if (std::holds_alternative<Reg>(bin.v)) {
        Reg src = std::get<Reg>(bin.v);
        switch (bin.op)
        {
            case Bin::Op::MOV: {
                auto it1 = m_registers.find(src.v);
                if (!it1) {
                    //std::cout << "type_error: assigning an unknown pointer or a number - r" << (int)src.v << "\n";
                    m_registers -= bin.dst.v;
                    break;
                }
                auto reg = reg_with_loc_t(bin.dst.v, loc);
                m_registers.insert(bin.dst.v, reg, it1.value());
                break;
            }
            case Bin::Op::ADD: {
                auto it1 = m_registers.find(src.v);
                if (it1) {
                    std::string s = std::to_string(static_cast<unsigned int>(src.v));
                    std::string s1 = std::to_string(static_cast<unsigned int>(bin.dst.v));
                    std::string desc = std::string("\taddition of two pointers, r") + s + " and r" + s1 + " not allowed\n";
                    report_type_error(desc, loc);
                    return;
                }
                else {
                    if (std::holds_alternative<ptr_with_off_t>(dst_reg)) {
                        /*
                        std::string s = std::to_string(static_cast<unsigned int>(bin.dst.v));
                        std::string desc = std::string("\toffset of the pointer r") + s + " unknown\n";
                        report_type_error(desc, loc);
                        return;
                        */
                        m_stack.set_to_top();
                        m_stack -= std::get<ptr_with_off_t>(dst_reg).get_offset();
                    }
                    else {
                        auto reg = reg_with_loc_t(bin.dst.v, loc);
                        m_registers.insert(bin.dst.v, reg, dst_reg);
                    }
                }
                break;
            }
            default:
                m_registers -= bin.dst.v;
                break;
        }
    }
    else {
        int imm = static_cast<int>(std::get<Imm>(bin.v).v);
        switch (bin.op)
        {
            case Bin::Op::ADD: {
                if (std::holds_alternative<ptr_with_off_t>(dst_reg)) {
                    auto ptr_with_off = std::get<ptr_with_off_t>(dst_reg);
                    ptr_with_off.set_offset(ptr_with_off.get_offset() + imm);
                    auto reg = reg_with_loc_t(bin.dst.v, loc);
                    m_registers.insert(bin.dst.v, reg, ptr_with_off);
                }
                else {
                    auto reg = reg_with_loc_t(bin.dst.v, loc);
                    m_registers.insert(bin.dst.v, reg, dst_reg);
                }
                break;
            }
            default: {
                m_registers -= bin.dst.v;
                break;
            }
        }
    }
}

void type_domain_t::do_load(const Mem& b, const Reg& target_reg, location_t loc, int print) {

    if (print > 0) {
        auto target_reg_loc = reg_with_loc_t(target_reg.v, loc);
        auto it = m_registers.find(target_reg_loc);
        if (it)
            print_annotated(b, it.value(), std::cout);
        else
            std::cout << "  " << b << ";\n";
        return;
    }
    int offset = b.access.offset;
    Reg basereg = b.access.basereg;

    auto it = m_registers.find(basereg.v);
    if (!it) {
        std::string s = std::to_string(static_cast<unsigned int>(basereg.v));
        std::string desc = std::string("\tloading from an unknown pointer, or from number - r") + s + "\n";
        report_type_error(desc, loc);
        return;
    }
    ptr_t type_basereg = it.value();

    if (std::holds_alternative<ptr_no_off_t>(type_basereg)) {
        //std::cout << "type_error: loading from either packet or shared region not allowed - r" << (int)basereg.v << "\n";
        m_registers -= target_reg.v;
        return;
    }

    ptr_with_off_t type_with_off = std::get<ptr_with_off_t>(type_basereg);
    int load_at = offset+type_with_off.get_offset();

    switch (type_with_off.get_region()) {
        case crab::region::T_STACK: {

            auto it = m_stack.find(load_at);

            if (!it) {
                //std::cout << "type_error: no field at loaded offset " << load_at << " in stack\n";
                m_registers -= target_reg.v;
                return;
            }
            ptr_t type_loaded = it.value();

            if (std::holds_alternative<ptr_with_off_t>(type_loaded)) {
                ptr_with_off_t type_loaded_with_off = std::get<ptr_with_off_t>(type_loaded);
                auto reg = reg_with_loc_t(target_reg.v, loc);
                m_registers.insert(target_reg.v, reg, type_loaded_with_off);
            }
            else {
                ptr_no_off_t type_loaded_no_off = std::get<ptr_no_off_t>(type_loaded);
                auto reg = reg_with_loc_t(target_reg.v, loc);
                m_registers.insert(target_reg.v, reg, type_loaded_no_off);
            }

            break;
        }
        case crab::region::T_CTX: {

            auto it = m_ctx->find(load_at);

            if (!it) {
                //std::cout << "type_error: no field at loaded offset " << load_at << " in context\n";
                m_registers -= target_reg.v;
                return;
            }
            ptr_no_off_t type_loaded = it.value();

            auto reg = reg_with_loc_t(target_reg.v, loc);
            m_registers.insert(target_reg.v, reg, type_loaded);
            break;
        }

        default: {
            assert(false);
        }
    }
}

void type_domain_t::do_mem_store(const Mem& b, const Reg& target_reg, location_t loc, int print) {

    if (print > 0) {
        std::cout << "  " << b << ";\n";
        return;
    }
    int offset = b.access.offset;
    Reg basereg = b.access.basereg;
    int width = b.access.width;

    auto it = m_registers.find(basereg.v);
    if (!it) {
        std::string s = std::to_string(static_cast<unsigned int>(basereg.v));
        std::string desc = std::string("\tstoring at an unknown pointer, or from number - r") + s + "\n";
        report_type_error(desc, loc);
        return;
    }
    ptr_t type_basereg = it.value();

    auto it2 = m_registers.find(target_reg.v);

    if (std::holds_alternative<ptr_with_off_t>(type_basereg)) {
        // base register is either CTX_P or STACK_P
        ptr_with_off_t type_basereg_with_off = std::get<ptr_with_off_t>(type_basereg);

        int store_at = offset+type_basereg_with_off.get_offset();
        if (type_basereg_with_off.get_region() == crab::region::T_STACK) {
            // type of basereg is STACK_P
            if (!it2) {
               // std::cout << "type_error: storing either a number or an unknown pointer - r" << (int)target_reg.v << "\n";
                m_stack -= store_at;
                return;
            }
            else {
                auto type_to_store = it2.value();
                /*
                if (std::holds_alternative<ptr_with_off_t>(type_to_store) &&
                        std::get<ptr_with_off_t>(type_to_store).r == crab::region::T_STACK) {
                    std::string s = std::to_string(static_cast<unsigned int>(target_reg.v));
                    std::string desc = std::string("\twe cannot store stack pointer, r") + s + ", into stack\n";
                    //report_type_error(desc, loc);
                    return;
                }
                else {
                */
                    for (auto i = store_at+1; i < store_at+width; i++) {
                        auto it3 = m_stack.find(i);
                        if (it3) {
                            std::string s = std::to_string(store_at);
                            std::string s1 = std::to_string(i);
                            std::string desc = std::string("\ttype being stored into stack at ") + s + " is overlapping with already stored at " + s1 + "\n";
                            report_type_error(desc, loc);
                            return;
                        }
                    }
                    m_stack.insert(store_at, type_to_store);
                    // revise: code below checks if there is already something stored at same location, the type should be the same -- it is very restricted and not required.
                    // However, when we support storing info like width of type, we need more checks
                    /*
                    auto it4 = m_stack.find(store_at);
                    if (it4) {
                        auto type_in_stack = it4.value();
                        if (type_to_store != type_in_stack) {
                            std::string s = std::to_string(store_at);
                            std::string desc = std::string("\ttype being stored at offset ") + s + " is not the same as already stored in stack\n";
                            report_type_error(desc, loc);
                            return;
                        }
                    }
                    else {
                        m_stack.insert(store_at, type_to_store);
                    }
                    */
                //}
            }
        }
        else if (type_basereg_with_off.get_region() == crab::region::T_CTX) {
            // type of basereg is CTX_P
            if (it2) {
                std::string s = std::to_string(static_cast<unsigned int>(target_reg.v));
                std::string desc = std::string("\twe cannot store a pointer, r") + s + ", into ctx\n";
                report_type_error(desc, loc);
                return;
            }
        }
        else
            assert(false);
    }
    else {
        // base register type is either PACKET_P, SHARED_P or STACK_P without known offset
        ptr_no_off_t type_basereg_no_off = std::get<ptr_no_off_t>(type_basereg);

        // if basereg is a stack_p with no offset, we do not store anything, and no type errors
        // if we later load with that pointer, we read nothing -- load is no-op
        if (it2 && type_basereg_no_off.get_region() != crab::region::T_STACK) {
            std::string s = std::to_string(static_cast<unsigned int>(target_reg.v));
            std::string desc = std::string("\twe cannot store a pointer, r") + s + ", into packet or shared\n";
            report_type_error(desc, loc);
            return;
        }
    }
}

void type_domain_t::operator()(const Mem& b, location_t loc, int print) {
    if (is_bottom()) return;

    if (std::holds_alternative<Reg>(b.value)) {
        if (b.is_load) {
            do_load(b, std::get<Reg>(b.value), loc, print);
        } else {
            do_mem_store(b, std::get<Reg>(b.value), loc, print);
        }
    } else {
        std::string s = std::to_string(static_cast<unsigned int>(std::get<Imm>(b.value).v));
        std::string desc = std::string("\tEither loading to a number (not allowed) or storing a number (not allowed yet) - ") + s + "\n";
        report_type_error(desc, loc);
        return;
    }
}

void type_domain_t::print_initial_types() {
    auto label = label_t::entry;
    location_t loc = location_t(std::make_pair(label, 0));
    std::cout << "\n" << *m_ctx << "\n";
    std::cout << m_stack << "\n";

    std::cout << "Initial register types:\n";
    auto r1_with_loc = reg_with_loc_t(R1_ARG, loc);
    auto it = m_registers.find(r1_with_loc);
    if (it) {
        std::cout << "  ";
        print_type(R1_ARG, it.value());
        std::cout << "\n";
    }
    auto r10_with_loc = reg_with_loc_t(R10_STACK_POINTER, loc);
    auto it2 = m_registers.find(r10_with_loc);
    if (it2) {
        std::cout << "  ";
        print_type(R10_STACK_POINTER, it2.value());
        std::cout << "\n";
    }
    std::cout << "\n";
}

void type_domain_t::operator()(const basic_block_t& bb, bool check_termination, int print) {
    auto label = bb.label();
    uint32_t curr_pos = 0;
    location_t loc;
    if (print > 0) {
        if (label == label_t::entry) {
            print_initial_types();
            m_is_bottom = false;
        }
        std::cout << label << ":\n";
    }

    for (const Instruction& statement : bb) {
        loc = location_t(std::make_pair(label, ++curr_pos));
        if (print > 0) std::cout << " " << curr_pos << ".";
        std::visit([this, loc, print](const auto& v) { std::apply(*this, std::make_tuple(v, loc, print)); }, statement);
        if (print > 0 && error_location->first == loc->first && error_location->second == loc->second) std::cout << "type_error\n";
    }

    if (print > 0) {
        auto [it, et] = bb.next_blocks();
        if (it != et) {
            std::cout << "  "
            << "goto ";
            for (; it != et;) {
                std::cout << *it;
                ++it;
                if (it == et) {
                    std::cout << ";";
                } else {
                    std::cout << ",";
                }
            }
        }
        std::cout << "\n\n";
    }
}

void type_domain_t::set_require_check(check_require_func_t f) {}
