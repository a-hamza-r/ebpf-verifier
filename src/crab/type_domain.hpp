// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>

#include "crab/cfg.hpp"

namespace crab {

enum class region {
	T_CTX,
	T_STACK,
	T_PACKET,
	T_SHARED,
};

struct ptr_no_off_t {
    region r;

    ptr_no_off_t(region _r) : r(_r) {}
    ptr_no_off_t(const ptr_no_off_t& p) : r(p.r) {}

    friend bool operator==(const ptr_no_off_t& p1, const ptr_no_off_t& p2) {
        return (p1.r == p2.r);
    }

    friend bool operator!=(const ptr_no_off_t& p1, const ptr_no_off_t& p2) {
        return !(p2 == p1);
    }
};

struct ptr_with_off_t {
    region r;
    int offset;

    ptr_with_off_t(region _r, int _off) : r(_r), offset(_off) {}
    ptr_with_off_t(const ptr_with_off_t& p) : r(p.r), offset(p.offset) {}

    friend bool operator==(const ptr_with_off_t& p1, const ptr_with_off_t& p2) {
        return (p1.r == p2.r);
    }

    friend bool operator!=(const ptr_with_off_t& p1, const ptr_with_off_t& p2) {
        return !(p1 == p2);
    }
};

using ptr_t = std::variant<ptr_no_off_t, ptr_with_off_t>;


struct reg_with_loc_t {
    int r;
    std::pair<label_t, uint32_t> loc;

    reg_with_loc_t() : r(-1), loc(std::make_pair(label_t::entry, 0)) {}
    reg_with_loc_t(int _r, const label_t& l, uint32_t loc_instr) : r(_r), loc(std::make_pair(l, loc_instr)) {}

    bool operator==(const reg_with_loc_t& other) const {
        return (r != -1 && r == other.r);
    }
};
}

// adapted from top answer in
// https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
// works for now but needs to be checked again
template <>
struct std::hash<crab::reg_with_loc_t>
{
    std::size_t operator()(const crab::reg_with_loc_t& reg) const
    {
        using std::size_t;
        using std::hash;
        using std::string;

        // Compute individual hash values for first,
        // second and third and combine them using XOR
        // and bit shifting:

        return ((hash<int>()(reg.r)
               ^ (hash<int>()(reg.loc.first.from) << 1)) >> 1)
               ^ (hash<int>()(reg.loc.second) << 1);
    }
};

namespace crab {


using offset_to_ptr_no_off_t = std::unordered_map<uint64_t, ptr_no_off_t>;
using offset_to_ptr_t = std::unordered_map<uint64_t, ptr_t>;

struct ctx_t {
    offset_to_ptr_no_off_t packet_ptrs;

    ctx_t(const ebpf_context_descriptor_t* desc)
    {
        packet_ptrs.insert(std::make_pair(desc->data, crab::ptr_no_off_t(crab::region::T_PACKET)));
        packet_ptrs.insert(std::make_pair(desc->end, crab::ptr_no_off_t(crab::region::T_PACKET)));
    }
};

struct stack_t {
    offset_to_ptr_t ptrs;

    stack_t operator|(const stack_t& other) const {
        stack_t st{};
        for (auto& e : ptrs) {
            auto it = other.ptrs.find(e.first);
            if (it == other.ptrs.end() || it->second == e.second) {
                st.ptrs.insert(e);
            }
        }

        for (auto& e : other.ptrs) {
            auto it = ptrs.find(e.first);
            if (it == ptrs.end()) {
                st.ptrs.insert(e);
            }
        }
        return st;
    }

    void set_to_bottom() {
        ptrs.clear();
    }

    bool is_bottom() const {
        return ptrs.empty();
    }
};

using all_types_t = std::unordered_map<reg_with_loc_t, ptr_t>;

struct types_t {
    std::array<reg_with_loc_t, 11> vars;
    std::shared_ptr<crab::all_types_t> all_types;

    types_t() {}
    types_t(const types_t& other) : vars(other.vars), all_types(other.all_types) {}

    types_t operator|(const types_t& other) const {
        types_t v{};
        for (int i = 0; i < vars.size(); i++) {
            auto it1 = all_types->find(vars[i]);
            auto it2 = other.all_types->find(other.vars[i]);
            if (it1 != all_types->end() && it2 != other.all_types->end()) {
                if (it1->second == it2->second) {
                    v.vars[i] = vars[i];
                }
            }
        }
        v.all_types = all_types;
        return v;
    }

    void set_to_bottom() {
       vars = {};
    }

    bool is_bottom() const {
        for (auto& it : vars) {
            if (it.r != -1) return false;
        }
        return true;
    }
};

}

class type_domain_t final {

    crab::stack_t stack;
    crab::types_t types;
    std::shared_ptr<crab::ctx_t> ctx;
    label_t label;
    uint32_t m_curr_pos = 0;

  public:

  type_domain_t(const label_t& _l) : label(_l) {}
  type_domain_t(const crab::types_t& _types, const crab::stack_t& _st, const label_t& _l, std::shared_ptr<crab::ctx_t> _ctx)
            : stack(_st), types(_types), ctx(_ctx), label(_l) {}
  // eBPF initialization: R1 points to ctx, R10 to stack, etc.
  static type_domain_t setup_entry(std::shared_ptr<crab::ctx_t>, std::shared_ptr<crab::all_types_t>);
  // bottom/top
  static type_domain_t bottom();
  void set_to_top();
  void set_to_bottom();
  bool is_bottom() const;
  bool is_top() const;
  // inclusion
  bool operator<=(const type_domain_t& other) const;
  // join
  void operator|=(type_domain_t& other) const;
  void operator|=(const type_domain_t& other) const;
  type_domain_t operator|(type_domain_t& other) const;
  type_domain_t operator|(const type_domain_t& other) const&;
  // meet
  type_domain_t operator&(const type_domain_t& other) const;
  // widening
  type_domain_t widen(const type_domain_t& other) const;
  // narrowing
  type_domain_t narrow(const type_domain_t& other) const;

  //// abstract transformers
  void operator()(const Undefined &);
  void operator()(const Bin &);
  void operator()(const Un &) ;
  void operator()(const LoadMapFd &);
  void operator()(const Call &);
  void operator()(const Exit &);
  void operator()(const Jmp &);
  void operator()(const Mem &);
  void operator()(const Packet &);
  void operator()(const LockAdd &);
  void operator()(const Assume &);
  void operator()(const Assert &);
  void operator()(const basic_block_t& bb) {
      m_curr_pos = 0;
      label = bb.label();
      for (const Instruction& statement : bb) {
        m_curr_pos++;
        std::visit(*this, statement);
    }
  }

  private:

  void do_load(const Mem&, const Reg&);
  void do_mem_store(const Mem&, const Reg&);

}; // end type_domain_t
