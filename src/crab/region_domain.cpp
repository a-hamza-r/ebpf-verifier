// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/region_domain.hpp"

namespace crab {

static inline std::vector<std::set<int>> join_shared_ptr_aliases(
        const std::vector<std::set<int>>& A, const std::vector<std::set<int>>& B) {
    auto flattenSet = [](const std::vector<std::set<int>>& v) {
        std::set<int> s;
        for (const auto& x : v) {
            s.insert(x.begin(), x.end());
        }
        return s;
    };

    std::set<int> a, b, intersect;
    a = flattenSet(A);
    b = flattenSet(B);

    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
            std::inserter(intersect, intersect.begin()));

    std::vector<std::set<int>> powerset;
    powerset.push_back({});
    for (int n : intersect) {
        auto size = (size_t)powerset.size();
        for (size_t i = 0; i < size; i++) {
            auto newSet = powerset[(int)i];
            newSet.insert(n);
            powerset.push_back(newSet);
        }
    }

    std::vector<std::set<int>> result;
    for (const auto& s : powerset) {
        auto foundInA = std::find(A.begin(), A.end(), s);
        auto foundInB = std::find(B.begin(), B.end(), s);
        if (foundInA != A.end() || foundInB != B.end()) {
            result.push_back(s);
        }
        auto flattened = flattenSet(result);
        if (flattened.size() == intersect.size()) {
            break;
        }
    }
    return result;
}

region_ctx_t::region_ctx_t(const ebpf_context_descriptor_t* desc) {
    if (desc == nullptr) return;
    if (desc->data >= 0) {
        m_keys.push_back(desc->data);
    }
    if (desc->end >= 0) {
        m_keys.push_back(desc->end);
    }
    if (desc->meta >= 0) {
        m_keys.push_back(desc->meta);
    }
    m_size = std::max(0, desc->size);
}

bool region_ctx_t::packet_ptr_at(uint64_t key) const {
    return std::find(m_keys.begin(), m_keys.end(), key) != m_keys.end();
}

void region_registers_t::scratch_caller_saved_registers() {
    for (uint8_t r = R1_ARG; r <= R5_ARG; r++) {
        operator-=(register_t{r});
    }
}

void region_registers_t::forget_packet_ptrs() {
    // skip the R12_PKT_BEGIN, as its region type will remain packet_ptr_t
    for (uint8_t r = R0_RETURN_VALUE; r < NUM_REGISTERS-2; r++) {
        if (is_packet_ptr(find(register_t{r}))) {
            operator-=(register_t{r});
        }
    }
}

region_registers_t region_registers_t::operator|(const region_registers_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    region_registers_t joined_reg_types;

    // a hack to store region information at the start of a joined basic block
    // in join, we do not know the label of the bb, hence we store the information
    // at a bb that is not used anywhere else in the program, and later when we know
    // the bb label, we can fix
    location_t loc = location_t::top();

    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (m_cur_register_def[i] == nullptr || other.m_cur_register_def[i] == nullptr) continue;
        auto maybe_ptr1 = find(register_t{i});
        auto maybe_ptr2 = other.find(register_t{i});
        if (maybe_ptr1 && maybe_ptr2) {
            ptr_or_mapfd_t ptr_or_mapfd1 = *maybe_ptr1, ptr_or_mapfd2 = *maybe_ptr2;
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd1)
                        && std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd2)) {
                ptr_with_off_t ptr_with_off1 = std::get<ptr_with_off_t>(ptr_or_mapfd1);
                ptr_with_off_t ptr_with_off2 = std::get<ptr_with_off_t>(ptr_or_mapfd2);
                if (ptr_with_off1.get_region() == ptr_with_off2.get_region()) {
                    auto joined_ptr = ptr_with_off1 | ptr_with_off2;
                    joined_reg_types.insert(register_t{i}, loc, joined_ptr);
                }
            }
            else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1)
                    && std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
                mapfd_t mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
                mapfd_t mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
                joined_reg_types.insert(register_t{i}, loc, mapfd1 | mapfd2);
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1)
                    && std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
                joined_reg_types.insert(register_t{i}, loc, packet_ptr_t());
            }
        }
    }
    return joined_reg_types;
}

bool region_registers_t::operator<=(const region_registers_t& other) const {
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (other.m_cur_register_def[i] == nullptr) continue;
        if (m_cur_register_def[i] == nullptr) return false;
        auto maybe_ptr1 = find(register_t{i});
        auto maybe_ptr2 = other.find(register_t{i});
        if (maybe_ptr1 && maybe_ptr2) {
            ptr_or_mapfd_t ptr_or_mapfd1 = *maybe_ptr1, ptr_or_mapfd2 = *maybe_ptr2;
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd1)
                    && std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd2)) {
                ptr_with_off_t ptr_with_off1 = std::get<ptr_with_off_t>(ptr_or_mapfd1);
                ptr_with_off_t ptr_with_off2 = std::get<ptr_with_off_t>(ptr_or_mapfd2);
                if (!(ptr_with_off1 <= ptr_with_off2)) return false;
            }
            else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1)
                    && std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
                mapfd_t mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
                mapfd_t mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
                if (!(mapfd1 <= mapfd2)) return false;
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1)
                    && std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
                continue;
            }
            else return false;
        }
    }
    return true;
}

region_registers_t region_registers_t::widen(const region_registers_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    region_registers_t joined_reg_types;

    location_t loc = location_t::top();
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (other.m_cur_register_def[i] == nullptr) continue;
        auto maybe_ptr1 = find(register_t{i});
        auto maybe_ptr2 = other.find(register_t{i});
        if (maybe_ptr1 && maybe_ptr2) {
            ptr_or_mapfd_t ptr_or_mapfd1 = *maybe_ptr1, ptr_or_mapfd2 = *maybe_ptr2;
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd1)
                        && std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd2)) {
                ptr_with_off_t ptr_with_off1 = std::get<ptr_with_off_t>(ptr_or_mapfd1);
                ptr_with_off_t ptr_with_off2 = std::get<ptr_with_off_t>(ptr_or_mapfd2);
                if (ptr_with_off1.get_region() == ptr_with_off2.get_region()) {
                    auto ptr_with_off = ptr_with_off1.widen(ptr_with_off2);
                    joined_reg_types.insert(register_t{i}, loc, ptr_with_off);
                }
            }
            else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1)
                    && std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
                mapfd_t mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
                mapfd_t mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
                auto map_fd = mapfd1.widen(mapfd2);
                joined_reg_types.insert(register_t{i}, loc, map_fd);
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1)
                    && std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
                joined_reg_types.insert(register_t{i}, loc, packet_ptr_t());
            }
        }
    }
    return joined_reg_types;
}

void region_registers_t::operator-=(register_t var) {
    if (is_bottom()) {
        return;
    }
    m_cur_register_def[var] = nullptr;
}

void region_registers_t::set_to_bottom() {
    m_is_bottom = true;
}

void region_registers_t::set_to_top() {
    m_registers_env = std::make_shared<global_env_region_registers_t>();
    m_cur_register_def = live_registers_t{nullptr};
    m_is_bottom = false;
}

bool region_registers_t::is_bottom() const { return m_is_bottom; }

bool region_registers_t::is_top() const {
    if (m_is_bottom) { return false; }
    if (m_registers_env == nullptr) return true;
    for (auto &it : m_cur_register_def) {
        if (it != nullptr) return false;
    }
    return true;
}

void region_registers_t::insert(register_t reg, const location_t& loc, const ptr_or_mapfd_t& type) {
    register_location_t register_location = register_location_t{reg, loc};
    m_registers_env->insert_or_assign(register_location, type);
    m_cur_register_def[reg] = std::make_shared<register_location_t>(register_location);
}

std::optional<ptr_or_mapfd_t> region_registers_t::find(register_location_t reg) const {
    auto it = m_registers_env->find(reg);
    if (it == m_registers_env->end()) return {};
    return it->second;
}

std::optional<ptr_or_mapfd_t> region_registers_t::find(register_t key) const {
    if (m_cur_register_def[key] == nullptr) return {};
    return find(*(m_cur_register_def[key]));
}

void region_registers_t::adjust_bb_for_registers(location_t loc) {
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (auto it = find(register_t{i})) {
            insert(register_t{i}, loc, *it);
        }
    }
}

bool region_stack_t::operator<=(const region_stack_t& other) const {
    size_t size1 = m_cells.size();
    size_t size2 = other.m_cells.size();
    if (size2 > size1) return false;
    for (auto const &kv : other.m_cells) {
        auto it = m_cells.find(kv.first);
        if (it == m_cells.end()) return false;
        auto ptr_or_mapfd_cells1 = it->second;
        auto ptr_or_mapfd_cells2 = kv.second;
        auto ptr_or_mapfd1 = ptr_or_mapfd_cells1.first;
        auto ptr_or_mapfd2 = ptr_or_mapfd_cells2.first;
        int width1 = ptr_or_mapfd_cells1.second;
        int width2 = ptr_or_mapfd_cells2.second;
        if (width1 != width2) return false;
        if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd1) &&
                std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd2)) {
            auto ptr_with_off1 = std::get<ptr_with_off_t>(ptr_or_mapfd1);
            auto ptr_with_off2 = std::get<ptr_with_off_t>(ptr_or_mapfd2);
            if (!(ptr_with_off1 <= ptr_with_off2)) return false;
        }
        else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1) &&
                std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
            auto mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
            auto mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
            if (!(mapfd1 <= mapfd2)) return false;
        }
        else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1) &&
                std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
            continue;
        }
        else return false;
    }
    return true;
}

region_stack_t region_stack_t::operator|(const region_stack_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    region_stack_t joined_stack;
    for (auto const&kv: m_cells) {
        auto maybe_ptr_or_mapfd_stack_cell = other.find(kv.first);
        if (maybe_ptr_or_mapfd_stack_cell) {
            auto ptr_or_mapfd_stack_cell1 = kv.second;
            auto ptr_or_mapfd_stack_cell2 = *maybe_ptr_or_mapfd_stack_cell;
            auto ptr_or_mapfd1 = ptr_or_mapfd_stack_cell1.first;
            auto ptr_or_mapfd2 = ptr_or_mapfd_stack_cell2.first;
            int width1 = ptr_or_mapfd_stack_cell1.second;
            int width2 = ptr_or_mapfd_stack_cell2.second;
            // this should be fixed in the future
            int width_joined = std::min(width1, width2);
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd2)) {
                auto ptr_with_off1 = std::get<ptr_with_off_t>(ptr_or_mapfd1);
                auto ptr_with_off2 = std::get<ptr_with_off_t>(ptr_or_mapfd2);
                if (ptr_with_off1.get_region() == ptr_with_off2.get_region()) {
                    joined_stack.store(kv.first, ptr_with_off1 | ptr_with_off2, width_joined);
                }
            }
            else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
                auto mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
                auto mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
                joined_stack.store(kv.first, mapfd1 | mapfd2, width_joined);
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
                joined_stack.store(kv.first, packet_ptr_t(), width_joined);
            }
        }
    }
    return joined_stack;
}

region_stack_t region_stack_t::widen(const region_stack_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    region_stack_t joined_stack;
    for (auto const&kv: m_cells) {
        auto maybe_ptr_or_mapfd_cells = other.find(kv.first);
        if (maybe_ptr_or_mapfd_cells) {
            auto ptr_or_mapfd_cells1 = kv.second;
            auto ptr_or_mapfd_cells2 = *maybe_ptr_or_mapfd_cells;
            auto ptr_or_mapfd1 = ptr_or_mapfd_cells1.first;
            auto ptr_or_mapfd2 = ptr_or_mapfd_cells2.first;
            int width1 = ptr_or_mapfd_cells1.second;
            int width2 = ptr_or_mapfd_cells2.second;
            // this should be fixed in the future
            int width_joined = std::min(width1, width2);
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd2)) {
                auto ptr_with_off1 = std::get<ptr_with_off_t>(ptr_or_mapfd1);
                auto ptr_with_off2 = std::get<ptr_with_off_t>(ptr_or_mapfd2);
                if (ptr_with_off1.get_region() == ptr_with_off2.get_region()) {
                    auto ptr_with_off = ptr_with_off1.widen(ptr_with_off2);
                    joined_stack.store(kv.first, ptr_with_off, width_joined);
                }
            }
            else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
                auto mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
                auto mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
                auto mapfd = mapfd1.widen(mapfd2);
                joined_stack.store(kv.first, mapfd, width_joined);
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
                joined_stack.store(kv.first, packet_ptr_t(), width_joined);
            }
        }
    }
    return joined_stack;
}

void region_stack_t::operator-=(uint64_t key) {
    auto it = find(key);
    if (it)
        m_cells.erase(key);
}

void region_stack_t::operator-=(const std::vector<uint64_t>& keys) {
    for (auto &key : keys) {
       *this -= key;
    }
}

void region_stack_t::set_to_bottom() {
    m_cells.clear();
    m_is_bottom = true;
}

void region_stack_t::set_to_top() {
    m_cells.clear();
    m_is_bottom = false;
}

region_stack_t region_stack_t::bottom() {
    region_stack_t stk;
    stk.set_to_bottom();
    return stk;
};

region_stack_t region_stack_t::top() { return region_stack_t(); }

bool region_stack_t::is_bottom() const { return m_is_bottom; }

bool region_stack_t::is_top() const {
    if (m_is_bottom)
        return false;
    return m_cells.empty();
}

void region_stack_t::store(uint64_t key, ptr_or_mapfd_t value, int width) {
    m_cells.insert_or_assign(key, std::make_pair(value, width));
}

size_t region_stack_t::size() const {
    return m_cells.size();
}

std::vector<uint64_t> region_stack_t::get_keys() const {
    std::vector<uint64_t> keys;
    keys.reserve(size());

    for (auto const&kv : m_cells) {
        keys.push_back(kv.first);
    }
    return keys;
}


std::optional<ptr_or_mapfd_stack_cell_t> region_stack_t::find(uint64_t key) const {
    auto it = m_cells.find(key);
    if (it == m_cells.end()) return {};
    return it->second;
}

std::vector<uint64_t> region_stack_t::find_overlapping_cells(uint64_t start, int width) const {
    std::vector<uint64_t> overlapping_cells;
    auto it = m_cells.begin();
    while (it != m_cells.end() && it->first < start) {
        it++;
    }
    if (it != m_cells.begin()) {
        it--;
        auto key = it->first;
        auto width_key = it->second.second;
        if (key < start && key+width_key > start) overlapping_cells.push_back(key);
    }

    for (; it != m_cells.end(); it++) {
        auto key = it->first;
        if (key >= start && key < start+width) overlapping_cells.push_back(key);
        if (key >= start+width) break;
    }
    return overlapping_cells;
}

size_t region_domain_t::get_ctx_size() const {
    return m_ctx->get_size();
}

std::optional<ptr_or_mapfd_t> region_domain_t::find_ptr_or_mapfd_type(register_t reg) const {
    return m_registers.find(reg);
}

void region_domain_t::insert_in_registers(register_t reg, location_t loc,
        const ptr_or_mapfd_t& ptr) {
    m_registers.insert(reg, loc, ptr);
}

void region_domain_t::store_in_stack(uint64_t key, ptr_or_mapfd_t value, int width) {
    m_stack.store(key, value, width);
}

bool region_domain_t::is_bottom() const {
    if (m_is_bottom) return true;
    return (m_stack.is_bottom() || m_registers.is_bottom());
}

bool region_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_stack.is_top() && m_registers.is_top());
}

region_domain_t region_domain_t::bottom() {
    region_domain_t typ;
    typ.set_to_bottom();
    return typ;
}

void region_domain_t::set_to_bottom() {
    m_is_bottom = true;
    m_stack.set_to_bottom();
    m_registers.set_to_bottom();
}

void region_domain_t::set_to_top() {
    m_is_bottom = false;
    m_stack.set_to_top();
    m_registers.set_to_top();
}

std::optional<ptr_or_mapfd_t> region_domain_t::find_ptr_or_mapfd_at_loc(const register_location_t& reg) const {
    return m_registers.find(reg);
}

void region_domain_t::set_registers_to_top() {
    m_registers.set_to_top();
}

const std::vector<uint64_t>& region_domain_t::get_ctx_keys() const {
    return m_ctx->get_keys();
}

std::vector<uint64_t> region_domain_t::get_stack_keys() const {
    return m_stack.get_keys();
}

std::optional<ptr_or_mapfd_stack_cell_t> region_domain_t::find_in_stack(uint64_t key) const {
    return m_stack.find(key);
}

bool region_domain_t::operator<=(const region_domain_t& abs) const {
    // TODO: check for shared_ptr_aliases
    return (m_registers <= abs.m_registers && m_stack <= abs.m_stack);
}

void region_domain_t::operator|=(const region_domain_t& abs) {
    region_domain_t tmp{abs};
    operator|=(std::move(tmp));
}

void region_domain_t::operator|=(region_domain_t&& abs) {
    if (is_bottom()) {
        *this = abs;
        return;
    }
    *this = *this | std::move(abs);
}

region_domain_t region_domain_t::operator|(const region_domain_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    auto aliases = join_shared_ptr_aliases(m_shared_ptr_aliases, other.m_shared_ptr_aliases);
    return region_domain_t(m_registers | other.m_registers, m_stack | other.m_stack, m_ctx,
            std::move(aliases));
}

region_domain_t region_domain_t::operator|(region_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    auto aliases = join_shared_ptr_aliases(m_shared_ptr_aliases, other.m_shared_ptr_aliases);
    return region_domain_t(m_registers | std::move(other.m_registers),
            m_stack | std::move(other.m_stack), std::move(m_ctx), std::move(aliases));
}

region_domain_t region_domain_t::operator&(const region_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

region_domain_t region_domain_t::widen(const region_domain_t& other, bool to_constants) {
    if (is_bottom() || other.is_top()) {
        return other;
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    auto aliases = join_shared_ptr_aliases(m_shared_ptr_aliases, other.m_shared_ptr_aliases);
    return region_domain_t(m_registers.widen(other.m_registers), m_stack.widen(other.m_stack),
            other.m_ctx, std::move(aliases));
}

region_domain_t region_domain_t::narrow(const region_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

crab::bound_t region_domain_t::get_loop_count_upper_bound() const {
    // WARNING: Not implemented yet.
    return crab::bound_t{crab::number_t{0}};
}

void region_domain_t::initialize_loop_counter(const label_t& label) {
    // WARNING: Not implemented yet.
}

string_invariant region_domain_t::to_set() {
    return string_invariant{};
}

void region_domain_t::operator()(const Undefined &u, location_t loc) {
    // nothing to do here
}

void region_domain_t::operator()(const Exit &u, location_t loc) {}

void region_domain_t::operator()(const Jmp &u, location_t loc) {}


void region_domain_t::assume_cst(Condition::Op op, ptr_with_off_t&& shared_ptr, int64_t imm,
        register_t left, location_t loc) {
    // we only reach here when the ptr is shared ptr
    auto nullness = shared_ptr.get_nullness();
    auto set_nullness = [this, loc, left](nullness_t n, int id) {
        if (id == -1) {
            for (size_t i = 0; i < m_shared_ptr_aliases.size(); i++) {
                if (m_shared_ptr_aliases[i].count(left)) {
                    id = i;
                    break;
                }
            }
        }
        if (id == -1) return;
        for (const auto& s : m_shared_ptr_aliases[id]) {
            if (s <= 10) {
                auto type = m_registers.find(register_t{(uint8_t)s});
                if (is_shared_ptr(type)) {
                    auto shared_ptr = std::get<ptr_with_off_t>(*type);
                    shared_ptr.set_nullness(n);
                    m_registers.insert(register_t{(uint8_t)s}, loc, shared_ptr);
                }
            }
            else {
                auto offset = s - 11;
                auto type_with_width = m_stack.find(offset);
                if (!type_with_width) continue;
                auto type = type_with_width->first;
                if (is_shared_ptr(type)) {
                    auto shared_ptr = std::get<ptr_with_off_t>(type);
                    shared_ptr.set_nullness(n);
                    m_stack.store(offset, shared_ptr, type_with_width->second);
                }
            }
        }
    };
    if (imm == 0) {
        if (op == Condition::Op::EQ) {
            if (nullness == nullness_t::_NULL) {
                //m_registers.set_to_top();
            }
            else if (nullness == nullness_t::NOT_NULL) {
                m_registers.set_to_bottom();
            }
            else {
                auto id = shared_ptr.get_id();
                set_nullness(nullness_t::_NULL, id);
            }
        }
        else if (op == Condition::Op::NE) {
            if (nullness == nullness_t::NOT_NULL) {
                //m_registers.set_to_top();
            }
            else if (nullness == nullness_t::_NULL) {
                m_registers.set_to_bottom();
            }
            else {
                auto id = shared_ptr.get_id();
                set_nullness(nullness_t::NOT_NULL, id);
            }
        }
    }
}

void region_domain_t::operator()(const Assume& u, location_t loc) {
    // nothing to do here
}

void region_domain_t::operator()(const Assert& u, location_t loc) {
    // nothing to do here
}

void region_domain_t::operator()(const ZeroCtxOffset& u, location_t loc) {
    auto maybe_ptr_or_mapfd = m_registers.find(u.reg.v);
    if (is_ctx_ptr(maybe_ptr_or_mapfd)) {
        auto ctx_ptr = std::get<ptr_with_off_t>(*maybe_ptr_or_mapfd);
        if (ctx_ptr.get_offset() == interval_t{number_t{0}}) return;
    }
    m_errors.push_back(loc.to_string() + ": Non-zero context offset");
}

void region_domain_t::operator()(const basic_block_t& bb) {
    // nothing to do here
}

void region_domain_t::operator()(const Un& u, location_t loc) {
    m_registers -= register_t{u.dst.v};
}

// Get the start and end of the range of possible map fd values.
// In the future, it would be cleaner to use a set rather than an interval
// for map fds.
bool region_domain_t::get_map_fd_range(const Reg& map_fd_reg, int32_t* start_fd, int32_t* end_fd) const {
    auto maybe_type = m_registers.find(map_fd_reg.v);
    if (!is_mapfd_type(maybe_type)) return false;
    auto mapfd_type = std::get<mapfd_t>(*maybe_type);
    const auto& mapfd_interval = mapfd_type.get_mapfd();
    auto lb = mapfd_interval.lb().number();
    auto ub = mapfd_interval.ub().number();
    if (!lb || !lb->fits<int32_t>() || !ub || !ub->fits<int32_t>())
        return false;
    *start_fd = lb.value().cast_to<int32_t>();
    *end_fd = ub.value().cast_to<int32_t>();

    // Cap the maximum range we'll check.
    const int max_range = 32;
    return (*mapfd_interval.finite_size() < max_range);
}

// All maps in the range must have the same type for us to use it.
std::optional<uint32_t> region_domain_t::get_map_type(const Reg& map_fd_reg) const {
    int32_t start_fd, end_fd;
    if (!get_map_fd_range(map_fd_reg, &start_fd, &end_fd))
        return std::optional<uint32_t>();

    std::optional<uint32_t> type;
    for (int32_t map_fd = start_fd; map_fd <= end_fd; map_fd++) {
        EbpfMapDescriptor* map = &global_program_info->platform->get_map_descriptor(map_fd);
        if (map == nullptr)
            return std::optional<uint32_t>();
        if (!type.has_value())
            type = map->type;
        else if (map->type != *type)
            return std::optional<uint32_t>();
    }
    return type;
}

// All maps in the range must have the same inner map fd for us to use it.
std::optional<uint32_t> region_domain_t::get_map_inner_map_fd(const Reg& map_fd_reg) const {
    int32_t start_fd, end_fd;
    if (!get_map_fd_range(map_fd_reg, &start_fd, &end_fd))
        return std::optional<uint32_t>();

    std::optional<uint32_t> inner_map_fd;
    for (int map_fd = start_fd; map_fd <= end_fd; map_fd++) {
        EbpfMapDescriptor* map = &global_program_info->platform->get_map_descriptor(map_fd);
        if (map == nullptr)
            return std::optional<uint32_t>();
        if (!inner_map_fd.has_value())
            inner_map_fd = map->inner_map_fd;
        else if (map->type != *inner_map_fd)
            return std::optional<uint32_t>();
    }
    return inner_map_fd;
}

// We can deal with a range of key sizes.
interval_t region_domain_t::get_map_key_size(const Reg& map_fd_reg) const {
    int start_fd, end_fd;
    if (!get_map_fd_range(map_fd_reg, &start_fd, &end_fd))
        return interval_t::top();

    interval_t result = interval_t::bottom();
    for (int map_fd = start_fd; map_fd <= end_fd; map_fd++) {
        if (EbpfMapDescriptor* map = &global_program_info->platform->get_map_descriptor(map_fd)) {
            result = result | interval_t(map->key_size);
        } else {
            return interval_t::top();
        }
    }
    return result;
}

// We can deal with a range of value sizes.
interval_t region_domain_t::get_map_value_size(const Reg& map_fd_reg) const {
    int start_fd, end_fd;
    if (!get_map_fd_range(map_fd_reg, &start_fd, &end_fd))
        return interval_t::top();

    interval_t result = crab::interval_t::bottom();
    for (int map_fd = start_fd; map_fd <= end_fd; map_fd++) {
        if (EbpfMapDescriptor* map = &global_program_info->platform->get_map_descriptor(map_fd)) {
            result = result | interval_t(map->value_size);
        } else {
            return interval_t::top();
        }
    }
    return result;
}

// We can deal with a range of max_entries values.
interval_t region_domain_t::get_map_max_entries(const Reg& map_fd_reg) const {
    int start_fd, end_fd;
    if (!get_map_fd_range(map_fd_reg, &start_fd, &end_fd))
        return interval_t::top();

    interval_t result = interval_t::bottom();
    for (int map_fd = start_fd; map_fd <= end_fd; map_fd++) {
        if (EbpfMapDescriptor* map = &global_program_info->platform->get_map_descriptor(map_fd)) {
            result = result | interval_t(map->max_entries);
        }
        else {
            return interval_t::top();
        }
    }
    return result;
}

void region_domain_t::do_load_mapfd(register_t dst_reg, int mapfd, location_t loc) {
    const auto& platform = global_program_info->platform;
    const EbpfMapDescriptor& desc = platform->get_map_descriptor(mapfd);
    const EbpfMapValueType& map_value_type = platform->get_map_type(desc.type).value_type;
    auto mapfd_interval = interval_t{number_t{mapfd}};
    mapfd_t type{mapfd_interval, map_value_type};
    m_registers.insert(dst_reg, loc, type);
}

void region_domain_t::operator()(const LoadMapFd &u, location_t loc) {
    do_load_mapfd(register_t{u.dst.v}, u.mapfd, loc);
}

static EbpfRelocationDescriptor* find_relocation_descriptor(const int relocation_fd) {
    for (EbpfRelocationDescriptor& relocation : global_program_info->relocation_descriptors) {
        if (relocation.original_fd == relocation_fd) {
            return &relocation;
        }
    }
    return nullptr;
}


void region_domain_t::operator()(const LoadVariable& u, location_t loc) {
    const EbpfRelocationDescriptor* desc = find_relocation_descriptor(u.varfd);
    if (desc == nullptr) {
        throw std::runtime_error(std::string("relocation_fd not found"));
        m_registers -= register_t{u.dst.v};
        return;
    } else {
        auto type = ptr_with_off_t::shared_region_ptr(0, -1, nullness_t::MAYBE_NULL,
                              interval_t{number_t{desc->value_size}});
        m_registers.insert(u.dst.v, loc, type);
    }
}

void region_domain_t::set_aliases(int v, ptr_with_off_t& ptr) {
    size_t i = 0;
    for (; i < m_shared_ptr_aliases.size(); i++) {
        if (m_shared_ptr_aliases[(int)i].count(v) > 0) {
            break;
        }
    }
    if (i < m_shared_ptr_aliases.size()) m_shared_ptr_aliases[(int)i].erase(v);
    auto id = ptr.get_id();
    /* NOTE: this check is supposed to be for newly generated pointers, but it could also be the
     * case that an existing pointer has id == -1. This is because, at join of two ptr_with_off_t,
     * we set id == -1 for simplicity. Some code for correcting this is in the assume_cst
     * function, although it should work here as well but doesn't (ideally, may be not).
     * In future, check this again.
     */
    if (id == -1) {
        m_shared_ptr_aliases.push_back({v});
        ptr.set_id(m_shared_ptr_aliases.size() - 1);
    }
    else {
        m_shared_ptr_aliases[id].insert(v);
        ptr.set_id(id);
    }
}

void region_domain_t::do_call(const Call& u, const stack_cells_t& cells, location_t loc) {
    for (const auto& [offset, width] : cells) {
        m_stack -= m_stack.find_overlapping_cells(offset, width);
    }
    std::optional<Reg> maybe_fd_reg{};
    for (ArgSingle param : u.singles) {
        if (param.kind == ArgSingle::Kind::MAP_FD) maybe_fd_reg = param.reg;
        break;
    }
    register_t r0{R0_RETURN_VALUE};
    if (u.is_map_lookup) {
        if (maybe_fd_reg) {
            if (auto map_type = get_map_type(*maybe_fd_reg)) {
                if (global_program_info->platform->get_map_type(*map_type).value_type
                        == EbpfMapValueType::MAP) {
                    if (auto inner_map_fd = get_map_inner_map_fd(*maybe_fd_reg)) {
                        do_load_mapfd(r0, to_signed(*inner_map_fd), loc);
                        goto out;
                    }
                } else {
                    auto type = ptr_with_off_t::shared_region_ptr(0, -1, nullness_t::NOT_NULL,
                                               get_map_value_size(*maybe_fd_reg));
                    set_aliases((int)r0, type);
                    m_registers.insert(r0, loc, type);
                }
            }
        }
        else {
            // here, we only know that the return value is a pointer to a shared region
            auto type = ptr_with_off_t::shared_region_ptr(0, -1, nullness_t::MAYBE_NULL);
            set_aliases((int)r0, type);
            m_registers.insert(r0, loc, type);
        }
    }
    else {
        m_registers -= r0;
    }
out:
    m_registers.scratch_caller_saved_registers();
    if (u.reallocate_packet) {
        m_registers.forget_packet_ptrs();
    }
}

void region_domain_t::operator()(const Call& u, location_t loc) {
    // nothing to do here
}

void region_domain_t::operator()(const IncrementLoopCounter &u, location_t loc) {
    // WARNING: Not implemented yet.
}

void region_domain_t::operator()(const Atomic &u, location_t loc) {
    // WARNING: Not implemented yet.
}

void region_domain_t::operator()(const Packet& u, location_t loc) {
    m_registers -= register_t{R0_RETURN_VALUE};
    m_registers.scratch_caller_saved_registers();
}

void region_domain_t::check_valid_access(const ValidAccess &s, int width, location_t loc) {
    bool is_comparison_check = s.width == (Value)Imm{0};
    std::string loc_str = loc.to_string();
    // we reach here only if the register is a pointer or mapfd
    auto reg_ptr_or_mapfd_type = *(m_registers.find(s.reg.v));
    if (std::holds_alternative<ptr_with_off_t>(reg_ptr_or_mapfd_type)) {
        auto ptr_with_off_type = std::get<ptr_with_off_t>(reg_ptr_or_mapfd_type);
        auto offset = ptr_with_off_type.get_offset();
        auto offset_to_check = offset+interval_t{s.offset};
        auto offset_lb = offset_to_check.lb();
        auto offset_plus_width_ub = offset_to_check.ub()+bound_t{width};
        if (ptr_with_off_type.get_region() == region_t::R_STACK) {
            if (!(bound_t{STACK_BEGIN} <= offset_lb)) {
                m_errors.push_back(loc_str + ": Lower bound must be at least " +
                                   std::to_string(STACK_BEGIN));
                return;
            }
            if (!(offset_plus_width_ub <= bound_t{EBPF_STACK_SIZE})) {
                m_errors.push_back(loc_str + ": Upper bound must be at most " +
                                   std::to_string(EBPF_STACK_SIZE));
                return;
            }
        }
        else if (ptr_with_off_type.get_region() == region_t::R_CTX) {
            if (!(bound_t{CTX_BEGIN} <= offset_lb)) {
                m_errors.push_back(loc_str + ": Lower bound must be at least " +
                                   std::to_string(CTX_BEGIN));
                return;
            }
            if (!(offset_plus_width_ub <= bound_t{get_ctx_size()})) {
                m_errors.push_back(loc_str + ": Upper bound must be at most " +
                                   std::to_string(get_ctx_size()));
                return;
            }
        }
        else { // shared
            if (!(bound_t{SHARED_BEGIN} <= offset_lb)) {
                m_errors.push_back(loc_str + ": Lower bound must be at least " +
                                   std::to_string(SHARED_BEGIN));
                return;
            }
            if (!(offset_plus_width_ub <= ptr_with_off_type.get_region_size().lb())) {
                auto number_lb = ptr_with_off_type.get_region_size().lb().number();
                int lb_value = number_lb->cast_to<int>();
                m_errors.push_back(loc_str + ": Upper bound must be at most " +
                                   std::to_string(lb_value));
                return;
            }
            if (!is_comparison_check && !s.or_null) {
                auto nullness = ptr_with_off_type.get_nullness();
                if (nullness == nullness_t::NOT_NULL) {
                    m_errors.push_back(loc_str + ": Possible null access");
                }
            }
        }
    }
    else if (std::holds_alternative<packet_ptr_t>(reg_ptr_or_mapfd_type)) {
        // We do not handle packet ptr access in region domain
        return;
    }
    else {
        // mapfd
        if (!is_comparison_check) {
            m_errors.push_back("FDs cannot be dereferenced directly");
        }
    }
}

void region_domain_t::operator()(const ValidAccess &s, location_t loc) {
    // nothing to do here
}

region_domain_t region_domain_t::setup_entry(bool init_r1) {
    location_t loc{label_t::entry, 0};
    region_registers_t typ;
    if (init_r1) {
        const auto reg_r1 = register_t{R1_ARG};
        const auto ctx_begin_ptr = ptr_with_off_t::ctx_region_ptr(CTX_BEGIN);
        typ.insert(reg_r1, loc, ctx_begin_ptr);
    }
    const auto reg_r10 = register_t{R10_STACK_POINTER};
    const auto stack_end_ptr = ptr_with_off_t::stack_region_ptr(STACK_END);
    typ.insert(reg_r10, loc, stack_end_ptr);

    // Initialize R12 to point to the pkt_begin pointer, which is mainly needed in offset domain;
    // however, for consistency, the type is stored in the region domain.
    const auto reg_r12 = register_t{R12_PKT_BEGIN};
    const auto packet_begin_ptr = packet_ptr_t();
    typ.insert(reg_r12, loc, packet_begin_ptr);

    return region_domain_t{std::move(typ), region_stack_t::top(),
        std::make_shared<region_ctx_t>(global_program_info->type.context_descriptor)};
}

void region_domain_t::operator()(const TypeConstraint& s, location_t loc) {
    // nothing to do here
}

void region_domain_t::check_type(const TypeConstraint& s, bool is_numeric, location_t loc) {
    auto ptr_or_mapfd_opt = m_registers.find(s.reg.v);
    if (ptr_or_mapfd_opt) {
        // it is a pointer or mapfd
        auto ptr_or_mapfd_type = ptr_or_mapfd_opt.value();
        if (std::holds_alternative<mapfd_t>(ptr_or_mapfd_type)) {
            auto map_fd = std::get<mapfd_t>(ptr_or_mapfd_type);
            if (map_fd.has_type_map_programs()) {
                if (s.types == TypeGroup::map_fd_programs) return;
            } else {
                if (s.types == TypeGroup::map_fd) return;
            }
        }
        else {
            if (s.types == TypeGroup::pointer || s.types == TypeGroup::ptr_or_num) return;
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd_type)) {
                ptr_with_off_t ptr_with_off = std::get<ptr_with_off_t>(ptr_or_mapfd_type);
                if (ptr_with_off.get_region() == crab::region_t::R_CTX) {
                    if (s.types == TypeGroup::singleton_ptr) return;
                    if (s.types == TypeGroup::ctx) return;
                }
                else {
                    if (s.types == TypeGroup::mem || s.types == TypeGroup::mem_or_num) return;
                    if (ptr_with_off.get_region() == crab::region_t::R_SHARED) {
                        if (s.types == TypeGroup::shared) return;
                    }
                    else {
                        if (s.types == TypeGroup::singleton_ptr) return;
                        if (s.types == TypeGroup::stack || s.types == TypeGroup::stack_or_packet)
                            return;
                    }
                }
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd_type)) {
                if (s.types == TypeGroup::singleton_ptr) return;
                if (s.types == TypeGroup::mem || s.types == TypeGroup::mem_or_num) return;
                if (s.types == TypeGroup::packet || s.types == TypeGroup::stack_or_packet) return;
            }
        }
    }
    else if (is_numeric) {
        if (s.types == TypeGroup::number || s.types == TypeGroup::ptr_or_num
                || s.types == TypeGroup::mem_or_num)
            return;
    }
    m_errors.push_back(loc.to_string() + ": Invalid type");
}

void region_domain_t::update_ptr_or_mapfd(const ptr_or_mapfd_t& ptr_or_mapfd, const interval_t& change,
        location_t loc, register_t reg) {
    if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd)) {
        auto ptr_or_mapfd_with_off = std::get<ptr_with_off_t>(ptr_or_mapfd);
        auto updated_offset = ptr_or_mapfd_with_off.get_offset() + change;
        ptr_or_mapfd_with_off.set_offset(updated_offset);
        m_registers.insert(reg, loc, ptr_or_mapfd_with_off);
    }
    else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd)) {
        m_registers.insert(reg, loc, ptr_or_mapfd);
    }
    else {
        m_errors.push_back(loc.to_string() + ": Cannot update mapfd type through arithmetic");
        m_registers -= reg;
    }
}

void region_domain_t::operator()(const Bin& b, location_t loc) {
    // nothing to do here
}

void region_domain_t::do_bin(const Bin& bin,
                             std::optional<interval_t> dst_signed_interval_opt,
                             std::optional<interval_t> src_signed_interval_opt,
                             location_t loc) {

    auto dst_register = register_t{bin.dst.v};
    auto dst_ptr_or_mapfd_opt = m_registers.find(dst_register);
    bool is_numeric_dst = dst_signed_interval_opt.has_value();
    bool is_numeric_src = std::holds_alternative<Imm>(bin.v) || src_signed_interval_opt.has_value();

    if (is_numeric_dst && is_numeric_src) {
        m_registers -= dst_register;
        return;
    }

    using Op = Bin::Op;

    if (auto pimm = std::get_if<Imm>(&bin.v)) {
        int64_t imm;
        if (bin.is64) {
            // Use the full signed value.
            imm = to_signed(pimm->v);
        } else {
            // Use only the low 32 bits of the value.
            imm = gsl::narrow_cast<int64_t>(pimm->v);
        }
        auto imm_interval = interval_t{imm};
        switch (bin.op) {
            case Op::MOV: {
                // ra = imm, we forget the type in the region domain
                m_registers -= dst_register;
                break;
            }
            case Op::ADD: {
                // ra += imm
                if (imm == 0) break;
                if (!is_numeric_dst) {
                    update_ptr_or_mapfd(*dst_ptr_or_mapfd_opt, imm_interval, loc, dst_register);
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::SUB: {
                // ra -= imm
                if (imm == 0) break;
                if (!is_numeric_dst) {
                    update_ptr_or_mapfd(*dst_ptr_or_mapfd_opt, -imm_interval, loc, dst_register);
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            default: {
                // no other operations supported for region domain
                m_registers -= dst_register;
                break;
            }
        }
    }
    else {
        auto src_register = register_t{std::get<Reg>(bin.v).v};
        auto src_ptr_or_mapfd_opt = m_registers.find(src_register);
        switch (bin.op) {
            case Op::MOV: {
                // ra = rb
                if (!is_numeric_src) {
                    if (is_shared_ptr(*src_ptr_or_mapfd_opt)) {
                        auto shared_ptr = std::get<ptr_with_off_t>(*src_ptr_or_mapfd_opt);
                        set_aliases(dst_register, shared_ptr);
                        m_registers.insert(dst_register, loc, shared_ptr);
                    }
                    else {
                        m_registers.insert(dst_register, loc, *src_ptr_or_mapfd_opt);
                    }
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::ADD: {
                // ra += rb
                if (!is_numeric_dst && is_numeric_src) {
                    update_ptr_or_mapfd(*dst_ptr_or_mapfd_opt, *src_signed_interval_opt, loc,
                                        dst_register);
                }
                else if (is_numeric_dst && !is_numeric_src) {
                    update_ptr_or_mapfd(*src_ptr_or_mapfd_opt, *dst_signed_interval_opt, loc,
                                        dst_register);
                }
                else if (!is_numeric_dst && !is_numeric_src) {
                    // possibly adding two pointers
                    set_to_bottom();
                }
                break;
            }
            case Op::SUB: {
                // ra -= rb
                if (!is_numeric_dst && is_numeric_src) {
                    update_ptr_or_mapfd(*dst_ptr_or_mapfd_opt, -(*src_signed_interval_opt),
                                        loc, dst_register);
                }
                else if (is_numeric_dst && !is_numeric_src) {
                    m_registers -= dst_register;
                }
                else {
                    // ptr -= ptr
                    // this case already handled in the inference domain
                }
                break;
            }
            default: {
                // no other operations supported for region domain
                m_registers -= dst_register;
                break;
            }
        }
    }
    return;
}

bool region_domain_t::do_load(const Mem& b, const register_t& target_register, location_t loc) {

    int width = b.access.width;
    int offset = b.access.offset;
    Reg basereg = b.access.basereg;

    auto ptr_or_mapfd_opt = m_registers.find(basereg.v);
    bool is_stack_p = is_stack_ptr(ptr_or_mapfd_opt);
    bool is_ctx_p = is_ctx_ptr(ptr_or_mapfd_opt);
    if (!is_ctx_p && !is_stack_p) {
        // loading from either packet or shared region or mapfd does not happen in region domain
        m_registers -= target_register;
        return false;
    }

    auto type_with_off = std::get<ptr_with_off_t>(*ptr_or_mapfd_opt);
    auto p_offset = type_with_off.get_offset();
    auto offset_singleton = p_offset.singleton();

    std::string loc_str = loc.to_string();
    if (!offset_singleton) {
        m_errors.push_back(loc_str + ": Load at an unknown offset");
        m_registers -= target_register;
        return false;
    }

    auto ptr_offset = offset_singleton.value();
    auto load_at = (ptr_offset + offset).cast_to<uint64_t>();
    if (is_stack_p) {
        if (width != 1 && width != 2 && width != 4 && width != 8) {
            m_registers -= target_register;
            m_errors.push_back(loc_str + ": Invalid width for stack load");
            return false;
        }
        if (width != 8) {
            // we do not support loading pointers from stack with width != 8
            // can be relaxed, if needed
            m_registers -= target_register;
            m_errors.push_back(loc_str + ": Only 8-byte stack loads for pointers are supported");
            return false;
        }
        auto loaded = m_stack.find(load_at);
        if (!loaded) {
            // no field at loaded offset in stack, but possibly a number is there
            m_registers -= target_register;
            return false;
        }
        auto ptr_or_mapfd = loaded->first;
        if (is_shared_ptr(ptr_or_mapfd)) {
            auto shared_ptr = std::get<ptr_with_off_t>(ptr_or_mapfd);
            set_aliases((int)target_register, shared_ptr);
            m_registers.insert(target_register, loc, shared_ptr);
        }
        else {
            m_registers.insert(target_register, loc, ptr_or_mapfd);
        }
        return true;
    }
    else {
        if (m_ctx->packet_ptr_at(load_at)) {
            if (width != 4) {
                m_registers -= target_register;
                // Special case for packet pointers, that these are stored into ctx as 4-byte, and
                // also loaded as 4-byte.
                // This probably needs to be relaxed in the future.
                m_errors.push_back(loc_str + ": Only 4-byte ctx loads for pkt pointers are supported");
                return false;
            }
            m_registers.insert(target_register, loc, packet_ptr_t{});
            return true;
        }
        else {
            m_registers -= target_register;
        }
    }
    return false;
}

void region_domain_t::operator()(const Mem& m, location_t loc) {
    // nothing to do here
}

void region_domain_t::do_mem_store(const Mem& b, location_t loc) {

    std::optional<ptr_or_mapfd_t> targetreg_type = {};
    if (std::holds_alternative<Reg>(b.value)) {
        auto target_reg = std::get<Reg>(b.value);
        targetreg_type = m_registers.find(target_reg.v);
    }
    int offset = b.access.offset;
    Reg basereg = b.access.basereg;
    int width = b.access.width;

    auto maybe_basereg_type = m_registers.find(basereg.v);

    bool is_ctx_p = is_ctx_ptr(maybe_basereg_type);
    bool is_shared_p = is_shared_ptr(maybe_basereg_type);
    bool is_packet_p = is_packet_ptr(maybe_basereg_type);
    bool is_mapfd = is_mapfd_type(maybe_basereg_type);

    std::string loc_str = loc.to_string();
    if (is_mapfd) {
        m_errors.push_back(loc_str + ": Cannot store to a mapfd");
        return;
    }
    if (is_shared_p || is_packet_p || is_ctx_p) {
        if (targetreg_type) {
            m_errors.push_back(loc_str + ": Cannot store a pointer into regions other than stack");
            return;
        }
        else {
            // storing a number into a region does not affect the region
            return;
        }
    }

    // if the code reaches here, we are storing into a stack pointer
    auto basereg_type_with_off = std::get<ptr_with_off_t>(*maybe_basereg_type);
    auto offset_reg = basereg_type_with_off.get_offset();
    if (auto finite = offset_reg.finite_size()) {
        int finite_size = finite->cast_to<int>();
        const number_t lb = offset_reg.lb().number().value();
        uint64_t lb_n = lb.cast_to<uint64_t>();
        uint64_t store_at = lb_n + offset;
        m_stack -= m_stack.find_overlapping_cells(store_at, width + finite_size);

        auto offset_singleton = offset_reg.singleton();
        if (!offset_singleton) {
            m_errors.push_back(loc_str + ": Storing at an unknown offset into stack");
            return;
        }
        // if targetreg_type is empty, we are storing a number
        if (!targetreg_type) return;
        auto type = *targetreg_type;
        if (is_shared_ptr(type)) {
            auto shared_ptr = std::get<ptr_with_off_t>(type);
            set_aliases(store_at+11, shared_ptr);
            m_stack.store(store_at, shared_ptr, width);
        }
        else {
            m_stack.store(store_at, type, width);
        }
    }
    else {
        m_errors.push_back(loc_str + ": Storing at an unknown offset into stack");
    }
}

void region_domain_t::adjust_bb_for_types(location_t loc) {
    m_registers.adjust_bb_for_registers(loc);
}

} // namespace crab
