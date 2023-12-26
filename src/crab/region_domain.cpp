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

ctx_t::ctx_t(const ebpf_context_descriptor_t* desc)
{
    if (desc->data >= 0) {
        m_packet_ptrs[desc->data] = std::move(crab::packet_ptr_t{});
    }
    if (desc->end >= 0) {
        m_packet_ptrs[desc->end] = std::move(crab::packet_ptr_t{});
    }
    if (desc->meta >= 0) {
        m_packet_ptrs[desc->meta] = std::move(crab::packet_ptr_t{});
    }
    if (desc->size >= 0) {
        size = desc->size;
    }
}

std::vector<uint64_t> ctx_t::get_keys() const {
    std::vector<uint64_t> keys;
    keys.reserve(size);

    for (auto const&kv : m_packet_ptrs) {
        keys.push_back(kv.first);
    }
    return keys;
}

std::optional<packet_ptr_t> ctx_t::find(uint64_t key) const {
    auto it = m_packet_ptrs.find(key);
    if (it == m_packet_ptrs.end()) return {};
    return it->second;
}

void register_types_t::scratch_caller_saved_registers() {
    for (uint8_t r = R1_ARG; r <= R5_ARG; r++) {
        operator-=(register_t{r});
    }
}

void register_types_t::forget_packet_ptrs() {
    for (uint8_t r = R0_RETURN_VALUE; r < NUM_REGISTERS; r++) {
        if (is_packet_ptr(find(register_t{r}))) {
            operator-=(register_t{r});
        }
    }
}

register_types_t register_types_t::operator|(const register_types_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    auto region_env = std::make_shared<global_region_env_t>();
    register_types_t joined_reg_types(region_env);

    // a hack to store region information at the start of a joined basic block
    // in join, we do not know the label of the bb, hence we store the information
    // at a bb that is not used anywhere else in the program, and later when we know
    // the bb label, we can fix
    location_t loc = location_t{std::make_pair(label_t(-2, -2), 0)};

    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (m_cur_def[i] == nullptr || other.m_cur_def[i] == nullptr) continue;
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
                    joined_reg_types.insert(register_t{i}, loc, std::move(joined_ptr));
                }
            }
            else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1)
                    && std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
                mapfd_t mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
                mapfd_t mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
                joined_reg_types.insert(register_t{i}, loc, std::move(mapfd1 | mapfd2));
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1)
                    && std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
                joined_reg_types.insert(register_t{i}, loc, std::move(packet_ptr_t()));
            }
        }
    }
    return joined_reg_types;
}

void register_types_t::operator-=(register_t var) {
    if (is_bottom()) {
        return;
    }
    m_cur_def[var] = nullptr;
}

void register_types_t::set_to_bottom() {
    m_is_bottom = true;
}

void register_types_t::set_to_top() {
    m_region_env = std::make_shared<global_region_env_t>();
    m_cur_def = live_registers_t{nullptr};
    m_is_bottom = false;
}

bool register_types_t::is_bottom() const { return m_is_bottom; }

bool register_types_t::is_top() const {
    if (m_is_bottom) { return false; }
    if (m_region_env == nullptr) return true;
    for (auto &it : m_cur_def) {
        if (it != nullptr) return false;
    }
    return true;
}

void register_types_t::insert(register_t reg, const location_t& loc, const ptr_or_mapfd_t& type) {
    reg_with_loc_t reg_with_loc = reg_with_loc_t{reg, loc};
    (*m_region_env)[reg_with_loc] = type;
    m_cur_def[reg] = std::make_shared<reg_with_loc_t>(reg_with_loc);
}

std::optional<ptr_or_mapfd_t> register_types_t::find(reg_with_loc_t reg) const {
    auto it = m_region_env->find(reg);
    if (it == m_region_env->end()) return {};
    return it->second;
}

std::optional<ptr_or_mapfd_t> register_types_t::find(register_t key) const {
    if (m_cur_def[key] == nullptr) return {};
    return find(*(m_cur_def[key]));
}

void register_types_t::adjust_bb_for_registers(location_t loc) {
    for (uint8_t i = 0; i < NUM_REGISTERS; i++) {
        if (auto it = find(register_t{i})) {
            insert(register_t{i}, loc, *it);
        }
    }
}

stack_t stack_t::operator|(const stack_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    stack_t joined_stack;
    for (auto const&kv: m_ptrs) {
        auto maybe_ptr_or_mapfd_cells = other.find(kv.first);
        if (maybe_ptr_or_mapfd_cells) {
            auto ptr_or_mapfd_cells1 = kv.second;
            auto ptr_or_mapfd_cells2 = *maybe_ptr_or_mapfd_cells;
            auto ptr_or_mapfd1 = ptr_or_mapfd_cells1.first;
            auto ptr_or_mapfd2 = ptr_or_mapfd_cells2.first;
            int width1 = ptr_or_mapfd_cells1.second;
            int width2 = ptr_or_mapfd_cells2.second;
            int width_joined = std::min(width1, width2);
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd2)) {
                auto ptr_with_off1 = std::get<ptr_with_off_t>(ptr_or_mapfd1);
                auto ptr_with_off2 = std::get<ptr_with_off_t>(ptr_or_mapfd2);
                if (ptr_with_off1.get_region() == ptr_with_off2.get_region()) {
                    joined_stack.store(kv.first, std::move(ptr_with_off1 | ptr_with_off2),
                            width_joined);
                }
            }
            else if (std::holds_alternative<mapfd_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<mapfd_t>(ptr_or_mapfd2)) {
                auto mapfd1 = std::get<mapfd_t>(ptr_or_mapfd1);
                auto mapfd2 = std::get<mapfd_t>(ptr_or_mapfd2);
                joined_stack.store(kv.first, std::move(mapfd1 | mapfd2), width_joined);
            }
            else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd1) &&
                    std::holds_alternative<packet_ptr_t>(ptr_or_mapfd2)) {
                joined_stack.store(kv.first, std::move(packet_ptr_t()), width_joined);
            }
        }
    }
    return joined_stack;
}

void stack_t::operator-=(uint64_t key) {
    auto it = find(key);
    if (it)
        m_ptrs.erase(key);
}

void stack_t::operator-=(const std::vector<uint64_t>& keys) {
    for (auto &key : keys) {
       *this -= key;
    }
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

void stack_t::store(uint64_t key, ptr_or_mapfd_t value, int width) {
    m_ptrs[key] = std::make_pair(value, width);
}

size_t stack_t::size() const {
    return m_ptrs.size();
}

std::vector<uint64_t> stack_t::get_keys() const {
    std::vector<uint64_t> keys;
    keys.reserve(size());

    for (auto const&kv : m_ptrs) {
        keys.push_back(kv.first);
    }
    return keys;
}


std::optional<ptr_or_mapfd_cells_t> stack_t::find(uint64_t key) const {
    auto it = m_ptrs.find(key);
    if (it == m_ptrs.end()) return {};
    return it->second;
}

std::vector<uint64_t> stack_t::find_overlapping_cells(uint64_t start, int width) const {
    std::vector<uint64_t> overlapping_cells;
    auto it = m_ptrs.begin();
    while (it != m_ptrs.end() && it->first < start) {
        it++;
    }
    if (it != m_ptrs.begin()) {
        it--;
        auto key = it->first;
        auto width_key = it->second.second;
        if (key < start && key+width_key > start) overlapping_cells.push_back(key);
    }

    for (; it != m_ptrs.end(); it++) {
        auto key = it->first;
        if (key >= start && key < start+width) overlapping_cells.push_back(key);
        if (key >= start+width) break;
    }
    return overlapping_cells;
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

std::optional<ptr_or_mapfd_t> region_domain_t::find_ptr_or_mapfd_at_loc(const reg_with_loc_t& reg) const {
    return m_registers.find(reg);
}

void region_domain_t::set_registers_to_top() {
    m_registers.set_to_top();
}

size_t region_domain_t::ctx_size() const {
    return m_ctx->get_size();
}

std::vector<uint64_t> region_domain_t::get_ctx_keys() const {
    return m_ctx->get_keys();
}

std::vector<uint64_t> region_domain_t::get_stack_keys() const {
    return m_stack.get_keys();
}

std::optional<packet_ptr_t> region_domain_t::find_in_ctx(uint64_t key) const {
    return m_ctx->find(key);
}

std::optional<ptr_or_mapfd_cells_t> region_domain_t::find_in_stack(uint64_t key) const {
    return m_stack.find(key);
}

bool region_domain_t::operator<=(const region_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return true;
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
    return region_domain_t(m_registers | other.m_registers, m_stack | other.m_stack, other.m_ctx,
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
            m_stack | std::move(other.m_stack), other.m_ctx, std::move(aliases));
}

region_domain_t region_domain_t::operator&(const region_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

region_domain_t region_domain_t::widen(const region_domain_t& abs, bool to_constants) {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

region_domain_t region_domain_t::narrow(const region_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

crab::bound_t region_domain_t::get_loop_count_upper_bound() {
    // WARNING: Not implemented yet.
    return crab::bound_t{crab::number_t{0}};
}

void region_domain_t::initialize_loop_counter(const label_t& label) {
    // WARNING: Not implemented yet.
}

string_invariant region_domain_t::to_set() {
    return string_invariant{};
}

void region_domain_t::operator()(const Undefined &u, location_t loc) {}

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
                m_registers.set_to_top();
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
                m_registers.set_to_top();
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
        if (ctx_ptr.get_offset() == interval_t{crab::number_t{0}}) return;
    }
    //std::cout << "type error: Zero Offset assertion fail\n";
    m_errors.push_back("Zero Ctx Offset assertion fail");
}

void region_domain_t::operator()(const basic_block_t& bb, int print) {
    // nothing to do here
}

void region_domain_t::operator()(const Un& u, location_t loc) {
    m_registers -= u.dst.v;
}

// Get the start and end of the range of possible map fd values.
// In the future, it would be cleaner to use a set rather than an interval
// for map fds.
bool region_domain_t::get_map_fd_range(const Reg& map_fd_reg, int32_t* start_fd, int32_t* end_fd) const {
    auto maybe_type = m_registers.find(map_fd_reg.v);
    if (!is_mapfd_type(maybe_type)) return false;
    auto mapfd_type = std::get<mapfd_t>(*maybe_type);
    const auto& mapfd_interval = mapfd_type.get_mapfd().to_interval();
    auto lb = mapfd_interval.lb().number();
    auto ub = mapfd_interval.ub().number();
    if (!lb || !lb->fits_sint32() || !ub || !ub->fits_sint32())
        return false;
    *start_fd = (int32_t)lb.value();
    *end_fd = (int32_t)ub.value();

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
        if (EbpfMapDescriptor* map = &global_program_info->platform->get_map_descriptor(map_fd))
            result = result | interval_t(number_t(map->key_size));
        else
            return interval_t::top();
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
        if (EbpfMapDescriptor* map = &global_program_info->platform->get_map_descriptor(map_fd))
            result = result | crab::interval_t(number_t(map->value_size));
        else
            return interval_t::top();
    }
    return result;
}

void region_domain_t::do_load_mapfd(const register_t& dst_reg, int mapfd, location_t loc) {
    const auto& platform = global_program_info->platform;
    const EbpfMapDescriptor& desc = platform->get_map_descriptor(mapfd);
    const EbpfMapValueType& map_value_type = platform->get_map_type(desc.type).value_type;
    auto mapfd_interval = interval_t{number_t{mapfd}};
    auto type = mapfd_t(mapfd_interval, map_value_type);
    m_registers.insert(dst_reg, loc, type);
}

void region_domain_t::operator()(const LoadMapFd &u, location_t loc) {
    do_load_mapfd((register_t)u.dst.v, u.mapfd, loc);
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
    for (const auto& kv : cells) {
        auto offset = kv.first;
        auto width = kv.second;
        auto overlapping_cells
            = m_stack.find_overlapping_cells(offset, width);
        m_stack -= overlapping_cells;
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
                        do_load_mapfd(r0, (int)*inner_map_fd, loc);
                        goto out;
                    }
                } else {
                    auto type = ptr_with_off_t(region_t::T_SHARED, -1,
                            interval_t{number_t{0}}, nullness_t::MAYBE_NULL,
                            get_map_value_size(*maybe_fd_reg));
                    set_aliases((int)r0, type);
                    m_registers.insert(r0, loc, type);
                }
            }
        }
        else {
            auto type = ptr_with_off_t(region_t::T_SHARED, -1, interval_t{number_t{0}},
                    nullness_t::MAYBE_NULL);
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

void region_domain_t::operator()(const Callx &u, location_t loc) {
    // WARNING: Not implemented yet.
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

void region_domain_t::check_valid_access(const ValidAccess &s, int width) {
    bool is_comparison_check = s.width == (Value)Imm{0};

    auto maybe_ptr_or_mapfd_type = m_registers.find(s.reg.v);
    if (maybe_ptr_or_mapfd_type) {
        auto reg_ptr_or_mapfd_type = *maybe_ptr_or_mapfd_type;
        if (std::holds_alternative<ptr_with_off_t>(reg_ptr_or_mapfd_type)) {
            auto ptr_with_off_type = std::get<ptr_with_off_t>(reg_ptr_or_mapfd_type);
            auto offset = ptr_with_off_type.get_offset();
            auto offset_to_check = offset.to_interval()+interval_t{s.offset};
            auto offset_lb = offset_to_check.lb();
            auto offset_plus_width_ub = offset_to_check.ub()+bound_t{width};
            if (ptr_with_off_type.get_region() == region_t::T_STACK) {
                if (bound_t{STACK_BEGIN} <= offset_lb
                        && offset_plus_width_ub <= bound_t{EBPF_STACK_SIZE})
                    return;
            }
            else if (ptr_with_off_type.get_region() == region_t::T_CTX) {
                if (bound_t{CTX_BEGIN} <= offset_lb
                        && offset_plus_width_ub <= bound_t{ctx_size()})
                    return;
            }
            else { // shared
                if (crab::bound_t{SHARED_BEGIN} <= offset_lb &&
                        offset_plus_width_ub <= ptr_with_off_type.get_region_size().lb()) {
                    if (!is_comparison_check && !s.or_null) {
                        auto nullness = ptr_with_off_type.get_nullness();
                        if (nullness != nullness_t::NOT_NULL) {
                            m_errors.push_back("possible null access");
                        }
                        return;
                    }
                    return;
                }
            }
        }
        else if (std::holds_alternative<packet_ptr_t>(reg_ptr_or_mapfd_type)) {
            // We do not handle packet ptr access in region domain
            return;
        }
        else {
            // mapfd
            if (is_comparison_check) return;
            //std::cout << "type error: FDs cannot be dereferenced directly\n";
            m_errors.push_back("FDs cannot be dereferenced directly");
        }
        //std::cout << "type error: valid access assert fail\n";
        m_errors.push_back("valid access assert fail");
    }

}

void region_domain_t::operator()(const ValidAccess &s, location_t loc) {
    // nothing to do here
}

region_domain_t&& region_domain_t::setup_entry(bool init_r1) {

    std::shared_ptr<ctx_t> ctx
        = std::make_shared<ctx_t>(global_program_info.get().type.context_descriptor);

    register_types_t typ(std::make_shared<global_region_env_t>());

    auto loc = std::make_pair(label_t::entry, (unsigned int)0);
    if (init_r1) {
        auto ctx_ptr_r1 = ptr_with_off_t(region_t::T_CTX, -1, mock_interval_t{number_t{0}});
        typ.insert(register_t{R1_ARG}, loc, ctx_ptr_r1);
    }
    auto stack_ptr_r10 = ptr_with_off_t(region_t::T_STACK, -1,  mock_interval_t{number_t{512}});
    typ.insert(register_t{R10_STACK_POINTER}, loc, stack_ptr_r10);

    static region_domain_t inv(std::move(typ), stack_t::top(), ctx);
    return std::move(inv);
}

void region_domain_t::operator()(const TypeConstraint& s, location_t loc) {
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
            if (s.types == TypeGroup::non_map_fd) return;
            if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd_type)) {
                ptr_with_off_t ptr_with_off = std::get<ptr_with_off_t>(ptr_or_mapfd_type);
                if (ptr_with_off.get_region() == crab::region_t::T_CTX) {
                    if (s.types == TypeGroup::singleton_ptr) return;
                    if (s.types == TypeGroup::ctx) return;
                }
                else {
                    if (s.types == TypeGroup::mem || s.types == TypeGroup::mem_or_num) return;
                    if (ptr_with_off.get_region() == crab::region_t::T_SHARED) {
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
    else {
        // if we don't know the type, we assume it is a number
        if (s.types == TypeGroup::number || s.types == TypeGroup::ptr_or_num
                || s.types == TypeGroup::non_map_fd || s.types == TypeGroup::mem_or_num)
            return;
    }
    //std::cout << "type error: type constraint assert fail\n";
    m_errors.push_back("type constraint assert fail");
}

void region_domain_t::update_ptr_or_mapfd(ptr_or_mapfd_t&& ptr_or_mapfd, interval_t&& change,
        const location_t& loc, register_t reg) {
    if (std::holds_alternative<ptr_with_off_t>(ptr_or_mapfd)) {
        auto ptr_or_mapfd_with_off = std::get<ptr_with_off_t>(ptr_or_mapfd);
        auto offset = ptr_or_mapfd_with_off.get_offset();
        auto updated_offset = offset.to_interval() + change;
        ptr_or_mapfd_with_off.set_offset(updated_offset);
        m_registers.insert(reg, loc, ptr_or_mapfd_with_off);
    }
    else if (std::holds_alternative<packet_ptr_t>(ptr_or_mapfd)) {
        m_registers.insert(reg, loc, ptr_or_mapfd);
    }
    else {
        m_errors.push_back("mapfd register cannot be incremented/decremented");
        m_registers -= reg;
    }
}

void region_domain_t::operator()(const Bin& b, location_t loc) {
    // nothing to do here
}

interval_t region_domain_t::do_bin(const Bin& bin,
        const std::optional<interval_t>& src_signed_interval_opt,
        const std::optional<ptr_or_mapfd_t>& src_ptr_or_mapfd_opt,
        const std::optional<interval_t>& dst_signed_interval_opt,
        const std::optional<ptr_or_mapfd_t>& dst_ptr_or_mapfd_opt, location_t loc) {

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
        auto imm_interval = interval_t{number_t{imm}};
        switch (bin.op) {
            case Op::MOV: {
                // ra = imm, we forget the type in the region domain
                m_registers -= dst_register;
                break;
            }
            case Op::ADD: {
                // ra += imm
                if (imm == 0) break;
                if (dst_ptr_or_mapfd_opt) {
                    auto dst_ptr_or_mapfd = std::move(*dst_ptr_or_mapfd_opt);
                    update_ptr_or_mapfd(std::move(dst_ptr_or_mapfd),
                            std::move(imm_interval), loc, dst_register);
                }
                else {
                    m_registers -= dst_register;
                }
                break;
            }
            case Op::SUB: {
                // ra -= imm
                if (imm == 0) break;
                if (dst_ptr_or_mapfd_opt) {
                    auto dst_ptr_or_mapfd = std::move(*dst_ptr_or_mapfd_opt);
                    update_ptr_or_mapfd(std::move(dst_ptr_or_mapfd),
                            -std::move(imm_interval), loc, dst_register);
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
        switch (bin.op) {
            case Op::MOV: {
                // ra = rb
                if (src_ptr_or_mapfd_opt) {
                    auto src_ptr_or_mapfd = std::move(*src_ptr_or_mapfd_opt);
                    if (is_shared_ptr(src_ptr_or_mapfd)) {
                        auto shared_ptr = std::get<ptr_with_off_t>(src_ptr_or_mapfd);
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
                if (dst_ptr_or_mapfd_opt && src_signed_interval_opt) {
                    auto dst_ptr_or_mapfd = std::move(*dst_ptr_or_mapfd_opt);
                    auto src_signed = std::move(*src_signed_interval_opt);
                    update_ptr_or_mapfd(std::move(dst_ptr_or_mapfd),
                            std::move(src_signed), loc, dst_register);
                }
                else if (src_ptr_or_mapfd_opt && dst_signed_interval_opt) {
                    auto src_ptr_or_mapfd = std::move(*src_ptr_or_mapfd_opt);
                    auto dst_signed = std::move(*dst_signed_interval_opt);
                    update_ptr_or_mapfd(std::move(src_ptr_or_mapfd),
                            std::move(dst_signed), loc, dst_register);
                }
                else if (dst_signed_interval_opt && src_signed_interval_opt) {
                    // we do not deal with numbers in region domain
                    m_registers -= dst_register;
                }
                else {
                    // possibly adding two pointers
                    set_to_bottom();
                }
                break;
            }
            case Op::SUB: {
                // ra -= rb
                if (dst_ptr_or_mapfd_opt && src_signed_interval_opt) {
                    auto dst_ptr_or_mapfd = std::move(*dst_ptr_or_mapfd_opt);
                    auto src_signed = std::move(*src_signed_interval_opt);
                    update_ptr_or_mapfd(std::move(dst_ptr_or_mapfd),
                            -std::move(src_signed), loc, dst_register);
                }
                else if (src_ptr_or_mapfd_opt && dst_signed_interval_opt) {
                    m_registers -= dst_register;
                }
                else if (dst_signed_interval_opt && src_signed_interval_opt) {
                    // we do not deal with numbers in region domain
                    m_registers -= dst_register;
                }
                else {
                    // ptr -= ptr
                    if (std::holds_alternative<mapfd_t>(*dst_ptr_or_mapfd_opt) &&
                            std::holds_alternative<mapfd_t>(*src_ptr_or_mapfd_opt)) {
                        m_errors.push_back("mapfd registers subtraction not defined");
                    }
                    else if (same_region(*dst_ptr_or_mapfd_opt, *src_ptr_or_mapfd_opt)) {
                        if (std::holds_alternative<ptr_with_off_t>(*dst_ptr_or_mapfd_opt) &&
                                std::holds_alternative<ptr_with_off_t>(*src_ptr_or_mapfd_opt)) {
                            auto dst_ptr_with_off = std::get<ptr_with_off_t>(*dst_ptr_or_mapfd_opt);
                            auto src_ptr_with_off = std::get<ptr_with_off_t>(*src_ptr_or_mapfd_opt);
                            return (dst_ptr_with_off.get_offset().to_interval() -
                                src_ptr_with_off.get_offset().to_interval());
                        }
                    }
                    else {
                        // Assertions should make sure we only perform this on
                        // non-shared pointers, hence this should not happen
                        m_errors.push_back("subtraction between pointers of different region");
                    }
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
    return interval_t::bottom();
}

void region_domain_t::do_load(const Mem& b, const register_t& target_register, bool unknown_ptr,
        location_t loc) {

    if (unknown_ptr) {
        m_registers -= target_register;
        return;
    }

    int width = b.access.width;
    int offset = b.access.offset;
    Reg basereg = b.access.basereg;

    auto ptr_or_mapfd_opt = m_registers.find(basereg.v);
    bool is_stack_p = is_stack_ptr(ptr_or_mapfd_opt);
    bool is_ctx_p = is_ctx_ptr(ptr_or_mapfd_opt);
    if (!is_ctx_p && !is_stack_p) {
        // loading from either packet or shared region or mapfd does not happen in region domain
        m_registers -= target_register;
        return;
    }

    auto type_with_off = std::get<ptr_with_off_t>(*ptr_or_mapfd_opt);
    auto p_offset = type_with_off.get_offset();
    auto offset_singleton = p_offset.to_interval().singleton();

    if (is_stack_p) {
        if (!offset_singleton) {
            for (auto const& k : m_stack.get_keys()) {
                auto start = p_offset.lb();
                auto end = p_offset.ub()+number_t{offset+width-1};
                interval_t range{start, end};
                if (range[number_t{(int)k}]) {
                    //std::cout << "stack load at unknown offset, and offset range contains pointers\n";
                    m_errors.push_back("stack load at unknown offset, and offset range contains pointers");
                    break;
                }
            }
            m_registers -= target_register;
        }
        else {
            if (width != 1 && width != 2 && width != 4 && width != 8) {
                m_registers -= target_register;
                return;
            }
            auto ptr_offset = offset_singleton.value();
            auto load_at = (uint64_t)(ptr_offset + offset);

            auto loaded = m_stack.find(load_at);
            if (!loaded) {
                // no field at loaded offset in stack
                m_registers -= target_register;
                return;
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
        }
    }
    else {
        if (!offset_singleton) {
            for (auto const& k : m_ctx->get_keys()) {
                auto start = p_offset.lb();
                auto end = p_offset.ub()+crab::bound_t{offset+width-1};
                interval_t range{start, end};
                if (range[number_t{(int)k}]) {
                    //std::cout << "ctx load at unknown offset, and offset range contains pointers\n";
                    m_errors.push_back("ctx load at unknown offset, and offset range contains pointers");
                    break;
                }
            }
            m_registers -= target_register;
        }
        else {
            auto ptr_offset = offset_singleton.value();
            auto load_at = (uint64_t)(ptr_offset + offset);

            auto loaded = m_ctx->find(load_at);
            if (!loaded) {
                // no field at loaded offset in ctx
                m_registers -= target_register;
                return;
            }
            m_registers.insert(target_register, loc, *loaded);
        }
    }
}

void region_domain_t::operator()(const Mem& m, location_t loc) {
    // nothing to do here
}

void region_domain_t::do_mem_store(const Mem& b, location_t loc) {

    std::optional<ptr_or_mapfd_t> targetreg_type = {};
    bool target_is_reg = std::holds_alternative<Reg>(b.value);
    if (target_is_reg) {
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

    if (is_mapfd) {
        m_errors.push_back("storing into a mapfd register is not defined");
        return;
    }
    if (is_shared_p || is_packet_p || is_ctx_p) {
        if (targetreg_type) {
            m_errors.push_back("storing a pointer into a shared, packet or ctx pointer");
            return;
        }
        else {
            // storing a number into a region does not affect the region
            return;
        }
    }

    // if the code reaches here, we are storing into a stack pointer
    auto basereg_type_with_off = std::get<ptr_with_off_t>(*maybe_basereg_type);
    auto offset_singleton = basereg_type_with_off.get_offset().to_interval().singleton();
    if (!offset_singleton) {
        //std::cout << "type error: storing to a pointer with unknown offset\n";
        m_errors.push_back("storing to a pointer with unknown offset");
        return;
    }
    auto store_at = (uint64_t)offset+(uint64_t)offset_singleton.value();
    auto overlapping_cells = m_stack.find_overlapping_cells(store_at, width);
    m_stack -= overlapping_cells;

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

void region_domain_t::adjust_bb_for_types(location_t loc) {
    m_registers.adjust_bb_for_registers(loc);
}

} // namespace crab
