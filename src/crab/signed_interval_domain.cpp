// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "crab/signed_interval_domain.hpp"
#include "boost/endian/conversion.hpp"

namespace crab {

bool registers_signed_state_t::is_bottom() const {
    return m_is_bottom;
}

bool registers_signed_state_t::is_top() const {
    if (m_is_bottom) return false;
    if (m_interval_env == nullptr) return true;
    for (auto it : m_cur_def) {
        if (it != nullptr) return false;
    }
    return true;
}

void registers_signed_state_t::set_to_top() {
    m_interval_env = std::make_shared<global_interval_env_t>();
    m_cur_def = live_registers_t{nullptr};
    m_is_bottom = false;
}

void registers_signed_state_t::set_to_bottom() {
    m_is_bottom = true;
}

void registers_signed_state_t::insert(register_t reg, const location_t& loc, interval_t interval) {
    auto reg_with_loc = reg_with_loc_t{reg, loc};
    (*m_interval_env)[reg_with_loc] = mock_interval_t{interval};
    m_cur_def[reg] = std::make_shared<reg_with_loc_t>(reg_with_loc);
}

std::optional<mock_interval_t> registers_signed_state_t::find(reg_with_loc_t reg) const {
    auto it = m_interval_env->find(reg);
    if (it == m_interval_env->end()) return {};
    return it->second;
}

std::optional<mock_interval_t> registers_signed_state_t::find(register_t key) const {
    if (m_cur_def[key] == nullptr) return {};
    const reg_with_loc_t& reg = *(m_cur_def[key]);
    return find(reg);
}

registers_signed_state_t registers_signed_state_t::operator|(const registers_signed_state_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    registers_signed_state_t intervals_joined(m_interval_env);
    location_t loc = location_t(std::make_pair(label_t(-2, -2), 0));
    for (uint8_t i = 0; i < NUM_REGISTERS-1; i++) {
        if (m_cur_def[i] == nullptr || other.m_cur_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_def[i]));
        auto it2 = other.find(*(other.m_cur_def[i]));
        if (it1 && it2) {
            auto interval1 = it1->to_interval(), interval2 = it2->to_interval();
            intervals_joined.insert(register_t{i}, loc,
                    std::move(interval1 | interval2));
        }
    }
    return intervals_joined;
}

void registers_signed_state_t::adjust_bb_for_registers(location_t loc) {
    for (uint8_t i = 0; i < NUM_REGISTERS-1; i++) {
        if (auto it = find(register_t{i})) {
            insert(register_t{i}, loc, it->to_interval());
        }
    }
}

void registers_signed_state_t::operator-=(register_t var) {
    if (is_bottom()) return;
    m_cur_def[var] = nullptr;
}

bool stack_slots_signed_state_t::is_bottom() const {
    return m_is_bottom;
}

bool stack_slots_signed_state_t::is_top() const {
    if (m_is_bottom) return false;
    return m_interval_values.empty();
}

void stack_slots_signed_state_t::set_to_top() {
    m_interval_values.clear();
    m_is_bottom = false;
}

void stack_slots_signed_state_t::set_to_bottom() {
    m_is_bottom = true;
}

stack_slots_signed_state_t stack_slots_signed_state_t::top() {
    return stack_slots_signed_state_t(false);
}

std::optional<interval_cells_t> stack_slots_signed_state_t::find(uint64_t key) const {
    auto it = m_interval_values.find(key);
    if (it == m_interval_values.end()) return {};
    return it->second;
}

void stack_slots_signed_state_t::store(uint64_t key, mock_interval_t val, int width) {
    m_interval_values[key] = std::make_pair(val, width);
}

void stack_slots_signed_state_t::operator-=(uint64_t key) {
    auto it = find(key);
    if (it)
        m_interval_values.erase(key);
}

bool stack_slots_signed_state_t::all_numeric(uint64_t start_loc, int width) const {
    auto overlapping_cells = find_overlapping_cells(start_loc, width);
    if (overlapping_cells.empty()) return false;
    for (std::size_t i = 0; i < overlapping_cells.size()-1; i++) {
        int width_i = find(overlapping_cells[i]).value().second;
        if (overlapping_cells[i]+width_i != overlapping_cells[i+1]) return false;
    }
    return true;
}

void stack_slots_signed_state_t::fill_values(const std::vector<uint64_t>& keys,
        uint64_t start, int width) {
    if (keys[0] < start) {
        auto type = find(keys[0]);
        auto width_key = type.value().second;
        store(keys[0], interval_t::top(), width_key);
    }
    if (keys[0] > start) {
        store(start, interval_t::top(), keys[0]-start);
    }
    for (size_t i = 0; i < keys.size()-1; i++) {
        auto type = find(keys[i]);
        auto width_key = type.value().second;
        if (keys[i]+width_key != keys[i+1]) {
            store(keys[i]+width_key, interval_t::top(), keys[i+1]-(keys[i]+width_key));
        }
    }
    auto type = find(keys[keys.size()-1]);
    auto width_key = type.value().second;
    if (keys[keys.size()-1]+width_key < start+width) {
        store(keys[keys.size()-1]+width_key, interval_t::top(),
                start+width-(keys[keys.size()-1]+width_key));
    }
    if (keys[keys.size()-1]+width_key > start+width) {
        store(keys[keys.size()-1], interval_t::top(), width_key);
    }
}

void stack_slots_signed_state_t::remove_overlap(const std::vector<uint64_t>& keys, uint64_t start,
        int width) {
    for (auto& key : keys) {
        auto type = find(key);
        auto width_key = type.value().second;
        if (key < start) {
            int new_width = start-key;
            store(key, interval_t::top(), new_width);
        }
        if (key+width_key > start+width) {
            int new_width = key+width_key-(start+width);
            store(start+width, interval_t::top(), new_width);
        }
        if (key >= start) *this -= key;
    }
}

std::vector<uint64_t> stack_slots_signed_state_t::find_overlapping_cells(uint64_t start, int width) const {
    std::vector<uint64_t> overlapping_cells;
    auto it = m_interval_values.begin();
    while (it != m_interval_values.end() && it->first < start) {
        it++;
    }
    if (it != m_interval_values.begin()) {
        it--;
        auto key = it->first;
        auto width_key = it->second.second;
        if (key < start && key+width_key > start) overlapping_cells.push_back(key);
    }

    for (; it != m_interval_values.end(); it++) {
        auto key = it->first;
        if (key >= start && key < start+width) overlapping_cells.push_back(key);
        if (key >= start+width) break;
    }
    return overlapping_cells;
}

static inline void join_stack(const stack_slots_signed_state_t& stack1, uint64_t key1, int& loc1,
        const stack_slots_signed_state_t& stack2, uint64_t key2, int& loc2,
        interval_values_stack_t& interval_values_joined) {
    auto type1 = stack1.find(key1);    auto type2 = stack2.find(key2);
    auto& cells1 = type1.value();   auto& cells2 = type2.value();
    int width1 = cells1.second; int width2 = cells2.second;
    auto interval1 = cells1.first.to_interval();
    auto interval2 = cells2.first.to_interval();
    auto top_interval_mock = mock_interval_t(interval_t::top());
    if (key1 == key2) {
        if (width1 == width2) {
            interval_values_joined[key1] = std::make_pair(
                    mock_interval_t(interval1 | interval2), width1);
            loc1++; loc2++;
        }
        else if (width1 < width2) {
            interval_values_joined[key1] = std::make_pair(top_interval_mock, width1);
            loc1++;
        }
        else {
            interval_values_joined[key1] = std::make_pair(top_interval_mock, width2);
            loc2++;
        }
    }
    else if (key1 > key2) {
        if (key2+width2 > key1+width1) {
            interval_values_joined[key1] = std::make_pair(top_interval_mock, width1);
            loc1++;
        }
        else if (key2+width2 > key1) {
            interval_values_joined[key1] = std::make_pair(top_interval_mock, key2+width2-key1);
            loc2++;
        }
        else loc2++;
    }
    else {
        join_stack(stack2, key2, loc2, stack1, key1, loc1, interval_values_joined);
    }
}

stack_slots_signed_state_t stack_slots_signed_state_t::operator|(const stack_slots_signed_state_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    interval_values_stack_t interval_values_joined;
    auto stack1_keys = get_keys();
    auto stack2_keys = other.get_keys();
    int i = 0, j = 0;
    while (i < static_cast<int>(stack1_keys.size()) && j < static_cast<int>(stack2_keys.size())) {
        int key1 = stack1_keys[i], key2 = stack2_keys[j];
        join_stack(*this, key1, i, other, key2, j, interval_values_joined);
    }
    return stack_slots_signed_state_t(std::move(interval_values_joined));
}

size_t stack_slots_signed_state_t::size() const {
    return m_interval_values.size();
}

std::vector<uint64_t> stack_slots_signed_state_t::get_keys() const {
    std::vector<uint64_t> keys;
    keys.reserve(size());

    for (auto const&kv : m_interval_values) {
        keys.push_back(kv.first);
    }
    return keys;
}

bool signed_interval_domain_t::is_bottom() const {
    if (m_is_bottom) return true;
    return (m_registers_values.is_bottom() || m_stack_slots_values.is_bottom());
}

bool signed_interval_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_registers_values.is_top() && m_stack_slots_values.is_top());
}

signed_interval_domain_t signed_interval_domain_t::bottom() {
    signed_interval_domain_t interval;
    interval.set_to_bottom();
    return interval;
}

void signed_interval_domain_t::set_to_bottom() {
    m_is_bottom = true;
    m_registers_values.set_to_bottom();
    m_stack_slots_values.set_to_bottom();
}

void signed_interval_domain_t::set_registers_to_bottom() {
    m_registers_values.set_to_bottom();
}

void signed_interval_domain_t::set_registers_to_top() {
    m_registers_values.set_to_top();
}

void signed_interval_domain_t::set_to_top() {
    m_is_bottom = false;
    m_registers_values.set_to_top();
    m_stack_slots_values.set_to_top();
}

std::optional<interval_cells_t> signed_interval_domain_t::find_in_stack(uint64_t key) const {
    return m_stack_slots_values.find(key);
}

void signed_interval_domain_t::adjust_bb_for_types(location_t loc) {
    m_registers_values.adjust_bb_for_registers(loc);
}

std::vector<uint64_t> signed_interval_domain_t::get_stack_keys() const {
    return m_stack_slots_values.get_keys();
}

bool signed_interval_domain_t::all_numeric_in_stack(uint64_t start_loc, int width) const {
    return m_stack_slots_values.all_numeric(start_loc, width);
}

std::vector<uint64_t> signed_interval_domain_t::find_overlapping_cells_in_stack(uint64_t start_loc,
        int width) const {
    return m_stack_slots_values.find_overlapping_cells(start_loc, width);
}

void signed_interval_domain_t::remove_overlap_in_stack(const std::vector<uint64_t>& overlap,
        uint64_t start_loc, int width) {
    m_stack_slots_values.remove_overlap(overlap, start_loc, width);
}

void signed_interval_domain_t::fill_values_in_stack(const std::vector<uint64_t>& overlap,
        uint64_t start_loc, int width) {
    m_stack_slots_values.fill_values(overlap, start_loc, width);
}

std::optional<mock_interval_t> signed_interval_domain_t::find_interval_value(register_t reg) const {
    return m_registers_values.find(reg);
}

std::optional<mock_interval_t> signed_interval_domain_t::find_interval_at_loc(
        const reg_with_loc_t reg) const {
    return m_registers_values.find(reg);
}

void signed_interval_domain_t::insert_in_registers(register_t reg, location_t loc,
        interval_t interval) {
    m_registers_values.insert(reg, loc, interval);
}

void signed_interval_domain_t::store_in_stack(uint64_t key, mock_interval_t interval, int width) {
    m_stack_slots_values.store(key, interval, width);
}

bool signed_interval_domain_t::operator<=(const signed_interval_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return true;
}

void signed_interval_domain_t::operator|=(const signed_interval_domain_t& abs) {
    signed_interval_domain_t tmp{abs};
    operator|=(std::move(tmp));
}

void signed_interval_domain_t::operator|=(signed_interval_domain_t&& abs) {
    if (is_bottom()) {
        *this = abs;
        return;
    }
    *this = *this | std::move(abs);
}

signed_interval_domain_t signed_interval_domain_t::operator|(const signed_interval_domain_t& other) const {
    if (is_bottom() || other.is_top()) {
        return other;
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return signed_interval_domain_t(m_registers_values | other.m_registers_values,
            m_stack_slots_values | other.m_stack_slots_values);
}

signed_interval_domain_t signed_interval_domain_t::operator|(signed_interval_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return signed_interval_domain_t(
            m_registers_values | std::move(other.m_registers_values),
            m_stack_slots_values | std::move(other.m_stack_slots_values));
}

signed_interval_domain_t signed_interval_domain_t::operator&(const signed_interval_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

signed_interval_domain_t signed_interval_domain_t::widen(const signed_interval_domain_t& abs, bool to_constants) {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

signed_interval_domain_t signed_interval_domain_t::narrow(const signed_interval_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

crab::bound_t signed_interval_domain_t::get_loop_count_upper_bound() {
    /* WARNING: The operation is not implemented yet.*/
    return crab::bound_t{crab::number_t{0}};
}

void signed_interval_domain_t::initialize_loop_counter(const label_t& label) {
    /* WARNING: The operation is not implemented yet.*/
}

string_invariant signed_interval_domain_t::to_set() {
    return string_invariant{};
}

signed_interval_domain_t signed_interval_domain_t::setup_entry() {
    registers_signed_state_t registers(std::make_shared<global_interval_env_t>());
    signed_interval_domain_t interval(std::move(registers), stack_slots_signed_state_t::top());
    return interval;
}

void signed_interval_domain_t::operator()(const Un& u, location_t loc) {
    auto swap_endianness = [&](interval_t&& v, auto input, const auto& be_or_le) {
        if (std::optional<number_t> n = v.singleton()) {
            if (n->fits_cast_to_int64()) {
                input = (decltype(input))n.value().cast_to_sint64();
                decltype(input) output = be_or_le(input);
                m_registers_values.insert(u.dst.v, loc, interval_t{number_t{output}});
                return;
            }
        }
        m_registers_values.insert(u.dst.v, loc, interval_t::top());
    };

    auto mock_interval = m_registers_values.find(u.dst.v);
    if (!mock_interval) return;
    interval_t interval = mock_interval->to_interval();

    // Swap bytes.  For 64-bit types we need the weights to fit in a
    // signed int64, but for smaller types we don't want sign extension,
    // so we use unsigned which still fits in a signed int64.
    switch (u.op) {
    case Un::Op::BE16:
        swap_endianness(std::move(interval), uint16_t(0),
                boost::endian::native_to_big<uint16_t>);
        break;
    case Un::Op::BE32:
        swap_endianness(std::move(interval), uint32_t(0),
                boost::endian::native_to_big<uint32_t>);
        break;
    case Un::Op::BE64:
        swap_endianness(std::move(interval), int64_t(0),
                boost::endian::native_to_big<int64_t>);
        break;
    case Un::Op::LE16:
        swap_endianness(std::move(interval), uint16_t(0),
                boost::endian::native_to_little<uint16_t>);
        break;
    case Un::Op::LE32:
        swap_endianness(std::move(interval), uint32_t(0),
                boost::endian::native_to_little<uint32_t>);
        break;
    case Un::Op::LE64:
        swap_endianness(std::move(interval), int64_t(0),
                boost::endian::native_to_little<int64_t>);
        break;
    case Un::Op::NEG:
        auto reg_with_loc = reg_with_loc_t(u.dst.v, loc);
        m_registers_values.insert(u.dst.v, loc, -interval);
        break;
    }
}

void signed_interval_domain_t::operator()(const LoadMapFd& u, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Packet& u, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Assume& s, location_t loc) {
    // nothing to do here
}

bool signed_interval_domain_t::load_from_stack(register_t reg, uint64_t offset, location_t loc) {
    if (auto loaded = m_stack_slots_values.find(offset)) {
        m_registers_values.insert(reg, loc, (*loaded).first.to_interval());
        return true;
    }
    return false;
}

void signed_interval_domain_t::store_in_stack(const Mem& b, uint64_t offset, int width) {
    if (std::holds_alternative<Reg>(b.value)) {
        auto target_reg = std::get<Reg>(b.value);
        if (auto targetreg_mock_interval = m_registers_values.find(target_reg.v)) {
            auto targetreg_type = targetreg_mock_interval->to_interval();
            m_stack_slots_values.store(offset, targetreg_type, width);
        }
    }
    else {
        auto imm = (int64_t)std::get<Imm>(b.value).v;
        auto targetreg_type = interval_t{number_t{imm}};
        m_stack_slots_values.store(offset, targetreg_type, width);
    }
}

void signed_interval_domain_t::check_valid_access(const ValidAccess& s, interval_t&& interval,
        int width, bool check_stack_all_numeric) {

    if (check_stack_all_numeric) {
        auto start_interval = interval + interval_t{number_t{s.offset}};
        if (auto finite_size = start_interval.finite_size()) {
            if (auto start_interval_lb = start_interval.lb().number()) {
                auto start_offset = (uint64_t)(*start_interval_lb);
                int width_from_start = (int)(*finite_size) + width;
                if (!m_stack_slots_values.all_numeric(start_offset, width_from_start)) {
                    m_errors.push_back("Stack access not numeric");
                }
            }
            else {
                m_errors.push_back("Offset information not available");
            }
        }
        else {
            m_errors.push_back("Register interval not finite for stack access");
        }
    }
    else {
        bool is_comparison_check = s.width == (Value)Imm{0};
        if (!is_comparison_check) {
            if (s.or_null) {
                if (auto singleton = interval.singleton()) {
                    if (*singleton == number_t{0}) return;
                }
                m_errors.push_back("Non-null number");
            }
            else {
                m_errors.push_back("Only pointers can be dereferenced");
            }
        }
    }
}

void signed_interval_domain_t::operator()(const Undefined& u, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Bin& b, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Call&, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Exit&, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Jmp&, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Mem&, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Assert&, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const basic_block_t& bb, int print) {
    // nothing to do here
}

void signed_interval_domain_t::set_require_check(check_require_func_t f) {}

} // namespace crab
