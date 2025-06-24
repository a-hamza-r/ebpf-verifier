// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#include "signed_interval_domain.hpp"
#include "boost/endian/conversion.hpp"
#include "config.hpp"

namespace crab {

bool signed_interval_registers_t::is_bottom() const {
    return m_is_bottom;
}

bool signed_interval_registers_t::is_top() const {
    if (m_is_bottom) return false;
    if (m_registers_env == nullptr) return true;
    for (auto it : m_cur_register_def) {
        if (it != nullptr) return false;
    }
    return true;
}

void signed_interval_registers_t::set_to_top() {
    m_registers_env = std::make_shared<global_env_signed_registers_t>();
    m_cur_register_def = live_registers_t{nullptr};
    m_is_bottom = false;
}

void signed_interval_registers_t::set_to_bottom() {
    m_is_bottom = true;
}

void signed_interval_registers_t::insert(register_t reg, const location_t& loc, refinement_t rf) {
    auto register_location = register_location_t{reg, loc};
    (*m_registers_env)[register_location] = rf;
    m_cur_register_def[reg] = std::make_shared<register_location_t>(register_location);
}

std::optional<refinement_t> signed_interval_registers_t::find(register_location_t reg) const {
    auto it = m_registers_env->find(reg);
    if (it == m_registers_env->end()) return {};
    return it->second;
}

std::optional<refinement_t> signed_interval_registers_t::find(register_t key) const {
    if (m_cur_register_def[key] == nullptr) return {};
    const register_location_t& reg = *(m_cur_register_def[key]);
    return find(reg);
}

bool signed_interval_registers_t::inclusion(const signed_interval_registers_t& other,
                                            std::shared_ptr<slacks_t> slacks) const {
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (other.m_cur_register_def[i] == nullptr) continue;
        if (m_cur_register_def[i] == nullptr) return false;
        auto it1 = find(*(m_cur_register_def[i]));
        auto it2 = other.find(*(m_cur_register_def[i]));
        if (it1 && it2) {
            refinement_t rf1 = it1.value(), rf2 = it2.value();
            if (!(rf1.inclusion(rf2, slacks))) {
                return false;
            }
        }
    }
    return true;
}

signed_interval_registers_t signed_interval_registers_t::join(
    const signed_interval_registers_t& other, std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    signed_interval_registers_t refinements_joined;
    location_t loc = location_t::top();
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (m_cur_register_def[i] == nullptr || other.m_cur_register_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_register_def[i]));
        auto it2 = other.find(*(other.m_cur_register_def[i]));
        if (it1 && it2) {
            refinement_t rf1 = it1.value(), rf2 = it2.value();
            refinements_joined.insert(register_t{i}, loc, rf1.join(rf2, slacks));
        }
    }
    return refinements_joined;
}

signed_interval_registers_t signed_interval_registers_t::widen(
    const signed_interval_registers_t& other, std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    signed_interval_registers_t refinements_joined;
    location_t loc = location_t::top();
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (other.m_cur_register_def[i] == nullptr) continue;
        auto it1 = find(*(m_cur_register_def[i]));
        auto it2 = other.find(*(other.m_cur_register_def[i]));
        if (it1 && it2) {
            refinement_t rf1 = it1.value(), rf2 = it2.value();
            refinements_joined.insert(register_t{i}, loc, rf1.widen(rf2, slacks));
        }
    }
    return refinements_joined;
}

void signed_interval_registers_t::adjust_bb_for_registers(location_t loc) {
    for (uint8_t i = 0; i < NUM_REGISTERS-2; i++) {
        if (auto it = find(register_t{i})) {
            insert(register_t{i}, loc, it.value());
        }
    }
}

void signed_interval_registers_t::operator-=(register_t var) {
    if (is_bottom()) return;
    m_cur_register_def[var] = nullptr;
}

bool signed_interval_stack_t::is_bottom() const {
    return m_is_bottom;
}

bool signed_interval_stack_t::is_top() const {
    if (m_is_bottom) return false;
    return m_cells.empty();
}

void signed_interval_stack_t::set_to_top() {
    m_cells.clear();
    m_is_bottom = false;
}

void signed_interval_stack_t::set_to_bottom() {
    m_is_bottom = true;
}

signed_interval_stack_t signed_interval_stack_t::top() {
    return signed_interval_stack_t();
}

std::optional<signed_interval_stack_cell_t> signed_interval_stack_t::find(uint64_t key) const {
    auto it = m_cells.find(key);
    if (it == m_cells.end()) return {};
    return it->second;
}

void signed_interval_stack_t::store(uint64_t key, refinement_t val, int width) {
    m_cells[key] = std::make_pair(val, width);
}

void signed_interval_stack_t::operator-=(uint64_t key) {
    auto it = find(key);
    if (it)
        m_cells.erase(key);
}

bool signed_interval_stack_t::all_numeric(uint64_t start_loc, int width) const {
    auto overlapping_cells = find_overlapping_cells(start_loc, width);
    if (overlapping_cells.empty()) return false;
    for (std::size_t i = 0; i < overlapping_cells.size()-1; i++) {
        int width_i = find(overlapping_cells[i]).value().second;
        if (overlapping_cells[i]+width_i != overlapping_cells[i+1]) return false;
    }
    return true;
}

void signed_interval_stack_t::fill_values(const std::vector<uint64_t>& keys,
        uint64_t start, int width) {
    // numeric_refinement_top() creates interval_t::top()
    refinement_t top_rf = refinement_t::numeric_refinement_top();
    if (keys[0] < start) {
        auto type = find(keys[0]);
        auto width_key = type.value().second;
        store(keys[0], top_rf, width_key);
    }
    if (keys[0] > start) {
        store(start, top_rf, keys[0]-start);
    }
    for (size_t i = 0; i < keys.size()-1; i++) {
        auto type = find(keys[i]);
        auto width_key = type.value().second;
        if (keys[i]+width_key != keys[i+1]) {
            store(keys[i]+width_key, top_rf, keys[i+1]-(keys[i]+width_key));
        }
    }
    auto type = find(keys[keys.size()-1]);
    auto width_key = type.value().second;
    if (keys[keys.size()-1]+width_key < start+width) {
        store(keys[keys.size()-1]+width_key, top_rf,
                start+width-(keys[keys.size()-1]+width_key));
    }
    if (keys[keys.size()-1]+width_key > start+width) {
        store(keys[keys.size()-1], top_rf, width_key);
    }
}

void signed_interval_stack_t::remove_overlap(const std::vector<uint64_t>& keys, uint64_t start,
        int width) {
    // numeric_refinement_top() creates interval_t::top()
    refinement_t top_rf = refinement_t::numeric_refinement_top();
    for (auto& key : keys) {
        auto type = find(key);
        auto width_key = type.value().second;
        if (key < start) {
            int new_width = start-key;
            store(key, top_rf, new_width);
        }
        if (key+width_key > start+width) {
            int new_width = key+width_key-(start+width);
            store(start+width, top_rf, new_width);
        }
        if (key >= start) *this -= key;
    }
}

std::vector<uint64_t> signed_interval_stack_t::find_overlapping_cells(uint64_t start, int width) const {
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

bool signed_interval_stack_t::inclusion(const signed_interval_stack_t& other,
                                        std::shared_ptr<slacks_t> slacks) const {
    size_t size1 = size(), size2 = other.size();
    if (size2 > size1) return false;
    for (auto const &kv : other.m_cells) {
        auto key = kv.first;
        auto it = m_cells.find(key);
        if (it == m_cells.end()) return false;
        auto& cells1 = it->second; auto& cells2 = kv.second;
        auto rf1 = cells1.first;
        auto rf2 = cells2.first;
        auto width1 = cells1.second; auto width2 = cells2.second;
        if (width1 != width2) return false;
        if (!(rf1.inclusion(rf2, slacks))) {
            return false;
        }
    }
    return true;
}

using RefinementJoin = std::function<refinement_t(const refinement_t&, const refinement_t&)>;

static inline void join_stack(const signed_interval_stack_t& stack1, uint64_t key1, int& loc1,
        const signed_interval_stack_t& stack2, uint64_t key2, int& loc2,
        signed_interval_stack_cells_t& refinements_joined, RefinementJoin joinFunc) {
    auto type1 = stack1.find(key1);    auto type2 = stack2.find(key2);
    auto& cells1 = type1.value();   auto& cells2 = type2.value();
    int width1 = cells1.second; int width2 = cells2.second;
    refinement_t rf1 = cells1.first;
    refinement_t rf2 = cells2.first;
    // numeric_refinement_top() creates interval_t::top()
    refinement_t top_rf = refinement_t::numeric_refinement_top();
    if (key1 == key2) {
        if (width1 == width2) {
            refinements_joined[key1] = std::make_pair(joinFunc(rf1, rf2), width1);
            loc1++; loc2++;
        }
        else if (width1 < width2) {
            refinements_joined[key1] = std::make_pair(top_rf, width1);
            loc1++;
        }
        else {
            refinements_joined[key1] = std::make_pair(top_rf, width2);
            loc2++;
        }
    }
    else if (key1 > key2) {
        if (key2+width2 > key1+width1) {
            refinements_joined[key1] = std::make_pair(top_rf, width1);
            loc1++;
        }
        else if (key2+width2 > key1) {
            refinements_joined[key1] = std::make_pair(top_rf, key2+width2-key1);
            loc2++;
        }
        else loc2++;
    }
    else {
        join_stack(stack2, key2, loc2, stack1, key1, loc1, refinements_joined, joinFunc);
    }
}

signed_interval_stack_t signed_interval_stack_t::join(const signed_interval_stack_t& other,
                                                      std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    signed_interval_stack_cells_t refinements_joined;
    auto stack1_keys = get_keys();
    auto stack2_keys = other.get_keys();
    int i = 0, j = 0;
    while (i < static_cast<int>(stack1_keys.size()) && j < static_cast<int>(stack2_keys.size())) {
        int key1 = stack1_keys[i], key2 = stack2_keys[j];
        join_stack(*this, key1, i, other, key2, j, refinements_joined,
                   [slacks](const refinement_t& a, const refinement_t& b) {
                    return a.join(b, slacks);
                });
    }
    return signed_interval_stack_t(std::move(refinements_joined));
}

signed_interval_stack_t signed_interval_stack_t::widen(const signed_interval_stack_t& other,
                                                        std::shared_ptr<slacks_t> slacks) const {
    if (is_bottom() || other.is_top()) {
        return other;
    } else if (other.is_bottom() || is_top()) {
        return *this;
    }
    signed_interval_stack_cells_t refinements_joined;
    auto stack1_keys = get_keys();
    auto stack2_keys = other.get_keys();
    int i = 0, j = 0;
    while (i < static_cast<int>(stack1_keys.size()) && j < static_cast<int>(stack2_keys.size())) {
        int key1 = stack1_keys[i], key2 = stack2_keys[j];
        join_stack(*this, key1, i, other, key2, j, refinements_joined,
                [slacks](const refinement_t& a, const refinement_t& b) {
                   return a.widen(b, slacks);
                });
    }
    return signed_interval_stack_t(std::move(refinements_joined));
}

size_t signed_interval_stack_t::size() const {
    return m_cells.size();
}

std::vector<uint64_t> signed_interval_stack_t::get_keys() const {
    std::vector<uint64_t> keys;
    keys.reserve(size());

    for (auto const&kv : m_cells) {
        keys.push_back(kv.first);
    }
    return keys;
}

bool signed_interval_domain_t::is_bottom() const {
    if (m_is_bottom) return true;
    return (m_registers.is_bottom() || m_stack.is_bottom());
}

bool signed_interval_domain_t::is_top() const {
    if (m_is_bottom) return false;
    return (m_registers.is_top() && m_stack.is_top());
}

signed_interval_domain_t signed_interval_domain_t::bottom() {
    signed_interval_domain_t interval;
    interval.set_to_bottom();
    return interval;
}

void signed_interval_domain_t::set_to_bottom() {
    m_is_bottom = true;
    m_registers.set_to_bottom();
    m_stack.set_to_bottom();
}

void signed_interval_domain_t::set_registers_to_bottom() {
    m_registers.set_to_bottom();
}

void signed_interval_domain_t::set_registers_to_top() {
    m_registers.set_to_top();
}

void signed_interval_domain_t::set_to_top() {
    m_is_bottom = false;
    m_registers.set_to_top();
    m_stack.set_to_top();
}

std::optional<signed_interval_stack_cell_t> signed_interval_domain_t::find_in_stack(uint64_t key) const {
    return m_stack.find(key);
}

void signed_interval_domain_t::adjust_bb_for_types(location_t loc) {
    m_registers.adjust_bb_for_registers(loc);
}

std::vector<uint64_t> signed_interval_domain_t::get_stack_keys() const {
    return m_stack.get_keys();
}

bool signed_interval_domain_t::all_numeric_in_stack(uint64_t start_loc, int width) const {
    return m_stack.all_numeric(start_loc, width);
}

std::vector<uint64_t> signed_interval_domain_t::find_overlapping_cells_in_stack(uint64_t start_loc,
        int width) const {
    return m_stack.find_overlapping_cells(start_loc, width);
}

void signed_interval_domain_t::remove_overlap_in_stack(const std::vector<uint64_t>& overlap,
        uint64_t start_loc, int width) {
    m_stack.remove_overlap(overlap, start_loc, width);
}

void signed_interval_domain_t::fill_values_in_stack(const std::vector<uint64_t>& overlap,
        uint64_t start_loc, int width) {
    m_stack.fill_values(overlap, start_loc, width);
}

std::optional<refinement_t> signed_interval_domain_t::find_interval_value(register_t reg) const {
    return m_registers.find(reg);
}

std::optional<refinement_t> signed_interval_domain_t::find_interval_at_loc(
        const register_location_t reg) const {
    return m_registers.find(reg);
}

void signed_interval_domain_t::insert_in_registers(register_t reg, location_t loc,
        interval_t interval) {
    refinement_t rf = refinement_t::numeric_refinement(interval, m_slacks);
    m_registers.insert(reg, loc, rf);
}

void signed_interval_domain_t::insert_in_registers(register_t reg, location_t loc,
        refinement_t rf) {
    m_registers.insert(reg, loc, rf);
}

void signed_interval_domain_t::store_in_stack(uint64_t key, refinement_t interval, int width) {
    m_stack.store(key, interval, width);
}

void signed_interval_domain_t::store_in_stack(uint64_t key, interval_t interval, int width) {
    refinement_t rf = refinement_t::numeric_refinement(interval, m_slacks);
    m_stack.store(key, rf, width);
}

bool signed_interval_domain_t::operator<=(const signed_interval_domain_t& abs) const {
    return (m_registers.inclusion(abs.m_registers, m_slacks) &&
            m_stack.inclusion(abs.m_stack, m_slacks));
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
    return signed_interval_domain_t(
        m_registers.join(other.m_registers, m_slacks),
        m_stack.join(other.m_stack, m_slacks),
        m_slacks);
}

signed_interval_domain_t signed_interval_domain_t::operator|(signed_interval_domain_t&& other) const {
    if (is_bottom() || other.is_top()) {
        return std::move(other);
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    // TODO: check if we can use move semantics for m_slacks
    return signed_interval_domain_t(
        m_registers.join(std::move(other.m_registers), m_slacks),
        m_stack.join(std::move(other.m_stack), m_slacks),
        std::move(other.m_slacks));
}

signed_interval_domain_t signed_interval_domain_t::operator&(const signed_interval_domain_t& abs) const {
    /* WARNING: The operation is not implemented yet.*/
    return abs;
}

signed_interval_domain_t signed_interval_domain_t::widen(const signed_interval_domain_t& other, bool to_constants) {
    if (is_bottom() || other.is_top()) {
        return other;
    }
    else if (other.is_bottom() || is_top()) {
        return *this;
    }
    return signed_interval_domain_t(m_registers.widen(other.m_registers, m_slacks),
            m_stack.widen(other.m_stack, m_slacks), m_slacks);
}

signed_interval_domain_t signed_interval_domain_t::narrow(const signed_interval_domain_t& other) const {
    /* WARNING: The operation is not implemented yet.*/
    return other;
}

crab::bound_t signed_interval_domain_t::get_loop_count_upper_bound() const {
    /* WARNING: The operation is not implemented yet.*/
    return crab::bound_t{crab::number_t{0}};
}

void signed_interval_domain_t::initialize_loop_counter(const label_t& label) {
    /* WARNING: The operation is not implemented yet.*/
}

string_invariant signed_interval_domain_t::to_set() {
    return string_invariant{};
}

signed_interval_domain_t signed_interval_domain_t::setup_entry(std::shared_ptr<slacks_t> slacks) {
    return signed_interval_domain_t{signed_interval_registers_t{},
        signed_interval_stack_t::top(), slacks};
}

// Simple truncation function usable with swap_endianness().
template <class T>
constexpr T truncate(T x) noexcept {
    return x;
}

void signed_interval_domain_t::operator()(const Un& u, location_t loc) {
    // numeric_refinement_top() creates interval_t::top()
    auto top_rf = refinement_t::numeric_refinement_top(m_slacks);
    auto swap_endianness = [&](interval_t& v, auto be_or_le) {
        if (const auto n = v.singleton()) {
            if (n->fits_cast_to<int64_t>()) {
                auto interval = interval_t{be_or_le(n->cast_to<int64_t>())};
                refinement_t rf = refinement_t::numeric_refinement(interval, m_slacks);
                m_registers.insert(u.dst.v, loc, rf);
                return;
            }
        }
        m_registers.insert(u.dst.v, loc, top_rf);
    };

    auto rf_opt = m_registers.find(u.dst.v);
    if (!rf_opt) return;
    auto interval = rf_opt->get_interval_value(m_slacks);
    if (interval.is_bottom()) {
        m_registers.insert(u.dst.v, loc, top_rf);
        return;
    }

    // Swap bytes.  For 64-bit types we need the weights to fit in a
    // signed int64, but for smaller types we don't want sign extension,
    // so we use unsigned which still fits in a signed int64.
    switch (u.op) {
    case Un::Op::BE16:
        if (!thread_local_options.big_endian) {
            swap_endianness(interval, boost::endian::endian_reverse<uint16_t>);
        } else {
            swap_endianness(interval, truncate<uint16_t>);
        }
        break;
    case Un::Op::BE32:
        if (!thread_local_options.big_endian) {
            swap_endianness(interval, boost::endian::endian_reverse<uint32_t>);
        } else {
            swap_endianness(interval, truncate<uint32_t>);
        }
        break;
    case Un::Op::BE64:
        if (!thread_local_options.big_endian) {
            swap_endianness(interval, boost::endian::endian_reverse<uint64_t>);
        }
        break;
    case Un::Op::LE16:
        if (thread_local_options.big_endian) {
            swap_endianness(interval, boost::endian::endian_reverse<uint16_t>);
        } else {
            swap_endianness(interval, truncate<uint16_t>);
        }
        break;
    case Un::Op::LE32:
        if (thread_local_options.big_endian) {
            swap_endianness(interval, boost::endian::endian_reverse<uint32_t>);
        } else {
            swap_endianness(interval, truncate<uint32_t>);
        }
        break;
    case Un::Op::LE64:
        if (thread_local_options.big_endian) {
            swap_endianness(interval, boost::endian::endian_reverse<uint64_t>);
        }
        break;
    case Un::Op::SWAP16:
        swap_endianness(interval, boost::endian::endian_reverse<uint16_t>);
        break;
    case Un::Op::SWAP32:
        swap_endianness(interval, boost::endian::endian_reverse<uint32_t>);
        break;
    case Un::Op::SWAP64:
        swap_endianness(interval, boost::endian::endian_reverse<uint64_t>);
        break;
    default: // Un::Op::NEG
        break;
    }
}

void signed_interval_domain_t::operator()(const LoadMapFd& u, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const LoadVariable& u, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Packet& u, location_t loc) {
    // nothing to do here
}

void signed_interval_domain_t::operator()(const Assume& s, location_t loc) {
    // nothing to do here
}

bool signed_interval_domain_t::load_from_stack(register_t reg, uint64_t offset, location_t loc) {
    if (auto loaded = m_stack.find(offset)) {
        m_registers.insert(reg, loc, loaded->first);
        return true;
    }
    return false;
}

void signed_interval_domain_t::store_in_stack(const Mem& b, uint64_t offset, int width) {
    if (auto target_reg = std::get_if<Reg>(&b.value)) {
        if (auto rf_opt = m_registers.find(target_reg->v)) {
            m_stack.store(offset, *rf_opt, width);
        }
    }
    else {
        auto imm = to_signed(std::get<Imm>(b.value).v);
        auto rf = refinement_t::numeric_refinement(interval_t{imm}, m_slacks);
        m_stack.store(offset, rf, width);
    }
}

void signed_interval_domain_t::check_valid_access(const ValidAccess& s, interval_t interval,
        location_t loc) {
    std::string loc_str = loc.to_string();
    bool is_comparison_check = s.width == (Value)Imm{0};
    if (!is_comparison_check) {
        if (s.or_null) {
            if (auto singleton = interval.singleton()) {
                if (*singleton == number_t{0}) return;
            }
            m_errors.push_back(loc_str + ": Non-null number");
        }
        else {
            m_errors.push_back(loc_str + ": Only pointers can be dereferenced");
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

void signed_interval_domain_t::operator()(const basic_block_t& bb) {
    // nothing to do here
}

void signed_interval_domain_t::set_require_check(check_require_func_t f) {}

} // namespace crab
