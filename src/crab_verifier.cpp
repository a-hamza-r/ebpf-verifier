// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
/**
 *  This module is about selecting the numerical and memory domains, initiating
 *  the verification process and returning the results.
 **/
#include <cinttypes>

#include <ctime>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "crab/abstract_domain.hpp"
#include "crab/ebpf_domain.hpp"
#include "crab/type_domain.hpp"
#include "crab/region_domain.hpp"
#include "crab/fwd_analyzer.hpp"
#include "crab_utils/lazy_allocator.hpp"

#include "asm_syntax.hpp"
#include "crab_verifier.hpp"
#include "string_constraints.hpp"

using std::string;
using crab::ebpf_domain_t;
using crab::linear_constraint_t;

thread_local crab::lazy_allocator<program_info> global_program_info;
thread_local ebpf_verifier_options_t thread_local_options;

static checks_db generate_report(cfg_t& cfg,
                                 crab::invariant_table_t& pre_invariants,
                                 crab::invariant_table_t& post_invariants) {
    checks_db m_db;
    for (const label_t& label : cfg.sorted_labels()) {
        basic_block_t& bb = cfg.get_node(label);
        abstract_domain_t from_inv(pre_invariants.at(label));

        from_inv.set_require_check([&m_db, label](auto& inv, const linear_constraint_t& cst, const std::string& s) {
            if (inv.is_bottom())
                return true;
            if (cst.is_contradiction()) {
                m_db.add_warning(label, s);
                return false;
            }

            if (inv.entail(cst)) {
                // add_redundant(s);
                return true;
            } else if (inv.intersect(cst)) {
                // TODO: add_error() if imply negation
                m_db.add_warning(label, s);
                return false;
            } else {
                m_db.add_warning(label, s);
                return false;
            }
        });

        if (thread_local_options.check_termination) {
            // Pinpoint the places where divergence might occur.
            crab::bound_t min_instruction_count_upper_bound{crab::bound_t::plus_infinity()};
            for (const label_t& prev_label : bb.prev_blocks_set()) {
                crab::bound_t instruction_count = pre_invariants.at(prev_label).get_instruction_count_upper_bound();
                min_instruction_count_upper_bound = std::min(min_instruction_count_upper_bound, instruction_count);
            }

            crab::bound_t max_instructions{100000};
            crab::bound_t instruction_count_upper_bound = from_inv.get_instruction_count_upper_bound();
            if ((min_instruction_count_upper_bound < max_instructions) &&
                (instruction_count_upper_bound >= max_instructions))
                m_db.add_nontermination(label);

            m_db.max_instruction_count = std::max(m_db.max_instruction_count, instruction_count_upper_bound);
        }

        bool pre_bot = from_inv.is_bottom();

        from_inv(bb, thread_local_options.check_termination);

        if (!pre_bot && from_inv.is_bottom()) {
            m_db.add_unreachable(label, std::string("Code is unreachable after ") + to_string(bb.label()));
        }
    }
    return m_db;
}

static auto get_line_info(const InstructionSeq& insts) {
    std::map<int, btf_line_info_t> label_to_line_info;
    for (auto& [label, inst, line_info] : insts) {
        if (line_info.has_value())
            label_to_line_info.emplace(label.from, line_info.value());
    }
    return label_to_line_info;
}

static void print_report(std::ostream& os, const checks_db& db, const InstructionSeq& prog, bool print_line_info) {
    auto label_to_line_info = get_line_info(prog);
    os << "\n";
    for (auto [label, messages] : db.m_db) {
        for (const auto& msg : messages) {
            if (print_line_info) {
                auto line_info = label_to_line_info.find(label.from);
                if (line_info != label_to_line_info.end())
                    os << line_info->second;
            }
            os << label << ": " << msg << "\n";
        }
    }
    os << "\n";
    if (!db.maybe_nonterminating.empty()) {
        os << "Could not prove termination on join into: ";
        for (const label_t& label : db.maybe_nonterminating) {
            os << label << ", ";
        }
        os << "\n";
    }
}

static checks_db get_analysis_report(std::ostream& s, cfg_t& cfg, crab::invariant_table_t& pre_invariants,
                                     crab::invariant_table_t& post_invariants) {
    // Analyze the control-flow graph.
    checks_db db = generate_report(cfg, pre_invariants, post_invariants);
    if (thread_local_options.abstract_domain == abstract_domain_kind::TYPE_DOMAIN
            || thread_local_options.abstract_domain == abstract_domain_kind::REGION_DOMAIN) {
        auto state = post_invariants.at(label_t::exit);
        for (const label_t& label : cfg.sorted_labels()) {
            state(cfg.get_node(label), 0, thread_local_options.print_invariants ? 2 : 1);
        }
    }
    else if (thread_local_options.print_invariants) {
        for (const label_t& label : cfg.sorted_labels()) {
            s << "\nPre-invariant : " << pre_invariants.at(label) << "\n";
            s << cfg.get_node(label);
            s << "\nPost-invariant: " << post_invariants.at(label) << "\n";
        }
    }
    return db;
}

/* EXTEND FOR NEW DOMAINS */
static abstract_domain_t make_initial(const ebpf_verifier_options_t* options) {
    switch (options->abstract_domain) {
    case abstract_domain_kind::EBPF_DOMAIN: {
        ebpf_domain_t entry_inv = ebpf_domain_t::setup_entry(options->check_termination, true);
        return abstract_domain_t(entry_inv);
    }
    case abstract_domain_kind::REGION_DOMAIN: {
        region_domain_t entry_inv = region_domain_t::setup_entry();
        return abstract_domain_t(entry_inv);
    }
    case abstract_domain_kind::TYPE_DOMAIN: {
        type_domain_t entry_inv = type_domain_t::setup_entry();
        return abstract_domain_t(entry_inv);
    }
    default:
        // FIXME: supported abstract domains should be checked in check.cpp
        std::cerr << "error: unsupported abstract domain\n";
        std::exit(1);
    }
}

/* EXTEND FOR NEW DOMAINS */
static abstract_domain_t make_initial(abstract_domain_kind abstract_domain, const string_invariant& entry_invariant, bool setup_constraints) {

    switch (abstract_domain) {
    case abstract_domain_kind::EBPF_DOMAIN: {
        ebpf_domain_t entry_inv = entry_invariant.is_bottom()
                                      ? ebpf_domain_t::from_constraints({"false"}, setup_constraints)
                                      : ebpf_domain_t::from_constraints(entry_invariant.value(), setup_constraints);
        return abstract_domain_t(entry_inv);
    }
    case abstract_domain_kind::TYPE_DOMAIN: {
        // TODO
    }
    default:
        // FIXME: supported abstract domains should be checked in check.cpp
        std::cerr << "error: unsupported abstract domain\n";
        std::exit(1);
    }
}

crab_results get_ebpf_report(std::ostream& s, cfg_t& cfg, program_info info, const ebpf_verifier_options_t* options) {
    global_program_info = std::move(info);
    crab::domains::clear_global_state();
    crab::variable_t::clear_thread_local_state();
    thread_local_options = *options;

    try {

        abstract_domain_t entry_dom = make_initial(options);
        // Get dictionaries of pre-invariants and post-invariants for each basic block.
        auto [pre_invariants, post_invariants] =
            crab::run_forward_analyzer(cfg, std::move(entry_dom), options->check_termination);
        return crab_results(std::move(cfg),
			    std::move(pre_invariants), std::move(post_invariants),
			    std::move(get_analysis_report(s, cfg, pre_invariants, post_invariants)));
    } catch (std::runtime_error& e) {
        // Convert verifier runtime_error exceptions to failure.
        checks_db db;
        db.add_warning(label_t::exit, e.what());
        crab::invariant_table_t pre_invariants, post_invariants;
        return crab_results(std::move(cfg),
			    std::move(pre_invariants), std::move(post_invariants),
			    std::move(db));
    }
}

/// Returned value is true if the program passes verification.
bool run_ebpf_analysis(std::ostream& s, cfg_t& cfg, const program_info& info, const ebpf_verifier_options_t* options,
                       ebpf_verifier_stats_t* stats) {
    if (options == nullptr)
        options = &ebpf_verifier_default_options;
    checks_db report = get_ebpf_report(s, cfg, info, options).db;
    if (stats) {
        stats->total_unreachable = report.total_unreachable;
        stats->total_warnings = report.total_warnings;
        stats->max_instruction_count = report.get_max_instruction_count();
    }
    return (report.total_warnings == 0);
}

static string_invariant_map to_string_invariant_map(crab::invariant_table_t& inv_table) {
    string_invariant_map res;
    for (auto& [label, inv] : inv_table) {
        res.insert_or_assign(label, inv.to_set());
    }
    return res;
}

std::tuple<string_invariant, bool>
ebpf_analyze_program_for_test(std::ostream& os, const InstructionSeq& prog, const string_invariant& entry_invariant,
                              const program_info& info, const ebpf_verifier_options_t& options) {
    thread_local_options = options;
    global_program_info = info;
    assert(!entry_invariant.is_bottom());
    abstract_domain_t entry_inv = make_initial(abstract_domain_kind::EBPF_DOMAIN, entry_invariant, options.setup_constraints);
    if (entry_inv.is_bottom())
        throw std::runtime_error("Entry invariant is inconsistent");
    cfg_t cfg = prepare_cfg(prog, info, !options.no_simplify, false);
    auto [pre_invariants, post_invariants] = crab::run_forward_analyzer(cfg, entry_inv, options.check_termination);
    checks_db report = get_analysis_report(std::cerr, cfg, pre_invariants, post_invariants);
    print_report(os, report, prog, false);

    auto pre_invariant_map = to_string_invariant_map(pre_invariants);

    return {
        pre_invariant_map.at(label_t::exit),
        (report.total_warnings == 0)
    };
}

/// Returned value is true if the program passes verification.
crab_results ebpf_verify_program(std::ostream& os, const InstructionSeq& prog, const program_info& info,
                                 const ebpf_verifier_options_t* options, ebpf_verifier_stats_t* stats) {
    if (options == nullptr)
        options = &ebpf_verifier_default_options;

    // Convert the instruction sequence to a control-flow graph
    // in a "passive", non-deterministic form.
    cfg_t cfg = prepare_cfg(prog, info, !options->no_simplify);

    crab_results results = get_ebpf_report(os, cfg, info, options);
    checks_db& report = results.db;
    if (options->print_failures) {
        print_report(os, report, prog, options->print_line_info);
    }
    if (stats) {
        stats->total_unreachable = report.total_unreachable;
        stats->total_warnings = report.total_warnings;
        stats->max_instruction_count = report.get_max_instruction_count();
    }
    // return (report.total_warnings == 0);
    return results;
}

void ebpf_verifier_clear_thread_local_state() {
    crab::variable_t::clear_thread_local_state();
    crab::CrabStats::clear_thread_local_state();
    global_program_info.clear();
    crab::domains::clear_thread_local_state();
    crab::domains::SplitDBM::clear_thread_local_state();
}
