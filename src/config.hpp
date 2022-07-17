// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#pragma once

enum class abstract_domain_kind { EBPF_DOMAIN, TYPE_DOMAIN, OFFSET_DOMAIN, REGION_DOMAIN, CONST_PROP_DOMAIN };

struct ebpf_verifier_options_t {
    bool check_termination;
    bool assume_assertions;
    bool print_invariants;
    bool print_failures;
    bool no_simplify;
    abstract_domain_kind abstract_domain;

    // False to use actual map fd's, true to use mock fd's.
    bool mock_map_fds;

    // True to do additional checks for some things that would fail at runtime.
    bool strict;

    bool print_line_info;
    bool allow_division_by_zero;
    bool setup_constraints;

    bool dump_btf_types_json;
};

struct ebpf_verifier_stats_t {
    int total_unreachable;
    int total_warnings;
    int max_instruction_count;
};

extern const ebpf_verifier_options_t ebpf_verifier_default_options;
extern thread_local ebpf_verifier_options_t thread_local_options;
