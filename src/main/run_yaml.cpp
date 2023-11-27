// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#include <iostream>

#include "CLI11.hpp"

#include "ebpf_verifier.hpp"
#include "ebpf_yaml.hpp"

int main(int argc, char** argv) {
    CLI::App app{"Run YAML test cases"};

    std::string filename;
    app.add_option("path", filename, "YAML file.")->required()->type_name("FILE");

    std::string pattern;
    app.add_option("pattern", pattern, "Pattern for test cases to run (substring)")->type_name("PATTERN");

    bool verbose = false;
    app.add_flag("-v", verbose, "Verbose");

    bool quiet = false;
    app.add_flag("-q", quiet, "Never print code");

    CLI11_PARSE(app, argc, argv);
    bool res = true;
    bool tests_for_type_d = filename.find("types") != filename.npos;
    auto domain = tests_for_type_d ? abstract_domain_kind::TYPE_DOMAIN :
        abstract_domain_kind::EBPF_DOMAIN;
    foreach_suite(filename, [&](const TestCase& test_case) {
        if (!pattern.empty() && test_case.name.find(pattern) == test_case.name.npos) {
            return;
        }
        std::cout << test_case.name << ": " << std::flush;
        const auto& maybe_failure = run_yaml_test_case(domain, test_case, verbose);
        if (!quiet && (verbose || maybe_failure)) {
            std::cout << "\n";
            std::cout << "Pre-invariant:" << test_case.assumed_pre_invariant << "\n";
            print(test_case.instruction_seq, std::cout, {});
            std::cout << "Expected post-invariant: " << test_case.expected_post_invariant << "\n";
        }
        if (maybe_failure) {
            std::cout << "failed:\n";
            res = false;
            std::cout << "------\n";
            print_failure(*maybe_failure, std::cout);
            std::cout << "------\n";
        } else {
            std::cout << "pass\n";
        }
    });
    return !res;
}
