// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT
#include <catch2/catch_all.hpp>

#include "ebpf_verifier.hpp"
#include "ebpf_yaml.hpp"

// TODO: move out of this framework

#define YAML_CASE(path, domain) \
    TEST_CASE("YAML suite: " path, "[yaml]") { \
        foreach_suite(path, [&](const TestCase& test_case){ \
            std::optional<Failure> failure = run_yaml_test_case(domain, test_case); \
            if (failure) { \
                std::cout << "test case: " << test_case.name << "\n"; \
                print_failure(*failure, std::cout); \
            } \
            REQUIRE(!failure); \
        }); \
    }

YAML_CASE("test-data/add.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/assign.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/atomic.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/bitop.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/call.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/callx.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/udivmod.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/sdivmod.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/full64.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/jump.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/loop.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/movsx.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/packet.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/parse.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/sext.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/shift.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/stack.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/subtract.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/unop.yaml", abstract_domain_kind::EBPF_DOMAIN)
YAML_CASE("test-data/unsigned.yaml", abstract_domain_kind::EBPF_DOMAIN)
