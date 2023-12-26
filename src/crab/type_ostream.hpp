// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "crab/common.hpp"
#include "crab/offset_domain.hpp"

void print_numeric_register(std::ostream&, Reg, crab::interval_t);
void print_numeric_memory_cell(std::ostream&, int, int, crab::interval_t);
void print_non_numeric_register(std::ostream&, Reg, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::dist_t> = std::nullopt);
void print_non_numeric_memory_cell(std::ostream&, int, int, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::dist_t> = std::nullopt);
void print_register(std::ostream&, Reg, const std::optional<crab::ptr_or_mapfd_t>&,
        std::optional<crab::dist_t>, std::optional<crab::mock_interval_t>);
void print_memory_cell(std::ostream&, int, int, const std::optional<crab::ptr_or_mapfd_t>&,
        std::optional<crab::dist_t>, std::optional<crab::mock_interval_t>);

// Print select transformers
void print_annotated(std::ostream&, const Call&, std::optional<crab::ptr_or_mapfd_t>&,
        std::optional<crab::mock_interval_t>&);
void print_annotated(std::ostream&, const Bin&, std::optional<crab::ptr_or_mapfd_t>&,
        std::optional<crab::dist_t>&, std::optional<crab::mock_interval_t>&);
void print_annotated(std::ostream&, const LoadMapFd&, std::optional<crab::ptr_or_mapfd_t>&);
void print_annotated(std::ostream&, const Mem&, std::optional<crab::ptr_or_mapfd_t>&,
        std::optional<crab::dist_t>&, std::optional<crab::mock_interval_t>&);
void print_annotated(std::ostream&, const Un&, std::optional<crab::mock_interval_t>&);
