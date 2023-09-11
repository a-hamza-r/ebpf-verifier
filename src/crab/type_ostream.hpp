// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "string_constraints.hpp"
#include "asm_syntax.hpp"
#include "asm_ostream.hpp"
#include "crab/common.hpp"
#include "crab/offset_domain.hpp"

void print_ptr_or_mapfd_type(std::ostream&, const crab::ptr_or_mapfd_t&, std::optional<crab::dist_t>);
void print_number(std::ostream&, crab::interval_t);
void print_ptr_type(std::ostream&, const crab::ptr_or_mapfd_t& ptr, std::optional<crab::dist_t>);
void print_register(std::ostream& o, const Reg& r, const std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::dist_t>, std::optional<crab::mock_interval_t>);
void print_annotated(std::ostream& o, const Call& call, std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::mock_interval_t>&);
void print_annotated(std::ostream& o, const Bin& b, std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::dist_t>&, std::optional<crab::mock_interval_t>&);
void print_annotated(std::ostream& o, const LoadMapFd& u, std::optional<crab::ptr_or_mapfd_t>& p);
void print_annotated(std::ostream& o, const Mem& b, std::optional<crab::ptr_or_mapfd_t>& p, std::optional<crab::dist_t>&, std::optional<crab::mock_interval_t>&);
