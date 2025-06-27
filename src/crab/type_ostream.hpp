// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: MIT

#pragma once

#include "crab/offset_domain.hpp"

void print_numeric_register(std::ostream&, crab::register_t, crab::refinement_t, bool, 
                            std::shared_ptr<crab::slacks_t>);
void print_numeric_memory_cell(std::ostream&, int, int, crab::refinement_t, bool,
                               std::shared_ptr<crab::slacks_t>, crab::region_t);
void print_non_numeric_register(std::ostream&, crab::register_t, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::refinement_t>, std::shared_ptr<crab::slacks_t>);
void print_non_numeric_memory_cell(std::ostream&, int, int, const crab::ptr_or_mapfd_t& ptr,
        std::optional<crab::refinement_t>, std::shared_ptr<crab::slacks_t>, crab::region_t);
void print_register(std::ostream&, Reg, std::optional<crab::ptr_or_mapfd_t>,
                    std::optional<crab::refinement_t>,
                    std::optional<crab::refinement_t>,
                    bool,
                    std::shared_ptr<crab::slacks_t>);
void print_memory_cell(std::ostream&, int, int, std::optional<crab::ptr_or_mapfd_t>,
                       std::optional<crab::refinement_t>, std::optional<crab::refinement_t>,
                       std::shared_ptr<crab::slacks_t>, crab::region_t);

// Print select transformers
void print_annotated(std::ostream&, const Bin&,
                     std::optional<crab::ptr_or_mapfd_t>,
                     std::optional<crab::refinement_t>,
                     std::optional<crab::refinement_t>,
                     std::optional<crab::refinement_t>,
                     std::shared_ptr<crab::slacks_t>);
void print_annotated(std::ostream&, const LoadMapFd&,
                     std::optional<crab::ptr_or_mapfd_t>);
void print_annotated(std::ostream&, const Mem&,
                     std::optional<crab::ptr_or_mapfd_t>,
                     std::optional<crab::refinement_t>,
                     std::optional<crab::refinement_t>,
                     std::optional<crab::refinement_t>,
                     std::shared_ptr<crab::slacks_t>);
void print_annotated(std::ostream&, const Un&,
                     std::optional<crab::refinement_t>,
                     std::optional<crab::refinement_t>,
                     std::shared_ptr<crab::slacks_t>);
void print_bb(std::ostream&, const basic_block_t&);
void print_instr(std::ostream&, const Instruction&);
