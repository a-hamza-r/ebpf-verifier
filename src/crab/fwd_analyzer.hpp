// Copyright (c) Prevail Verifier contributors.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <map>

#include "config.hpp"
#include "crab/abstract_domain.hpp"

namespace crab {

using invariant_table_t = std::map<label_t, abstract_domain_t>;

std::pair<invariant_table_t, invariant_table_t> run_forward_analyzer(cfg_t& cfg, abstract_domain_t entry_inv);

} // namespace crab
