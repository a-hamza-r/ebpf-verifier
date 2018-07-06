/*
 * Copyright 2018 VMware, Inc
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <inttypes.h>
#include <sys/mman.h>
#include <assert.h>

#include "ubpf_int.h"
#include "abs_dom.h"
#include "abs_state.h"

void 
abs_initialize_entry(struct abs_state *state)
{
    for (int i = 0; i < 11; i++) {
        state->reg[i] = abs_dom_top;
    }
    state->reg[1] = abs_dom_ctx;
    state->reg[10] = abs_dom_stack;
    state->bot = false;
}

void 
abs_initialize_unreached(struct abs_state *state)
{
    state->bot = true;
}

void
abs_join(struct abs_state *state, struct abs_state other)
{
    if (state->bot) {
        *state = other;
    } else {
        for (int r = 0; r < 11; r++) {
            state->reg[r] = abs_dom_join(state->reg[r], other.reg[r]);
        } 
    }
}

static void
print_state(const char* header, struct abs_state *state, uint16_t pc)
{
    FILE *f = fopen("log.txt", "a");
    fprintf(f, "(%15s) %2u: ", header, pc);
    if (state->bot) {
        fprintf(f, "BOT");
    } else for (int r = 0; r < 11; r++) {
	    fprintf(f, "r%d=", r);
        abs_dom_print(f, state->reg[r]);
	    fprintf(f, "; ");
    }
    fprintf(f, "\n");
    fclose(f);
}

static int
access_width(uint8_t opcode) {
    switch (opcode) {
    case EBPF_OP_LDXB:
    case EBPF_OP_STB:
    case EBPF_OP_STXB: return 1;
    case EBPF_OP_LDXH:
    case EBPF_OP_STH:
    case EBPF_OP_STXH: return 2;
    case EBPF_OP_LDXW:
    case EBPF_OP_STW:
    case EBPF_OP_STXW: return 4;
    case EBPF_OP_LDXDW:
    case EBPF_OP_STDW:
    case EBPF_OP_STXDW: return 8;
    default: return -1;
    }
}

static bool
is_alu(uint8_t opcode)
{
    return (opcode & EBPF_CLS_MASK) == EBPF_CLS_ALU
        || (opcode & EBPF_CLS_MASK) == EBPF_CLS_ALU64;
}

void
abs_execute(struct abs_state *to, struct abs_state *from, struct ebpf_inst inst,
            int32_t imm, bool taken, uint16_t pc, char** errmsg)
{
    print_state("before command", from, pc);
    // (it can be an optimization for a specific domain, with default implementation as execute then join)
    struct abs_state state = *from;
    if (inst.opcode == EBPF_OP_EXIT) {
        if (!abs_dom_is_initialized(state.reg[0])) {
            *errmsg = ubpf_error("r0 must be initialized to a numerical value before exit");
	    return;
        }
    }
    if (inst.opcode == EBPF_OP_CALL) {
        state.reg[0] = abs_dom_call(inst, state.reg[1], state.reg[2], state.reg[3], state.reg[4], state.reg[5]);
        for (int r = 1; r <= 5; r++) {
            state.reg[r] = abs_dom_top;
        }
    } else if ((inst.opcode & EBPF_CLS_MASK) == EBPF_CLS_JMP) {
        struct abs_dom_value dummy = abs_dom_fromconst(inst.imm);
        // TODO: check feasibility; this might cause problems with pending.
        struct abs_dom_value *v2 = (inst.opcode | EBPF_SRC_IMM) ? &dummy : &state.reg[inst.src];
        abs_dom_assume(inst.opcode, taken, &state.reg[inst.dst], v2);
        if (abs_dom_is_bot(state.reg[inst.dst]) || abs_dom_is_bot(*v2)) {
            state.bot = true;
        }
    } else if (inst.opcode == EBPF_OP_LDDW) {
        state.reg[inst.dst] = abs_dom_fromconst((uint32_t)inst.imm | ((uint64_t)imm << 32));
    } else if (is_alu(inst.opcode)) {
        bool div = (inst.opcode == EBPF_OP_DIV_REG || inst.opcode == EBPF_OP_DIV64_REG);
        bool mod = (inst.opcode == EBPF_OP_MOD_REG || inst.opcode == EBPF_OP_MOD64_REG);
        if (div || mod) {
            bool is64 = (inst.opcode & EBPF_CLS_MASK) == EBPF_CLS_ALU64;
            if (abs_dom_maybe_zero(state.reg[inst.src], is64)) {
                *errmsg = ubpf_error("division by zero at PC %d", pc);
                return;
            }
        }
        state.reg[inst.dst] = abs_dom_alu(inst.opcode, inst.imm, state.reg[inst.dst], state.reg[inst.src]);
    } else {
        int width = access_width(inst.opcode);
        assert(width > 0);
        bool is_load = ((inst.opcode & EBPF_CLS_MASK) == EBPF_CLS_LD)
                    || ((inst.opcode & EBPF_CLS_MASK) == EBPF_CLS_LDX);

        uint8_t r = is_load ? inst.src : inst.dst;
        if (abs_dom_out_of_bounds(state.reg[r], inst.offset, width)) {
            *errmsg = ubpf_error("out of bounds memory %s at PC %d [r%d%+d]",
                                is_load ? "load" : "store", pc, is_load ? inst.src : inst.dst, inst.offset);
            return;
        }
        if (is_load) {
            state.reg[inst.dst] = abs_dom_unknown;
        }
    }

    print_state("before join", &state, pc);
    print_state("joining with", to, pc);
    bool bot = to->bot;

    abs_join(to, state);

    if (!bot)
        print_state("after join", to, pc);
    int _ = system("echo >> log.txt"); taken = (bool)_;
}
