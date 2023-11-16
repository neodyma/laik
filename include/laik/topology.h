/*
 * This file is part of the LAIK library.
 * Copyright (c) 2023 Lukas Heine <lu.heine@tum.de>
 *
 * LAIK is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, version 3 or later.
 *
 * LAIK is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LAIK_TOPOLOGY_H
#define LAIK_TOPOLOGY_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define top_mat_elm(r, c, mat) ((mat->matrix)[r * mat->nodecount + c])
#define top_mat_row(r, mat) (&(mat->matrix)[r * mat->nodecount])

// =============================================================================
// LAIK_TOPOLOGY BASE
// =============================================================================

typedef struct _Laik_CommMatrix {
    Laik_Instance* inst;
    size_t         nodecount;
    uint64_t*      matrix;
} Laik_CommMatrix;

typedef struct _Laik_Topology_Matrix {
    Laik_Instance* inst;
    size_t         nodecount;
    uint64_t*      matrix;
} Laik_Topology_Matrix;

typedef struct _Laik_Topology_Graph {
    uint64_t size;
} Laik_Topology_Graph;

enum _Laik_Topology_Which {
    LAIK_TOP_IS_MAT,
    LAIK_TOP_IS_GRAPH,
};

typedef struct _Laik_Topology {
    uint8_t which;
    union {
        Laik_Topology_Matrix* mat;
        Laik_Topology_Graph*  graph;
    } data;
} Laik_Topology;

enum _Laik_Reorder_Mapped {
    LAIK_RO_UNMAPPED = 0,
    LAIK_RO_OFFSET = 1,
};

typedef struct _Laik_Reordering_File {
    uint32_t nodecount;
    int      reordering[];  // use a VLA
} __attribute__((packed)) Laik_Reordering_File;

Laik_CommMatrix* laik_top_CommMatrix_from_SwitchStat(Laik_SwitchStat* ss);

Laik_CommMatrix* laik_top_CommMatrix_init(Laik_Instance* li);
void             laik_top_CommMatrix_free(Laik_CommMatrix* cm);
Laik_CommMatrix* laik_top_CommMatrix_update(Laik_CommMatrix* cm, size_t from, size_t to, int64_t amt);
void             laik_top_CommMatrix_sync(Laik_CommMatrix* cm);
Laik_CommMatrix* laik_top_CommMatrix_swapnodes(Laik_CommMatrix* cm, size_t from, size_t to);

int*        laik_top_reordering(Laik_Instance* li);
int*        laik_top_reordering_get(Laik_Instance* li);
Laik_Group* laik_allow_reordering(Laik_Instance* li, int phase);
int*        laik_top_do_reorder(Laik_CommMatrix* cm, Laik_Topology* top);

// =============================================================================
// LAIK_TOPOLOGY tauQAP
// =============================================================================

int* laik_top_do_reorder_QAP(Laik_CommMatrix* cm, Laik_Topology* top);

// todo virtual reordering func

#endif
