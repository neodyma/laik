
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

#include <laik-internal.h>
#include <laik/topology.h>

#define mat_elm(r, c, n, arr) (arr[r * n + c])
#define mat_row(r, n, arr)    (&arr[r * n])

// generate a CommMatrix from SwitchStat information
Laik_CommMatrix* laik_top_CommMatrix_from_SwitchStat(Laik_SwitchStat* ss) {}

// allocate a new CommMatrix
Laik_CommMatrix* laik_top_CommMatrix_init(size_t nodecount)
{
    Laik_CommMatrix* cm = calloc(1, sizeof(Laik_CommMatrix));

    if (!cm)
        laik_panic("Out of memory allocating CommMatrix object");

    cm->matrix = malloc(nodecount * nodecount);

    if (!cm->matrix)
        laik_panic("Out of memory allocating CommMatrix matrix");

    cm->nodecount = nodecount;
    return cm;
}

// free a CommMatrix
void laik_top_CommMatrix_free(Laik_CommMatrix* cm)
{
    free(cm->matrix);
    free(cm);
}

// add a new transfer to the CM
Laik_CommMatrix* laik_top_CommMatrix_update(Laik_CommMatrix* cm, size_t from, size_t to, int64_t amt)
{
    mat_elm(from, to, cm->nodecount, cm->matrix) += amt;
    mat_elm(to, from, cm->nodecount, cm->matrix) += amt;
    return cm;
}

// swap two nodes of the CM
Laik_CommMatrix* laik_top_CommMatrix_swapnodes(Laik_CommMatrix* cm, size_t from, size_t to)
{
    uint64_t tmp[cm->nodecount];
    memcpy(tmp, mat_row(to, cm->nodecount, cm->matrix), cm->nodecount);
    memcpy(mat_row(to, cm->nodecount, cm->matrix), mat_row(from, cm->nodecount, cm->matrix), cm->nodecount);
    memcpy(mat_row(from, cm->nodecount, cm->matrix), tmp, cm->nodecount);
    return cm;
}
