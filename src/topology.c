
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

// generate a CommMatrix from SwitchStat information
Laik_CommMatrix* laik_top_CommMatrix_from_SwitchStat(Laik_SwitchStat* ss) {}

// allocate a new CommMatrix
Laik_CommMatrix* laik_top_CommMatrix_init(Laik_Instance* li)
{
    Laik_CommMatrix* cm = calloc(1, sizeof(Laik_CommMatrix));

    if (!cm)
        laik_panic("Out of memory allocating CommMatrix object");

    size_t nodecount = li->locations;
    cm->matrix = calloc((nodecount * nodecount), sizeof(*cm->matrix));

    if (!cm->matrix)
        laik_panic("Out of memory allocating CommMatrix matrix");

    cm->nodecount = nodecount;
    cm->inst      = li;
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
    top_mat_elm(from, to, cm) += amt;
    top_mat_elm(to, from, cm) += amt;
    return cm;
}

void laik_top_CommMatrix_sync(Laik_CommMatrix* cm)
{
    const Laik_Backend* b = cm->inst->backend;
    if(!b || !b->matsync)
        laik_panic("backend or matrix sync unavailable");

    cm->in_sync = true;
    (b->matsync)(cm);
    cm->in_sync = false;
}

// swap two nodes of the CM
Laik_CommMatrix* laik_top_CommMatrix_swapnodes(Laik_CommMatrix* cm, size_t from, size_t to)
{
    uint64_t tmp[cm->nodecount];
    memcpy(tmp, top_mat_row(to, cm), cm->nodecount);
    memcpy(top_mat_row(to, cm), top_mat_row(from, cm), cm->nodecount);
    memcpy(top_mat_row(from, cm), tmp, cm->nodecount);
    return cm;
}
