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
#include <stdlib.h>
#include <string.h>

typedef struct _Laik_CommMatrix
{
    uint64_t* matrix;
    size_t    nodecount;
} Laik_CommMatrix;

Laik_CommMatrix* laik_top_CommMatrix_from_SwitchStat(Laik_SwitchStat* ss);

Laik_CommMatrix* laik_top_CommMatrix_init(size_t nodecount);
void             laik_top_CommMatrix_free(Laik_CommMatrix* cm);
Laik_CommMatrix* laik_top_CommMatrix_update(Laik_CommMatrix* cm, size_t from, size_t to, int64_t amt);
Laik_CommMatrix* laik_top_CommMatrix_swapnodes(Laik_CommMatrix* cm, size_t from, size_t to);

#endif
