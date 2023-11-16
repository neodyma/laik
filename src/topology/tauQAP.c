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

int*     top_QAP_construction(Laik_CommMatrix*, Laik_Topology_Matrix*);
int*     top_QAP_improvement(Laik_CommMatrix*, Laik_Topology_Matrix*);
int*     top_QAP_cyclicSearch(Laik_CommMatrix*, Laik_Topology_Matrix*, int*);
int*     top_QAP_pairxchg(int*, size_t, size_t);
uint64_t top_QAP_totalCost(Laik_CommMatrix*, Laik_Topology_Matrix*, int*);
uint64_t top_QAP_commLoad(Laik_CommMatrix*, size_t, size_t*);
uint64_t top_QAP_coreDist(Laik_Topology_Matrix*, size_t, size_t*);

// optimize given pattern for topology using QAP
// @param cm communication matrix @param top topology struct
// @returns reordering
int* laik_top_do_reorder_QAP(Laik_CommMatrix* cm, Laik_Topology* top)
{
    if (top->which != LAIK_TOP_IS_MAT) return NULL;
    return top_QAP_improvement(cm, top->data.mat);
}

// do the QAP construction method
int* top_QAP_construction(Laik_CommMatrix* cm, Laik_Topology_Matrix* top) {}

// do the QAP improvement runs
int* top_QAP_improvement(Laik_CommMatrix* cm, Laik_Topology_Matrix* top)
{
    return top_QAP_cyclicSearch(cm, top, top_QAP_construction(cm, top));
}

// do cyclic search around given order
// @param cm communication matrix @param top topology matrix
// @param initial initial ordering
// @returns initial, overwritten with best solution
int* top_QAP_cyclicSearch(Laik_CommMatrix* cm, Laik_Topology_Matrix* top, int* initial)
{
    const size_t n = cm->nodecount;

    int*     best_sol = initial;
    uint64_t best_cost = top_QAP_totalCost(cm, top, initial);
    int      current_sol[n];
    uint64_t current_cost = best_cost;
    memcpy(current_sol, best_sol, n * sizeof(int));

    size_t i = 0, j = 1;
    for (size_t k = 0; k < n * n; k++) {
        top_QAP_pairxchg(current_sol, i, j);
        current_cost = top_QAP_totalCost(cm, top, current_sol);
        if (current_cost < best_cost) {
            // always have the current best solution in memory for early exit
            memcpy(best_sol, current_sol, n * sizeof(int));
            best_cost = current_cost;
        }
        if (j < n - 1)
            j += 1;
        else if (j == n && i < n - 2) {
            i += 1;
            j = i + 2;
        }
        else {
            i = 1;
            j = 2;
        }
    }

    return best_sol;
}

// exchange two list elements
int* top_QAP_pairxchg(int* order, size_t i, size_t j)
{
    int tmp = order[i];
    order[i] = order[j];
    order[j] = tmp;

    return order;
}

// get total weight for given reordering
uint64_t top_QAP_totalCost(Laik_CommMatrix* cm, Laik_Topology_Matrix* top, int* order)
{
    uint64_t cost = 0;
    for (size_t i = 0; i < cm->nodecount; i++)
        for (size_t j = 0; j < cm->nodecount; j++) cost += top_mat_elm(order[i], order[j], cm) * top_mat_elm(i, j, top);
    return cost;
}

int cmpfunc(const void* a, const void* b)
{
    const size_t* a_ptr = a;
    const size_t* b_ptr = b;
    return *a_ptr == *b_ptr;
}

// calculate total communication load between process and already assigned processes
// @param assigned should always be sorted
uint64_t top_QAP_commLoad(Laik_CommMatrix* mat, size_t process, size_t* assigned)
{
    uint64_t load = 0;
    for (size_t i = 0; i < mat->nodecount; i++)
        if (i != process && bsearch(&i, assigned, mat->nodecount, sizeof(size_t), cmpfunc) != NULL) {
            load += top_mat_elm(process, i, mat);
            load += top_mat_elm(i, process, mat);
        }
    return load;
}

// calculate total distance between node and already assigned nodes
// ensure @param assigned is always sorted
uint64_t top_QAP_coreDist(Laik_Topology_Matrix* top, size_t node, size_t* assigned)
{  // maybe use https://stackoverflow.com/a/25689059
    uint64_t dist = 0;
    for (size_t i = 0; i < top->nodecount; i++)
        if (i != node && bsearch(&i, assigned, top->nodecount, sizeof(size_t), cmpfunc) != NULL)
            dist += top_mat_elm(node, i, top);
    return dist;
}
