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
#include <sys/queue.h>

int*     top_QAP_construction(Laik_CommMatrix*, Laik_TopologyMatrix*);
int*     top_QAP_improvement(Laik_CommMatrix*, Laik_TopologyMatrix*);
int*     top_QAP_cyclicSearch(Laik_CommMatrix*, Laik_TopologyMatrix*, int*);
int*     top_QAP_pairxchg(int*, size_t, size_t);
uint64_t top_QAP_totalCost(Laik_CommMatrix*, Laik_TopologyMatrix*, int*);
uint64_t top_QAP_commLoad(Laik_CommMatrix*, size_t, size_t, int*);
uint64_t top_QAP_coreDist(Laik_TopologyMatrix*, size_t, size_t, int*);
int*     top_QAP_list_assign(int*, size_t, size_t*, size_t);

// optimize given pattern for topology using QAP
// @param cm communication matrix @param top topology struct
// @returns reordering
int* laik_top_do_reorder_QAP(Laik_CommMatrix* cm, Laik_Topology* top)
{
    if (top->which != LAIK_TOP_IS_MAT) return NULL;
    return top_QAP_improvement(cm, top->data.mat);
}

// sort Laik_Top_IndexedElement
int top_func_sortIndexed(const void* a, const void* b)
{
    const Laik_Topology_IndexedElement* a_ptr = a;
    const Laik_Topology_IndexedElement* b_ptr = b;
    return a_ptr->val - b_ptr->val;
}

Laik_Topology_IndexedElement top_QAP_findIndexed(Laik_Topology_IndexedElement* arr, const size_t len, const int minmax)
{
    Laik_Topology_IndexedElement* best = arr;
    if (minmax == LAIK_TOP_MIN)
        for (size_t i = 1; i < len; i++) best = ((arr+i)->val < best->val) ? (arr+i) : best;
    else
        for (size_t i = 1; i < len; i++) best = ((arr+i)->val > best->val) ? (arr+i) : best;
    return *best;
}

// do the QAP construction method
// @returns heap pointer to reordering
int* top_QAP_construction(Laik_CommMatrix* cm, Laik_TopologyMatrix* top)
{
    size_t n          = cm->nodecount;
    int*   reordering = calloc(n, sizeof(int));
    if (!reordering) return NULL;

    int* identity = calloc(n, sizeof(int));
    if (!identity) return NULL;
    for (size_t i = 0; i < n; i++) identity[i] = (int)i;

    // save BOTH assigned and unassigned in one array:
    // [.., (assigned), .., | .., (unassigned), ..]
    // index keeps track of separator (== #assigned)
    // subarrays always sorted
    Laik_Topology_IndexedElement* loads = calloc(n, sizeof(Laik_Topology_IndexedElement));
    Laik_Topology_IndexedElement* dists = calloc(n, sizeof(Laik_Topology_IndexedElement));
    if (!loads || !dists) return NULL;

    int* procs = calloc(n, sizeof(int));
    int* cores = calloc(n, sizeof(int));
    if (!procs || !cores) return NULL;

    for (size_t i = 0; i < n; i++) {
        procs[i] = (int)i;
        cores[i] = (int)i;
    }

    size_t   assigned_procs = 0, assigned_cores = 0;
    uint64_t maxload_index = 0, mindist_index = 0;

    for (size_t i = 0; i < n; i++) {
        uint64_t maxload = 0, mindist = -1;
        uint64_t cur_load = top_QAP_commLoad(cm, i, n, identity);
        uint64_t cur_dist = top_QAP_coreDist(top, i, n, identity);
        if (cur_load > maxload) {
            maxload       = cur_load;
            maxload_index = i;
        }
        if (cur_dist < mindist) {
            mindist       = cur_dist;
            mindist_index = i;
        }
    }

    reordering[mindist_index] = maxload_index;
    top_QAP_list_assign(procs, n, &assigned_procs, maxload_index);
    top_QAP_list_assign(cores, n, &assigned_cores, mindist_index);

    // construct the reorder map rank by rank
    for (size_t i = 1; i < n; i++) {
        size_t num_unassigned = n - assigned_procs;
        for (size_t j = 0; j < num_unassigned; j++) {
            loads[j].index = procs[assigned_procs + j];
            loads[j].val   = top_QAP_commLoad(cm, procs[assigned_procs + j], assigned_procs, procs);
            dists[j].index = cores[assigned_cores + j];
            dists[j].val   = top_QAP_coreDist(top, cores[assigned_cores + j], assigned_cores, cores);
        }
            
        Laik_Topology_IndexedElement maxload_elm = top_QAP_findIndexed(loads, num_unassigned, LAIK_TOP_MAX);
        Laik_Topology_IndexedElement mindist_elm = top_QAP_findIndexed(dists, num_unassigned, LAIK_TOP_MIN);

        // assign the matching and add processes to assigned list
        reordering[mindist_elm.index] = maxload_elm.index;
        top_QAP_list_assign(procs, n, &assigned_procs, maxload_elm.index);
        top_QAP_list_assign(cores, n, &assigned_cores, mindist_elm.index);
    }

    free(identity);
    free(procs);
    free(cores);

    return reordering;
}

// do the QAP improvement runs
int* top_QAP_improvement(Laik_CommMatrix* cm, Laik_TopologyMatrix* top)
{
    return top_QAP_cyclicSearch(cm, top, top_QAP_construction(cm, top));
}

// do cyclic search around given order
// @param cm communication matrix @param top topology matrix
// @param initial initial ordering
// @returns initial, overwritten with best solution
int* top_QAP_cyclicSearch(Laik_CommMatrix* cm, Laik_TopologyMatrix* top, int* initial)
{
    const size_t n = cm->nodecount;

    int*     best_sol  = initial;
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
        else  // undo the exchange to reset to best solution
            top_QAP_pairxchg(current_sol, i, j);
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
    int tmp  = order[i];
    order[i] = order[j];
    order[j] = tmp;

    return order;
}

// get total weight for given reordering
uint64_t top_QAP_totalCost(Laik_CommMatrix* cm, Laik_TopologyMatrix* top, int* order)
{
    uint64_t cost = 0;
    for (size_t i = 0; i < cm->nodecount; i++)
        for (size_t j = 0; j < cm->nodecount; j++) cost += top_mat_elm(order[i], order[j], cm) * top_mat_elm(i, j, top);
    return cost;
}

int top_func_findEqual(const void* a, const void* b)
{
    const size_t* a_ptr = a;
    const size_t* b_ptr = b;
    return *a_ptr == *b_ptr;
}

// calculate total communication load between process and already assigned processes
// @param assigned should always be sorted
uint64_t top_QAP_commLoad(Laik_CommMatrix* mat, size_t process, size_t a_len, int* assigned)
{
    uint64_t load = 0;
    for (size_t i = 0; i < mat->nodecount; i++)
        if (i != process && bsearch(&i, assigned, a_len, sizeof(size_t), top_func_findEqual) != NULL) {
            load += top_mat_elm(process, i, mat);
            load += top_mat_elm(i, process, mat);
        }
    return load;
}

// calculate total distance between node and already assigned nodes
// ensure @param assigned is always sorted
uint64_t top_QAP_coreDist(Laik_TopologyMatrix* top, size_t node, size_t a_len, int* assigned)
{  // maybe use https://stackoverflow.com/a/25689059
    uint64_t dist = 0;
    for (size_t i = 0; i < top->nodecount; i++)
        if (i != node && bsearch(&i, assigned, a_len, sizeof(size_t), top_func_findEqual) != NULL)
            dist += top_mat_elm(node, i, top);
    return dist;
}

int top_func_sort(const void* a, const void* b)
{
    const int* a_ptr = a;
    const int* b_ptr = b;
    return *a_ptr - *b_ptr;
}

// move an element from unassigned to assigned part of the virtual list
// [.., (assigned), .., | .., (unassigned), ..]
// and update separator accordingly
// @param list list @param sep separator @param elm element to move
// @returns head of sorted lists
int* top_QAP_list_assign(int* list, size_t len, size_t* sep, size_t elm)
{
    int* where = bsearch(&elm, list + *sep, len - *sep, sizeof(int), top_func_sort);
    // element is already assigned
    if (where < list + *sep) return list;

    top_QAP_pairxchg(list, (*sep)++, (where - list));
    qsort(list, *sep, sizeof(int), top_func_sort);
    qsort(list + *sep, len - *sep, sizeof(int), top_func_sort);

    return list;
}
