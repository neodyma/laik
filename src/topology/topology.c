
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

// allocate a new CommMatrix
Laik_CommMatrix* laik_top_CommMatrix_init(Laik_Instance* li)
{
    Laik_CommMatrix* cm = calloc(1, sizeof(Laik_CommMatrix));

    if (!cm) laik_panic("Out of memory allocating CommMatrix object");

    size_t nodecount = li->locations;
    cm->matrix       = calloc((nodecount * nodecount), sizeof(*cm->matrix));

    if (!cm->matrix) laik_panic("Out of memory allocating CommMatrix matrix");

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
    // we could add a factor here to weigh certain paths
    top_mat_elm(from, to, cm) += amt;
    // top_mat_elm(to, from, cm) += amt;
    return cm;
}

// zero out given matrix
Laik_CommMatrix* laik_top_CommMatrix_reset(Laik_CommMatrix* cm)
{
    memset(cm->matrix, 0, cm->nodecount * sizeof(uint64_t));

    return cm;
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

// generate CommMatrix from Transition object
// TODO: for now work on external matrix
Laik_CommMatrix* laik_top_CommMatrix_add_Transition(Laik_CommMatrix* cm, Laik_Transition* tr)
{
    // Laik_CommMatrix* cm = laik_top_CommMatrix_init(tr->group->inst);

    // from -> to
    for (int i = 0; i < tr->sendCount; i++)
        laik_top_CommMatrix_update(cm, laik_myid(tr->group), tr->send[i].toTask, laik_range_size(&(tr->send[i].range)));

    // find all targets for reduction

    // for (int i = 0; i < tr->redCount; i++)
    // laik_top_CommMatrix_update(cm, laik_myid(tr->group), tr->red[i]);

    return cm;
}

// get reordered indices or NULL if no reordering is set
int* laik_top_reordering(Laik_Instance* li)
{  // this really should return u64*.. but locationid is int so we use that
    if (li->locationmap) return li->locationmap;

    char* reorderfile = getenv("LAIK_REORDER_FILE");
    char* reorderstr  = getenv("LAIK_REORDERING");
    if (!reorderfile && !reorderstr)  // no reordering set
        return NULL;

    char* reorderdyn = getenv("LAIK_REORDER_LIVE");
    if (reorderdyn) {
        li->locationmap = laik_top_do_reorder(laik_world(li)->comm_matrix, laik_top_Topology_from_sng(li));
        return li->locationmap;
    }

    // parse the env string and set reorderings
    if (reorderstr) {  // e.g. LAIK_REORDERING=2.3,0.4,5.1
        // TODO reverse lookup if not set?
        // e.g. 5 -> 1 but 1 is not mapped to another rank
        // in this case we should do a bidirectional map
        // idea, always map bidirectional except when already set?
        const size_t sz = li->locations * sizeof(int);
        li->locationmap = calloc(sz, 1);
        if (!li->locationmap) laik_panic("Reordering map could not be allocated!");

        char *kv, *save;
        laik_log(2, "Creating reorder map");
        // TODO: strtoul error handling
        while ((kv = strtok_r(reorderstr, ".", &save)) && !(reorderstr = NULL)) {
            if (!kv) break;
            size_t k = strtoul(kv, NULL, 0);

            kv = strtok_r(reorderstr, ",", &save);
            if (!kv) break;
            int v = strtol(kv, NULL, 0);

            if ((int)k >= li->locations) continue;

            li->locationmap[k] = v + LAIK_RO_OFFSET;
            // laik_log(2, "rank %lu -> rank %d\n", k, v);
        }

        if (reorderfile && (li->mylocationid == 0)) {  // both envs set -> dump reordering to file
            FILE* file = fopen(reorderfile, "wb");
            laik_log(2, "writing map to file %s\n", reorderfile);
            if (!file) laik_panic("Reordering file could not be opened!");

            if (fwrite(li->locationmap, 1, sz, file) != sz) laik_panic("Error writing to reordering file!");

            fclose(file);
        }
    }
    else if (reorderfile) {  // TODO: actually read and error detect
        FILE* file = fopen(reorderfile, "rb");
        if (!file) laik_panic("Reordering file could not be opened!");

        struct stat filestat;
        if (fstat(fileno(file), &filestat)) laik_panic("Reordering file stats could not be loaded!");

        if (!S_ISREG(filestat.st_mode) || filestat.st_size <= 0) laik_panic("Invalid reordering file!");

        // ...
    }
    return li->locationmap;
}

int* laik_top_reordering_get(Laik_Instance* li) { return li->locationmap; }

// allow LAIK to reorder processes
Laik_Group* laik_allow_reordering(Laik_Instance* li)
{
    if ((laik_epoch(li) | laik_phase(li)) == 0) {
        // initial reordering
        Laik_Group* g = laik_clone_group(li->world);
        laik_set_world(li, g);

        // either use LAIK_REORDERING or calculate on the fly
        int* reordermap = laik_top_reordering(li);
        if (reordermap && (reordermap[li->mylocationid] != LAIK_RO_UNMAPPED)) {
            int newid = reordermap[g->myid] - LAIK_RO_OFFSET;
            laik_log(2, "%s: mylocation %3d mapped to %3d", li->mylocation, g->myid, newid);
            g->myid = newid;
        }

        // propagate to backend
        li->backend->updateGroup(li->world);
    }
    else {
        // reordering after init, needs data movement
    }

    return li->world;
}

// return index of first differing char
size_t strcmp_index(const char* a, const char* b)
{
    size_t index = 0;
    for (; index < strlen(a) && index < strlen(b); index++)
        if (a[index] != b[index]) break;
    return index;
}

// generate the topology for running on SuperMUC-NG
Laik_Topology* laik_top_Topology_from_sng(Laik_Instance* li)
{
    // check if hostname conforms to iXXrXXcXXsXX
    char* loc = laik_mylocation(li);
    if (strlen(loc) < 10 || loc[0] != 'i' || loc[3] != 'r' || loc[6] != 'c' || loc[9] != 's') return NULL;

    Laik_Topology* top = laik_top_Topology_init(li, LAIK_TOP_IS_MAT);
    laik_sync_location(li);  // hostnames in li->location

    // uint64_t             weights[]     = {1, 4, 0, 0, 15};
    // int                  myloc         = li->mylocationid;
    uint64_t             hop_weights[] = {2, 10, 10, 10, 40};
    Laik_TopologyMatrix* mat           = top->data.mat;

    // local matrix
    // for (size_t i = 0; i < li->locations; i++) {
    //     if (i == myloc) continue;
    //     uint64_t weight = 0;
    //     // weight = total sum of weights across hops
    //     if (strcmp_index(li->mylocation, li->location[i]) < 3)
    //         weight = hop_weights[4];
    //     else if (strcmp_index(li->mylocation, li->location[i]) < 6)
    //         weight = hop_weights[3];
    //     else if (strcmp_index(li->mylocation, li->location[i]) < 9)
    //         weight = hop_weights[2];
    //     else if (strcmp_index(li->mylocation, li->location[i]) < 12)
    //         weight = hop_weights[1];
    //     else if (strcmp_index(li->mylocation, li->location[i]) == 12)
    //         weight = hop_weights[0];

    //     top_mat_elm(i, myloc, mat) = weight;
    //     top_mat_elm(myloc, i, mat) = weight;
    // }

    // global matrix
    for (size_t i = 0; i < (size_t)li->locations; i++) {
        for (size_t j = i; j < (size_t)li->locations; j++) {
            if (i == j) continue;
            uint64_t weight = 0;

            if (strcmp_index(li->location[i], li->location[j]) < 3)
                weight = hop_weights[4];
            else if (strcmp_index(li->location[i], li->location[j]) < 6)
                weight = hop_weights[3];
            else if (strcmp_index(li->location[i], li->location[j]) < 9)
                weight = hop_weights[2];
            else if (strcmp_index(li->location[i], li->location[j]) < 12)
                weight = hop_weights[1];
            else if (strcmp_index(li->location[i], li->location[j]) == 12)
                weight = hop_weights[0];

            top_mat_elm(i, j, mat) = weight;
            top_mat_elm(j, i, mat) = weight;
        }
    }

    return top;
}

Laik_Topology* laik_top_Topology_init(Laik_Instance* li, int which)
{
    // only matrices for now
    // if (which != LAIK_TOP_IS_MAT && which != LAIK_TOP_IS_GRAPH) return NULL;
    Laik_Topology* top = NULL;
    if (which == LAIK_TOP_IS_MAT) {
        top = calloc(sizeof(Laik_Topology), 1);
        if (!top) return NULL;

        // TODO: change this to TopologyMatrix_init
        top->data.mat = laik_top_Topology_TopopologyMatrix_init(li);
        if (!top->data.mat) return NULL;

        top->which = LAIK_TOP_IS_MAT;
    }

    return top;
}

void laik_top_Topology_free(Laik_Topology* top)
{
    if (top->which == LAIK_TOP_IS_MAT) {
        free(top->data.mat->matrix);
        free(top->data.mat);
    }
    else
        free(top->data.graph);
}

Laik_TopologyMatrix* laik_top_Topology_TopopologyMatrix_init(Laik_Instance* li)
{
    Laik_TopologyMatrix* topmat = calloc(sizeof(Laik_TopologyMatrix), 1);
    if (!topmat) return NULL;

    topmat->inst      = li;
    topmat->nodecount = li->locations;
    topmat->matrix    = calloc(sizeof(uint64_t), topmat->nodecount * topmat->nodecount);
    if (!topmat->matrix) {
        free(topmat);
        return NULL;
    }

    return topmat;
}

int* laik_top_do_reorder(Laik_CommMatrix* cm, Laik_Topology* top)
{
    if (!cm || !top) return NULL;
    return laik_top_do_reorder_QAP(cm, top);
}
