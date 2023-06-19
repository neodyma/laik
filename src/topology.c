
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
    cm->matrix       = calloc((nodecount * nodecount), sizeof(*cm->matrix));

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
    if (!b || !b->matsync)
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

// get reordered indices or NULL if no reordering is set
int* laik_top_reordering(Laik_Instance* li)
{ // this really should return u64*.. but locationid is int so we use that
    if (li->locationmap)
        return li->locationmap;

    char* reorderfile = getenv("LAIK_REORDER_FILE");
    char* reorderstr  = getenv("LAIK_REORDERING");
    if (!reorderfile && !reorderstr) // no reordering set
        return NULL;

    // parse the env string and set reorderings
    if (reorderstr) { // e.g. LAIK_REORDERING=2.3,0.4,5.1
        // TODO reverse lookup if not set?
        // e.g. 5 -> 1 but 1 is not mapped to another rank
        // in this case we should do a bidirectional map
        // idea, always map bidirectional except when already set?
        const size_t sz = li->locations * sizeof(int);
        li->locationmap = calloc(sz, 1);
        if (!li->locationmap)
            laik_panic("Reordering map could not be allocated!");

        char *kv, *save;
        laik_log(2, "Creating reorder map");
        // TODO: strtoul error handling
        while ((kv = strtok_r(reorderstr, ".", &save)) && !(reorderstr = NULL)) {
            if (!kv)
                break;
            size_t k = strtoul(kv, NULL, 0);

            kv = strtok_r(reorderstr, ",", &save);
            if (!kv)
                break;
            int v = strtol(kv, NULL, 0);

            if ((int)k >= li->locations)
                continue;

            li->locationmap[k] = v + LAIK_RO_OFFSET;
            laik_log(2, "rank %lu -> rank %d\n", k, v);
        }

        if (reorderfile && (li->mylocationid == 0)) { // both envs set -> dump reordering to file
            FILE* file = fopen(reorderfile, "wb");
            laik_log(2, "writing map to file %s\n", reorderfile);
            if (!file)
                laik_panic("Reordering file could not be opened!");

            if (fwrite(li->locationmap, 1, sz, file) != sz)
                laik_panic("Error writing to reordering file!");

            fclose(file);
        }
    }
    else if (reorderfile) { // TODO: actually read and error detect
        FILE* file = fopen(reorderfile, "rb");
        if (!file)
            laik_panic("Reordering file could not be opened!");

        struct stat filestat;
        if (fstat(fileno(file), &filestat))
            laik_panic("Reordering file stats could not be loaded!");

        if (!S_ISREG(filestat.st_mode) || filestat.st_size <= 0)
            laik_panic("Invalid reordering file!");
    }
    return li->locationmap;
}

int* laik_top_reordering_get(Laik_Instance* li)
{
    return li->locationmap;
}
