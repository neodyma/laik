/* 
 * This file is part of the LAIK parallel container library.
 * Copyright (c) 2017 Josef Weidendorfer
 */

#include "laik-internal.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// counter for space ID, just for debugging
static int space_id = 0;

// counter for partitioning ID, just for debugging
static int part_id = 0;

// helpers

static
int getSpaceStr(char* s, Laik_Space* spc)
{
    switch(spc->dims) {
    case 1:
        return sprintf(s, "[0-%lu]", spc->size[0]-1);
    case 2:
        return sprintf(s, "[0-%lu/0-%lu]",
                       spc->size[0]-1, spc->size[1]-1);
    case 3:
        return sprintf(s, "[0-%lu/0-%lu/0-%lu]",
                       spc->size[0]-1, spc->size[1]-1, spc->size[2]-1);
    }
    return 0;
}


static
int getIndexStr(char* s, int dims, Laik_Index* idx, bool minus1)
{
    uint64_t i1 = idx->i[0];
    uint64_t i2 = idx->i[1];
    uint64_t i3 = idx->i[2];
    if (minus1) {
        i1--;
        i2--;
        i3--;
    }

    switch(dims) {
    case 1:
        return sprintf(s, "%lu", i1);
    case 2:
        return sprintf(s, "%lu/%lu", i1, i2);
    case 3:
        return sprintf(s, "%lu/%lu/%lu", i1, i2, i3);
    }
    return 0;
}

// is the given slice empty?
bool laik_slice_isEmpty(int dims, Laik_Slice* slc)
{
    if (slc->from.i[0] >= slc->to.i[0])
        return true;

    if (dims>1) {
        if (slc->from.i[1] >= slc->to.i[1])
            return true;

        if (dims>2) {
            if (slc->from.i[2] >= slc->to.i[2])
                return true;
        }
    }
    return false;
}


// returns false if intersection of ranges is empty
static
bool intersectRange(uint64_t from1, uint64_t to1, uint64_t from2, uint64_t to2,
                    uint64_t* resFrom, uint64_t* resTo)
{
    if (from1 >= to2) return false;
    if (from2 >= to1) return false;
    *resFrom = (from1 > from2) ? from1 : from2;
    *resTo = (to1 > to2) ? to2 : to1;
    return true;
}

// get the intersection of 2 slices; return 0 if intersection is empty
Laik_Slice* laik_slice_intersect(int dims, Laik_Slice* s1, Laik_Slice* s2)
{
    static Laik_Slice s;

    if (!intersectRange(s1->from.i[0], s1->to.i[0],
                        s2->from.i[0], s2->to.i[0],
                        &(s.from.i[0]), &(s.to.i[0])) ) return 0;
    if (dims>1) {
        if (!intersectRange(s1->from.i[1], s1->to.i[1],
                            s2->from.i[1], s2->to.i[1],
                            &(s.from.i[1]), &(s.to.i[1])) ) return 0;
        if (dims>2) {
            if (!intersectRange(s1->from.i[2], s1->to.i[2],
                                s2->from.i[2], s2->to.i[2],
                                &(s.from.i[2]), &(s.to.i[2])) ) return 0;
        }
    }
    return &s;
}

static
Laik_Slice* sliceFromSpace(Laik_Space* s)
{
    static Laik_Slice slc;

    slc.from.i[0] = 0;
    slc.from.i[1] = 0;
    slc.from.i[2] = 0;
    slc.to.i[0] = s->size[0];
    slc.to.i[1] = s->size[1];
    slc.to.i[2] = s->size[2];

    return &slc;
}

static
int getSliceStr(char* s, int dims, Laik_Slice* slc)
{
    if (laik_slice_isEmpty(dims, slc))
        return sprintf(s, "(empty)");

    int off;
    off  = sprintf(s, "[");
    off += getIndexStr(s+off, dims, &(slc->from), false);
    off += sprintf(s+off, "-");
    off += getIndexStr(s+off, dims, &(slc->to), true);
    off += sprintf(s+off, "]");
    return off;
}


static
int getPartitioningTypeStr(char* s, Laik_PartitionType type)
{
    switch(type) {
    case LAIK_PT_All:    return sprintf(s, "all");
    case LAIK_PT_Stripe: return sprintf(s, "stripe");
    case LAIK_PT_Master: return sprintf(s, "master");
    }
    return 0;
}

static
int getAccessPermissionStr(char* s, Laik_AccessPermission ap)
{
    switch(ap) {
    case LAIK_AP_ReadOnly:  return sprintf(s, "readonly");
    case LAIK_AP_WriteOnly: return sprintf(s, "writeonly");
    case LAIK_AP_ReadWrite: return sprintf(s, "readwrite");
    case LAIK_AP_Plus:      return sprintf(s, "plus-red");
    case LAIK_AP_Times:     return sprintf(s, "times-red");
    case LAIK_AP_Min:       return sprintf(s, "min-red");
    case LAIK_AP_Max:       return sprintf(s, "max-red");
    }
    return 0;
}


static
int getTransitionStr(char* s, Laik_Transition* t)
{
    int off = 0;

    if (t->redCount>0) {
        off += sprintf(s+off, "  %d reds: ", t->redCount);
        for(int i=0; i<t->redCount; i++) {
            if (i>0) off += sprintf(s+off, ", ");
            off += getAccessPermissionStr(s+off, t->redOp[i]);
            off += getSliceStr(s+off, t->dims, &(t->red[i]));
            off += sprintf(s+off, " => %s (%d)",
                           (t->redRoot[i] == -1) ? "all":"master",
                           t->redRoot[i]);
        }
        off += sprintf(s+off, "\n");
    }

    if (t->sendCount>0) {
        off += sprintf(s+off, "  %d sends: ", t->sendCount);
        for(int i=0; i<t->sendCount; i++) {
            if (i>0) off += sprintf(s+off, ", ");
            off += getSliceStr(s+off, t->dims, &(t->send[i]));
            off += sprintf(s+off, " => %d", t->sendTo[i]);
        }
        off += sprintf(s+off, "\n");
    }

    if (t->recvCount>0) {
        off += sprintf(s+off, "  %d recvs: ", t->recvCount);
        for(int i=0; i<t->recvCount; i++) {
            if (i>0) off += sprintf(s+off, ", ");
            off += sprintf(s+off, "%d => ", t->recvFrom[i]);
            off += getSliceStr(s+off, t->dims, &(t->recv[i]));
        }
        off += sprintf(s+off, "\n");
    }
    return off;
}



// create a new index space object (initially invalid)
Laik_Space* laik_new_space(Laik_Instance* i)
{
    Laik_Space* space = (Laik_Space*) malloc(sizeof(Laik_Space));

    space->id = space_id++;
    space->name = strdup("space-0     ");
    sprintf(space->name, "space-%d", space->id);

    space->inst = i;
    space->dims = 0; // invalid
    space->first_partitioning = 0;

    // append this space to list of spaces used by LAIK instance
    space->next = i->firstspace;
    i->firstspace = space;

    return space;
}

// create a new index space object with an initial size
Laik_Space* laik_new_space_1d(Laik_Instance* i, uint64_t s1)
{
    Laik_Space* space = laik_new_space(i);
    space->dims = 1;
    space->size[0] = s1;

#ifdef LAIK_DEBUG
    char s[100];
    getSpaceStr(s, space);
    printf("LAIK %d/%d - new 1d space '%s': %s\n",
           space->inst->myid, space->inst->size, space->name, s);
#endif

    return space;
}

Laik_Space* laik_new_space_2d(Laik_Instance* i,
                              uint64_t s1, uint64_t s2)
{
    Laik_Space* space = laik_new_space(i);
    space->dims = 2;
    space->size[0] = s1;
    space->size[1] = s2;

#ifdef LAIK_DEBUG
    char s[100];
    getSpaceStr(s, space);
    printf("LAIK %d/%d - new 2d space '%s': %s\n",
           space->inst->myid, space->inst->size, space->name, s);
#endif

    return space;
}

Laik_Space* laik_new_space_3d(Laik_Instance* i,
                              uint64_t s1, uint64_t s2, uint64_t s3)
{
    Laik_Space* space = laik_new_space(i);
    space->dims = 3;
    space->size[0] = s1;
    space->size[1] = s2;
    space->size[2] = s3;

#ifdef LAIK_DEBUG
    char s[100];
    getSpaceStr(s, space);
    printf("LAIK %d/%d - new 3d space '%s': %s\n",
           space->inst->myid, space->inst->size, space->name, s);
#endif

    return space;
}

// free a space with all resources depending on it (e.g. paritionings)
void laik_free_space(Laik_Space* s)
{
    free(s->name);

    // TODO
}

// set a space a name, for debug output
void laik_set_space_name(Laik_Space* s, char* n)
{
    s->name = strdup(n);
}

// change the size of an index space, eventually triggering a repartitiong
void laik_change_space_1d(Laik_Space* s, uint64_t s1)
{
    assert(s->dims == 1);
    if (s->size[0] == s1) return;

    s->size[0] = s1;

    // TODO: notify partitionings about space change
}

void laik_change_space_2d(Laik_Space* s,
                          uint64_t s1, uint64_t s2)
{
    assert(0); // TODO
}

void laik_change_space_3d(Laik_Space* s,
                          uint64_t s1, uint64_t s2, uint64_t s3)
{
    assert(0); // TODO
}


// create a new partitioning on a space
Laik_Partitioning*
laik_new_partitioning(Laik_Space* s,
                      Laik_PartitionType pt,
                      Laik_AccessPermission ap)
{
    Laik_Partitioning* p;
    p = (Laik_Partitioning*) malloc(sizeof(Laik_Partitioning));

    p->id = part_id++;
    p->name = strdup("partng-0     ");
    sprintf(p->name, "partng-%d", p->id);

    p->space = s;
    p->next = s->first_partitioning;
    s->first_partitioning = p;

    p->permission = ap;
    p->type = pt;
    p->group = laik_world(s->inst);
    p->pdim = 0;

    p->base = 0;
    p->haloWidth = 0;

    p->bordersValid = false;
    p->borders = 0;

    return p;
}

Laik_Partitioning*
laik_new_base_partitioning(Laik_Space* space,
                           Laik_PartitionType pt,
                           Laik_AccessPermission ap)
{
    Laik_Partitioning* p;
    p = laik_new_partitioning(space, pt, ap);

#ifdef LAIK_DEBUG
    char s[100];
    getPartitioningTypeStr(s, p->type);
    getAccessPermissionStr(s+50, p->permission);
    printf("LAIK %d/%d - new partitioning '%s': type %s, access %s, group %d\n",
           space->inst->myid, space->inst->size, p->name,
           s, s+50, p->group->gid);
#endif

    return p;
}

// for multiple-dimensional spaces, set dimension to partition (default is 0)
void laik_set_partitioning_dimension(Laik_Partitioning* p, int d)
{
    assert((d >= 0) && (d < p->space->dims));
    p->pdim = d;
}


// create a new partitioning based on another one on the same space
Laik_Partitioning*
laik_new_coupled_partitioning(Laik_Partitioning* p,
                              Laik_PartitionType pt,
                              Laik_AccessPermission ap)
{
    Laik_Partitioning* partitioning;
    partitioning = laik_new_partitioning(p->space, pt, ap);

    assert(0); // TODO

    return partitioning;
}

// create a new partitioning based on another one on a different space
// this also needs to know which dimensions should be coupled
Laik_Partitioning*
laik_new_spacecoupled_partitioning(Laik_Partitioning* p,
                                   Laik_Space* s, int from, int to,
                                   Laik_PartitionType pt,
                                   Laik_AccessPermission ap)
{
    Laik_Partitioning* partitioning;
    partitioning = laik_new_partitioning(p->space, pt, ap);

    assert(0); // TODO

    return partitioning;
}

// free a partitioning with related resources
void laik_free_partitioning(Laik_Partitioning* p)
{
    free(p->name);

    // TODO
}


// give a partitioning a name, for debug output
void laik_set_partitioning_name(Laik_Partitioning* p, char* n)
{
    p->name = strdup(n);
}


// make sure partitioning borders are up to date
void laik_update_partitioning(Laik_Partitioning* p)
{
    Laik_Slice* baseBorders = 0;
    int pdim = p->pdim;
    int basepdim;

    if (p->base) {
        laik_update_partitioning(p->base);
        baseBorders = p->base->borders;
        basepdim = p->base->pdim;
        // sizes of coupled dimensions should be equal
        assert(p->space->size[pdim] == p->base->space->size[basepdim]);
    }

    if (p->bordersValid) return;

    int count = p->group->size;
    if (!p->borders)
        p->borders = (Laik_Slice*) malloc(count * sizeof(Laik_Slice));

    // partition according to dimension 0
    uint64_t size = p->space->size[pdim];
    uint64_t idx = 0;
    uint64_t inc = size / count;
    if (inc * count < size) inc++;

    for(int task = 0; task < count; task++) {
        Laik_Slice* b = &(p->borders[task]);
        switch(p->type) {
        case LAIK_PT_All:
            b->from.i[0] = 0;
            b->from.i[1] = 0;
            b->from.i[2] = 0;
            b->to.i[0] = p->space->size[0];
            b->to.i[1] = p->space->size[1];
            b->to.i[2] = p->space->size[2];
            break;

        case LAIK_PT_Stripe:
            b->from.i[0] = 0;
            b->from.i[1] = 0;
            b->from.i[2] = 0;
            b->to.i[0] = p->space->size[0];
            b->to.i[1] = p->space->size[1];
            b->to.i[2] = p->space->size[2];

            b->from.i[pdim] = idx;
            idx += inc;
            if (idx > size) idx = size;
            b->to.i[pdim] = idx;
            break;

        case LAIK_PT_Master:
            b->from.i[0] = 0;
            b->from.i[1] = 0;
            b->from.i[2] = 0;
            b->to.i[0] = (task == 0) ? p->space->size[0] : 0;
            b->to.i[1] = (task == 0) ? p->space->size[1] : 0;
            b->to.i[2] = (task == 0) ? p->space->size[2] : 0;
            break;

        case LAIK_PT_Copy:
            assert(baseBorders);
            b->from.i[0] = 0;
            b->from.i[1] = 0;
            b->from.i[2] = 0;
            b->to.i[0] = p->space->size[0];
            b->to.i[1] = p->space->size[1];
            b->to.i[2] = p->space->size[2];

            b->from.i[pdim] = baseBorders[task].from.i[basepdim];
            b->to.i[pdim] = baseBorders[task].to.i[basepdim];
            break;

        default:
            assert(0); // TODO
            break;
        }
    }
    p->bordersValid = true;

#ifdef LAIK_DEBUG
    char str[1000];
    int off;
    off = sprintf(str, "partitioning '%s' (group %d) updated: ",
                  p->name, p->group->gid);
    for(int task = 0; task < count; task++) {
        if (task>0)
            off += sprintf(str+off, ", ");
        off += sprintf(str+off, "%d:", task);
        off += getSliceStr(str+off, p->space->dims, &(p->borders[task]));
    }
    printf("LAIK %d/%d - %s\n",
           p->group->inst->myid, p->group->inst->size, str);
#endif
}



// append a partitioning to a partioning group whose consistency should
// be enforced at the same point in time
void laik_append_partitioning(Laik_PartGroup* g, Laik_Partitioning* p)
{
    assert(0); // TODO
}

// Calculate communication required for transitioning between partitionings
Laik_Transition* laik_calc_transitionP(Laik_Partitioning* from,
                                       Laik_Partitioning* to)
{
    Laik_Transition* t = (Laik_Transition*) malloc(sizeof(Laik_Transition));
    t->sendCount = 0;
    t->recvCount = 0;
    t->redCount = 0;

    assert(from->space == to->space);
    int dims = from->space->dims;
    t->dims = dims;

    int myid = from->group->inst->myid;
    int count = from->group->size;

    // something to send?
    switch(from->permission) {
    case LAIK_AP_Max:
    case LAIK_AP_Min:
    case LAIK_AP_Plus:
    case LAIK_AP_Times:
        // reductions always should involve everyone
        assert(from->type == LAIK_PT_All);
        if ((to->permission == LAIK_AP_ReadOnly) ||
            (to->permission == LAIK_AP_ReadWrite)) {
                assert(t->redCount < COMMSLICES_MAX);
                assert((to->type == LAIK_PT_Master) ||
                       (to->type == LAIK_PT_All));
                t->red[t->redCount] = *sliceFromSpace(from->space);
                t->redOp[t->redCount] = from->permission;
                t->redRoot[t->redCount] = (to->type == LAIK_PT_All) ? -1 : 0;
                t->redCount++;
        }
        break;
    case LAIK_AP_ReadWrite:
    case LAIK_AP_WriteOnly:
        switch(from->type) {
        case LAIK_PT_Master:
        case LAIK_PT_All:
        case LAIK_PT_Stripe:
            if (!laik_slice_isEmpty(dims, &(from->borders[myid]))) {
                for(int task = 0; task < count; task++) {
                    Laik_Slice* s;
                    s = laik_slice_intersect(dims,
                                             &(from->borders[myid]),
                                             &(to->borders[task]));
                    if (s == 0) continue;

                    assert(t->sendCount < COMMSLICES_MAX);
                    t->send[t->sendCount] = *s;
                    t->sendTo[t->sendCount] = task;
                    t->sendCount++;
                }
            }
            break;
        default:
            break;
        }
        break;
    default:
        break;
    }

    switch(from->permission) {
    // aggregation do not receive...
    case LAIK_AP_Max:
    case LAIK_AP_Min:
    case LAIK_AP_Plus:
    case LAIK_AP_Times:
        break;
    default:
        // something to receive?
        switch(to->permission) {
        case LAIK_AP_ReadWrite:
        case LAIK_AP_ReadOnly:
            switch(to->type) {
            case LAIK_PT_Master:
            case LAIK_PT_All:
            case LAIK_PT_Stripe:
                if (!laik_slice_isEmpty(dims, &(to->borders[myid]))) {
                    for(int task = 0; task < count; task++) {
                        Laik_Slice* s;
                        s = laik_slice_intersect(dims,
                                                 &(to->borders[myid]),
                                                 &(from->borders[task]));
                        if (s == 0) continue;

                        assert(t->recvCount < COMMSLICES_MAX);
                        t->recv[t->recvCount] = *s;
                        t->recvFrom[t->recvCount] = task;
                        t->recvCount++;
                    }
                }
                break;
            default:
                break;
            }
            break;
        default:
            break;
        }
    }

#ifdef LAIK_DEBUG
        char s[1000];
        getTransitionStr(s, t);
        printf("LAIK %d/%d - transition %s => %s:\n%s",
               from->space->inst->myid, from->space->inst->size,
               from->name, to->name, s);
#endif


    return t;
}

// Calculate communication for transitioning between partitioning groups
Laik_Transition* laik_calc_transitionG(Laik_PartGroup* from,
                                       Laik_PartGroup* to)
{
    Laik_Transition* t;

    assert(0); // TODO
}

// enforce consistency for the partitioning group, depending on previous
void laik_enforce_consistency(Laik_Instance* i, Laik_PartGroup* g)
{
    assert(0); // TODO
}

// set a weight for each participating task in a partitioning, to be
//  used when a repartitioning is requested
void laik_set_partition_weights(Laik_Partitioning* p, int* w)
{
    assert(0); // TODO
}


// change an existing base partitioning
void laik_repartition(Laik_Partitioning* p, Laik_PartitionType pt)
{
    assert(0); // TODO
}


// couple different LAIK instances via spaces:
// one partition of calling task in outer space is mapped to inner space
void laik_couple_nested(Laik_Space* outer, Laik_Space* inner)
{
    assert(0); // TODO
}



