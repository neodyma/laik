/*
 * This file is part of the LAIK library.
 * Copyright (c) 2017, 2018 Josef Weidendorfer <Josef.Weidendorfer@gmx.de>
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


#ifndef LAIK_DEBUG_H
#define LAIK_DEBUG_H

#include <stdint.h>  // for uint64_t
#include "data.h"    // for Laik_SwitchStat
#include "space.h"   // for Laik_DataFlow, Laik_Index, Laik_Partitioning
#include "action.h"  // for Laik_Action

char* laik_at_str(Laik_ActionType t);

void laik_log_IntList(int len, int* list);
void laik_log_PrettyInt(uint64_t v);
void laik_log_Space(Laik_Space* spc);
void laik_log_Index(int dims, const Laik_Index* idx);
void laik_log_Range(Laik_Range* range);
void laik_log_Reduction(Laik_ReductionOperation op);
void laik_log_DataFlow(Laik_DataFlow flow);
void laik_log_Transition(Laik_Transition* t, bool showActions);
void laik_log_RangeList(Laik_RangeList* list);
void laik_log_RangeFilter(Laik_RangeFilter*);
void laik_log_Partitioning(Laik_Partitioning* p);
void laik_log_SwitchStat(Laik_SwitchStat* ss);
void laik_log_Action(Laik_Action* a, Laik_ActionSeq* as);
void laik_log_ActionSeq(Laik_ActionSeq* as, bool showActions);
void laik_log_Checksum(char* buf, int count, Laik_Type* t);
void laik_log_CommMatrix(Laik_CommMatrix* cm);
void laik_print_CommMatrix_to_file(Laik_CommMatrix* cm);

// write action sequence at level 1 if <changed> is true, prepend with title
void laik_log_ActionSeqIfChanged(bool changed, Laik_ActionSeq* as, char* title);


#endif // LAIK_DEBUG_H
