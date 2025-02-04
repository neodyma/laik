import itertools
import igraph
import numpy as np
from toptypes import HostGraph, getNodeChildren, UINT64_MAX
from typing import Iterable


class TauQAP:
    def __init__(self, comm_mat: list, top_graph: HostGraph, hostnames: list[str]) -> None:
        self.comm_mat = comm_mat
        self.top_graph = top_graph
        self.hostnames = hostnames

    def __str__(self) -> str:
        return (
            "TauQAP\n"
            + "mat     {}\n".format(self.comm_mat)
            + "graph   {}\n".format(self.top_graph)
            + "hosts   {}".format(self.hostnames)
        )

    def solve(self):
        # return [self.hostnames.index(x) for x in self.doQAP()]
        return [x for x in self.doQAP()]

    def doQAP(self):
        # return self.doConstructionMethod()
        return self.doImprovementMethod()

    # construct a reordering based on communication load and distances
    def doConstructionMethod(self):
        reordering = [0 for _ in range(len(self.hostnames))]
        initial_loads = [self.calcCommLoad(x, list(range(len(self.hostnames)))) for x in range(len(self.hostnames))]
        initial_dists = [self.calcCoreDist(x, list(range(len(self.hostnames)))) for x in range(len(self.hostnames))]

        # print("initial process loads:  ", initial_loads)
        # print("initial core distances: ", initial_dists)

        max_load = initial_loads.index(max(initial_loads))
        min_dist = initial_dists.index(min(initial_dists))

        # the core with the lowest total distance receives the rank with highest comm load
        # print("assigned {} to {}".format(max_load, min_dist))

        reordering[min_dist] = max_load
        unassigned_procs = list(range(len(self.hostnames)))
        unassigned_procs.remove(max_load)
        assigned_procs = [max_load]
        unassigned_cores = list(range(len(self.hostnames)))
        unassigned_cores.remove(min_dist)
        assigned_cores = [min_dist]

        for n in range(1, len(self.hostnames)):
            # calculate all loads and dists for every still unassigned process
            # print("unassigned procs for n {}: ".format(n), unassigned_procs)
            # print("unassigned cores for n {}: ".format(n), unassigned_cores)

            loads = [self.calcCommLoad(proc, assigned_procs) for proc in unassigned_procs]
            dists = [self.calcCoreDist(core, assigned_cores) for core in unassigned_cores]

            # print("loads for n {}: loads, assigned".format(n), loads, assigned_procs)
            # print("dists for n {}: dists, assigned".format(n), dists, assigned_cores)
            # print("index of max load: ", loads[::-1].index(max(loads)), "-> proc", unassigned_procs[loads.index(max(loads))])
            # print("index of min dist: ", dists[::-1].index(min(dists)), "-> core", unassigned_cores[dists.index(min(dists))])


            # get element from unassigned list which matches minmax load/dist
            max_load = unassigned_procs[loads.index(max(loads))]
            min_dist = unassigned_cores[dists.index(min(dists))]
            # print(max_load, min_dist)            
            reordering[min_dist] = max_load
            # print("assigned {} to {}".format(max_load, min_dist))

            unassigned_procs.remove(max_load)
            assigned_procs.append(max_load)
            unassigned_cores.remove(min_dist)
            assigned_cores.append(min_dist)

        # print("QAP Construction: ", reordering, self.totalCost(reordering))
        return reordering

    # iteratively improve initial reordering (identity)
    def doImprovementMethod(self):
        # return self.cyclicSearch(list(range(len(self.hostnames))))[0]
        res = self.cyclicSearch(self.doConstructionMethod())
        # print("total QAP cost: ", res[1])
        return res[0]

    # calculate total communication load between process and already assigned processes
    # initial process: treat every process as assigned
    def calcCommLoad(self, process: int, assigned: Iterable[int]) -> int:
        load = 0
        for i in range(len(self.hostnames)):
            if i == process or i not in assigned:
                continue
            load += self.comm_mat[process][i]
            load += self.comm_mat[i][process]
        return load

    # calculate total distance between node and already assigned nodes
    def calcCoreDist(self, node: int, assigned: Iterable[int]):
        dist = 0
        for i in range(len(self.hostnames)):
            if i == node or i not in assigned:
                continue
            dist += self.top_graph.topMatrix[node][i]
            # dist += self.top_graph.topMatrix[i][node]
        return dist

    def cyclicSearch(self, initial: Iterable):
        # if self.totalCost(initial) > self.totalCost(range(len(self.hostnames))):
            # best_sol, best_cost = list(range(len(self.hostnames))), self.totalCost(range(len(self.hostnames)))
        # else:
        best_sol, best_cost = initial, self.totalCost(initial)
        current_sol, current_cost = best_sol, best_cost

        i, j, n = 0, 1, len(self.hostnames)
        for _ in range(n ** 2): # for i, for j
            current_sol = self.pairExchange(best_sol, i, j)
            current_cost = self.totalCost(current_sol)
            if current_cost < best_cost:
                best_sol = current_sol
                best_cost = current_cost
            if j < n - 1:
                j += 1
            elif j == n and i < n - 2:
                i += 1
                j = i + 2
            else:
                i = 1
                j = 2

        return best_sol, best_cost
    
    def pairExchange(self, order, i, j) -> list:
        reorder = order.copy()
        reorder[i], reorder[j] = reorder[j], reorder[i]
        return reorder

    def totalCost(self, order) -> int:
        cost = 0
        for i in range(len(self.hostnames)):
            for j in range(len(self.hostnames)):
                cost += self.comm_mat[order[i]][order[j]] * self.top_graph.topMatrix[i][j]
        return cost

if __name__ == "__main__":
    qap = TauQAP([[1,0,0,2],[2,0,0,0],[2,0,0,0],[2,0,0,0]], HostGraph(igraph.Graph(), [[1,10,10,1],[10,1,1,1],[10,1,1,1],[1,1,1,1]], [], []), ["0","1","2","3"])
    print(qap.cyclicSearch([0,1,2,3]))
