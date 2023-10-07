import itertools
import igraph
from toptypes import HostGraph, getNodeChildren, UINT64_MAX


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
        return [self.hostnames.index(x) for x in self.doQAP()]

    def doQAP(self):
        return

    # calculate total communication load between process and already assigned processes
    # initial process: treat every process as assigned
    def calcCommLoad(self, process: int, assigned: list) -> int:
        load = 0
        for i in range(len(self.hostnames)):
            if i == process or i not in assigned:
                continue
            load += self.comm_mat[process][i]
            load += self.comm_mat[i][process]
        return load

    # calculate total distance  between node and already assigned nodes
    def calcCoreDistance(self, node: int, assigned: list):
        dist = 0
        for i in range(len(self.hostnames)):
            if i == node or i not in assigned:
                continue
            dist += self.top_graph.topMatrix[node][i]
            dist += self.top_graph.topMatrix[i][node]
        return dist
