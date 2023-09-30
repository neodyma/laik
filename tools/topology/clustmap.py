import itertools
import igraph
from toptypes import HostGraph, getNodeChildren, UINT64_MAX

class ClustMap:
    def __init__(self, comm_mat: list, top_graph: HostGraph, hostnames: list[str]) -> None:
        self.comm_mat = comm_mat
        self.mod_mat = comm_mat
        self.top_graph = top_graph
        self.hostnames = hostnames
        self.M = []

    def __str__(self) -> str:
        return (
            "ClustMap\n"
            + "mat     {}\n".format(self.comm_mat)
            + "graph   {}\n".format(self.top_graph)
            + "hosts   {}".format(self.hostnames)
        )

    def doClustMap(self, S, A):
        if len(S) == len(A):
            self.MapClusters(S, A)
        elif len(S) < len(A):
            self.reMap(S, A)
        elif len(S) > len(A):
            for subcluster in S:
                self.doClustMap(subcluster, A)

    def MapClusters(self, S, A):
        j = 0
        for leaf in self.getLeaves(A):
            self.M[0][j] = S[j]
            self.M[1][j] = A[j]
        return
    
    def reMap(self, S, A):
        i = 0
        j = 0
        for process, core in zip(A, S):
            self.M[0][i] = process[i]
            self.M[1][i] = core[j]
        return

    def getLeaves(self, A) -> list:
        return []
