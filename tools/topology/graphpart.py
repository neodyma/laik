import itertools
import igraph
from toptypes import HostGraph, getNodeChildren, UINT64_MAX


class NECGraphPart:
    def __init__(self, comm_mat: list, top_graph: HostGraph, hostnames: list[str]) -> None:
        self.comm_mat = comm_mat
        self.mod_mat = comm_mat
        self.top_graph = top_graph
        self.hostnames = hostnames

    def __str__(self) -> str:
        return (
            "NECGraphPart\n"
            + "mat     {}\n".format(self.comm_mat)
            + "graph   {}\n".format(self.top_graph)
            + "hosts   {}".format(self.hostnames)
        )
