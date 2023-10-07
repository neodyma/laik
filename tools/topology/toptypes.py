import igraph
from dataclasses import dataclass

UINT64_MAX = 0xFFFFFFFFFFFFFFFF


@dataclass
class CommStats:
    commGraph: igraph.Graph
    commMatrix: list
    hostnames: list


@dataclass
class HostGraph:
    graph: igraph.Graph
    topMatrix: list
    layers: list
    weights: list


def getNodeChildren(graph: HostGraph, node: str, nodelayer: int) -> list:
    if nodelayer + 1 >= len(graph.layers):
        return []
    connected_nodes = set(graph.graph.vs[graph.graph.neighbors(node)]["name"]).intersection(graph.layers[nodelayer + 1])
    return list(connected_nodes)
