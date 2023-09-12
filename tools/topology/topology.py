#!/usr/bin/python3
import argparse
from functools import reduce
import igraph
import itertools
import math
import re
from dataclasses import dataclass


@dataclass
class CommStats:
    commGraph: igraph.Graph
    commMatrix: list
    hostnames: list


# LAIK_LOG_FILE -> commGraph, commMatrix, hostnames
def parseCommStats(logfiles) -> CommStats:
    commMatrix = list

    marker = "Communication Matrix"
    backend = "MPI backend initialized"
    R_BE = re.compile(r"MPI backend initialized \(at (.+), rank (\d+)/(\d+)\).*$")
    # == LAIK-(0000)-L(03) (0006).(01)  (0):(00).(004) | (Communication Matrix:) (\n)
    R_EX = re.compile(r".+LAIK-(\d+)-L(\d+) (\d+).(\d+)  (\d+):(\d+).(\d+) \| (.+)")

    ranks = 0
    matrices = []
    hostnames = list

    for file in logfiles:
        with open(file, "r") as logfile:
            for line in logfile:
                if backend in line:  # set number of ranks
                    output = R_BE.search(line)
                    if ranks == 0:
                        ranks = int(output.group(3))
                        commMatrix = [[0 for _ in range(ranks)] for _ in range(ranks)]
                        hostnames = ["" for _ in range(ranks)]
                    hostnames[int(output.group(2))] = output.group(1)
                    continue
                rex = R_EX.search(line)
                if rex is not None:
                    if marker in line:
                        matrices.append((rex.group(2), rex.group(3)))  # (rank, sequence)
                    elif (rex.group(2), rex.group(3)) in matrices and int(rex.group(4)) > 2:
                        # matched a matrix row
                        content = R_EX.split(line)[-2].split()
                        if content is None:
                            raise ValueError
                        if len(content) != ranks + 2:  # misformed matrix output
                            raise IndexError
                        values = [int(x) for x in content[2:]]
                        commMatrix[int(content[0])] = [sum(x) for x in zip(commMatrix[int(content[0])], values)]

    # undirected: double transfer values!
    # commGraph = igraph.Graph.Weighted_Adjacency(commMatrix, mode="undirected")
    print(f"Parsed matrices from {len(logfiles)} files.")
    commGraph = igraph.Graph.Weighted_Adjacency(commMatrix)
    return CommStats(commGraph, commMatrix, hostnames)


@dataclass
class HostGraph:
    graph: igraph.Graph
    layers: list


# generate a host graph based on the supermuc-ng node naming scheme
# i01r01c01s01
def generateHostGraph(hostnames: list[str]) -> HostGraph:
    topGraph = igraph.Graph()
    R_LRZ = re.compile(r"(i\d+)(r\d+)(c\d+)(s\d+):(\d+)")  # too complicated

    isls = []
    racks = []
    cabs = []
    srvs = []

    topGraph.add_vertices(hostnames)
    # we want to build a "tree" of the used topology
    # need to insert nodes from the top down
    for host in hostnames:
        if host[0:3] not in isls:
            isls.append(host[0:3])
        if host[0:6] not in racks:
            racks.append(host[0:6])
        if host[0:9] not in cabs:
            cabs.append(host[0:9])
        if host[0:12] not in srvs:
            srvs.append(host[0:12])

    topGraph.add_vertices(isls)
    topGraph.add_vertices(racks)
    topGraph.add_vertices(cabs)
    topGraph.add_vertices(srvs)

    isls_edges = list(itertools.combinations(isls, r=2))
    racks_edges = list(filter(lambda xs: xs[0] == xs[1][:3], itertools.combinations(isls + racks, r=2)))
    cabs_edges = list(filter(lambda xs: xs[0] == xs[1][:6], itertools.combinations(isls + racks + cabs, r=2)))
    srvs_edges = list(filter(lambda xs: xs[0] == xs[1][:9], itertools.combinations(isls + racks + cabs + srvs, r=2)))
    node_edges = list(
        filter(lambda xs: xs[0] == xs[1][:12], itertools.combinations(isls + racks + cabs + srvs + hostnames, r=2))
    )

    # weights: (https://doku.lrz.de/hoechstleistungsrechner-10333235.html)
    # intra node: set to 1
    # inter node, on island: full omnipath
    # intra island: 3.75:1 compared to intra island
    topGraph.add_edges(node_edges, dict(weight=[1 for _ in range(len(node_edges))]))
    topGraph.add_edges(srvs_edges, dict(weight=[4 for _ in range(len(srvs_edges))]))
    topGraph.add_edges(cabs_edges, dict(weight=[0 for _ in range(len(cabs_edges))]))
    topGraph.add_edges(racks_edges, dict(weight=[0 for _ in range(len(racks_edges))]))
    topGraph.add_edges(isls_edges, dict(weight=[15 for _ in range(len(isls_edges))]))

    layers = []
    # filter root nodes with only one leaf to minimize tree
    if len(isls) < 2 and len(racks) < 2:
        topGraph.delete_vertices(isls)
    else:
        layers.append(isls)
    if len(racks) < 2 and len(cabs) < 2:
        topGraph.delete_vertices(racks)
    else:
        layers.append(racks)
    if len(cabs) < 2 and len(srvs) < 2:
        topGraph.delete_vertices(cabs)
    else:
        layers.append(cabs)

    layers.append(srvs)
    layers.append(hostnames)

    return HostGraph(topGraph, layers)


def getNodeChildren(graph: HostGraph, node: str, nodelayer: int) -> list:
    if nodelayer + 1 >= len(graph.layers):
        return []
    connected_nodes = set(graph.graph.vs[graph.graph.neighbors(node)]["name"]).intersection(graph.layers[nodelayer + 1])
    return list(connected_nodes)


# solve the embedding problem for given graphs and return acceptable reordering
# todo we probably need the hostnames sooner or later
#   for now on test we always get successive hosts
def reorder(comm_graph, top_graph) -> list:
    return [1, 2, 3]


def treeMatch(comm_mat, top_graph, hostnames) -> list:
    solver = TreeMatch(comm_mat, top_graph, hostnames)
    return solver.solve()


class TreeMatch:
    def __init__(self, comm_mat: list, top_graph: HostGraph, hostnames: list[str]):
        self.comm_mat = comm_mat
        self.mod_mat = comm_mat
        self.top_graph = top_graph
        self.hostnames = hostnames

    def arity(self, depth):
        return self.top_graph.graph.degree(self.top_graph.layers[depth][0]) - 1

    def solve(self):
        # return self.extendCommMatrix(1)
        # return self.aggregateCommMatrix((list(range(8)), list(range(8,16))))
        # return self.groupProcesses(2)
        return [self.hostnames.index(x) for x in self.doTreeMatch()]

    def doTreeMatch(self):
        groups = [[] for _ in range(len(self.top_graph.layers))]

        print("Starting TreeMatch for {} groups.".format(len(groups)))

        for d in range(len(self.top_graph.layers) - 1, 1, -1):
            p = len(self.mod_mat)
            if p % self.arity(d - 1) != 0:
                self.mod_mat = self.extendCommMatrix(d)
            groups[d] = self.groupProcesses(self.mod_mat, d)
            self.mod_mat = self.aggregateCommMatrix(groups[d])
        return [x for xs in groups[-1] for x in xs]

    def groupProcesses(self, matrix: list, cur_depth: int):
        # TODO: this only works for balanced trees for now
        UINT64_MAX = 0xFFFFFFFFFFFFFFFF
        # print("processes         ", len(self.hostnames))
        # print("groupsize         ", self.arity(cur_depth - 1))
        # TODO: this should only work with connected nodes instead of the entire layer
        l = list(itertools.combinations(self.top_graph.layers[cur_depth], self.arity(cur_depth - 1)))

        print("Grouping process combinations with arity {} at depth {}.".format(self.arity(cur_depth - 1), cur_depth))

        # print("combinations      ", l)
        # print("num_combinations  ", len(l))
        d = dict(list(map(lambda x: (hash(x) & UINT64_MAX, x), l)))
        G = igraph.Graph()
        G.add_vertices(list(map(lambda x: str(hash(x) & UINT64_MAX), l)))
        for group1 in l:
            for group2 in l:
                if group1 == group2:
                    continue
                for node in group1:
                    if node in group2:
                        G.add_edge(str(hash(group1) & UINT64_MAX), str(hash(group2) & UINT64_MAX))
                        # print("incompatible groups {} ~ {} added to set.".format(str(group1), str(group2)))

        # print("graph_groups      ", len(G.vs))
        # print("groups to find    ", len(self.comm_mat) // self.arity(cur_depth - 1))
        # TODO: improve independent set algorithm
        independent_set = G.independent_vertex_sets(
            len(self.comm_mat) // self.arity(cur_depth - 1), len(self.comm_mat) // self.arity(cur_depth - 1)
        )

        print("found {} possible combinations out of {} groups.".format(len(independent_set), len(G.vs)))

        print(
            "node children of {}: {}".format(
                self.top_graph.layers[cur_depth][0],
                getNodeChildren(self.top_graph, self.top_graph.layers[cur_depth][0], cur_depth),
            )
        )

        # print("independentset[0] ", independent_set[0], l[independent_set[0][0]], l[independent_set[0][1]])
        # TODO: filter the sets greedy
        min_group = []
        min_weight = 0xFFFFFFFF
        for matching in independent_set:
            cur_weight = self.groupWeight(matrix, [l[g] for g in matching])
            if cur_weight < min_weight:
                min_group = [l[g] for g in matching]
                min_weight = cur_weight

        return min_group

    # return the number of cross-node groups
    # TODO factor in the actual transfer cost in this
    # weight should just be the total communication cost of this group
    # TODO this only works on leaf level like this
    def groupWeight(self, matrix, group):
        weight = 0
        if len(group) == 0:
            return weight
        # for match in group:
            # TODO: calculate ranking based on communication matrix
            # use the outbound group communication to weigh
            # matchset = set(list(map(lambda x: x[:12], match)))
            # weight += len(matchset) - 1

        # print("calculating outbound communication for group {}".format(group))

        for match in group:
            weight -= self.elementSum([match], 0, 0)

        # exit()


        print("calculated weight {} for group {}".format(weight, group))

        return weight

    # aggregate submatrices according to groups
    def aggregateCommMatrix(self, groups):
        n = len(groups)
        r = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    r[i][j] = self.elementSum(groups, i, j)
        return r

    # element wise sum across group submatrices
    def elementSum(self, groups, i, j):
        acc = 0
        for a in groups[i]:
            for b in groups[j]:
                acc += self.comm_mat[self.hostnames.index(a)][self.hostnames.index(b)]

        # print("element wise sum for group {}: {}".format(groups, acc))

        return acc

    # extend matrix by arity of lower depth lines and columns
    # if cur_depth points to single root, ariry = degree
    # else vertex degree-1 due to minimal graph
    # treematch as described wants balanced nodes, does this change anything? ..
    def extendCommMatrix(self, cur_depth) -> list:
        p = len(self.mod_mat)
        k = self.arity(cur_depth)
        for l in self.mod_mat:
            l.extend([0] * k)
        for _ in range(k):
            self.mod_mat.append([0 for _ in range(p)])
        return self.mod_mat


def optimize(optimizer, *args) -> list:
    return globals()[optimizer](*args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="LAIK topology optimizer",
        description="Parse communication metadata, create structured representations, optimize communication paths and generate reordering parameters.",
    )
    parser.add_argument("-i", "--ilog", nargs="+", help="Input LAIK_LOG_FILEs")
    parser.add_argument("-c", "--cg", help="Communication graph")
    parser.add_argument("-t", "--tg", help="Topology graph")
    parser.add_argument("-r", help=f"Reorder using ansatz [treeMatch, ..]")
    args = parser.parse_args()

    # we have an input logfile, let's convert it to a usable Graph
    if args.ilog is not None:
        cg = parseCommStats(args.ilog)
        if args.cg is not None:
            # dump to file
            exit
        else:
            igraph.plot(
                cg.commGraph,
                target="graph.svg",
                edge_label=[edge for edge in cg.commGraph.es["weight"]],
                vertex_label=[x for x in range(len(cg.commGraph.vs))],
            )

            # for line in cg.commMatrix:
            # print(line)

            hostgraph = generateHostGraph(list(map(lambda s: s.strip("'"), cg.hostnames)))
            hostgraphlayout = hostgraph.graph.layout_reingold_tilford(mode="in", root=hostgraph.layers[0])

            igraph.plot(
                hostgraph.graph,
                target="hostgraph.svg",
                margin=50,
                layout=hostgraphlayout,
                edge_label=[edge for edge in hostgraph.graph.es["weight"]],
                vertex_label=[n[-3:] for n in hostgraph.graph.vs["name"]],
                vertex_size=20,
                vertex_label_dist=1.5,
                vertex_label_angle=math.pi / 4,
            )

            print(optimize("treeMatch", cg.commMatrix, hostgraph, list(map(lambda s: s.strip("'"), cg.hostnames))))

            for line in cg.commMatrix:
                print(line)

    # we have an input matrix and topology, optimize and output reordering
    elif args.cg is not None and args.tg is not None:
        exit
