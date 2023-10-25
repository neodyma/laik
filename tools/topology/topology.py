#!/usr/bin/python3
import argparse
from functools import reduce
import random
from toptypes import CommStats, HostGraph, getNodeChildren
from treematch import TreeMatch
from qap import TauQAP
import igraph
import itertools
import math
import numpy as np
import re


# LAIK_LOG_FILE -> commGraph, commMatrix, hostnames
def parseCommStats(logfiles: list) -> CommStats:
    commMatrix = list

    str_marker = "Communication Matrix"
    str_backend = "MPI backend initialized"
    R_BE = re.compile(r"MPI backend initialized \(at (.+), rank (\d+)/(\d+)\).*$")
    # == LAIK-(0000)-L(03) (0006).(01)  (0):(00).(004) | (Communication Matrix:) (\n)
    R_EX = re.compile(r".+LAIK-(\d+)-L(\d+) (\d+).(\d+)  (\d+):(\d+).(\d+) \| (.+)")

    ranks = 0
    matrices = []
    hostnames = list

    for file in logfiles:
        with open(file, "r") as logfile:
            for line in logfile:
                if str_backend in line:  # set number of ranks
                    output = R_BE.search(line)
                    if output is None:
                        continue
                    if ranks == 0:
                        ranks = int(output.group(3))
                        commMatrix = [[0 for _ in range(ranks)] for _ in range(ranks)]
                        hostnames = ["" for _ in range(ranks)]
                    hostnames[int(output.group(2))] = output.group(1)  # type: ignore
                    continue
                rex = R_EX.search(line)
                if rex is not None:
                    if str_marker in line:
                        matrices.append((rex.group(2), rex.group(3)))  # (rank, sequence)
                    elif (rex.group(2), rex.group(3)) in matrices and int(rex.group(4)) > 2:
                        # matched a matrix row
                        content = R_EX.split(line)[-2].split()
                        if content is None:
                            raise ValueError
                        if len(content) != ranks + 2:  # misformed matrix output
                            raise IndexError
                        values = [int(x) for x in content[2:]]
                        commMatrix[int(content[0])] = [sum(x) for x in zip(commMatrix[int(content[0])], values)]  # type: ignore

    # undirected: double transfer values!
    # commGraph = igraph.Graph.Weighted_Adjacency(commMatrix, mode="undirected")
    print(f"Parsed matrices from {len(logfiles)} files.")
    commGraph = igraph.Graph.Weighted_Adjacency(commMatrix)
    return CommStats(commGraph, commMatrix, hostnames)  # type: ignore


# generate a host graph based on the supermuc-ng node naming scheme
# i01r01c01s01
def generateHostTopology(hostnames: list[str]) -> HostGraph:
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
    weights = [1, 4, 0, 0, 15]
    topGraph.add_edges(node_edges, dict(weight=[weights[0] for _ in range(len(node_edges))]))
    topGraph.add_edges(srvs_edges, dict(weight=[weights[1] for _ in range(len(srvs_edges))]))
    topGraph.add_edges(cabs_edges, dict(weight=[weights[2] for _ in range(len(cabs_edges))]))
    topGraph.add_edges(racks_edges, dict(weight=[weights[3] for _ in range(len(racks_edges))]))
    topGraph.add_edges(isls_edges, dict(weight=[weights[4] for _ in range(len(isls_edges))]))

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

    topArray = np.array(topGraph.distances(weights="weight"), dtype=int)
    return HostGraph(topGraph, topArray[0 : len(hostnames), 0 : len(hostnames)].tolist(), layers, weights)


# solve the embedding problem for given graphs and return acceptable reordering
# todo we probably need the hostnames sooner or later
#   for now on test we always get successive hosts
def generate_LAIK_REORDERING(reordering: list) -> str:
    reorderstr = "LAIK_REORDERING="
    for index, node in enumerate(reordering):
        reorderstr += "{}.{},".format(index, node)

    return reorderstr[:-1]


def treeMatch(comm_mat, top_graph, hostnames) -> list:
    solver = TreeMatch(comm_mat, top_graph, hostnames)
    return solver.solve()


def tauQAP(comm_mat, top_graph, hostnames) -> list:
    solver = TauQAP(comm_mat, top_graph, hostnames)
    return solver.solve()


def optimize(optimizer, *args) -> list:
    return globals()[optimizer](*args)


# create artificial communication matrix with known ideal reordering
def generateGroupedComms(num_procs: int, procs_per_cluster: int) -> tuple:
    # ideal_order = list(itertools.permutations(range(num_procs)))[random.randrange(math.factorial(num_procs))]
    ideal_order = list(range(num_procs))
    random.shuffle(ideal_order)

    synth_matrix = [[0 for _ in range(num_procs)] for _ in range(num_procs)]
    for chunk in np.array_split(ideal_order, len(ideal_order) // procs_per_cluster):
        comm_weight = random.randrange(2**14)
        for n1 in chunk:
            for n2 in chunk:
                if n1 != n2:
                    synth_matrix[n1][n2] = comm_weight
    return synth_matrix, ideal_order

    # create clusters with high-communicating processes corresponding to procs_per_cpu
    # e.g. num_procs=8, procs_per_cpu=2 -> pairs of communicating processes that should get mapped to same cpu


def matchedReorderGroups(order1: list, order2: list, num_procs: int, procs_per_cluster: int):
    set1, set2 = set(), set()
    for chunk in np.array_split(order1, len(order1) // procs_per_cluster):
        set1.add(frozenset(list(chunk)))
    for chunk in np.array_split(order2, len(order2) // procs_per_cluster):
        set2.add(frozenset(list(chunk)))
    # set2 = frozenset(np.split(order1, len(order1) // procs_per_cluster))
    return True if len(set1.difference(set2)) == 0 else False


def parserSetup() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="LAIK topology optimizer",
        description="Parse communication metadata, create structured representations, optimize communication paths and generate reordering parameters.",
    )
    parser.add_argument("-i", "--ilog", nargs="+", help="Input LAIK_LOG_FILEs")
    parser.add_argument("-c", "--cg", help="Communication graph")
    parser.add_argument("-t", "--tg", help="Topology graph")
    parser.add_argument("-o", "--out", help="Output matrix file")
    parser.add_argument("-r", help=f"Reorder using ansatz [treeMatch, ..]")
    return parser


if __name__ == "__main__":
    parser = parserSetup()
    args = parser.parse_args()

    # we have an input logfile, let's convert it to a usable Graph
    if args.ilog is not None:
        comm_stats = parseCommStats(args.ilog)
        if args.out is not None:
            np.savetxt(args.out, np.array(comm_stats.commMatrix, dtype=int), "%10d")
        if args.cg is not None:
            # dump to file
            exit(1)
        else:
            igraph.plot(
                comm_stats.commGraph,
                target="graph.svg",
                edge_label=[edge for edge in comm_stats.commGraph.es["weight"]],
                vertex_label=[x for x in range(len(comm_stats.commGraph.vs))],
            )

            # for line in cg.commMatrix:
            # print(line)

            hostgraph = generateHostTopology(list(map(lambda s: s.strip("'"), comm_stats.hostnames)))
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

            np.set_printoptions(linewidth=168, edgeitems=4)
            print(np.array(comm_stats.commMatrix, dtype=int))

            # tm = optimize("treeMatch", comm_stats.commMatrix, hostgraph, list(map(lambda s: s.strip("'"), comm_stats.hostnames)))
            # print("treeMatch: ", tm)
            # print(generate_LAIK_REORDERING(tm))

            qap = optimize("tauQAP", comm_stats.commMatrix, hostgraph, list(map(lambda s: s.strip("'"), comm_stats.hostnames)))
            print("tauQAP:    ", qap)
            print(generate_LAIK_REORDERING(qap))

            # print(np.array(hostgraph.topMatrix, dtype=int))

            # test_comms = generateGroupedComms(16, 2)
            # res = optimize("tauQAP", test_comms[0], hostgraph, list(map(lambda s: s.strip("'"), comm_stats.hostnames)))
            # print(res)
            # print(optimize("treeMatch", test_comms[0], hostgraph, list(map(lambda s: s.strip("'"), comm_stats.hostnames))))

            # print(matchedReorderGroups(res, test_comms[1], 16, 2))

    # we have an input matrix and topology, optimize and output reordering
    elif args.cg is not None and args.tg is not None:
        exit(1)
