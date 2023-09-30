#!/usr/bin/python3
import argparse
from functools import reduce
from toptypes import CommStats, HostGraph, getNodeChildren
from treematch import TreeMatch
import igraph
import itertools
import math
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
                    hostnames[int(output.group(2))] = output.group(1) # type: ignore
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
                        commMatrix[int(content[0])] = [sum(x) for x in zip(commMatrix[int(content[0])], values)] # type: ignore

    # undirected: double transfer values!
    # commGraph = igraph.Graph.Weighted_Adjacency(commMatrix, mode="undirected")
    print(f"Parsed matrices from {len(logfiles)} files.")
    commGraph = igraph.Graph.Weighted_Adjacency(commMatrix)
    return CommStats(commGraph, commMatrix, hostnames) # type: ignore


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

    return HostGraph(topGraph, layers, weights)


# solve the embedding problem for given graphs and return acceptable reordering
# todo we probably need the hostnames sooner or later
#   for now on test we always get successive hosts
def reorder(comm_graph, top_graph) -> list:
    return [1, 2, 3]


def treeMatch(comm_mat, top_graph, hostnames) -> list:
    solver = TreeMatch(comm_mat, top_graph, hostnames)
    return solver.solve()


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
            exit(1)
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

            # for line in cg.commMatrix:
                # print(line)

    # we have an input matrix and topology, optimize and output reordering
    elif args.cg is not None and args.tg is not None:
        exit(1)
