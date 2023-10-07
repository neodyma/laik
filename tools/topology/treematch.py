import itertools
import igraph
from toptypes import HostGraph, getNodeChildren, UINT64_MAX


class TreeMatch:
    def __init__(self, comm_mat: list, top_graph: HostGraph, hostnames: list[str]) -> None:
        self.comm_mat = comm_mat
        self.mod_mat = comm_mat
        self.top_graph = top_graph
        self.hostnames = hostnames

    def __str__(self) -> str:
        return (
            "TreeMatch\n"
            + "mat     {}\n".format(self.comm_mat)
            + "graph   {}\n".format(self.top_graph)
            + "hosts   {}".format(self.hostnames)
        )

    def arity(self, depth):
        # TODO: support different arities of unbalanced trees by adding node parameter
        node = 0
        return len(getNodeChildren(self.top_graph, self.top_graph.layers[depth][node], depth))

    def solve(self):
        # return self.extendCommMatrix(1)
        # return self.aggregateCommMatrix((list(range(8)), list(range(8,16))))
        # return self.groupProcesses(2)
        return [self.hostnames.index(x) for x in self.doTreeMatch()]

    def doTreeMatch(self):
        groups = [[] for _ in range(len(self.top_graph.layers))]
        print("Starting TreeMatch for {} groups.".format(len(groups)))

        # process the tree bottom-up
        for d in range(len(self.top_graph.layers) - 1, 0, -1):
            p = len(self.mod_mat)
            # extend communication matrix if needed
            if p % self.arity(d - 1) != 0:
                self.mod_mat = self.extendCommMatrix(d)
            # find suitable grouping of processes at current tree level
            groups[d] = self.groupProcesses(d, 0.5)
            # aggregate matrix by groups for next round of optimization
            self.mod_mat = self.aggregateCommMatrix(groups[d])
        return [x for xs in groups[-1] for x in xs]

    def groupProcesses(self, cur_depth: int, percentile: float):
        # TODO: this only works for balanced trees for now
        # print("processes         ", len(self.hostnames))
        # print("groupsize         ", self.arity(cur_depth - 1))
        # TODO: this should only work with connected nodes instead of the entire layer

        # all possible combinations of elements at the current layer
        l = list(itertools.combinations(self.top_graph.layers[cur_depth], self.arity(cur_depth - 1)))
        print("grouping process combinations with arity {} at depth {}.".format(self.arity(cur_depth - 1), cur_depth))

        # print("combinations      ", l)
        # print("num_combinations  ", len(l))
        # d = dict(list(map(lambda x: (hash(x) & UINT64_MAX, x), l)))
        G = igraph.Graph()
        G.add_vertices(list(map(lambda x: str(hash(x) & UINT64_MAX), l)))
        for group1 in l:
            for group2 in l:
                if group1 == group2:
                    continue
                for node in group1:
                    if node in group2:
                        G.add_edge(str(hash(group1) & UINT64_MAX), str(hash(group2) & UINT64_MAX))

        # print("graph_groups      ", len(G.vs))
        # print("groups to find    ", len(self.comm_mat) // self.arity(cur_depth - 1))
        # TODO: improve independent set algorithm

        # find independent sets of groups for possible solutions
        independent_set = G.independent_vertex_sets(
            len(self.comm_mat) // self.arity(cur_depth - 1), len(self.comm_mat) // self.arity(cur_depth - 1)
        )

        print(
            "found {} possible combinations out of {} groups at depth {}.".format(
                len(independent_set), len(G.vs), cur_depth
            )
        )

        # print("node children of {}: {}".format(self.top_graph.layers[cur_depth][0],
        # getNodeChildren(self.top_graph, self.top_graph.layers[cur_depth][0], cur_depth),))

        # print("independentset[0] ", independent_set[0], l[independent_set[0][0]], l[independent_set[0][1]])
        # TODO: filter the sets greedy
        return self.optimizeGroups(l, independent_set, percentile)

    def optimizeGroups(self, full_set, independent_set, percentile):
        if len(independent_set) < 1:
            return []
        
        min_group = []
        min_weight = 0xFFFFFFFF

        optimal_weight = 0  # TODO: this should be the maximum achievable weight over all groups
        cutoff = optimal_weight * percentile

        iterations = 0
        for matching in independent_set:
            iterations += 1
            cur_weight = self.groupWeight([full_set[g] for g in matching])
            if cur_weight < min_weight:
                min_group = [full_set[g] for g in matching]
                min_weight = cur_weight
            # elif cur_weight < cutoff:
                # min_group = [full_set[g] for g in matching]
                # break

        # print("found matching within {}-percentile after {} iterations.".format(int(100*percentile), iterations))
        return min_group
    
    def findPossibleWeight(self, tuple_size):

        return 0

    # return the number of cross-node groups
    # TODO factor in the actual transfer cost in this
    # weight should just be the total communication cost of this group
    # TODO this only works on leaf level like this
    def groupWeight(self, group):
        weight = 0
        if len(group) == 0:
            return weight

        for match in group:
            # print("finding weight for match {}".format(match))
            weight -= self.elementSum([match], 0, 0)  # * weight(link)
            # print("calculated weight {} for set {}".format(self.elementSum([match], 0, 0), match))

        # print("calculated weight {} for group {}".format(weight, group))
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
