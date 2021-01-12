import itertools
from collections import deque
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc

from .._util import CapitalData


class Tree_Alignment:
    def __init__(self):
        self.__successors1 = None
        self.__postorder1 = None
        self.__tree1 = None
        self.__successors2 = None
        self.__postorder2 = None
        self.__tree2 = None
        self.__forestdistance = None
        self.__traceforest = None
        self.__treedistance = None
        self.__tracetree = None
        self.__alignmentcost = None

    def tree_alignment(
        self,
        adata1,
        adata2,
        cost=1.0,
        N_1=2000,
        N_2=2000
    ):

        COST = cost
        gene_list = self.sort_data(
            adata1, adata2, N_1, N_2)

        adata1.uns["capital"]["intersection_genes"] = np.array(
            gene_list, dtype=object)
        adata2.uns["capital"]["intersection_genes"] = np.array(
            gene_list, dtype=object)

        self._dp(adata1, adata2, gene_list, COST)
        alignedtree = self._traceback()

        path_cluster_list = []
        source_node = list(nx.topological_sort(alignedtree))[0]
        for node in list(alignedtree.nodes):
            if alignedtree.out_degree(node) == 0:
                cluster_list = nx.shortest_path(
                    alignedtree, source=source_node, target=node)
                route1 = [i[0] for i in cluster_list]
                route2 = [i[1] for i in cluster_list]
                path_cluster_list.append([route1, route2])

        alignmentdict = {"alignment{:03d}".format(i):
                            {"data1": clusters[0],
                            "data2": clusters[1]}
                            for i, clusters in enumerate(path_cluster_list)}

        aligned_data = CapitalData(
            adata1.copy(),
            adata2.copy(),
            alignedtree,
            np.array([self.__alignmentcost], dtype=int),
            np.array(gene_list, dtype=object),
            alignmentdict,
        )

        return aligned_data

    def _set_initial_condition(
        self,
        data1,
        data2,
        cost=1.0
    ):

        self.__successors1 = data1.uns["capital"]["tree"]["successors"]
        self.__postorder1 = data1.uns["capital"]["tree"]["postorder"]
        self.__tree1 = nx.convert_matrix.from_pandas_adjacency(
            data1.uns["capital"]["tree"]["tree"], create_using=nx.DiGraph)
        self.__successors2 = data2.uns["capital"]["tree"]["successors"]
        self.__postorder2 = data2.uns["capital"]["tree"]["postorder"]
        self.__tree2 = nx.convert_matrix.from_pandas_adjacency(
            data2.uns["capital"]["tree"]["tree"],create_using=nx.DiGraph)

        # get combination of children
        # D(F1[i],F2[j]) is stored in forestdistance.loc[i,j]
        # D(F1[i1,i2],F2[j]) is stored in forestdistance.loc["(i1,i2)",j]
        # D({T1[i]},F2[j]) is stored in forestdistance.loc["(i,)", j]
        forest1_combinations = []
        for child in self.__successors1.values():
            if child.size == 1:
                children = list(itertools.combinations(child, 1))
                forest1_combinations.extend(children)
            elif child.size >= 1:
                for k in range(1, child.size):
                    children = list(itertools.combinations(child, k))
                    forest1_combinations.extend(children)

        forest2_combinations = []
        for child in self.__successors2.values():
            if child.size == 1:
                children = list(itertools.combinations(child, 1))
                forest2_combinations.extend(children)
            elif child.size >= 1:
                for k in range(1, child.size):
                    children = list(itertools.combinations(child, k))
                    forest2_combinations.extend(children)

        forest1 = [i for i in list(self.__tree1.nodes)] + \
            forest1_combinations + ["#"]
        forest2 = [j for j in list(self.__tree2.nodes)] + \
            forest2_combinations + ["#"]
        forest1 = list(map(str, forest1))
        forest2 = list(map(str, forest2))
        forest = pd.DataFrame(index=forest1, columns=forest2)
        forest.loc["#", "#"] = 0

        tree = pd.DataFrame(
            index=list(map(str, list(self.__tree1))) + ["#"],
            columns=list(map(str, list(self.__tree2))) + ["#"])
        tree.loc["#", "#"] = 0

        self.__forestdistance = forest
        self.__traceforest = pd.DataFrame(index=forest1, columns=forest2)
        self.__treedistance = tree
        self.__tracetree = pd.DataFrame(
            index=list(map(str, list(self.__tree1))) + ["#"],
            columns=list(map(str, list(self.__tree2))) + ["#"])

        COST = cost

        for i in self.__postorder1:
            size, successors = self._get_successors(self.__successors1, i)
            if size == 1:
                # D(F1[i],θ) = Σ D(T1[ik],θ)
                self._setF(i, "#", self._getT(successors[0], "#"))
                self._setT(i, "#", self._getF(i, "#") + COST)
            else:
                # D(F1[i],θ) = Σ D(T1[ik],θ)
                tmp = 0
                for ichild in successors:
                    tmp += self._getT(ichild, "#")
                self._setF(i, "#", tmp)

                # D({T1[ip],...,T1[iq]},θ) = D(T1[ip],θ) + ... + D(T1[iq],θ)
                for k in range(1, size):
                    children = list(itertools.combinations(successors, k))

                    for ichild in children:
                        tmp = 0
                        for k in ichild:
                            tmp += self._getT(k, "#")
                        self._setF(ichild, "#", tmp)
                self._setT(i, "#", self._getF(i, "#") + COST)

        for j in self.__postorder2:
            size, successors = self._get_successors(self.__successors2, j)
            if size == 1:
                self._setF("#", j, self._getT("#", successors[0]))
                self._setT("#", j, self._getF("#", j) + COST)
            else:
                tmp = 0
                for jchild in successors:
                    tmp += self._getT("#", jchild)
                self._setF("#", j, tmp)

                for k in range(1, size):
                    children = list(itertools.combinations(successors, k))

                    for jchild in children:
                        tmp = 0
                        for k in jchild:
                            tmp += self._getT("#", k)
                        self._setF("#", jchild, tmp)
                self._setT("#", j, self._getF("#", j) + COST)

    @property
    def forestdistance(self):
        return self.__forestdistance

    @property
    def traceforest(self):
        return self.__traceforest

    @property
    def tracetree(self):
        return self.__tracetree

    @property
    def treedistance(self):
        return self.__treedistance

    def sort_data(
        self,
        adata1,
        adata2,
        N_1=None,
        N_2=None
    ):
        if N_1 is not None:
            adata1 = adata1.raw.to_adata()
            sc.pp.highly_variable_genes(adata1, n_top_genes=N_1)
            adata1 = adata1[:, adata1.var['highly_variable']]
        elif N_1 is None:
            pass

        if N_2 is not None:
            adata2 = adata2.raw.to_adata()
            sc.pp.highly_variable_genes(adata2, n_top_genes=N_2)
            adata2 = adata2[:, adata2.var['highly_variable']]
        elif N_2 is None:
            pass

        s1 = set(adata1.var.index)
        s2 = set(adata2.var.index)
        intersection_list = list(s1.intersection(s2))

        if len(intersection_list) < 2:
            raise ValueError("highly variable genes of intersection of data1 and data2 are not enough "\
                                "to calculate the cost of a tree alignment. \n"\
                                "Specify num_genes1 and num_genes2 carefully.")

        print("{} genes are used to calculate cost of tree alignment.\n".format(
            len(intersection_list)))

        return intersection_list

    # cluster_centroid: pd.DataFrame
    # index is cluster name, columns is gene name, X is gene expression level
    def _calculate_cluster_centroid_for_genes(
        self,
        adata,
        gene_list,
    ):
        groupby = adata.uns["capital"]["tree"]["annotation"]
        filtered_data = adata.raw.to_adata()[:, gene_list]
        cluster_centroid_data = np.empty((0, filtered_data.n_vars))
        clustername = filtered_data.obs[groupby].unique().tolist()

        for i in clustername:
            a_cluster_data = filtered_data[filtered_data.obs[groupby] == "{}".format(
                i)].to_df()
            a_cluster_median = a_cluster_data.median(axis=0).values
            cluster_centroid_data = np.vstack(
                (cluster_centroid_data, a_cluster_median)
            )
        cluster_centroid = pd.DataFrame(
            cluster_centroid_data,
            index=clustername,
            columns=filtered_data.var_names
        )
        return cluster_centroid
    # return length of i's children and tuple of children

    def _get_successors(self, successors, i):
        size = successors[i].size
        successor = tuple([str(k) for k in successors[i]])
        if len(successor) == 0:
            successor = ("#")

        return size, successor

    def _setF(self, i, j, distance):
        if isinstance(i, tuple):
            if len(i) == 0:
                i = "#"
            i = str(i)

        if isinstance(j, tuple):
            if len(j) == 0:
                j = "#"
            j = str(j)

        if i not in self.__forestdistance.index:
            print("Error: {} does not exist in forestdistance index.".format(i))

        if j not in self.__forestdistance.columns:
            print("Error: {} does not exist in forestdistance columns.".format(j))

        self.__forestdistance.loc[i, j] = distance

    def _settraceF(self, i, j, trace):
        if isinstance(i, tuple):
            if len(i) == 0:
                i = "#"
            i = str(i)

        if isinstance(j, tuple):
            if len(j) == 0:
                j = "#"
            j = str(j)
        if i not in self.__traceforest.index:
            print("Error: {} does not exist in traceforest index.".format(i))

        if j not in self.__traceforest.columns:
            print("Error: {} does not exist in traceforest columns.".format(j))
        self.__traceforest.loc[i, j] = trace

    def _settraceT(self, i, j, trace):
        if isinstance(i, tuple):
            if len(i) == 0:
                i = "#"
        i = str(i)

        if isinstance(j, tuple):
            if len(j) == 0:
                j = "#"
        j = str(j)

        if i not in self.__tracetree.index:
            print("Error: {} does not exist in tracetree index.".format(i))

        if j not in self.__tracetree.columns:
            print("Error: {} does not exist in tracetree columns.".format(j))

        self.__tracetree.loc[i, j] = trace

    def _setT(self, i, j, distance):
        if isinstance(i, tuple):
            if len(i) == 0:
                i = "#"
        i = str(i)

        if isinstance(j, tuple):
            if len(j) == 0:
                j = "#"
        j = str(j)

        if i not in self.__treedistance.index:
            print("Error: {} does not exist in treedistance index.".format(i))

        if j not in self.__treedistance.columns:
            print("Error: {} does not exist in treedistance columns.".format(j))

        self.__treedistance.loc[i, j] = distance

    def _getF(self, i, j, parent1="Nan", parent2="Nan"):
        if isinstance(i, tuple):
            if len(i) == 0:
                i = "#"
        i = str(i)

        if isinstance(j, tuple):
            if len(j) == 0:
                j = "#"
        j = str(j)

        if i not in self.__forestdistance.index:
            i = str(parent1)

        if j not in self.__forestdistance.columns:
            j = str(parent2)

        F = self.__forestdistance.loc[i, j]

        return F

    def _getT(self, i, j):
        i = str(i)
        j = str(j)

        if i not in self.__treedistance.index:
            print("Error: TreeDistance index called does not exist.")

        if j not in self.__treedistance.columns:
            print("Error: TreeDistance columns called does not exist.")

        T = self.__treedistance.loc[i, j]

        return T

    def _cal1(self, A, B):
        if not isinstance(A, tuple):
            _, A = self._get_successors(self.__successors1, A)

        if not isinstance(B, tuple):
            _, B = self._get_successors(self.__successors2, B)

        mintemp = 1000
        trace = "Nan"
        temp = 0

        for k in A:
            for l in B:
                Asub = tuple([i for i in A if i != k])
                Bsub = tuple([j for j in B if j != l])

                temp = self._getF(Asub, Bsub) + self._getT(k, l)

                if mintemp > temp:
                    mintemp = temp
                    if Bsub == ():
                        Bsub = "#"
                    if l == ():
                        l = "#"
                    if Asub == ():
                        Asub = "#"
                    if k == ():
                        k = "#"

                    trace = [1, [Asub, Bsub], [k, l]]

        return mintemp, trace

    def _cal2(self, A, B, cost):
        COST = cost
        parentA = "Nan"
        parentB = "Nan"
        if not isinstance(A, tuple):
            parentA = A
            _, A = self._get_successors(self.__successors1, A)

        if not isinstance(B, tuple):
            parentB = B
            _, B = self._get_successors(self.__successors2, B)

        Bprime = [()]
        for m in range(1, len(B)+1):
            for Bp in list(itertools.combinations(B, m)):
                Bprime.extend([Bp])

        mintemp = 1000
        trace = "Nan"
        temp = 0

        for k in A:
            for Bp in Bprime:
                Asub = tuple([i for i in A if i != k])
                Bsub = tuple([j for j in B if not j in Bp])
                if k == "#":
                    temp = self._getF(Asub, Bsub, parent2=parentB) + \
                        self._getF("#", Bp, parent2=parentB)
                else:
                    temp = self._getF(Asub, Bsub, parent2=parentB) + \
                        self._getF(k, Bp, parent2=parentB) + COST

                if mintemp > temp:
                    mintemp = temp
                    if Bsub == ():
                        Bsub = "#"
                    if Bp == ():
                        Bp = "#"
                    if Asub == ():
                        Asub = "#"
                    if k == ():
                        k = "#"
                    if not "{}".format(Bsub) in self.__forestdistance.columns:
                        Bsub = parentB
                    if not "{}".format(Bp) in self.__forestdistance.columns:
                        Bp = parentB

                    trace = [2, [[Asub, Bsub], [k, Bp]], (k, "#")]

        return mintemp, trace

    # min D(A-Aprime,B-{T2[jq]}) + D(Aprime,F2[jq]) + cost
    def _cal3(self, A, B, cost):
        COST = cost
        parentA = "Nan"
        parentB = "Nan"

        if not isinstance(A, tuple):
            parentA = A
            _, A = self._get_successors(self.__successors1, A)

        if not isinstance(B, tuple):
            parentB = B
            _, B = self._get_successors(self.__successors2, B)

        Aprime = [()]
        for m in range(1, len(A)+1):
            for Ap in list(itertools.combinations(A, m)):
                Aprime.extend([Ap])

        mintemp = 1000
        trace = "Nan"
        temp = 0

        for l in B:
            for Ap in Aprime:

                Asub = tuple([i for i in A if i not in Ap])
                Bsub = tuple([j for j in B if j != l])
                temp = self._getF(Asub, Bsub, parent1=parentA) + \
                    self._getF(Ap, l, parent1=parentA) + COST

                if mintemp > temp:
                    mintemp = temp
                    if Bsub == ():
                        Bsub = "#"
                    if l == ():
                        l = "#"
                    if Asub == ():
                        Asub = "#"
                    if Ap == ():
                        Ap = "#"

                    if not "{}".format(Asub) in self.__forestdistance.index:
                        Asub = parentA
                    if not "{}".format(Ap) in self.__forestdistance.index:
                        Ap = parentA
                    trace = [3, [[Asub, Bsub], [Ap, l]], ("#", l)]

        return mintemp, trace

    # trace is a list that has 3 pairs like [a,(b),(c)]
    # a is a record that show which calculation was done
    # (b) is a pair or 2 pairs that goes to the search in traceforest for the next traceback
    # (c) is a pair that goes to stack or foreststack

    def _calculateForest(
        self,
        A,
        B,
        cost
    ):
        COST = cost
        min1, trace1 = self._cal1(A, B)
        min2, trace2 = self._cal2(A, B, COST)
        min3, trace3 = self._cal3(A, B, COST)

        forestdistance = np.array([min1, min2, min3], dtype=object).min()
        trace = [trace1, trace2, trace3][np.array([min1, min2, min3], dtype=object).argmin()]

        return forestdistance, trace

    # trace is a list that has 3 pairs like [(a),(b),(c)]
    # (a) is a record that show which calculation was done
    # (b) is a pair of the match result for D(T[i],T[j])
    # (c) is a pair that goes to the stack
    def _calculateTreedistance(
        self,
        i,
        j,
        mincost
    ):
        MINCOST = mincost
        _, successor1 = self._get_successors(self.__successors1, i)
        _, successor2 = self._get_successors(self.__successors2, j)

        distancelist = []
        for l in successor2:
            distancelist.append(
                [self._getT("#", j) + self._getT(i, l) - self._getT("#", l), [0, ("#", j), (i, l)]])
        for k in successor1:
            distancelist.append(
                [self._getT(i, "#") + self._getT(k, j) - self._getT(k, "#"), [0, (i, "#"), (k, j)]])

        distancelist.append([self._getF(i, j) + MINCOST, [1, (i, j), (i, j)]])

        array = np.array(distancelist, dtype=object)

        treedistance = array[:, 0].min()
        trace = array[array[:, 0].argmin(), 1]

        return treedistance, trace

    def _dp(
        self,
        adata1,
        adata2,
        gene_list,
        cost
    ):
        # np.warning cause warning below
        # VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences is deprecated.
        # If you meant to do this, you must specify 'dtype=object' when creating the ndarray.
        # however it causes inside pandas and to recover this wanrnig, we need to rewrite all the code above
        # and it still get the result we need. we will get this right soon.
        np.warnings.filterwarnings(
            'ignore', category=np.VisibleDeprecationWarning)

        COST = cost

        self._set_initial_condition(adata1, adata2, cost=COST)
        cluster_centroid1 = self._calculate_cluster_centroid_for_genes(
            adata1, gene_list)
        cluster_centroid2 = self._calculate_cluster_centroid_for_genes(
            adata2, gene_list)

        for i in self.__postorder1:
            for j in self.__postorder2:
                df = pd.DataFrame(
                    {"A": cluster_centroid1.loc[i], "B": cluster_centroid2.loc[j]})
                mincost = 1 - df.corr(method="spearman").iloc[0, 1]

                size1, successor1 = self._get_successors(self.__successors1, i)
                size2, successor2 = self._get_successors(self.__successors2, j)

                Alist = []
                if size1 == 0:
                    pass
                elif size1 == 1:
                    Alist.extend([(successor1)])
                else:
                    for m in range(1, len(successor1)):
                        for l in list(itertools.combinations(successor1, m)):
                            Alist.extend([l])

                Alist.extend([i])
                if size2 == 0:
                    pass
                elif size2 == 1:
                    for A in Alist:
                        fdistance, ftrace = self._calculateForest(
                            A, (successor2), COST)
                        self._setF(A, (successor2), fdistance)
                        self._settraceF(A, (successor2), ftrace)
                else:
                    for m in range(1, len(successor2)):
                        for B in list(itertools.combinations(successor2, m)):
                            for A in Alist:
                                fdistance, ftrace = self._calculateForest(
                                    A, B, COST)
                                self._setF(A, B, fdistance)
                                self._settraceF(A, B, ftrace)

                for A in Alist:
                    fdistance, ftrace = self._calculateForest(A, j, COST)
                    self._setF(A, j, fdistance)
                    self._settraceF(A, j, ftrace)

                tdistance, ttrace = self._calculateTreedistance(i, j, mincost)
                self._setT(i, j, tdistance)
                self._settraceT(i, j, ttrace)

    def _traceback(self):
        G = nx.DiGraph()
        G.add_node("tempnode")
        parent = "tempnode"
        stack = deque()
        stack.append(
            (self.__postorder1[-1], self.__postorder2[-1], parent))

        while len(stack) != 0:
            i, j, parent = stack.pop()

            if i == "#" and j == "#":
                continue

            elif i == "#":
                if j != "#":
                    H = nx.dfs_tree(self.__tree2, j)
                    G = nx.compose(G, H)
                    G.add_edge(parent, j)
                    G = nx.relabel_nodes(
                        G, dict(zip(list(H.nodes), [("#", k) for k in H.nodes])))
                continue

            elif j == "#":
                if i != "#":
                    H = nx.dfs_tree(self.__tree1, i)
                    G = nx.compose(G, H)
                    G.add_edge(parent, i)
                    G = nx.relabel_nodes(
                        G, dict(zip(list(H.nodes), [(k, "#") for k in H.nodes])))
                continue

            elif i != "#" and j != "#":
                tree_result = self.__tracetree.loc[i, j]
                forest_result = self.__traceforest.loc[i, j]

                if tree_result[0] == 0:
                    if tree_result[1][0] == "#" and tree_result[1][1] != "#":
                        H = nx.dfs_tree(self.__tree2,
                                        source=tree_result[1][1])
                        Hprime = nx.dfs_tree(
                            self.__tree2, source=tree_result[2][1])
                        H.remove_nodes_from(list(Hprime.nodes))
                        G = nx.compose(G, H)
                        G.add_edge(parent, tree_result[1][1])
                        G = nx.relabel_nodes(
                            G, dict(zip(list(H.nodes), [("#", k) for k in H.nodes])))

                    elif tree_result[1][0] != "#" and tree_result[1][1] == "#":
                        H = nx.dfs_tree(self.__tree1,
                                        source=tree_result[1][0])
                        Hprime = nx.dfs_tree(
                            self.__tree1, source=tree_result[2][0])
                        H.remove_nodes_from(list(Hprime.nodes))
                        G = nx.compose(G, H)
                        G.add_edge(parent, tree_result[1][0])
                        G = nx.relabel_nodes(
                            G, dict(zip(list(H.nodes), [(k, "#") for k in H.nodes])))

                    stack.append(
                        [tree_result[2][0], tree_result[2][1], tree_result[1]])

                elif tree_result[0] == 1:
                    G.add_node(tree_result[1])
                    G.add_edge(parent, tree_result[1])

                    foreststack = deque()
                    foreststack.append([i, j, (i, j)])

                    while len(foreststack) != 0:
                        i_f, j_f, p_f = foreststack.pop()

                        if i_f == "#" and j_f == "#":
                            continue

                        elif i_f == "#":
                            if j_f != "#":
                                for j_tmp in j_f:
                                    H = nx.dfs_tree(
                                        self.__tree2, source=j_tmp)
                                    G = nx.compose(G, H)
                                    G.add_edge(p_f, j_tmp)
                                    G = nx.relabel_nodes(
                                        G, dict(zip(list(H.nodes), [("#", k) for k in H.nodes])))

                        elif j_f == "#":
                            if i_f != "#":
                                for i_tmp in i_f:
                                    H = nx.dfs_tree(
                                        self.__tree1, source=i_tmp)
                                    G = nx.compose(G, H)
                                    G.add_edge(p_f, i_tmp)
                                    G = nx.relabel_nodes(
                                        G, dict(zip(list(H.nodes), [(k, "#") for k in H.nodes])))

                        elif i_f != "#" and j_f != "#":
                            i_f = "{}".format(i_f)
                            j_f = "{}".format(j_f)

                            forest_result = self.__traceforest.loc[i_f, j_f]

                            if forest_result[0] == 1:
                                stack.append(
                                    [forest_result[2][0], forest_result[2][1], p_f])
                                foreststack.append(
                                    [forest_result[1][0], forest_result[1][1], p_f])

                            elif forest_result[0] == 2:
                                foreststack.append(
                                    [forest_result[1][0][0], forest_result[1][0][1], p_f])
                                if forest_result[1][1][1] != "#":
                                    G.add_node(forest_result[2])
                                    G.add_edge(p_f, forest_result[2])
                                    foreststack.append(
                                        [forest_result[1][1][0], forest_result[1][1][1], forest_result[2]])

                                elif forest_result[1][1][1] == "#":
                                    H = nx.dfs_tree(
                                        self.__tree1, forest_result[1][1][0])
                                    G = nx.compose(G, H)
                                    G.add_edge(p_f, forest_result[1][1][0])
                                    G = nx.relabel_nodes(
                                        G, dict(zip(list(H.nodes), [(k, "#") for k in H.nodes])))

                            elif forest_result[0] == 3:
                                foreststack.append(
                                    [forest_result[1][0][0], forest_result[1][0][1], p_f])

                                if forest_result[1][1][0] != "#":
                                    G.add_node(forest_result[2])
                                    G.add_edge(p_f, forest_result[2])
                                    foreststack.append(
                                        [forest_result[1][1][0], forest_result[1][1][1], forest_result[2]])

                                elif forest_result[1][1][0] == "#":
                                    H = nx.dfs_tree(
                                        self.__tree2, forest_result[1][1][1])
                                    G = nx.compose(G, H)
                                    G.add_edge(p_f, forest_result[1][1][1])
                                    G = nx.relabel_nodes(
                                        G, dict(zip(list(H.nodes), [("#", k) for k in H.nodes])))

        G.remove_node("tempnode")
        alignmentcost = self.__treedistance.loc[self.__postorder1[-1],
                                                self.__postorder2[-1]]

        self.__alignmentcost = alignmentcost/len(G)

        return G
