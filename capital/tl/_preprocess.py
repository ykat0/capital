import math
import os
import sys
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata
from scipy import stats

def read_file(filename, transpose=False):
    adata = None
    if os.path.exists(filename):
        if os.path.isdir(filename):
            adata = sc.read_10x_mtx(filename)

        elif os.path.isfile(filename):
            name, filetype = os.path.splitext(filename)
            if filetype == ".txt":
                print()
                adata = sc.read_text(filename)

            if filetype == ".csv":
                adata = sc.read_csv(filename)

            if filetype == ".h5ad":
                adata = sc.read(filename)

        else:
            print(
                "ERROR: the format must be [H5AD|CSV|TXT] for file or 10x-MTX for directory.")
            sys.exit()

        if transpose:
            adata = adata.transpose()
    elif not os.path.exists(filename):
        sys.exit("ERROR: no such file or directory.")

    if not isinstance(adata.X, np.ndarray):
        X = adata.X.toarray()
        adata = anndata.AnnData(X, obs=adata.obs, var=adata.var)
    return adata

class Preprocessing:
    def __init__(self):
        pass

    def preprocessing_rawdata(
        self,
        adata,
        Min_Genes=200,
        Min_Cells=3,
        Min_Mean=0.0125,
        Max_Mean=3,
        Min_Disp=0.5,
        N_pcs=50,
        n_Top_genes=2000,
        K=10,
        magic_imputation=False,
        copy=False
    ):

        # adata = adata.copy() if copy else adata
        sc.pp.filter_cells(adata, min_genes=Min_Genes)
        sc.pp.filter_genes(adata, min_cells=Min_Cells)
        sc.pp.normalize_total(adata, exclude_highly_expressed=True)
        sc.pp.log1p(adata)

        if magic_imputation is True:
            #     set np.random.seed for magic to create same imputation
            np.random.seed(1)
            sce.pp.magic(
                adata,
                name_list='all_genes',
                knn=5
            )

        sc.pp.highly_variable_genes(
            adata,
            min_mean=Min_Mean,
            max_mean=Max_Mean,
            min_disp=Min_Disp,
            n_top_genes=n_Top_genes
        )
        adata.raw = adata
        # same as adata = adata[:,adata.var['highly_variable']] but inplace
        adata._inplace_subset_var(adata.var['highly_variable'])
        sc.tl.pca(
            adata,
            n_comps=N_pcs
        )
        sc.pp.neighbors(adata,
                        n_neighbors=K,
                        n_pcs=N_pcs
                        )
        sc.tl.diffmap(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata)
        sc.tl.paga(adata, groups='leiden')

        return adata if copy else None

    # get source node of tree
    # calcurating variance of gene expression level and
    # set node which has max variance as source node
    def calculate_root_node(
        self,
        adata,
        threshold,
        groupby="leiden"
    ):

        cluster_centroid_data = self.calculate_cluster_centroid(
            adata,
            dimension='raw',
            groupby=groupby
        )

        cluster_centroid = pd.DataFrame(
            cluster_centroid_data, index=adata.obs[groupby].cat.categories)

        cluster_centroid = cluster_centroid.loc[:,
                                                (cluster_centroid.sum(axis=0) != 0)]
        cor, pval = stats.spearmanr(cluster_centroid)
        gene_adjacency_matrix = np.where(abs(cor) > threshold, True, False)
        scEnergy = []
        for i in range(len(cluster_centroid.index)):
            tmp = 0
            for j in range(len(cluster_centroid.columns)):
                if cluster_centroid.iloc[i, j] != 0:
                    tmp += -cluster_centroid.iloc[i, j]*math.log(
                        cluster_centroid.iloc[i, j]/cluster_centroid.iloc[i, :][gene_adjacency_matrix[j]].sum())
            scEnergy.append(tmp)
        root_node = cluster_centroid.index[scEnergy.index(max(scEnergy))]

        return root_node

    # cluster_centroid: pd.DataFrame
    # index is cluster name, columns is gene name, X is gene expression level
    # argument "dimension" is either pca or diffmap

    def calculate_cluster_centroid(
        self,
        adata,
        dimension="pca",
        groupby="leiden"
    ):

        if dimension == "pca":
            X_dimension = "X_pca"
        elif dimension == "diffmap":
            X_dimension = "X_diffmap"
        elif dimension == "raw":
            X_dimension = "raw"
        else:
            raise ValueError(
                "Argument 'dimension' must be 'pca' or 'diffmap'.")

        if dimension in ["pca", "diffmap"]:
            clustername = adata.obs[groupby].cat.categories
            cluster_centroid_data = np.empty(
                (0, adata.obsm[X_dimension].shape[1]))
            for i in clustername:
                a_cluster_data = pd.DataFrame(
                    adata[adata.obs[groupby] == "{}".format(i)].obsm[X_dimension])
                a_cluster_median = a_cluster_data.median(axis=0).values
                cluster_centroid_data = np.vstack(
                    (cluster_centroid_data, a_cluster_median))
        else:
            clustername = adata.obs[groupby].cat.categories
            cluster_centroid_data = np.empty(
                (0, adata.X.shape[1]))
            for i in clustername:
                a_cluster_data = pd.DataFrame(
                    adata[adata.obs[groupby] == "{}".format(i)].X)
                a_cluster_median = a_cluster_data.median(axis=0).values
                cluster_centroid_data = np.vstack(
                    (cluster_centroid_data, a_cluster_median))

        return cluster_centroid_data

    def _convert_to_tree(
        self,
        adata,
        root_node,
        cluster_centroid_data,
        method,
        groupby="leiden"
    ):
        if root_node in adata.obs[groupby].cat.categories:
            pass
        else:
            raise ValueError(
                "Root {} does not exist in the clusters. Choose the root again.".format(root_node))

        if method == "paga":
            Adjacency_matrix = np.array(
                adata.uns['paga']["connectivities_tree"].todense()
            )

            G = nx.from_numpy_matrix(Adjacency_matrix, create_using=nx.Graph())

            # get names of cluster from adata.uns['paga']['groups'] as string
            groups_key = adata.uns['paga']['groups']
            groups_names = adata.obs[groups_key].cat.categories

            # relabel nodes name as names set in adata.obs of adata.uns['paga']['groups']
            G = nx.relabel_nodes(G, dict(zip(G, groups_names)))

        elif method == "gauss":
            from scipy.spatial.distance import pdist, squareform
            ka = 6
            y = pdist(cluster_centroid_data, 'euclid')
            Y = squareform(y)
            sigma = np.sort(Y, axis=1)[:, ka]
            affinity_matrix = np.exp(np.square(Y/sigma) * -1)
            matrix = 1 - affinity_matrix
            G = nx.from_numpy_matrix(matrix, create_using=nx.Graph())
            G = nx.relabel_nodes(
                G, dict(zip(G, adata.obs[groupby].cat.categories)))

        elif method == "euclid":
            from scipy.spatial.distance import pdist, squareform
            y = pdist(cluster_centroid_data, 'euclid')
            Y = squareform(y)
            G = nx.from_numpy_matrix(Y, create_using=nx.Graph())
            G = nx.relabel_nodes(
                G, dict(zip(G, adata.obs[groupby].cat.categories)))

        T = nx.minimum_spanning_tree(G)
        sorted(T.edges(data=True))
        tree = nx.dfs_tree(T, root_node)

        return tree

    def _set_postorder_successors(
        self,
        tree,
        root_node
    ):

        postorder_tmp = list(nx.dfs_postorder_nodes(tree, source=root_node))
        postorder = np.array([str(i) for i in postorder_tmp], dtype=object)
        successors = dict(
                        zip(
                            map(str, list(tree.nodes)),
                            [np.array(list(tree.successors(i)), dtype=object) for i in tree]
                            )
                        )

        return postorder, successors

    # combining methods above
    def trajectory_tree(
        self,
        adata,
        root_node=None,
        tree=None,
        groupby="leiden",
        method="euclid",
        dimension="pca",
        copy=False
    ):

        # if adata.X is scipy.sparsematrix, adata.X is converted to numpy array.
        if not isinstance(adata.X, np.ndarray):
            adata_tmp = adata.copy()
            adata_tmp.X = adata.X.toarray()

        else:
            adata_tmp = adata.copy()
        # add numpy array:cluster centroids matrix to adata.uns["cluster_centroid"]
        # index of cluster_centorid is adata.obs["leiden"].cat.categories
        cluster_centroid_data = self.calculate_cluster_centroid(
                                adata_tmp,
                                dimension=dimension,
                                groupby=groupby
                            )
        adata.uns["cluster_centroid"] = cluster_centroid_data

        # add networkx.classes.digraph.DiGraph of trajectory to adata.uns["tree"]["tree"]
        # add postorder list to adata.uns["tree"]["postorder"]
        # add successors dictionary to adata.uns["tree"]["successors"]
        # add root node str to adata.uns["tree"]["root_node"]
        adata.uns["capital"] = {}
        adata.uns["capital"]["tree"] = {}

        tree_dict = adata.uns["capital"]["tree"]
        tree_dict["annotation"] = groupby

        if tree is None:
            if root_node is not None:
                tree_dict["root_node"] = root_node
            else:
                root_node = self.calculate_root_node(
                    adata_tmp, threshold=0.5, groupby=groupby)
                tree_dict["root_node"] = root_node
                print("{} is set as a root node".format(root_node))

            tree = self._convert_to_tree(
                    adata_tmp,
                    root_node,
                    cluster_centroid_data,
                    method,
                    groupby=groupby
                    )

            postorder, successors = self._set_postorder_successors(tree, root_node)

            tree_dict["tree"] = nx.convert_matrix.to_pandas_adjacency(
                tree, nodelist=tree.nodes())
            tree_dict["postorder"] = postorder
            tree_dict["successors"] = successors

        else:
            if type(tree) is nx.DiGraph:
                root_node_tmp = nx.topological_sort(tree)[0]
                if root_node_tmp != root_node:
                    raise ValueError("root node of tree and passed argument 'root_node' did not match.")
                else:
                    root_node = root_node_tmp
                postorder, successors = self._set_postorder_successors(tree, root_node)
            elif type(tree) is nx.Graph and root_node is not None:
                tree = nx.dfs_tree(tree, root_node)
                postorder, successors = self._set_postorder_successors(tree, root_node)
            else:
                raise ValueError("Argument 'tree' must be nx.DiGraph or nx.Graph with 'root_node' argument passed")

            tree_dict["root_node"] = root_node
            tree_dict["tree"] = nx.convert_matrix.to_pandas_adjacency(
                tree, nodelist=tree.nodes())
            tree_dict["postorder"] = postorder
            tree_dict["successors"] = successors

        return adata if copy else None
