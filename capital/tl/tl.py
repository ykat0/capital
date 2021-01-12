import os
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
from anndata import AnnData
from pandas.api.types import is_categorical_dtype
from typing import Union, Optional, Literal

from ._preprocess import Preprocessing
from ._tree_alignment import Tree_Alignment
from ._dynamic_time_warping import DPT, DynamicTimeWarping
from .._util import CapitalData


def preprocessing(
    adata: AnnData,
    Min_Genes: int = 200,
    Min_Cells: int = 3,
    Min_Mean: float = 0.0125,
    Max_Mean: float = 3,
    Min_Disp: float = 0.5,
    N_pcs: int = 50,
    n_Top_genes: int = 2000,
    K: int = 10,
    magic_imputation: bool = False,
):
    """\
    The recipe for preprocessing raw count data.

    In adata.raw, all genes are stored, so that those genes can be used in later calculation in CAPITAL.

    The recipe runs the following steps:
    ::

        import scanpy as sc
        sc.pp.filter_cells(adata, min_genes=Min_Genes)
        sc.pp.filter_genes(adata, min_cells=Min_Cells)
        sc.pp.normalize_total(adata, exclude_highly_expressed=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=Min_Mean, max_mean=Max_Mean, min_disp=Min_Disp, n_top_genes=n_Top_genes)
        adata.raw = adata
        adata = adata[:,adata.var['highly_variable']]
        sc.tl.pca(adata, n_comps=N_pcs)
        sc.pp.neighbors(adata, n_neighbors=K, n_pcs=N_pcs)
        sc.tl.diffmap(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata)
        sc.tl.paga(adata, groups='leiden')

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
        to cells and columns to genes.
    Min_Genes : int
        The number of genes to filter in scanpy.pp.filter_cells(),
        by default 200.
    Min_Cells : int
        The number of cells to filter in scanpy.pp.filter_genes(),
        by default 3.
    Min_Mean : float
        The minimum mean that is filtered to calculate highly variable genes.
        Look scanpy.pp.highly_variable_genes(),
        by default 0.0125.
    Max_Mean : int
        The maxmum mean that is filtered to calculate highly variable genes.
        Look scanpy.pp.highly_variable_genes(),
        by default 3.
    Min_Disp : float
        The minimum dispersion that is filtered to calculate highly variable genes.
        Look scanpy.pp.highly_variable_genes(),
        by default 0.5.
    N_pcs : int
        The number of principal components used,
        by default 50.
    n_Top_genes : int
        The number of highly variable genes,
        by default 2000.
    K : int
        The size of a local neighborhood used for manifold approximation,
        by default 10.
    magic_imputation : bool
        If `True`, MAGIC imputation is done,
        by default `False`.
    """

    if not isinstance(adata, AnnData):
        raise ValueError("preprocessing() expects AnnData argument")

    pp = Preprocessing()
    pp.preprocessing_rawdata(
        adata,
        Min_Genes=Min_Genes,
        Min_Cells=Min_Cells,
        Min_Mean=Min_Mean,
        Max_Mean=Max_Mean,
        Min_Disp=Min_Disp,
        N_pcs=N_pcs,
        n_Top_genes=n_Top_genes,
        K=K,
        magic_imputation=magic_imputation,
        copy=False
    )


def trajectory_tree(
    adata: AnnData,
    root_node: Optional[str] = None,
    method: Literal["euclid", "gauss", "paga"] = "euclid",
    dimension: Literal["pca", "diffmap"] ="pca",
    tree: Optional[nx.DiGraph] = None,
    groupby: str = "leiden"
):
    """\
    Calculate a trajectory tree.

    Specify AnnData and the root name in the trajectory.

    You can pass your metadata in adata.obs
    by passing a key in adata.obs to `groupby`.

    You can also define your own tree by drawing it as nx.DiGraph
    and pass it to `tree`.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    root_node : Optional[str]
        A given cluster is used as the root for the tree.
        If `None`, it calculates a most likely root of the tree,
        but this will not be reliable for the moment.
        Please set the root cluster, by default `None`.
    method : Literal["euclid", "gauss", "paga"]
        Method for calculating the tree, by default "euclid."
    dimension : Literal["pca", "diffmap"]
        The data for calculating the trajectory tree, by default "pca."
    tree : Optional[nx.DiGraph]
        If `None`, it calculates the trajectory tree. If nx.DiGraph is given,
        the tree is constructed from the nx.DiGraph, by default `None`.
    groupby : str
        Key for categorical in `adata.obs`. You can pass your
        metadata of clusters, by default "leiden."
    """
    if method not in ["euclid", "paga", "gauss"]:
        raise ValueError("Argument 'method' must be 'euclid','paga' or 'gauss'.")
    if dimension not in ['pca', 'diffmap']:
        raise ValueError("Argument 'dimension' must be 'pca' or 'diffmap'.")

    if not isinstance(adata, AnnData):
        raise ValueError("trajectory_tree() expects AnnData argument")
    else:
        if groupby not in adata.obs.keys():
            raise ValueError("Did not find adata.obs[{}]".format(groupby))
        if "X_{}".format(dimension) not in list(adata.obsm):
            raise ValueError(
                "Did not find 'X_{}' in adata.obsm. Run scanpy.tl to calculate first.".format(dimension))

    # when adata.obs[groupby] is metadata and not categorized, it cause errors.
    # maybe there are better ways to set adata.obs, but one code below does jobs enough.
    if not is_categorical_dtype(adata.obs[groupby]):
        adata.obs[groupby] = adata.obs[groupby].astype('category')

    pp = Preprocessing()
    pp.trajectory_tree(
        adata,
        root_node=root_node,
        tree=tree,
        method=method,
        dimension=dimension,
        groupby=groupby
    )


def _tree_align_dpt(
    adata1: AnnData,
    adata2: AnnData,
    num_genes1: int = 2000,
    num_genes2: int = 2000,
    cost: float = 1.0,
    use_raw: bool = True
):
    """\
    Just a combination of tree_alignment() and dpt().

    Parameters
    ----------
    adata1 : AnnData
        [description]
    adata2 : AnnData
        [description]
    num_genes1 : int, optional
        [description], by default 2000.
    num_genes2 : int, optional
        [description], by default 2000.
    cost : float, optional
        [description], by default 1.0.
    use_raw : bool, optional
        [description], by default `True`.

    Returns
    -------
    capital_data: [CapitalData]
        [description]

    """
    if not isinstance(adata1, AnnData):
        raise ValueError("tree_align() expect AnnData argument")
    if not isinstance(adata2, AnnData):
        raise ValueError("tree_align() expect AnnData argument")

    if 'capital' not in adata1.uns:
        raise ValueError(
            "'capital' is not found in adata1.uns. Run cp.trajectory_tree() first.")

    if 'capital' not in adata2.uns:
        raise ValueError(
            "'capital' is not found in adata2.uns. Run cp.trajectory_tree() first.")

    if use_raw:
        if adata1.raw == None:
            adata1.raw = adata1
        if adata2.raw == None:
            adata2.raw = adata2
    else:
        adata1.raw = adata1
        adata2.raw = adata2

    print("Calculating tree alignment")
    tree_align = Tree_Alignment()
    aligned_data = tree_align.tree_alignment(
        adata1,
        adata2,
        cost=cost,
        N_1=num_genes1,
        N_2=num_genes2
    )
    print("Calculating pseudotime for each alignment")
    dpt = DPT()
    dpt.dpt_for_alignments(
        aligned_data,
        no_prune=False
    )

    print("Calculation finished.")

    return aligned_data


def tree_alignment(
    adata1: AnnData,
    adata2: AnnData,
    num_genes1: int = 2000,
    num_genes2: int = 2000,
    cost: float = 1.0,
    use_raw: bool = True
):
    """\
    Calculate an alignment of the two trajectory trees.

    Use a dynamic programming algorithm in `[P Bille05] <https://www.sciencedirect.com/science/article/pii/S0304397505000174>`_ .
    Use `capital.pl.tree_alignment` to obtain the alignment tree.

    Parameters
    ----------
    adata1 : AnnData
        The annotated data matrix.
    adata2 : AnnData
        The annotated data matrix.
    num_genes1 : int
        Number of highly variable genes in adata1 to calculate intersection of the genes.
        If not `None`, it recalculates highly variable genes using genes in adata.raw,
        if `None`, it uses the genes in adata.obs,
        by default 2000.
    num_genes2 : int
        Same as num_genes1, but for adata2, by default 2000.
    cost : float
        Gap cost for the tree alignment calculation, by default 1.0.
    use_raw : bool
        If `True`, it uses adata.raw to calculate highly variable genes.
        If `False` and adata.raw was `None`, it resets adata.raw as adata, by default `True`.

    Returns
    -------
    capital_data: [CapitalData]
        The data matrices containing the results of CAPTIAL.

    """
    if not isinstance(adata1, AnnData):
        raise ValueError("tree_align() expect AnnData argument")
    if not isinstance(adata2, AnnData):
        raise ValueError("tree_align() expect AnnData argument")

    if 'capital' not in adata1.uns:
        raise ValueError("'capital' is not found in adata1.uns. Run cp.trajectory_tree() first.")

    if 'capital' not in adata2.uns:
        raise ValueError("'capital' is not found in adata2.uns. Run cp.trajectory_tree() first.")

    if use_raw:
        if adata1.raw == None:
            adata1.raw = adata1
        if adata2.raw == None:
            adata2.raw = adata2
    else:
        adata1.raw = adata1
        adata2.raw = adata2

    print("Calculating tree alignment")
    tree_align = Tree_Alignment()
    aligned_data = tree_align.tree_alignment(
        adata1,
        adata2,
        cost=cost,
        N_1=num_genes1,
        N_2=num_genes2
    )
    print("Calculation finished.")

    return aligned_data

def dpt(
    aligned_data: CapitalData,
    alignment: Union[str, list, None] = None,
    no_prune: bool = False
):
    """\
    Calculate pseudotime for the alignments.

    Use `scanpy.tl.dpt() <https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.dpt.html#scanpy.tl.dpt>`_.

    If there are many alignments, you are recommended to specify which alignment is used to calculate pseudotime to save time.
    Use `cdata.alignmentlist` and `capital.pl.tree_alignment` to take a look at alignment data.

    Parameters
    ----------
    aligned_data : CapitalData
        The data matrices containing the results of CAPTIAL.
    alignment : Union[str, list, None]
        Specify which alignment is used to calculate pseudotime.
        If `None`, it calculates pseudotime for all alignments, by default `None`.
    no_prune : bool
        Specify global alignment or local alignment.
        If `True`, it calculates pseudotime using all the cells in the clusters in the alignment.
        If `False`, it uses all the cells except ones that have no corresponding cluster in the root or the leaf of the alignment,
        by default `False`.
    """

    if not isinstance(aligned_data, CapitalData):
        raise ValueError(
                        "dpt() expect CapitalData argument."
                        "Run tree_alignment() first and pass the return as the argument as below.\n"
                        "cdata = cp.tree_alignment(adata1, adata2)\n"
                        "cp.dpt(cdata)\n"
                        )
    dpt = DPT()
    dpt.dpt_for_alignments(
        aligned_data,
        alignment=alignment,
        no_prune=no_prune
    )


def dtw(
    aligned_data: CapitalData,
    gene: Union[str, list, np.ndarray],
    alignment: Union[str, list, None] = None
):
    """\
    Calculate dynamic time warping for genes.

    This requires that `dpt()` has been run for the specified alignments.

    Parameters
    ----------
    aligned_data : CapitalData
        The data matrices containing the results of CAPTIAL.
    gene : Union[str, list, np.ndarray]
        Genes for calculating dynamic time warping.
    alignment : Union[str, list, None]
        Specify which alignment is used to calculate dynamic time warping.
        If `None`, it calculates dynamic time warping for all the alignments, by default `None`.
    """

    if not isinstance(aligned_data, CapitalData):
        raise ValueError(
                        "dtw() expect CapitalData argument." \
                        "Run tree_align() first and pass the return as the argument as below.\n"\
                        "cdata = cp.tree_alignment(adata1, adata2)\n"\
                        "cp.dpt(cdata)\n"\
                        "cp.dtw(cdata, genelist)"
                        )
    dtw = DynamicTimeWarping()
    dtw.dtw_for_alignments(
        aligned_data,
        gene,
        alignment=alignment
    )


def genes_similarity_score(
    aligned_data: CapitalData,
    gene: Union[str, list, np.ndarray, None] = None,
    alignment: Union[str, list, None] = None,
    min_disp: float = 2.0,
):
    """\
    Calculate similarity scores using dynamic time warping `[H.Sakoe78] <https://ieeexplore.ieee.org/document/1163055>`_.

    Use `tslearn.metrics.dtw <https://tslearn.readthedocs.io/en/stable/gen_modules/metrics/tslearn.metrics.dtw.html#r34771c3a90b5-1>`_ \
    for each gene and calculate the similarity of gene expression dynamics.

    Use `cdata.similarity_score` to access the results. Higher on the list, more similar the gene expression dynamics are.

    Parameters
    ----------
    aligned_data : CapitalData
        The data matrices containing the results of CAPTIAL.
    gene : Union[str, list, np.ndarray]
        Genes for calculating similarity, by default `None`.
    alignment : Union[str, list, None]
        Specify alignments to calculate similarity. If `None`, it computes for all the alignments,
        by default `None`.
    min_disp : float
        If `gene` is `None`, the union of genes that have bigger dispersion than min_disp
        in the clusters of the alignment is used, by default 2.0.
    """
    if not isinstance(aligned_data, CapitalData):
        raise ValueError(
            "genes_similarity_score() expect CapitalData argument.")

    dtw = DynamicTimeWarping()
    dtw.get_genes_similarity_score(
        aligned_data,
        gene=gene,
        alignment=alignment,
        min_disp=min_disp,
    )


def read_capital_data(
    dirname: str,
    adata1_name: Optional[str] = None,
    adata2_name: Optional[str] = None
):
    """\
    Reading CAPITAL data.

    Three files (capital_data.npz, adata1.h5ad, adata2.hda5) must be in the same directory.

    Parameters
    ----------
    dirname : str
        Path to the CAPITAL data directory.
    adata1_name : Optional[str]
        Specify adata1 file name if the file name is changed since created, by default `None`.
    adata2_name : Optional[str]
        Specify adata2 file name if the file name is changed since created, by default `None`.
    """
    capital_data = np.load(os.path.join(dirname, "capital_data.npz"), allow_pickle=True)

    filedata = capital_data['filenamedata']
    adata1_tmp_name = filedata[0][0]
    adata2_tmp_name = filedata[1][0]

    adata1_filename = adata1_name if adata1_name is not None else adata1_tmp_name
    adata2_filename = adata2_name if adata2_name is not None else adata2_tmp_name

    adata1 = sc.read(os.path.join(dirname, adata1_filename))
    adata2 = sc.read(os.path.join(dirname, adata2_filename))

    alignedtree_adjacencymatrix = capital_data["alignedtree"]
    alignedtree_node = capital_data["alignedtree_node"]
    df_aligned_tree = pd.DataFrame(alignedtree_adjacencymatrix, index=alignedtree_node, columns=alignedtree_node)
    aligned_tree = nx.convert_matrix.from_pandas_adjacency(df_aligned_tree, create_using=nx.DiGraph)
    genes_for_tree_align = capital_data["genes_for_tree_align"]
    alignmentcost = capital_data["alignmentcost"]
    alignmentdict = capital_data["alignmentdict"][()]
    similarity_score = capital_data["similarity_score"][()]

    alignd_data = CapitalData(
        adata1,
        adata2,
        aligned_tree,
        alignmentcost,
        genes_for_tree_align,
        alignmentdict
        )

    alignd_data.similarity_score = similarity_score

    return alignd_data

