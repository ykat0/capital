import os
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.lines import Line2D
from networkx.drawing.nx_pydot import graphviz_layout
from anndata import AnnData
from typing import Union, Optional
from .._util import CapitalData


def dtw(
    aligned_data: CapitalData,
    gene: Union[str, list],
    alignment: Union[str, list, None] = None,
    data1_name: Optional[str] = "data1",
    data2_name: Optional[str] = "data2",
    ncols: int = 2,
    widthspace: float = 0.10,
    heightspace: float = 0.30,
    fontsize: float = 12,
    legend_fontsize: float = 12,
    dpi: int = 600,
    show: bool = True,
    save: Union[str, bool] = False,
):
    """\
    Plot the results of dynamic time warping.

    Parameters
    ----------
    aligned_data : CapitalData
        The data matrices containing the results of CAPTIAL.
    gene : Union[str, list]
        Keys for annotations of genes.
    alignment : Union[str, list, None], optional
        Keys for alignments to be plotted. If `None`, all alignments will be plotted, by default `None`.
    data1_name : Optional[str], optional
        Text of data1's legend, by default "data1".
    data2_name : Optional[str], optional
        Text of data2's legend, by default "data2".
    ncols : int
        Number of panels per row, by default 2.
    widthspace : float
        Width of space in the panels, by default 0.10.
    heightspace : float
        Height of space in the panels, by default 0.30.
    fontsize : float
        Font size of the title, by default 12.
    legend_fontsize : float, optional
        Font size of the legend, by default 12.
    show : bool,
        If `True`, show the figure. If `False`, return the figure, by default `True`.
    save : Union[str, bool]
        If `True` or `str`, save the figure. If a path is specified as `str`, the figure is saved in the path, by default `False`.
    """

    if not isinstance(aligned_data, CapitalData):
        ValueError("draw_dtw() expects an CapitalData argument.")

    plot_list = __set_plot_list(aligned_data, alignment, gene)

    data1 = aligned_data.adata1
    data2 = aligned_data.adata2
    groupby1 = data1.uns["capital"]["tree"]["annotation"]
    groupby2 = data2.uns["capital"]["tree"]["annotation"]
    groupby_colors1 = groupby1 + "_colors"
    groupby_colors2 = groupby2 + "_colors"

    if groupby_colors1 not in data1.uns:
        from scanpy.plotting._utils import _set_default_colors_for_categorical_obs
        _set_default_colors_for_categorical_obs(data1, groupby1)
    if groupby_colors2 not in data2.uns:
        from scanpy.plotting._utils import _set_default_colors_for_categorical_obs
        _set_default_colors_for_categorical_obs(data2, groupby2)

    fig, grid = __panel(ncols, len(plot_list), widthspace=widthspace, heightspace=heightspace)

    for count, (alignment, genename) in enumerate(plot_list):
        ax = fig.add_subplot(grid[count])
        dtw_dic = aligned_data.alignmentdict[alignment]
        path = dtw_dic[genename]["path"]
        ordered_cells1 = dtw_dic[genename]["ordered_cells1"]
        ordered_cells2 = dtw_dic[genename]["ordered_cells2"]

        ax.set_title("{}_{}".format(alignment, genename), fontsize=fontsize)
        ax.set_xlabel("Pseudotime", fontsize=fontsize)
        ax.set_aspect(0.8)
        ax.tick_params(
            labelbottom=True,
            labelleft=False,
            labelright=False,
            labeltop=False)
        ax.tick_params(
            bottom=False,
            left=False,
            right=False,
            top=False)
        ax.grid(False)

        y1 = data1[ordered_cells1, :].obs["{}_dpt_pseudotime".format(alignment)]
        y2 = data2[ordered_cells2, :].obs["{}_dpt_pseudotime".format(alignment)]

        clusters1 = data1[ordered_cells1, :].obs[groupby1]
        colors1 = data1[ordered_cells1, :].uns[groupby_colors1]
        clusters2 = data2[ordered_cells2, :].obs[groupby2]
        colors2 = data2[ordered_cells2, :].uns[groupby_colors2]

        if not type(colors1) == np.ndarray:
            colors1 = np.array([colors1])

        if not type(colors2) == np.ndarray:
            colors2 = np.array([colors2])

        dic1 = dict(
            zip(data1[ordered_cells1, :].obs[groupby1].cat.categories, colors1))
        dic2 = dict(
            zip(data2[ordered_cells2, :].obs[groupby2].cat.categories, colors2))

        colorlist1 = clusters1.replace(dic1)
        colorlist2 = clusters2.replace(dic2)

        ax.scatter(
            np.array(list(y1)), np.ones(len(y1)),
            color=colorlist1,
            zorder=-2,
        )

        ax.scatter(
            np.array(list(y2)), np.zeros(len(y2)),
            color=colorlist2,
            zorder=-2,
        )

        n_col1 = np.ceil(len(dic1.keys()) / 5).astype(int)
        ordered_clusters1 = [
            cluster for cluster in
            aligned_data.alignmentdict[alignment]['data1']
            if cluster != "#"
            ]
        patches1 = []
        for cluster in ordered_clusters1:
            patches1.append(Line2D(
                                range(1), range(1),
                                marker='o', color=dic1[cluster],
                                label=cluster, linewidth=0
                                )
                            )
        legend1 = ax.legend(handles=patches1,
                            labels=ordered_clusters1,
                            bbox_to_anchor=(1.05, 0.5),
                            loc='lower left',
                            ncol=n_col1,
                            title="{}".format(data1_name),
                            title_fontsize=legend_fontsize
                            )

        n_col2 = np.ceil(len(dic2.keys()) / 5).astype(int)
        ordered_clusters2 = [
            cluster for cluster in
            aligned_data.alignmentdict[alignment]['data2']
            if cluster != "#"
        ]
        patches2 = []
        for cluster in ordered_clusters2:
            patches2.append(Line2D(range(1), range(1),
                            marker='o', color=dic2[cluster], label=cluster, linewidth=0)
                            )
        legend2 = ax.legend(handles=patches2,
                            labels=ordered_clusters2,
                            bbox_to_anchor=(1.05, 0.5),
                            loc='upper left',
                            ncol=n_col2,
                            title="{}".format(data2_name),
                            title_fontsize=legend_fontsize
                            )

        ax.add_artist(legend1)
        ax.add_artist(legend2)

        for i, j in path:
            i = int(i)
            j = int(j)
            ax.plot((y1[i], y2[j]), (1, 0),
                    color='grey', alpha=0.5, zorder=-3)

        ax.set_rasterization_zorder(-1)

    if save:
        if isinstance(save, str):
            save_dir = save
            os.makedirs(os.path.dirname(save_dir), exist_ok=True)
        elif save is True:
            os.makedirs("./figures", exist_ok=True)
            save_dir = "./figures/dtw.png"
        fig.savefig(save_dir, bbox_inches='tight', pad_inches=0.1, dpi=dpi)
    if show:
        plt.show()
        plt.close()
    else:
        plt.close()
        return fig

def gene_expression_trend(
    aligned_data: CapitalData,
    gene: Union[str, list],
    alignment: Union[str, list, None] = None,
    outliers: list = [100, 0],
    polyfit_dimension: int = 3,
    switch_psedotime: bool = False,
    data1_name: Optional[str] = "data1",
    data2_name: Optional[str] = "data2",
    ncols: int = 2,
    widthspace: float = 0.5,
    heightspace: float = 0.30,
    dpi: int = 600,
    show: bool = True,
    save: Union[str, bool] = False,
):
    """\
    Plot gene expression trend.

    Parameters
    ----------
    aligned_data : CapitalData
        The data matrices containing the results of CAPTIAL.
    gene : Union[str, list]
        Keys for annotations of genes.
    alignment : Union[str, list, None], optional
        Keys for alignments to be plotted. If `None`, all alignments will be plotted, by default `None`.
    outliers : list
        Outliers for gene expression, by default [100, 0].
    polyfit_dimension : int
        Degree of the fitting polynomial, by default 3.
    switch_psedotime : bool, optional
        If `False`, data1's pseudotime is used to plot two expression trends. If `True`, switch it to data2's pseudotime, by default `False`.
    data1_name : Optional[str]
        Text of data1's legend, by default "data1".
    data2_name : Optional[str]
        Text of data2's legend, by default "data2".
    ncols : int
        Number of panels per row, by default 2.
    widthspace : float
        Width of space in the panels, by default 0.5.
    heightspace : float
        Height of space in the panels, by default 0.30.
    show : bool
        If `True`, show the figure. If `False`, return the figure, by default `True`.
    save : Union[str, bool]
        If `True` or `str`, save the figure. If a path is specified as `str`, the figure is saved in the path, by default `None`.
    """

    if not isinstance(aligned_data, CapitalData):
        ValueError(
            "draw_gene_expression_trend() expects an CapitalData argument.")

    plot_list = __set_plot_list(aligned_data, alignment, gene)

    data1 = aligned_data.adata1
    data2 = aligned_data.adata2

    fig, grid = __panel(ncols, len(plot_list), widthspace=widthspace, heightspace=heightspace)

    max_out, min_out = outliers[0], outliers[1]
    Num_quantiles = 10
    count = 0
    for count, (alignment, genename) in enumerate(plot_list):
        ax = fig.add_subplot(grid[count])
        dtw_dic = aligned_data.alignmentdict[alignment]
        path = dtw_dic[genename]["path"]
        ordered_cells1 = dtw_dic[genename]["ordered_cells1"]
        ordered_cells2 = dtw_dic[genename]["ordered_cells2"]

        ax.set_title("{}, {}".format(alignment, genename))
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("Expression level")

        pseudotime = []
        data1_expression_level = []
        data2_expression_level = []

        expression1 = data1.raw.to_adata(
        )[ordered_cells1, genename].X.T[0]
        expression2 = data2.raw.to_adata(
        )[ordered_cells2, genename].X.T[0]

        if switch_psedotime:
            pseudotime = data2[ordered_cells2, :].obs["{}_dpt_pseudotime".format(
                alignment)][[j for _, j in path]].values
        else:
            pseudotime = data1[ordered_cells1, :].obs["{}_dpt_pseudotime".format(
                alignment)][[i for i, _ in path]].values

        data1_expression_level = expression1[[i for i, _ in path]]
        data2_expression_level = expression2[[j for _, j in path]]

        array = np.array(
            [pseudotime, data1_expression_level, data2_expression_level])

        # when all or most of the gene expressions are 0, np.percentile causes error
        # catching error by try and except, but it needs to be fixed
        qcut_index = pd.qcut(
            array[0], Num_quantiles, labels=False, duplicates='drop')

        gene_expression1_list = []
        for i in range(Num_quantiles):
            qcut_array = array[:, qcut_index == i]
            try:
                min_outlier, max_outlier = np.percentile(
                    qcut_array[1], q=[min_out, max_out]
                )
                gene_expression1 = qcut_array[:, (qcut_array[1] <= max_outlier) & (
                    qcut_array[1] >= min_outlier)][[0, 1]]
            except IndexError as e:
                gene_expression1 = qcut_array[[0, 1]]
            gene_expression1_list.append(gene_expression1)

        array1_poly1d = np.concatenate(gene_expression1_list, axis=1)

        gene_expression2_list = []
        for i in range(Num_quantiles):
            qcut_array = array[:, qcut_index == i]
            try:
                min_outlier, max_outlier = np.percentile(
                    qcut_array[2], q=[min_out, max_out]
                )
                gene_expression2 = qcut_array[:, (qcut_array[2] <= max_outlier) & (
                    qcut_array[2] >= min_outlier)][[0, 2]]
            except IndexError as e:
                gene_expression2 = qcut_array[[0, 2]]
            gene_expression2_list.append(gene_expression2)

        array2_poly1d = np.concatenate(gene_expression2_list, axis=1)

        y1 = np.poly1d(np.polyfit(array1_poly1d[0], array1_poly1d[1], polyfit_dimension))(
            array1_poly1d[0])
        y2 = np.poly1d(np.polyfit(array2_poly1d[0], array2_poly1d[1], polyfit_dimension))(
            array2_poly1d[0])

        ax.scatter(array[0], array[1],
                    color="tomato", alpha=0.5, rasterized=True)
        ax.scatter(array[0], array[2],
                    color="lightskyblue", alpha=0.5, rasterized=True)

        if data1_name is None:
            ax.plot(array1_poly1d[0], y1, color="red")
        else:
            ax.plot(array1_poly1d[0], y1,
                    color="red", label=data1_name)

        if data2_name is None:
            ax.plot(array2_poly1d[0], y2, color="blue")
        else:
            ax.plot(array2_poly1d[0], y2,
                    color="blue", label=data2_name)

        if data1_name is None and data2_name is None:
            pass
        else:
            ax.legend(
                bbox_to_anchor=(1.05, 1),
                loc='upper left',
                borderaxespad=0,
                fontsize=12
                )

    if save:
        if isinstance(save, str):
            save_dir = save
            os.makedirs(os.path.dirname(save_dir), exist_ok=True)
        elif save is True:
            os.makedirs("./figures", exist_ok=True)
            save_dir = "./figures/expression_trend.png"
        fig.savefig(save_dir,  bbox_inches='tight', pad_inches=0.1, dpi=dpi)
    if show:
        plt.show()
        plt.close()
    else:
        plt.close()
        return fig



def trajectory_tree(
    adata: AnnData,
    figsize: tuple = (10, 8),
    node_size: int = 1200,
    font_size: int = 18,
    show: bool = True,
    dpi: int = 600,
    save: Union[str, bool] = False,
):
    """\
    Plot a trajectory tree.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix.
    figsize : tuple,
        Size of figure, by default (10, 8).
    node_size : int,
        Size of nodes, by default 1200.
    font_size : int,
        Font size, by default 18.
    show : bool,
        If `True`, show the figure. If `False`, return figure, by default `True`.
    save : Union[str, bool]
        If `True` or `str`, save the figure. If a path is specified as `str`, the figure is saved in the path, by default `False`.
    """

    if not isinstance(adata, AnnData):
        ValueError("draw_tree() expects an AnnData argument.")

    groupby = adata.uns["capital"]["tree"]["annotation"]
    groupby_colors = groupby + "_colors"
    if groupby_colors not in adata.uns:
        from scanpy.plotting._utils import _set_default_colors_for_categorical_obs
        _set_default_colors_for_categorical_obs(adata, groupby)

    tree = nx.convert_matrix.from_pandas_adjacency(
        adata.uns["capital"]["tree"]["tree"], create_using=nx.DiGraph)

    clusters = list(adata.obs[groupby].cat.categories)
    dic = dict(zip(clusters, adata.uns[groupby_colors]))
    colorlist = list(pd.Series(list(tree.nodes())).replace(dic))

    plt.figure(figsize=figsize)
    pos = graphviz_layout(tree, prog='dot')
    nx.draw(
        tree,
        pos,
        node_color=colorlist,
        font_weight='bold',
        font_size=font_size,
        node_size=node_size,
        with_labels=True,
        edge_color="gray",
        arrows=True
        )

    if save:
        if isinstance(save, str):
            save_dir = save
            os.makedirs(os.path.dirname(save_dir), exist_ok=True)
        elif save is True:
            os.makedirs("./figures", exist_ok=True)
            save_dir = "./figures/tree.png"
        plt.savefig(save_dir,  bbox_inches='tight', pad_inches=0.1, dpi=dpi)
    if show:
        plt.show()
        plt.close()
    else:
        plt.close()
        return fig


def tree_alignment(
    aligned_data: CapitalData,
    figsize: tuple = (10, 8),
    node_size: int = 1200,
    font_size: int = 18,
    dpi: int = 600,
    show: bool = True,
    save: Union[str, bool] = False,
):
    """\
    Plot the alignment of the trees.

    Parameters
    ----------
    aligned_data : CapitalData
        [description]
    figsize : tuple
        [description], by default (10, 8).
    node_size : int
        [description], by default 1200.
    font_size : int
        [description], by default 18.
    show : bool,
        If `True`, show the figure. If `False`, return figure, by default `True`.
    save : Union[str, bool]
        If `True` or `str`, save the figure. If a path is specified as `str`, the figure is saved in the path, by default `None`.
    """

    if not isinstance(aligned_data, CapitalData):
        ValueError("draw_alignmenttree() expects an CapitalData argument.")

    tree = aligned_data.alignedtree
    mapping = {i: str(i).replace("'", "").replace(
        "(", "").replace(")", "") for i in tree.nodes()}
    tree = nx.relabel_nodes(tree, mapping)
    plt.figure(figsize=figsize)
    pos = graphviz_layout(tree, prog='dot')
    nx.draw(tree, pos, node_color='lightskyblue', font_weight='bold', font_size=font_size,
            node_size=node_size, with_labels=True, edge_color="gray", arrows=True)

    if save:
        if isinstance(save, str):
            save_dir = save
            os.makedirs(os.path.dirname(save_dir), exist_ok=True)
        elif save is True:
            os.makedirs("./figures", exist_ok=True)
            save_dir = "./figures/alignmenttree.png"
        plt.savefig(save_dir,  bbox_inches='tight', pad_inches=0.1,dpi = dpi)
    if show:
        plt.show()
    else:
        plt.close()
    plt.close()

def __panel(
    ncols,
    num_panels,
    widthspace=0.15,
    heightspace=0.30
):
    from matplotlib import gridspec

    n_panels_x = min(ncols, num_panels)
    n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
    fig = plt.figure(
        figsize=(
            n_panels_x * rcParams['figure.figsize'][0] * (1.20),
            n_panels_y * rcParams['figure.figsize'][1],
        ),
    )
    grid = gridspec.GridSpec(
        nrows=n_panels_y,
        ncols=n_panels_x,
        hspace=heightspace,
        wspace=widthspace,
    )
    return fig, grid


def __set_plot_list(
    aligned_data,
    alignment,
    gene
):
    alignmentlist = []
    if alignment is None:
        alignmentlist = list(aligned_data.alignmentdict.keys())
    else:
        if isinstance(alignment, list):
            alignmentlist = alignment
        elif isinstance(alignment, str):
            alignmentlist = [alignment]
        else:
            raise ValueError("alignment must be list or str of alignment \
            e.g 'alignment000' or ['alignment000','alignment001', ...].")

    alignment_check = [
        alignment for alignment in alignmentlist
        if alignment not in aligned_data.alignmentdict.keys()
    ]
    if len(alignment_check) > 0:
        raise ValueError(
            "Alignments were not found. {}".format(alignment_check)
        )

    if isinstance(gene, list):
        genenamelist = gene
    elif isinstance(gene, str):
        genenamelist = [gene]
    elif isinstance(gene, np.ndarray):
        genenamelist = list(gene)
    else:
        raise ValueError("gene must be list, str or np.ndarray.")

    plot_list_tmp = list(itertools.product(alignmentlist, genenamelist))
    plot_list = [
        (alignment, gene) for alignment, gene in plot_list_tmp
        if gene in aligned_data.alignmentdict[alignment]
    ]

    if len(plot_list) == 0:
        raise ValueError(
            "Genes were not found in alignments. Run capital.dtw() first.")

    non_plot_list = [
        gene + " not in " + alignment for alignment, gene in plot_list_tmp
        if gene not in aligned_data.alignmentdict[alignment]
    ]

    if len(non_plot_list) > 0:
        print("The genes below were not found in the alignments. Run capital.dtw() or specify the alignments to draw.")
        if len(non_plot_list) < 5:
            print(non_plot_list)
        else:
            print("{} and more...".format(non_plot_list[:5]))

    return plot_list
