import numpy as np
import scanpy as sc
from scipy.spatial.distance import pdist, squareform
from tslearn.metrics import dtw_path, dtw
from sklearn.preprocessing import MinMaxScaler
from .._util import CapitalData


class DPT:
    def __init__(self):
        pass

    def _dpt_for_an_alignment(
        self,
        adata,
        cluster_list,
        alignment_id,
        copy=False
    ):

        adata = adata.copy() if copy else adata
        groupby = adata.uns["capital"]["tree"]["annotation"]
        root_cluster = cluster_list[0]
        # get  cells only in the clusters of cluster_list
        adata_dpt = adata[adata.obs[groupby].isin(
            cluster_list)].copy()

        if not isinstance(adata.X, np.ndarray):
            X = adata_dpt.X.toarray()
        else:
            X = adata_dpt.X

        Y = pdist(X, 'euclidean')
        distancearray = squareform(Y)

        distancearray = np.sum(distancearray, axis=1)
        loc = np.argsort(distancearray)
        count = 0
        while adata_dpt.obs[groupby][np.where(loc == count)[0][0]] != root_cluster:
            count += 1

        # "iroot" is a cell that are the source of
        root_cell = np.flatnonzero(adata_dpt.obs[groupby])[
            np.where(loc == count)[0][0]]
        adata_dpt.uns['iroot'] = root_cell

        # process diffusion maps and dpt, calculate "dpt_pseudotime"
        sc.tl.diffmap(adata_dpt)
        sc.tl.dpt(adata_dpt)

        adata.obs["{}_dpt_pseudotime".format(
            alignment_id)] = adata_dpt.obs["dpt_pseudotime"]
        adata.uns["capital"]["pseudotime"]["{}".format(alignment_id)] = {}

        each_pseudotime_dict = adata.uns["capital"]["pseudotime"]["{}".format(
            alignment_id)]
        each_pseudotime_dict["clusters"] = np.array(cluster_list, dtype=object)
        each_pseudotime_dict["iroot"] = root_cell
        # add dpt_pseudotime for the clusters in one alignment to adata.obs["alignment000_dpt_pseudotime"]
        # add clusters name used in the alignment  to adata.uns["capital"]["pseudotime"]["alignment000"]["clusters"]
        # add cell id of root cell of dpt to adata.uns["capital"]["pseudotime"]["alignment000"]["iroot"]
        return adata if copy else None

    def dpt_for_alignments(
        self,
        aligned_data: CapitalData,
        alignment=None,
        no_prune=False
    ):
        groupby1 = aligned_data.adata1.uns["capital"]["tree"]["annotation"]
        groupby2 = aligned_data.adata2.uns["capital"]["tree"]["annotation"]
        aligned_data.adata1.uns["capital"]["pseudotime"] = {}
        aligned_data.adata2.uns["capital"]["pseudotime"] = {}

        alignment_id_list = []
        if alignment is None:
            alignment_id_list = list(aligned_data.alignmentdict.keys())
        else:
            if isinstance(alignment, list):
                alignment_id_list = alignment
            elif isinstance(alignment, str):
                alignment_id_list = [alignment]
            else:
                raise ValueError(
                    "alignment must be list or str of alignment. "
                    "e.g. 'alignment000' or ['alignment000','alignment001', ...].")

        for alignment_id in alignment_id_list:
            route1 = aligned_data.alignmentdict[alignment_id]["data1"]
            route2 = aligned_data.alignmentdict[alignment_id]["data2"]
            if all([i == "#" for i in route1]):
                continue
            if all([i == "#" for i in route2]):
                continue

            if not no_prune:
                tmp = []
                for i in range(len(route1)):
                    if route1[i] == "#" or route2[i] == "#":
                        tmp.append(i)
                    if route1[i] != "#" and route2[i] != "#":
                        break
                for i in reversed(range(len(route1))):
                    if route1[i] == "#" or route2[i] == "#":
                        tmp.append(i)
                    if route1[i] != "#" and route2[i] != "#":
                        break
                cluster_list1 = [route1[i]
                                for i in list(range(len(route1))) if i not in tmp]
                cluster_list2 = [route2[i]
                                for i in list(range(len(route2))) if i not in tmp]

            cluster_list1 = [
                node for node in route1
                if node in aligned_data.adata1.obs[groupby1].values
            ]
            cluster_list2 = [
                node for node in route2
                if node in aligned_data.adata2.obs[groupby2].values
            ]

            if len(cluster_list1) == 0 or len(cluster_list2) == 0:
                break

            self._dpt_for_an_alignment(
                aligned_data.adata1, cluster_list1, alignment_id)
            self._dpt_for_an_alignment(
                aligned_data.adata2, cluster_list2, alignment_id)


class DynamicTimeWarping():
    def __init__(self):
        pass

    def dtw_for_alignments(
        self,
        aligned_data: CapitalData,
        gene,
        alignment=None,
    ):
        groupby1 = aligned_data.adata1.uns["capital"]["tree"]["annotation"]
        groupby2 = aligned_data.adata2.uns["capital"]["tree"]["annotation"]

        alignment_id_list = []
        if alignment is None:
            alignment_id_list = list(aligned_data.alignmentdict.keys())
        else:
            if isinstance(alignment, list):
                alignment_id_list = alignment
            elif isinstance(alignment, str):
                alignment_id_list = [alignment]
            else:
                raise ValueError(
                    "alignment must be list or str of alignment. "\
                    "e.g. 'alignment000' or ['alignment000','alignment001', ...].")

        if isinstance(gene, list):
            genenamelist = gene
        elif isinstance(gene, str):
            genenamelist = [gene]
        elif isinstance(gene, np.ndarray):
            genenamelist = list(gene)
        else:
            raise ValueError("gene must be list, str or np.ndarray.")

        for alignment_id in alignment_id_list:
            cluster_list1 = aligned_data.alignmentdict[alignment_id]["data1"]
            cluster_list2 = aligned_data.alignmentdict[alignment_id]["data2"]

            adata_dpt1 = aligned_data.adata1[aligned_data.adata1.obs[groupby1].isin(
                cluster_list1)].copy()
            adata_dpt1 = adata_dpt1[adata_dpt1.obs.sort_values(
                "{}_dpt_pseudotime".format(alignment_id)).index].copy()
            adata_dpt2 = aligned_data.adata2[aligned_data.adata2.obs[groupby2].isin(
                cluster_list2)].copy()
            adata_dpt2 = adata_dpt2[adata_dpt2.obs.sort_values(
                "{}_dpt_pseudotime".format(alignment_id)).index].copy()

            for genename in genenamelist:
                ordered_cells1, ordered_cells2, path, dist = self._applying_dtw_to_clusters(
                    adata_dpt1, adata_dpt2, genename)
                result = aligned_data.alignmentdict[alignment_id]
                result[genename] = {"ordered_cells1": ordered_cells1,
                                    "ordered_cells2": ordered_cells2,
                                    "path": path,
                                    }

    # data used to do dynamic time warping are stored in
    # file1_ordered_data, file2_ordered_data, paths
    # if nodes tha are compared are empty, data are stored as "#"
    # if the opponent side of the node is empty, the data is sorted and stored,
    # but the another node and path are stored as "#"

    # When using MAGIC there aren't many outlier or dropout,
    # so outlier aren't taken any
    # expecting Anndata with one gene

    def _applying_dtw_to_clusters(
        self,
        adata1,
        adata2,
        genename,
        min_percentile_outlier=0,
        max_percentile_outlier=100
    ):
        expression1 = adata1.raw.to_adata()[:, genename].X
        expression2 = adata2.raw.to_adata()[:, genename].X

        if not isinstance(expression1, np.ndarray):
            expression1 = expression1.toarray()
        if not isinstance(expression2, np.ndarray):
            expression2 = expression2.toarray()

        # when all or most of the gene expressions are 0, np.percentile causes error
        # catching error by try and except, but it needs to be fixed
        # excluding cells that have too low or too high gene expression
        try:
            min_outlier, max_outlier = np.percentile(
                expression1, q=[min_percentile_outlier, max_percentile_outlier]
            )
            gene_expression1 = expression1[(expression1 <= max_outlier) & (
                expression1 >= min_outlier)]
            ordered_cells1 = np.array(adata1[(
                expression1 <= max_outlier) & (expression1 >= min_outlier)].obs_names.to_list())
        except IndexError as e:
            gene_expression1 = expression1
            ordered_cells1 = np.array(adata1.obs_names.to_list())

        # excluding cells that have too low or too high gene expression
        try:
            min_outlier, max_outlier = np.percentile(
                expression2, q=[min_percentile_outlier, max_percentile_outlier]
            )

            gene_expression2 = expression2[(
                expression2 <= max_outlier) & (expression2 >= min_outlier)]

            ordered_cells2 = np.array(adata2[(
                expression2 <= max_outlier) & (expression2 >= min_outlier)].obs_names.to_list())
        except IndexError as e:
            gene_expression2 = expression2
            ordered_cells2 = np.array(adata2.obs_names.to_list())

        path, dist = dtw_path(
            gene_expression1, gene_expression2)

        return ordered_cells1, ordered_cells2, path, dist

    def get_genes_similarity_score(
        self,
        aligned_data: CapitalData,
        gene=None,
        alignment=None,
        min_disp=1.0,
    ):

        groupby1 = aligned_data.adata1.uns["capital"]["tree"]["annotation"]
        groupby2 = aligned_data.adata2.uns["capital"]["tree"]["annotation"]

        alignment_id_list = []
        if alignment is None:
            alignment_id_list = list(aligned_data.alignmentdict.keys())
        else:
            if isinstance(alignment, list):
                alignment_id_list = alignment
            elif isinstance(alignment, str):
                alignment_id_list = [alignment]
            else:
                raise ValueError(
                    "alignment must be list or str of alignment. "\
                    "e.g. 'alignment000' or ['alignment000','alignment001', ...].")

        genenamelist = []
        if gene is not None:
            if isinstance(gene, list):
                genenamelist = gene
            elif isinstance(gene, str):
                genenamelist = [gene]
            elif isinstance(gene, np.ndarray):
                genenamelist = list(gene)
            else:
                raise ValueError("gene must be list, str or np.ndarray.")


        if aligned_data.similarity_score is None:
            dic_similarity_score = {}
        else:
            dic_similarity_score = aligned_data.similarity_score

        for alignment_id in alignment_id_list:
            cluster_list1 = aligned_data.alignmentdict[alignment_id]["data1"]
            cluster_list2 = aligned_data.alignmentdict[alignment_id]["data2"]

            adata_dpt1 = aligned_data.adata1[aligned_data.adata1.obs[groupby1].isin(
                cluster_list1)].copy()
            adata_dpt1 = adata_dpt1[adata_dpt1.obs.sort_values(
                "{}_dpt_pseudotime".format(alignment_id)).index].raw.to_adata()
            adata_dpt2 = aligned_data.adata2[aligned_data.adata2.obs[groupby2].isin(
                cluster_list2)].copy()
            adata_dpt2 = adata_dpt2[adata_dpt2.obs.sort_values(
                "{}_dpt_pseudotime".format(alignment_id)).index].raw.to_adata()

            if gene is None:
                sc.pp.highly_variable_genes(adata_dpt1)
                sc.pp.highly_variable_genes(adata_dpt2)
                s1 = set(adata_dpt1.var.index)
                s2 = set(adata_dpt2.var.index)
                genenamelist = list(s1.intersection(s2))
                disp1 = adata_dpt1[:, genenamelist].var["dispersions_norm"]
                disp2 = adata_dpt2[:, genenamelist].var["dispersions_norm"]

                genenamelist = disp1[(
                    disp1 > min_disp) | (disp2 > min_disp)].index.values

            print("Calculating similarity score of {} genes in {}".format(
                len(genenamelist), alignment_id))

            score_list = self._get_dtw_score(
                adata_dpt1, adata_dpt2, genenamelist)
            ar = np.array([genenamelist, score_list], dtype=object)
            ar = ar[:, np.argsort(ar[1])][0]
            dic_similarity_score[alignment_id] = ar

        aligned_data.similarity_score = dic_similarity_score
        print("Calculating finished")

    def _get_dtw_score(
        self,
        adata1,
        adata2,
        genenamelist,
    ):
        dist_list = []
        expression1 = adata1[:, genenamelist].X
        expression2 = adata2[:, genenamelist].X

        if not isinstance(expression1, np.ndarray):
            expression1 = expression1.toarray()
        if not isinstance(expression2, np.ndarray):
            expression2 = expression2.toarray()

        mmscaler = MinMaxScaler()
        expression1 = mmscaler.fit_transform(expression1)
        expression2 = mmscaler.fit_transform(expression2)

        for i in range(len(genenamelist)):
            dist = dtw(
                expression1[:, i], expression2[:, i])
            dist_list.append(dist)

        return dist_list
