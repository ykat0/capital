import os
import networkx as nx
import numpy as np
from anndata import AnnData
import dataclasses
from typing import Final

@dataclasses.dataclass
class CapitalData:
    adata1: Final[AnnData] = None
    adata2: Final[AnnData] = None
    alignedtree: Final[nx.Graph] = None
    alignmentcost: Final[float] = None
    genes_for_tree_align: Final[list] = None
    alignmentdict: dict = dataclasses.field(
        default_factory=dict
    )
    alignmentlist: list = dataclasses.field(
        default_factory=list,
        init=False
    )
    similarity_score: dict = dataclasses.field(
        default_factory=dict,
        init=False
    )

    def __post_init__(self):
        dic = self.alignmentdict
        self.alignmentlist = [
            (i, dic[i]["data1"], dic[i]["data2"])
            for i in dic.keys()
            ]

    def write(
        self,
        dirname=None,
        adata1_name=None,
        adata2_name=None
    ):
        """\
        Write CapitalData.

        Save data in a file. In the file, capital_data.npz, adata1.h5ad and adata2.h5ad are saved.

        Parameters
        ----------
        dirname : str,
            Path to the directory, if `None`, saves in `./capital_data/`, by default None
        adata1_name : str
            File name of adata1, by default None
        adata2_name : str
            File name of adata2, by default None
        """
        if dirname is None:
            dirname = "./capital_data/"
        os.makedirs(dirname, exist_ok=True)

        adata1_filename = adata1_name if adata1_name is not None else "adata1.h5ad"
        adata2_filename = adata2_name if adata2_name is not None else "adata2.h5ad"

        filenamedata = np.array(
                    [[str(adata1_filename)],
                        [str(adata2_filename)],
                        ["capital_data.npz"]],
                    dtype=object
                )
        self.adata1.write(os.path.join(dirname, adata1_filename))
        self.adata2.write(os.path.join(dirname, adata2_filename))

        np.savez_compressed(os.path.join(dirname, "capital_data.npz"),
                            filenamedata=filenamedata,
                            genes_for_tree_align=self.genes_for_tree_align,
                            alignmentcost=self.alignmentcost,
                            alignmentdict=self.alignmentdict,
                            similarity_score=self.similarity_score,
                            alignedtree=nx.convert_matrix.to_numpy_array(self.alignedtree),
                            alignedtree_node=np.array(self.alignedtree.nodes())
                            )
