import os
import scanpy as sc


def setty19(
    fpath: str = None
):
    """
    A preprocessed dataset of `Setty2019 <https://www.nature.com/articles/s41587-019-0068-4>`_  \
    used in our work.
    This downloads the dataset from `here <http://www.med.osaka-u.ac.jp/pub/rna/ykato/project/capital/>`__.

    Parameters
    ----------
    fpath : str
        If a path is specified as `str`, the dataset is saved in the path.
        If `None`, the dataset is saved in `../data/capital_dataset`, by default `None`.

    Returns
    -------
    An :class:`~anndata.AnnData` object
    """

    URL = "http://www.med.osaka-u.ac.jp/pub/rna/ykato/project/capital/data/setty2019_magic_n2000_k40_r4.h5ad"

    if fpath:
        download_filepath = fpath
    else:
        download_filepath = "../data/capital_dataset/setty19_capital.h5ad"

    if os.path.exists(download_filepath):
        print(f"The dataset already exist in {download_filepath}.")
        print(f"Reading datasets from {download_filepath}.")
        adata = sc.read(download_filepath)

    else:
        os.makedirs(os.path.dirname(download_filepath), exist_ok=True)
        print("Downloading the dataset.")
        adata = sc.read(download_filepath, backup_url=URL)
        print("Download completed. The dataset is saved in {}".format(download_filepath))

    return adata


def paul15(
    fpath: str = None
):
    """
    A preprocessed dataset of `Paul15 <https://www.cell.com/cell/fulltext/S0092-8674(15)01493-2>`_  \
    used in our work.
    This downloads the dataset from `here <https://github.com/ykat0/capital/tree/master/data>`__.


    Parameters
    ----------
    fpath : str
        If a path is specified as `str`, the dataset is saved in the path.
        If `None`, the dataset is saved in `../data/capital_dataset`, by default `None`.

    Returns
    -------
    An :class:`~anndata.AnnData` object
    """

    URL = "https://github.com/ykat0/capital/raw/master/data/paul2015_magic_n500_k10_r7.h5ad"

    if fpath:
        download_filepath = fpath
    else:
        download_filepath = "../data/capital_dataset/paul15_capital.h5ad"

    if os.path.exists(download_filepath):
        print(f"The dataset already exist in {download_filepath}.")
        print(f"Reading datasets from {download_filepath}.")
        adata = sc.read(download_filepath)

    else:
        os.makedirs(os.path.dirname(download_filepath), exist_ok=True)
        print("Downloading the dataset.")
        adata = sc.read(download_filepath, backup_url=URL)
        print("Download completed. The dataset is saved in {}".format(download_filepath))

    return adata


def velten17(
    fpath: str = None
):
    """
    A preprocessed dataset of `Velten2017 <https://www.nature.com/articles/ncb3493>`_  \
    used in our work.
    This downloads the dataset from `here <https://github.com/ykat0/capital/tree/master/data>`__.

    Parameters
    ----------
    fpath : str
        If a path is specified as `str`, the dataset is saved in the path.
        If `None`, the dataset is saved in `../data/capital_dataset`, by default `None`.

    Returns
    -------
    An :class:`~anndata.AnnData` object
    """

    URL = "https://github.com/ykat0/capital/raw/master/data/velten2017_por_magic_n250_k6_r0.h5ad"

    if fpath:
        download_filepath = fpath
    else:
        download_filepath = "../data/capital_dataset/velten17_capital.h5ad"

    if os.path.exists(download_filepath):
        print(f"The dataset already exist in {download_filepath}.")
        print(f"Reading datasets from {download_filepath}.")
        adata = sc.read(download_filepath)

    else:
        os.makedirs(os.path.dirname(download_filepath), exist_ok=True)
        print("Downloading the dataset.")
        adata = sc.read(download_filepath, backup_url=URL)
        print("Download completed. The dataset is saved in {}".format(download_filepath))

    return adata
