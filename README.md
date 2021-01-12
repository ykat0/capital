# CAPITAL

### Comparative pseudotime analysis of single-cell RNA-seq data

Last updated: 2021-01-12

We present CAPITAL, a computational method for comparing pseudotime trajectories with tree alignment whereby trajectories including branchings can be automatically compared.

## Installation
* CAPITAL (ver. 1.0.0) in Python

### Requirements
* Python>=3.8 ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended)
* anndata>=0.7.4
* graphtools
* graphviz
* h5py<=2.10
* leidenalg
* magic-impute
* matplotlib>=3.1.2
* networkx>=2.3
* pandas>=0.21
* pydot
* scanpy>=1.6
* scikit-learn
* scipy>=1.4
* scprep
* tslearn

### Install on Linux, Windows (WSL) and macOS
Create a new environment for CAPITAL (recommended):
```
$ conda create -n capital python=3.8 graphviz
```
Next, activate the environment and pull CAPITAL from PyPI:
```
$ conda activate capital
$ pip install capital
```

## Usage
Read the documentation (available very soon). CAPITAL uses a single-cell analysis toolkit [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) in its implementation so that one can also use Scanpy's useful functions including preprocessing, plotting and datasets in the CAPITAL environment.

## Reference
Reiichi Sugihara, Yuki Kato, Tomoya Mori and Yukio Kawahara,
**Alignment of time-course single-cell RNA-seq data with CAPITAL**,
Preprint *bioRxiv* at [https://doi.org/10.1101/859751](https://doi.org/10.1101/859751), 2019.

Reiichi Sugihara, Yuki Kato, Tomoya Mori and Yukio Kawahara,
**Comparative pseudotime analysis of single-cell RNA-seq data with CAPITAL**,
*submitted*.

---
If you have any questions, please contact [Yuki Kato](http://www.med.osaka-u.ac.jp/pub/rna/ykato/en/).
*Graduate School of Medicine, Osaka University, Japan*