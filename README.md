[![PyPI version](https://badge.fury.io/py/capital.svg)](https://badge.fury.io/py/capital)
[![Documentation Status](https://readthedocs.org/projects/capital/badge/?version=latest)](https://capital.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/221850492.svg)](https://zenodo.org/badge/latestdoi/221850492)

# CAPITAL

## Alignment of single-cell trajectory trees

Last updated: 2022-12-16

We present CAPITAL, a computational method for comparing pseudotime trajectories with tree alignment whereby trajectories including branches can be automatically compared.

## Installation
* CAPITAL (ver. 1.0.14) in Python

### Requirements
* Python>=3.8 ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended)
* anndata>=0.7.4
* graphtools
* graphviz
* h5py>=2.10
* leidenalg
* magic-impute
* matplotlib==3.5.2
* networkx<=2.8.3
* pandas>=0.21
* pydot
* scanpy==1.9.1
* scikit-learn
* scipy>=1.4
* scprep
* tslearn

### Install on Linux, Windows (WSL) and macOS
Create a new environment for CAPITAL (recommended):
```
conda create -n capital python=3.9 graphviz
```
Next, activate the environment and pull CAPITAL from PyPI:
```
conda activate capital
```
```
pip install capital
```

## Usage
Read the [documentation](https://capital.readthedocs.io/en/latest/). CAPITAL uses a single-cell analysis toolkit [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) in its implementation so that one can also use Scanpy's useful functions including preprocessing, plotting and datasets in the CAPITAL environment.

## Code Ocean
We also provide a Code Ocean [compute capsule](https://codeocean.com/capsule/5673663/tree/v1) to reproduce our results.

## Reference
Reiichi Sugihara, Yuki Kato, Tomoya Mori and Yukio Kawahara,
**Alignment of single-cell trajectory trees with CAPITAL**,
*Nature Communications*, vol. 13, 5972, 2022. [[Link]](https://www.nature.com/articles/s41467-022-33681-3)

---
If you have any questions, please contact [Yuki Kato](https://www.med.osaka-u.ac.jp/pub/rna/ykato/en/index.html).
*Graduate School of Medicine, Osaka University, Japan*.