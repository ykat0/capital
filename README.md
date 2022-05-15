[![PyPI version](https://badge.fury.io/py/capital.svg)](https://badge.fury.io/py/capital)
[![Documentation Status](https://readthedocs.org/projects/capital/badge/?version=latest)](https://capital.readthedocs.io/en/latest/?badge=latest)

# CAPITAL

## Alignment of single-cell trajectory trees

Last updated: 2022-05-15

We present CAPITAL, a computational method for comparing pseudotime trajectories with tree alignment whereby trajectories including branches can be automatically compared.

## Installation
* CAPITAL (ver. 1.0.7) in Python

### Requirements
* Python>=3.8 ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended)
* anndata>=0.7.4
* graphtools
* graphviz
* h5py>=2.10
* leidenalg
* magic-impute
* matplotlib>=3.1.2
* networkx>=2.3
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
$ conda create -n capital python=3.8 graphviz
```
Next, activate the environment and pull CAPITAL from PyPI:
```
$ conda activate capital
$ pip install capital
```

## Usage
Read the [documentation](https://capital.readthedocs.io/en/latest/). CAPITAL uses a single-cell analysis toolkit [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) in its implementation so that one can also use Scanpy's useful functions including preprocessing, plotting and datasets in the CAPITAL environment.