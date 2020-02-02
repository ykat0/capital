# CAPITAL

### Alignment of time-course single-cell RNA-seq data

Last updated: 2020-02-02

We present CAPITAL, a method for comparing pseudotime trajectories with tree alignment whereby trajectories including branching can be compared without any knowledge of paths to be compared.

## Installation
* CAPITAL (ver. 0.1.3) (**capital-0.1.3.tar.gz**) in Python

### Requirements
* Python >= 3.6 ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended)
  
* leidenalg
* matplotlib
* networkx>=2.4
* numpy
* pandas
* scanpy
* scikit-learn
* scipy
* tslearn

### Install on Linux, Windows (WSL) and macOS
0. Create a new environment for CAPITAL if you want to keep your own Python environment built with conda:
```
$ conda create -n capital python=3.7
```
Then, activate the environment:
```
$ source activate capital
```

1. Install [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) with conda:
```
$ conda install -c bioconda scanpy
```
This command will install almost all of the required modules other than "leidenalg" and "tslearn," so you will have to install them manually by:
```
$ conda install leidenalg tslearn
```

2. Download the tarball, and type the followings in your terminal:
```
$ tar zxf capital-0.1.3.tar.gz
$ cd capital-0.1.3
```

## Reference
Reiichi Sugihara, Yuki Kato, Tomoya Mori and Yukio Kawahara,
**Alignment of time-course single-cell RNA-seq data with CAPITAL**,
Preprint *bioRxiv* at [https://doi.org/10.1101/859751](https://doi.org/10.1101/859751), 2019.

---
If you have any questions, please contact [Yuki Kato](http://www.med.osaka-u.ac.jp/pub/rna/ykato/en/)  
*Graduate School of Medicine, Osaka University, Japan*
