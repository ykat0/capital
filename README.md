# CAPITAL

### Alignment of time-course single-cell RNA-seq data

Last updated: 2019-12-05

We present CAPITAL, a method for comparing pseudotime trajectories with tree alignment whereby trajectories including branching can be compared without any knowledge of paths to be compared.

## Requirements
* Python 3

### Required modules
* matplotlib
* networkx
* numpy
* pandas
* scanpy
* scikit-learn
* scipy
* tslearn

## Recommended usage
Install [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) with Anaconda (Miniconda):
```
$ conda install scanpy
```
This will install almost all of the required modules other than "tslearn," so you will have to install it manually by:
```
$ conda instal tslearn
```
Download the tarball, and type the followings in your terminal:
```
$ tar zxf capital-0.0.2.tar.gz
$ cd capital-0.0.2
$ ./capital.py
```
**Note:**
The current version 0.0.2 of CAPITAL is only for test use.
More stable and user-friendly version 1.0.0 will be available soon.

## Reference
Reiichi Sugihara, Yuki Kato, Tomoya Mori and Yukio Kawahara,
**Alignment of time-course single-cell RNA-seq data with CAPITAL**,
Preprint *bioRxiv* at [https://doi.org/10.1101/859751](https://doi.org/10.1101/859751), 2019.

---
If you have any questions, please contact [Yuki Kato](http://www.med.osaka-u.ac.jp/pub/rna/ykato/)  
*Graduate School of Medicine, Osaka University, Japan*
