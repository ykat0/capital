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
## Pipeline
CAPITAL consists of three python codes. The standard use of CAPITAL is illustrated with the following steps:

### Step 1: run pre_capital.py to preprocess raw data of scRNA-seq gene expression for each experiment
```
$ ./pre_capital.py [option]* /path/to/data1
$ ./pre_capital.py [option]* /path/to/data2
```

#### Usage (pre_capital.py)
```
usage: ./pre_capital.py [option]* <data>

positional arguments:
  data <STR>            path to the raw data of scRNA-seq gene expression
                        profile

optional arguments:
  -h, --help            show this help message and exit
  -t, --transpose       transpose the data [off]
  --min-genes <INT>     minimum number of genes expressed to keep [200]
  --min-cells <INT>     minimum number of cells expressed to keep [3]
  -g <INT>, --top-n-genes <INT>
                        number of highly variable genes to keep [1000]
  -k <INT>, --neighbors <INT>
                        size k of local neighborhood used to compute a
                        k-nearest neighbor graph [10]
  --no-save             results are not saved [on: saved in ./processed_data]
  -n <STR>, --name <STR>
                        save data as <name>.h5ad
  --save-fig            save a UMAP PDF figure in ./figures [off]
```

## Reference
Reiichi Sugihara, Yuki Kato, Tomoya Mori and Yukio Kawahara,
**Alignment of time-course single-cell RNA-seq data with CAPITAL**,
Preprint *bioRxiv* at [https://doi.org/10.1101/859751](https://doi.org/10.1101/859751), 2019.

---
If you have any questions, please contact [Yuki Kato](http://www.med.osaka-u.ac.jp/pub/rna/ykato/en/)  
*Graduate School of Medicine, Osaka University, Japan*
