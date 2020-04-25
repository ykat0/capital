# CAPITAL

### Alignment of time-course single-cell RNA-seq data

Last updated: 2020-04-24

We present CAPITAL, a method for comparing pseudotime trajectories with tree alignment whereby trajectories including branching can be compared without any knowledge of paths to be compared.

## Installation
* CAPITAL (ver. 0.1.11) in Python

### Requirements
* Python>=3.6 ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended)
* leidenalg
* scanpy
* tslearn

### Install on Linux, Windows (WSL) and macOS
0. Create a new environment for CAPITAL if you want to keep your own Python environment built with conda:
```
$ conda create -n capital python=3.8
```
Then, activate the environment:
```
$ conda activate capital
```

1. Install [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) and others with conda:
```
$ conda install -c bioconda scanpy
$ conda install leidenalg tslearn
```

2. Download the tarball, and type the followings in your terminal:
```
$ tar zxf capital-0.1.11.tar.gz
$ cd capital-0.1.11
```

## Pipeline
CAPITAL consists of three python codes. The standard use of CAPITAL is illustrated with the following steps:

### Step 1: run pre_capital.py to preprocess raw data of scRNA-seq gene expression for each experiment
```
$ ./pre_capital.py [option]* <data1>
$ ./pre_capital.py [option]* <data2>
```

#### Usage: pre_capital.py
```
positional arguments:
  data <STR>            path to the raw [TXT|CSV|H5AD] data of scRNA-seq gene
                        expression profile

optional arguments:
  -h, --help            show this help message and exit
  -t, --transpose       transpose the data [off]
  --min-genes <INT>     minimum number of genes expressed to keep [200]
  --min-cells <INT>     minimum number of cells expressed to keep [3]
  -n <INT>, --top-n-genes <INT>
                        number of highly variable genes to keep [1000]
  -k <INT>, --neighbors <INT>
                        size k of local neighborhood used to compute a
                        k-nearest neighbor graph [10]
  -p <INT>, --n_pcs <INT>
                        number of PCs used to compute a k-nearest neighbor
                        graph and a tree [50]
  --no-save             results are not saved [on: saved in ./processed_data]
  -f <STR>, --filename <STR>
                        save data as <filename>.h5ad and umap_<filename>.pdf
  --save-fig            save a UMAP PDF figure in ./figures [off]
```

### Step 2: run capital.py to predict pseudotime trajectories for two preprocessed data, and compute a trajectory alignment along with associated dynamic time warping of each aligned path
```
$ ./capital.py [option]* <data1> <data2> <root1> <root2> <genes>
```

#### Usage: capital.py
```
positional arguments:
  data1 <STR>           path to the preprocessed expression H5AD data for
                        experiment 1 generated with pre_capital.py
  data2 <STR>           path to the preprocessed expression H5AD data for
                        experiment 2 generated with pre_capital.py
  root1 <STR>           root cluster in data1
  root2 <STR>           root cluster in data2
  genes <STR>           path to the file that contains gene names to be
                        analyzed (one gene per line)

optional arguments:
  -h, --help            show this help message and exit
  -c <INT>, --gapcost <INT>
                        gap cost used to calculate tree alignment [3]
  -m <INT>, --n-genes1 <INT>
                        number of highly variable genes in data1 [2000]
  -n <INT>, --n-genes2 <INT>
                        number of highly variable genes in data2 [2000]
  -M {euclid,gauss,paga}, --method {euclid,gauss,paga}
                        method used to calculate tree (PAGA adjacency matrix
                        is used by default) [paga]
  -l, --local-align     calculate dynamic time warping on local ailgnment
                        [off]
  -t, --tune            tuning mode, which affects naming of the result
                        directory and never saves H5AD data [off]
```

### Step 3: run draw_capital.py to draw figures on dynamic time warping and/or expression dynamics for a gene in aligned_data created in Step 2
```
$ ./draw_capital.py [option]* <alignment>
```

#### Usage: draw_capital.py
```
positional arguments:
  alignment <STR>     path to the directory for aligned data generated with
                      capital.py (e.g.
                      ./aligned_data/data1_data2/gene/alignment001)

optional arguments:
  -h, --help          show this help message and exit
  --dtw <STR>         path to the file for a (e.g. PDF) figure on dynamic time
                      warping
  --dyn <STR>         path to the file for a (e.g. PDF) figure on gene
                      expression dynamics
  --data1-name <STR>  data1 name on the plot
  --data2-name <STR>  data2 name on the plot
```

## Reference
Reiichi Sugihara, Yuki Kato, Tomoya Mori and Yukio Kawahara,
**Alignment of time-course single-cell RNA-seq data with CAPITAL**,
Preprint *bioRxiv* at [https://doi.org/10.1101/859751](https://doi.org/10.1101/859751), 2019.

---
If you have any questions, please contact [Yuki Kato](http://www.med.osaka-u.ac.jp/pub/rna/ykato/en/)  
*Graduate School of Medicine, Osaka University, Japan*
