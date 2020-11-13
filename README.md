# CAPITAL

### Comparative pseudotime analysis of single-cell RNA-seq data

Last updated: 2020-11-13

We present CAPITAL, a method for comparing pseudotime trajectories with tree alignment whereby trajectories including branching can be compared without any knowledge of paths to be compared.

## Installation
* CAPITAL (ver. 0.2.3) in Python

More user-friendly version (1.0.0) that can be used in an interactive development environment such as [JupyterLab](https://jupyter.org/) will be available soon.

### Requirements
* Python>=3.8 ([Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended)
* graphtools
* h5py<=2.10
* leidenalg
* magic-impute
* pygraphviz
* scanpy>=1.6
* scprep
* tslearn

### Option
* jupyterlab

### Install on Linux, Windows (WSL) and macOS
0. Create a new environment for CAPITAL if you want to keep your own Python environment built with conda:
```
$ conda create -n capital python=3.8
```
Then, activate the environment:
```
$ conda activate capital
```

1. Install [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) and others via conda:
```
$ conda install -c bioconda scanpy
$ conda install graphtools leidenalg pygraphviz scprep tslearn
```
Note: the above command may automatically install h5py>=3.0, which will cause the problem of having bytes labels in H5Ad data. To circumvent this, try:
```
$ conda install h5py=2
```

2. Install [MAGIC](https://magic.readthedocs.io/en/stable/) via pip:
```
$ pip install --no-deps magic-impute
```

3. Download the tarball, and type the followings in your terminal:
```
$ tar zxf capital-0.2.3.tar.gz
$ cd capital-0.2.3
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
  data <STR>            path to the raw [H5AD|CSV|TXT] file or 10x-MTX directory of scRNA-seq gene expression profile

optional arguments:
  -h, --help            show this help message and exit
  -t, --transpose       transpose the data [off]
  --min-genes <INT>     keep cells with at least <INT> genes expressed [200]
  --min-cells <INT>     keep genes that are expressed in at least <INT> cells [3]
  --magic               impute gene expression data with MAGIC [off]
  -n <INT>, --top-n-genes <INT>
                        number of highly variable genes to keep [2000]
  -p <INT>, --n-pcs <INT>
                        number of principal components for computing a k-nearest neighbor graph and a tree [50]
  -k <INT>, --neighbors <INT>
                        compute an <INT>-nearest neighbor graph [10]
  --no-save             results are not saved [off: saved in ./processed_data]
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
  data1 <STR>           path to the preprocessed expression H5AD data for experiment 1 generated with pre_capital.py
  data2 <STR>           path to the preprocessed expression H5AD data for experiment 2 generated with pre_capital.py
  root1 <STR>           root cluster in data1
  root2 <STR>           root cluster in data2
  genes <STR>           path to the file that contains gene names to be analyzed (one gene per line)

optional arguments:
  -h, --help            show this help message and exit
  -M {euclid,gauss,paga}, --method {euclid,gauss,paga}
                        method for calculating a tree [euclid]
  -d {pca,diffmap}, --dimension {pca,diffmap}
                        dimension metric for calculating a tree [pca]
  -c <FLOAT>, --gapcost <FLOAT>
                        gap cost for calculating a tree alignment [1.0]
  -m <INT>, --n-genes1 <INT>
                        number of highly variable genes in data1 [2000]
  -n <INT>, --n-genes2 <INT>
                        number of highly variable genes in data2 [2000]
  -t, --tune            tuning mode, which affects naming of the result directory and never saves H5AD data [off]
  --no-prune            disable pruning of space nodes on edges of each alignment path for dynamic time warping [off]
```

### Step 3: run draw_capital.py to draw figures on dynamic time warping and/or expression dynamics for a gene in aligned_data created in Step 2
```
$ ./draw_capital.py [option]* <alignment> <genes>
```

#### Usage: draw_capital.py
```
positional arguments:
  alignment <STR>     path to the directory for aligned data generated with capital.py (e.g. ./aligned_data/data1_data2/alignment000/)
  genes <STR>         path to the file that contains gene names to be analyzed (one gene per line)

optional arguments:
  -h, --help          show this help message and exit
  --dtw <STR>         path to the directory for a figure on dynamic time warping. Unless specified, the figure will be saved in each gene file in    
                      the alignment directory
  --dyn <STR>         path to the directory for a figure on gene expression dynamics. Unless specified, the figure will be saved in each gene file   
                      in the alignment directory
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
