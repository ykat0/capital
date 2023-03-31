# Changelog

## Ver. 1.1.2 (2023-03-22)
* Fix minor bugs

## Ver. 1.1.1 (2023-03-22)
* Add a new tutorial `Use user defined pseudotime`
* Fix minor bugs

## Ver. 1.1.0 (2023-03-17)
* Add a `pseudotime` parameter to `cp.tl.dtw()`, `cp.tl.genes_similarity_score()`, and `cp.pl` that enable users to use user-defined pseudotime in dynamic time warping
* `cp.tl.trajectory_tree()` now accepts a user-preprocessed dataset
* Fix minor bugs

## Ver. 1.0.14 (2022-12-16)
* Specify the matplotlib version as 3.5.2 due to errors in plotting figures

## Ver. 1.0.13 (2022-07-23)
* Fix minor bugs

## Ver. 1.0.12 (2022-07-23)
* Add `multi_gene` parameter to `capital.pl.gene_expression_trend()` so that users can plot a gene expression trend using results of multi genes calculated in `capital.tl.dtw()`
* Fix minor bugs

## Ver. 1.0.11 (2022-07-14)
* Change function `cp.tl.dpt()` to calculate the root cell so that it is derived from all cells in the datasets and used for all linear alignments
* Fix minor bugs

## Ver. 1.0.10 (2022-07-05)
* Fix minor bugs

## Ver. 1.0.9 (2022-07-05)
* Updated `capital.tl.dtw()`
* Fix minor bugs

## Ver. 1.0.8 (2022-06-28)
* Specify the Networkx version as 2.8.3 due to errors in plotting figures
* Fix minor bugs

## Ver. 1.0.7 (2022-05-15)
* Fix minor bugs

## Ver. 1.0.6 (2022-05-15)
* Fix compatibility with scanpy 1.9.1

## Ver. 1.0.5 (2022-05-14)
* Fix minor bugs
* Delete synthetic_dataset1 and synthetic_dataset2 from `capital.dataset` function

## Ver. 1.0.4 (2021-12-03)
* Add some new arguments in plotting
* Fix minor bugs

## Ver. 1.0.3 (2021-10-08)
* Adjust CAPITAL codes to the recent Scanpy update to 1.8.1
* Fix minor bugs

## Ver. 1.0.2 (2021-03-16)
* Fine-tune codes to the Code Ocean compute capsule, changing the locations of directories for saving results

## Ver. 1.0.1 (2021-03-03)
* Specify the Scanpy version as 1.6 for pip installation to guarantee the reproducibility of results

## Ver. 1.0.0 (2021-01-13)
* Provide a brand new system so that one can use CAPITAL interactively, which is useful for the interactive development environment JupyterLab
* Upload CAPITAL to PyPI to install it via 'pip' command
* Upload a tutorial and API using a theme provided by Read the Docs
* Upload datasets used in our work, which are available via CAPITAL's dataset functions
* Improve running time
* Fix minor bugs
* Move early versions of CAPITAL to "legacy" directory

## Ver. 0.2.3 (2020-11-13)
* Add '--magic' option to impute gene expression data with [MAGIC](https://magic.readthedocs.io/en/stable/) in pre_capital.py
* Fix a minor bug in draw_capital.py
* Improve calculation processes in all codes

## Ver. 0.2.2 (2020-08-05)
* Fix a bug for 10x Genomics data in capital.py
* Replace '--local-align' option with '--no-prune' in capital.py
* Add a new function for automatically saving figures on trajectory trees and their alignment in capital.py
* Add a new positional argument \<genes\> in draw_capital.py, so that a user can input multiple genes all together
* Add a new function for automatically saving figures on dynamic time warping and expression dynamics in the alignment directory generated with capital.py if a user does not specify where to save in draw_capital.py
* Remove '--save' option in draw_capital.py, meaning that a figure file will be always saved
* Remove '--showfig' option in draw_capital.py
* Add x-label (and y-label) on plots generated with draw_capital.py
* Adjust the layout of a figure legend to the figure area if --data1-name (--data2-name) is specified

## Ver. 0.2.1 (2020-06-21)
* Fix a bug on the process without cluster in capital.py
* Fix a bug on the key error of leiden_colors in draw_capital.py

## Ver. 0.2.0 (2020-06-01)
* Adjust CAPITAL codes to the recent Scanpy update to 1.5
* Add a function for 10x-MTX directory input in pre_capital.py

## Ver. 0.1.14 (2020-05-15)
* Improve the naming of the result directory with respect to gapcost in the tuning mode in capital.py

## Ver. 0.1.13 (2020-05-08)
* Change the data type of gapcost to float in capital.py
* Keep the value of a gapcost in the result directory when the tuning mode is on in capital.py

## Ver. 0.1.12 (2020-05-01)
* Fix a bug in the statement on concatenation of AnnData in capital.py
* Keep the name of a method for computing a tree in the result directory when the tuning mode is on in capital.py

## Ver. 0.1.11 (2020-04-24)
* Add '-p' option to specify the number of principal components to construct a nearest neighbor graph in pre_capital.py
* Add '-M' option to specify the way of computing a minimum spanning tree (Euclid distance, adaptive Gaussian kernel or PAGA connectivity) in capital.py
* Improve some processes for speed-up

## Ver. 0.1.10 (2020-04-10)
* Fix a bug on dynamic programing calculation in capital.py
* Fix a bug on displaying a figure for tree alignment in capital.py

## Ver. 0.1.9 (2020-02-15)
* Compress output H5AD data generated with pre_capital.py and capital.py to save disk space
* Disable output of H5AD data when '-t' option is used in capital.py
* Change the hierarchy of generated aligned_data directory to remove redundancy
* Fix a bug that occurs when no directory for saving figures is found in draw_capital.py
* Change some option symbols in draw_capital.py for usability

## Ver. 0.1.8 (2020-02-13)
* Show the number of predicted clusters in standard out in pre_capital.py
* Change the hierarchy of generated directories in capital.py
* Improve the naming of result directories in the tuning mode in capital.py

## Ver. 0.1.7 (2020-02-11)
* Change part of option symbols in pre_capital.py and capital.py to avoid inconsistency between scripts as much as possible
* Add a tuning mode in capital.py, which affects only naming of result directories

## Ver. 0.1.6 (2020-02-10)

* Improve a method for saving a PDF figure in pre_capital.py
* Add a function for showing a progress bar when calculating dynamic time warping for each gene in capital.py
* Change part of directories/files naming in capital.py

## Ver. 0.1.5 (2020-02-07)

* Fix a bug on print and change generated directories/files naming in capital.py

## Ver. 0.1.4 (2020-02-03)

* Fix a bug for the input of one aligned cluster pair in draw_capital.py

## Ver. 0.1.3 (2020-02-02)

* Add further explanations to help messages in pre_capital.py and draw_capital.py

## Ver. 0.1.2 (2020-01-31)

* Fix a bug on an attribute error associated with networkx in capital.py

## Ver. 0.1.1 (2020-01-31)

* Add a function for absolute paths input in capital.py
* Allow multiple genes input in capital.py
* Change the definition of a cluster centroid to a vector of the median of expression of each gene in the cluster
* Add a function for showing normalized alignment distance
* Fix a bug in traceback of tree alignment in capital.py

## Ver. 0.1.0 (2020-01-21)

* Add command line option parse
* Split codes into pre_capital.py, capital.py and draw_capital.py

## Ver. 0.0.2 (2019-12-05)

* Fix a minor bug in capital.py

## Ver. 0.0.1 (2019-11-29)

* Upload capital.py for test use