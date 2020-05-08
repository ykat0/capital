# Changelog

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