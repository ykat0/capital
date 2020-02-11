# Changelog

## Ver. 0.1.7 (2020-02-11)
* Change part of option symbols in pre_capital.py and capital.py to avoid inconsistency between scripts as much as possible
* Add a tuning mode in capital.py, which affects only naming of the result directory

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

* Upload capital.py