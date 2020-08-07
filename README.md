# Minimum Redundancy Maximum Relevance (MRMR) variable selection, version 1.0

The code implemented in R software performs supervised variable selection of a given (possibly high-dimensional) dataset by the MRMR method, i.e. Minimum Redundancy-Maximum Relevance. The computations, which were tested over three real datasets,
include automatic choices of all parameters and compare various measures of relevance and redundancy as well as various classifiers.

Feel free to use or modify the code.

## Requirements

You need to install these package of R software: MASS, glmnet, e1071, pamr, rda, rrlda.

Available at https://cran.r-project.org/web/packages/robustbase/index.html

## Usage

* DimReduction.R: main file with relevance and redundancy measures, variable selection (including the search for its optimal parameter, penalizing redundancy with respect to relevance)
* Classifiers.R: auxiliary file comparing various classifiers (suitable for high-dimensional data)

## Authors
  * Jan Kalina, The Czech Academy of Sciences, Institute of Computer Science
  * Anna Schlenker, The First Medical Faculty, Charles University & Faculty of biomedical engineering, Czech Technical University

## Contact

Do not hesitate to contact us (kalina@cs.cas.cz) or write an Issue.

## How to cite

Please consider citing the following:

Kalina J, Schlenker A (2015): A robust supervised variable selection for noisy high-dimensional data. BioMed Research International, 2015, Article 320385.

## Acknowledgement

This work was supported by the Czech Science Foundation grant 19-05704S.