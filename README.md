# TransmissionCluster
  TransmissionCluster is a GUI Python program that computes the optimal genetic distance threshold and identifies clusters in a phylogenetic tree.

  Given a Newick (.nwk) formatted tree, genetic distance threshold *d*, and an optional support threshold *s*, TransmissionCluster finds the minimum number of clusters of samples in the tree such that:

  1. The maximum pairwise distance between samples in the cluster is at most *d*.\*
  2. Samples cannot be connected by branches with support less than or equal to *s*.
  3. The samples in the cluster must define a clade in the tree.

  For a tree with *n* samples, the algorithm is O(*n*).

  >\*TransmissionCluster can compute the distance threshold that maximizes the number of non-singleton clusters over all thresholds from 0 to *d* in steps of 0.0001 (i.e., 0, 0.001, 0.002, ..., *d*).


## Usage

  1. Clone or download this repository.
  2. Ensure that all dependencies are installed (see below).
  3. Run using one of the following terminal commands:
  - `./TransmissionCluster.py`
  - `python3 TransmissionCluster.py`

  **Required Parameters:**
  - A Newick formatted phylogenetic tree (.nwk)
  - Output Filename
  - Genetic Distance Threshold
    - Note: When using the distance-free method this threshold is used as the upper bound of possible thresholds to compute.

  **Optional Parameters:**
  - 'Compute Best Distance Threshold' (when unchecked TransmissionCluster uses the specified threshold)
  - Posterior/Bootstrap Support Threshold
  - Pairwise Distance File (TransmissionCluster will generate an edge list to graph clusters when a distance file is provided.)
    - *See documentation for details and formatting requirements.*


## Dependencies
  * [Python 3](https://www.python.org/downloads/)
  * [TreeSwift](https://github.com/niemasd/TreeSwift)
    - Install via command line using: `pip install treeswift`
  * [PySimpleGUI](https://pypi.org/project/PySimpleGUI/)
    - Install via command line using: `pip install pysimplegui`


### Acknowledgements
TransmissionCluster builds upon methodology from [TreeCluster](https://github.com/niemasd/TreeCluster).
