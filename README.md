daee 0.1.3
====

daee - Data Analysis for Ecology and Evolution

Set of functions for analysis of ecological and evolutionary data. The functions are available as an R package only to facilitate installation.

List of functions:
- Calculate Phylogenetic Eigenvector Regression (PVR) with eigenvector selection;
- Extract residuals from Mantel test;
- Generate Linear Models (LM) with all possible combinations of variables included in the full model;
- Show information about a label in a phylogenetic tree;
- Makes node labels;
- Add species in a phylogenetic tree;
- Organize a list in a single matrix.


## Installation
  
To install the latest version of this package, use [`devtools`](https://github.com/hadley/devtools):

```r
require(devtools)
install_github("vanderleidebastiani/daee”)
```

Require last version of [`SYNCSA`](https://github.com/vanderleidebastiani/SYNCSA) and 
[`PCPS`](https://github.com/vanderleidebastiani/PCPS) packages.
