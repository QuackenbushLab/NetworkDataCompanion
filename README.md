# NetworkDataCompanion

An R library of utilities for performing analyses on TCGA and GTEx data using the Network Zoo (https://netzoo.github.io). This is the engine behind [this Nextflow workflow](https://github.com/QuackenbushLab/tcga-data-nf/tree/main).

# Installing

The package requires the "GenomicDataCommons", "edgeR", "recount", "recount3" bioconductor packages to be installed prior to the main installation.

To install and use the library you can use the following: 

```R
install.packages("devtools")
devtools::install()  # provided you are in the project folder
library(NetworkDataCompanion) # load the library in your code
```

# Usage

```R
## load libraries
library(NetworkDataCompanion)
library(recount)
library(recount3)
library(GenomicDataCommons)

## Obtain expression data from recount
TCGA_lung <- recount3::create_rse_manual(
  project = "LUAD",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

## Generate NetworkDataCompanion instance
obj <- CreateNetworkDataCompanionObject(project_name = "LUAD")

## use package functionality, for example
exp_TCGA_lung <- obj$logTPMNormalization(TCGA_lung)
```

A lot more details can be found in [this vignette](./vignettes/introduction.Rmd).

# Structure of the repo
- ```R``` contains the source code of our functions.
- ```vignettes``` contains a tutorial on how to use the package.
- ```tests``` contains extensive tests of the implemented functions.
- ```insts``` contains the data needed to run the analyses on TCGA and GTEx.
- ```man``` is used to generate documentation. 

# Building

We rely on roxygen2 for package building.

```R
install.packages("tinytex")
tinytex::install_tinytex()

library(roxygen2)
roxygen2::roxygenize()

install.packages("devtools")
library(devtools)
devtools::build()
devtools::document()
build_manual(path = ".")
```
