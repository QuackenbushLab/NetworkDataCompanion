# NetworkDataCompanion
**This package is part of the [Reproducible processing of TCGA regulatory networks](https://www.biorxiv.org/content/biorxiv/early/2024/11/07/2024.11.05.622163.full.pdf) paper available on BioRxiv**

- `NetworkDataCompanion` is an R library of utilities for performing analyses on TCGA and GTEx data using the [Network Zoo](https://netzoo.github.io)
- `NetworkDataCompanion` is also the engine behind the [tcga-data-nf workflow](https://workflowhub.eu/workflows/1306)
- The code for `tcga-data-nf` is available in the QuackenbushLab GitHub: https://github.com/QuackenbushLab/tcga-data-nf

<img align="center" width="90%" src="overview.png">

# Installing

You can install this package with `devtools` using the following code:

```{R}
library(devtools)
devtools::install_github("QuackenbushLab/NetworkDataCompanion")
library(NetworkDataCompanion)
```

If you want to work with the source, you can also clone the repository and install it from the root directory locally once you have done so. 

```R
library(devtools)
devtools::install()  # provided you are in the project folder
library(NetworkDataCompanion) 
```

We have noticed that sometimes the following packages may require separate installation: "GenomicDataCommons", "edgeR", "recount", "recount3". If you have problems with these, please raise an issue here: https://github.com/QuackenbushLab/NetworkDataCompanion/issues or directly and/contact the maintainer (Kate Shutta, kshutta at hsph.harvard.edu). 

# Usage

- A quickstart guide is located in the base directory (```quickstart.Rmd``` and ```quickstart.html```). 
- You can view the rendered html of the quickstart here: [NetworkDataCompanion Quickstart](https://htmlpreview.github.io/?https://github.com/QuackenbushLab/NetworkDataCompanion/blob/main/quickstart.html)
- Example usage of all functions is also available in the unit tests: https://github.com/QuackenbushLab/NetworkDataCompanion/tree/main/tests/testthat

# Structure of the repo
- ```R``` contains the source code of our functions.
- ```tests``` contains extensive tests of the implemented functions.
- ```inst``` contains external data needed to run the analyses on TCGA and GTEx (e.g., gene mapping files)
- ```man``` is used to generate documentation. 

# Citation
```bibtex
@article{fanfani2024reproducible,
  title={Reproducible processing of TCGA regulatory networks},
  author={Fanfani, Viola and Shutta, Katherine H and Mandros, Panagiotis and Fischer, Jonas and Saha, Enakshi and Micheletti, Soel and Chen, Chen and Ben Guebila, Marouen and Lopes-Ramos, Camila Miranda and Quackenbush, John},
  journal={bioRxiv},
  pages={2024--11},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

If you find this package useful, feel free to star this repository!
