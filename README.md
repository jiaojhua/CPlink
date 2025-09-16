# CPlink
Mapping clinical phenotype to single-cell and spatial omics profiles 

## Overview
We introduce CPlink, an interpretable, unified, and flexible computational framework that identifies clinical phenotype-associated cells or spots from single-cell or spatial omics data leveraging reference bulk omics data for enhanced disease prediction and biological discovery.
https://github.com/jiaojhua/CPlink/blob/main/vignettes/Figure%201.jpg

# Installation
To run ``CPlink`` R package, install from GitHub through ``devtools`` directly:
```R
install.packages('devtools')
library(devtools)
devtools::install_github("jiaojhua/CPlink")
```

# Tutorials

* For ST datasets, please see [here](https://github.com/jiaojhua/CPlink/blob/main/vignettes/Tutorial-ST.ipynb), datasets are available at this [link].

* For scRNA-seq datasets, please see [here](https://github.com/jiaojhua/CPlink/blob/main/vignettes/Tutorial-scRNA-seq.ipynb), datasets are available at this [link].

* For scATAC-seq datasets, please see [here](https://github.com/jiaojhua/CPlink/blob/main/vignettes/Tutorial-scATAC-seq.ipynb), datasets are available at this [link].

# Dependencies
- Seurat
- Signac
- CelliD
- DESeq2
- edgeR
- limma
- methods
- proxyC
- Matrix
- parallel
- doParallel
- survival
- foreach

# Contact
If you have any questions, please feel free to contact Jiao Hua (jhua@stu.hit.edu.cn).
