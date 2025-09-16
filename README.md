# CPlink
Mapping clinical phenotype to single-cell and spatial omics profiles 

## Overview
We introduce CPlink, an interpretable, unified, and flexible computational framework that identifies clinical phenotype-associated cells or spots from single-cell or spatial omics data leveraging reference bulk omics data for enhanced disease prediction and biological discovery.

<p align="center">
<img  src="vignettes/Figure%201.jpg" width="700" height=auto > 
</p>

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
- Seurat v4.4.0
- Signac v1.13.0
- CelliD v1.10.1
- DESeq2 v1.42.1
- survival v3.7.0
- proxyC v0.4.1
- Matrix v1.6.5
- parallel v4.3.3
- doParallel v1.0.17
- methods v4.2.0
- foreach v1.5.2
- edgeR v4.0.16
- limma 3.58.1

# Contact
If you have any questions, please feel free to contact Jiao Hua (jhua@stu.hit.edu.cn).
