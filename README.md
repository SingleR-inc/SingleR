# SingleR - Single-cell Recognition

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/SingleR.svg)](https://bioconductor.org/packages/release/bioc/html/SingleR.html "How long has SingleR been in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/release/SingleR.svg)](https://bioconductor.org/packages/stats/bioc/SingleR/ "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months")
[![Support posts](https://bioconductor.org/shields/posts/SingleR.svg)](https://support.bioconductor.org/t/SingleR/ "Support site activity for SingleR, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")

**Current build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/SingleR.svg)](https://bioconductor.org/packages/release/bioc/html/SingleR.html#archives "Whether SingleR release is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/SingleR.svg)](https://bioconductor.org/packages/release/bioc/html/SingleR.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/SingleR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/SingleR "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/SingleR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/SingleR/ "Bioconductor release build")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/SingleR.svg)](https://bioconductor.org/packages/devel/bioc/html/SingleR.html#archives "Whether SingleR devel is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/SingleR.svg)](https://bioconductor.org/packages/devel/bioc/html/SingleR.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/devel/bioc/SingleR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/SingleR "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/SingleR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/SingleR/ "Bioconductor devel build")

Recent advances in single cell RNA-seq (scRNA-seq) have enabled an unprecedented level of granularity in characterizing gene expression changes in disease models. 
Multiple single cell analysis methodologies have been developed to detect gene expression changes and to cluster cells by similarity of gene expression. 
However, the classification of clusters by cell type relies heavily on known marker genes, and the annotation of clusters is performed manually. 
This strategy suffers from subjectivity and limits adequate differentiation of closely related cell subsets. 
Here, we present SingleR, a novel computational method for unbiased cell type recognition of scRNA-seq. 
SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently. 

For more informations please refer to the manuscript: [Aran, Looney, Liu et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology (2019)](https://www.nature.com/articles/s41590-018-0276-y)

This repository contains a simplified, more performant version of `SingleR`. 
The original repository containing the legacy version can be found [here](https://github.com/dviraran/SingleR). 
**This version does not support the browser application that accompanied the original version.**

## Installation

This is the __development__ version of the R/Bioconductor package SingleR. It may contain unstable or untested new features. If you are looking for the __release__ version of this package please go to its official [Bioconductor landing page](https://bioconductor.org/packages/SingleR) and follow the instructions there to install it.

If you were really looking for this development version, then you can install it via:

```r
install.packages("BiocManager")
BiocManager::install("SingleR", version = "devel")
```

Alternatively, you can install it from GitHub using the [devtools](https://github.com/hadley/devtools "devtools") package.

```r
install.packages("devtools")
library(devtools)
install_github("LTLA/SingleR")
```

## Usage

The `SingleR()` function annotates each cell in a test dataset given a reference dataset with known labels. Documentation and basic examples can be accessed with `?SingleR`.

Both basic and advanced examples can be found in the [SingleR book](https://ltla.github.io/SingleRBook/).

### Usage with Seurat/SingleCellExperiment objects

`SingleR()` is made to be workflow/package agnostic - if you can get a matrix of normalized counts, you can use it.
`SingleCellExperiment` objects can be used directly. 
`Seurat` objects can be converted to `SingleCellExperiment` objects via Seurat's `as.SingleCellExperiment()` function or their normalized counts can be retrieved via `GetAssayData` or `FetchData`.

`SingleR` results labels can be easily added back to the metadata of these objects as well:

```R
seurat.obj[["SingleR.labels"]] <- singler.results$labels

# Or if `method="cluster"` was used:
seurat.obj[["SingleR.cluster.labels"]] <- 
        singler.results$labels[match(seurat.obj[[]][["my.input.clusters"]], rownames(singler.results))]
```

## Scalability

`SingleR` performs well on large numbers of cells - annotating 100k cells with fine-grain labels typically takes under an hour using a single processing core. 
Using broad labels can reduce the time to under 15 minutes, though run times will vary between datasets and the reference dataset used.

## Contributors

SingleR was originally developed by Dvir Aran. 
This refactor was initiated by Aaron Lun, with additional contributions from Daniel Bunis, Friederike Dündar, and Jared Andrews.

[Issues](https://github.com/LTLA/SingleR/issues) and [pull requests](https://github.com/LTLA/SingleR/pulls) are welcome.
