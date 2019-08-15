# SingleR - Single-cell Recognition

Recent advances in single cell RNA-seq (scRNA-seq) have enabled an unprecedented level of granularity in characterizing gene expression changes in disease models. 
Multiple single cell analysis methodologies have been developed to detect gene expression changes and to cluster cells by similarity of gene expression. 
However, the classification of clusters by cell type relies heavily on known marker genes, and the annotation of clusters is performed manually. 
This strategy suffers from subjectivity and limits adequate differentiation of closely related cell subsets. 
Here, we present SingleR, a novel computational method for unbiased cell type recognition of scRNA-seq. 
SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently. 

For more informations please refer to the manuscript: [Aran, Looney, Liu et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology (2019)](https://www.nature.com/articles/s41590-018-0276-y)

## Install

```R
devtools::install_github('LTLA/SingleR')
```

## Updates

**08.14.2019**
This repository contains a simplified, more performant version of `SingleR`. 
It is currently in the process of being added to Bioconductor. 
The original repository can be found [here](https://github.com/dviraran/SingleR). 
This version does not support the browser application that accompanied the original version.

## Usage

The `SingleR()` function annotates each cell in a test dataset given a reference dataset with known labels.
The package directly provides a number of reference datasets generated from bulk RNA-seq of pure cell types.
There are two data sets from human cells (Human Primary Cell Atlas and Blueprint/ENCODE) and two data sets from mouse cells (e.g. Immunological Genome Project).
More details can be found in the datasets vignette.

Each reference dataset is obtained with a specific function: `HumanPrimaryCellAtlasData()`, `BlueprintEncodeData()`, `ImmGenData()`, `MouseBulkData()`.
Here, we show an example with data from HPCA:


```R
library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()
```

The newly generated `SummarizedExperiment` object can then be used for the annotation of your scRNA-seq dataset.
To illustrate this, we will use a hESC dataset from the `scRNAseq` pacakge, subsetted for the sake of speed.

```R
library(scRNAseq)
library(scater)

hESCs <- LaMannoBrainData('human-es')
hESCs <- hESCs[,1:100] # for demo-purposes only!

# Grab the common genes between the sets.
common <- intersect(rownames(hESCs), rownames(hpca.se))
hpca.se <- hpca.se[common,]
hESCs <- hESCs[common,]

# Test and reference sets should always be log normalized. The included reference sets are already normalized. 
hESCs <- logNormCounts(hESCs)
```

The reference data sets all come with two sets of cell labels, `label.main` and `label.fine`. 
The "main" labels tend to be less fine-grained than the "fine" labels, e.g., instead of specifying different subsets of T cells, the main labels will lump them all together under the label "T cells". 
For our example, we will use the more specific fine-grained cell type labels.

```R
pred.hpca <- SingleR(test = hESCs, ref = hpca.se, 
    labels = hpca.se$label.fine, assay.type.ref = "normcounts")
table(pred.hpca$labels)
```

We can now visualize our results via `plotScoreHeatmap()`, which visualizes the scores for all cells across all reference labels.
This allows users to inspect the confidence of the predicted labels across the dataset.

```R
plotScoreHeatmap(pred.hpca)
```

This is the most basic use of `SingleR()`, more advanced examples can be found in the vignette.

### Usage with Seurat/SingleCellExperiment objects

`SingleR()` is made to be workflow/package agnostic - if you can get a matrix of normalized counts, you can use it.
`SingleCellExperiment` objects can be used directly. 
`Seurat` objects can be converted to `SingleCellExperiment` objects via Seurat's `as.SingleCellExperiment()` function or their normalized counts can be retrieved via `GetAssayData` or `FetchData`.

`SingleR` results labels can be easily added back to the metadata of these objects as well:

```R
seurat.obj[["SingleR.labels"]] <- singler.results$labels

# Or if `method="cluster"` was used:
seurat.obj[["SingleR.cluster.labels"]] <- 
        singler.results$labels[match(seurat.obj@meta.data[["my.input.clusters"]], singler.results$clusts)]
```

## Scalability

`SingleR` performs well on large numbers of cells - annotating 100k cells with fine-grain labels typically takes under an hour using a single processing core. 
Using broad labels can reduce the time to under 15 minutes, though run times will vary between datasets and the reference dataset used.

## Contributors

SingleR was originally developed by Dvir Aran. 
This refactor was initiated by Aaron Lun, with additional contributions from Daniel Bunis, Friederike Düendar, and Jared Andrews.

[Issues](https://github.com/LTLA/SingleR/issues) and [pull requests](https://github.com/LTLA/SingleR/pulls) are welcome.