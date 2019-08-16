df <- data.frame(
    Title = sprintf("Mouse bulk RNA-seq %s", c("logcounts", "colData")), # This can be the exact file name (if self-describing) or a more complete description.
    Description = sprintf("%s 358 bulk RNA-seq samples of sorted cell types collected from GEO", 
        c("Matrix of log-normalized expression values from", "Per-sample metadata containing cell type labels of")),
    RDataPath = file.path("SingleR", "mouse.rnaseq","1.0.0", c("logcounts.rds", "coldata.rds")),
    BiocVersion="3.10", # The first Bioconductor version the resource was made available for.
    Genome=NA, # Can be NA.
    SourceType="RDA", #  Format of original data, e.g., FASTA, BAM, BigWig, etc. ‘getValidSourceTypes()’ for currently acceptable values
    SourceUrl="https://github.com/dviraran/SingleR/tree/master/data",
    SourceVersion="mouse.rnaseq.rda",
    Species="Mus musculus", # getSpeciesList, validSpecies, or suggestSpecies(); can be NA
    TaxonomyId="10090",
    Coordinate_1_based=NA, #TRUE, FALSE, NA
    DataProvider="Benayoun Lab",
    Maintainer="Friederike Duendar <frd2007@med.cornell.edu>",
    RDataClass="character", # R / Bioconductor class the data are stored in, e.g., GRanges, SummarizedExperiment, ExpressionSet etc.
    DispatchClass="Rds", # Determines how data are loaded into R.
    stringsAsFactors = FALSE
)

write.csv(file="../extdata/metadata-mouse.rnaseq.csv", df, row.names=FALSE)
