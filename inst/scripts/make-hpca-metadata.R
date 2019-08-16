df <- data.frame(
    Title = sprintf("Human Primary Cell Atlas %s", c("logcounts", "colData")), # This can be the exact file name (if self-describing) or a more complete description.
    Description = sprintf("%s human primary cell atlas data with 38 main cell types (169 subtypes)", 
        c("Matrix of log-normalized expression values from the", "Per-sample metadata containing cell type labels of")),
    RDataPath = file.path("SingleR", "hpca","1.0.0", c("logcounts.rds", "coldata.rds")),
    BiocVersion="3.10", # The first Bioconductor version the resource was made available for.
    Genome=NA, # Can be NA.
    SourceType="RDA", #  Format of original data, e.g., FASTA, BAM, BigWig, etc. ‘getValidSourceTypes()’ for currently acceptable values
    SourceUrl="https://github.com/dviraran/SingleR/tree/master/data",
    SourceVersion="hpca.rda",
    Species="Homo sapiens", # getSpeciesList, validSpecies, or suggestSpecies(); can be NA
    TaxonomyId="9606",
    Coordinate_1_based=NA, #TRUE, FALSE, NA
    DataProvider="Human Primary Cell Atlas",
    Maintainer="Friederike Duendar <frd2007@med.cornell.edu>",
    RDataClass="character", # R / Bioconductor class the data are stored in, e.g., GRanges, SummarizedExperiment, ExpressionSet etc.
    DispatchClass="Rds", # Determines how data are loaded into R.
    stringsAsFactors = FALSE
)

write.csv(file="../extdata/metadata-hpca.csv", df, row.names=FALSE)

