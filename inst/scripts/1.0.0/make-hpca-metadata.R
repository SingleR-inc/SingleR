df <- data.frame(
    Title = sprintf("Human Primary Cell Atlas %s", c("logcounts", "colData")), 
    Description = sprintf("%s human primary cell atlas data with 38 main cell types (169 subtypes)", 
        c("Matrix of log-normalized expression values from the", "Per-sample metadata containing cell type labels of")),
    RDataPath = file.path("SingleR", "hpca","1.0.0", c("logcounts.rds", "coldata.rds")),
    BiocVersion="3.10", 
    Genome=NA,
    SourceType="RDA", 
    SourceUrl="https://github.com/dviraran/SingleR/tree/master/data",
    SourceVersion="hpca.rda",
    Species="Homo sapiens", # getSpeciesList, validSpecies, or suggestSpecies(); can be NA
    TaxonomyId="9606",
    Coordinate_1_based=NA, 
    DataProvider="Dvir Aran",
    Maintainer="Friederike Duendar <frd2007@med.cornell.edu>",
    RDataClass=c("matrix", "DFrame"),
    DispatchClass="Rds", 
    stringsAsFactors = FALSE
)

write.csv(file="../extdata/1.0.0/metadata-hpca.csv", df, row.names=FALSE)

