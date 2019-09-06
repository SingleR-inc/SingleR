df <- data.frame(
    Title = sprintf("Mouse bulk RNA-seq %s", c("logcounts", "colData")), 
    Description = sprintf("%s 358 bulk RNA-seq samples of sorted cell types collected from GEO", 
        c("Matrix of log-normalized expression values from", "Per-sample metadata containing cell type labels of")),
    RDataPath = file.path("SingleR", "mouse.rnaseq","1.0.0", c("logcounts.rds", "coldata.rds")),
    BiocVersion="3.10",
    Genome=NA, 
    SourceType="RDA", 
    SourceUrl="https://github.com/dviraran/SingleR/tree/master/data",
    SourceVersion="mouse.rnaseq.rda",
    Species="Mus musculus", 
    TaxonomyId="10090",
    Coordinate_1_based=NA,
    DataProvider="Dvir Aran",
    Maintainer="Friederike Duendar <frd2007@med.cornell.edu>",
    RDataClass=c("matrix", "DataFrame"),
    DispatchClass="Rds", 
    stringsAsFactors = FALSE
)

write.csv(file="../extdata/metadata-mouse.rnaseq.csv", df, row.names=FALSE)
