df <- data.frame(
    Title = sprintf("Monaco Immune Cell RNA-seq %s", c("logcounts", "colData")), 
    Description = sprintf("%s 114 RNA-seq samples of immune cells from GSE107011", 
        c("Matrix of log-normalized expression values from", "Per-sample metadata containing cell type labels of")),
    RDataPath = file.path("SingleR", "monaco_immune","1.0.0", c("logcounts.rds", "coldata.rds")),
    BiocVersion="3.10", 
    Genome=NA, 
    SourceType="TXT", 
    SourceUrl="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011",
    SourceVersion="GSE107011",
    Species="Homo sapiens", 
    TaxonomyId="9606",
    Coordinate_1_based=NA, 
    DataProvider="GEO",
    Maintainer="Jared Andrews <jared.andrews@wustl.edu>",
    RDataClass=c("matrix", "DataFrame"), 
    DispatchClass="Rds", 
    stringsAsFactors = FALSE
)

write.csv(file="../extdata/metadata-monaco_immune.csv", df, row.names=FALSE)
