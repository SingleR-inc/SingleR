write.csv(file="../extdata/metadata-hpca.csv", 
          data.frame(
              Title = sprintf("Human Primary Cell Atlas: %s", c("normcounts", "colData")), # This can be the exact file name (if self-describing) or a more complete description.
              Description = sprintf("%s human primary cell atlas data with 38 main cell types (169 subtypes)", 
                                    c("Matrix of norm. expression values from the", "Per-sample metadata containing the cell type labels of")),
              RDataPath = file.path("SingleR", "hpca","1.0.0"), 
                                    c("normcounts.rds", "coldata.rds"),
              BiocVersion="3.10", # The first Bioconductor version the resource was made available for.
              Genome=NA, # Can be NA.
              SourceType="RDA", #  Format of original data, e.g., FASTA, BAM, BigWig, etc. ‘getValidSourceTypes()’ for currently acceptable values
              SourceUrl=c(
                  "https://github.com/dviraran/SingleR/tree/master/data",
                  "https://github.com/dviraran/SingleR/tree/master/data"
              ),
              SourceVersion=c(
                  "hpca.rda",
                  "hpca.rda"),
              Species="Homo sapiens", # getSpeciesList, validSpecies, or suggestSpecies(); can be NA
              TaxonomyId="9606",
              Coordinate_1_based=NA, #TRUE, FALSE, NA
              DataProvider="Human Primary Cell Atlas",
              Maintainer="Friederike Duendar <frd2007@med.cornell.edu>",
              RDataClass="character", # R / Bioconductor class the data are stored in, e.g., GRanges, SummarizedExperiment, ExpressionSet etc.
              DispatchClass="Rds", # Determines how data are loaded into R.
              stringsAsFactors = FALSE
          ),
          row.names=FALSE)
