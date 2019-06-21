# This tests the trainSingleR function.
# library(testthat); library(SingleR); source("setup.R"); source("test-train.R")

test_that("trainSingleR works correctly for genes='de'", {
    out <- trainSingleR(training, training$label)
    expect_identical(out$search$mode, "de")

    # Checking that the original expression is correctly returned,
    # and that the NN indices are correctly constructed.
    for (u in unique(training$label)) {
        current <- u == training$label
        expect_identical(out$original.exprs[[u]], assay(training)[,current])
        expect_s4_class(out$nn.indices[[u]], "BiocNeighborIndex")

        # NOTE: not true in general, as zero-variance cells get discarded; see below.
        expect_identical(nrow(out$nn.indices[[u]]), sum(current))
    }

    # Checking the structure of the DE gene set.
    expect_identical(names(out$search$extra), unique(training$label))
    for (u in names(out$search$extra)) {
        expect_identical(names(out$search$extra[[u]]), names(out$search$extra))
        expect_identical(out$search$extra[[u]][[u]], character(0))

        # Genes in opposite directions should not intersect.
        for (j in names(out$search$extra)) {
            combined <- intersect(out$search$extra[[u]][[j]], out$search$extra[[j]][[u]])
            expect_identical(combined, character(0))
        }
    }

    expect_identical(training$common.genes, unlist(training$search$extra))
})

test_that("trainSingleR works correctly for genes='sd'", {
    out <- trainSingleR(training, training$label, genes='sd')
    expect_identical(out$search$mode, "sd")

    # Checking the structure of the extras (a median matrix).
    expect_identical(colnames(out$search$extra), unique(training$label))
    expect_identical(rownames(out$search$extra), rownames(training))

    # Checking the selected genes.
    expect_identical(out$common.genes, rownames(training)[rowSds(out$search$extra) > 1])
})

test_that("trainSingleR works correctly for genes='all'", {
    out <- trainSingleR(training, training$label, genes='all')
    expect_identical(out$common.genes, rownames(training))

    ref <- trainSingleR(training, training$label, genes='sd')
    expect_identical(ref$search, out$search)
})

test_that("trainSingleR works correctly for a list of lists of genes", {
    set.seed(92) # As NN index construction uses the random seed.
    ref <- trainSingleR(training, training$label, genes='de')

    set.seed(92)
    out <- trainSingleR(training, training$label, genes=SingleR:::.get_genes_by_de(assay(training), training$label))

    expect_identical(ref, out)
})

test_that("trainSingleR is robust to no-variance cells", {
    sce <- training
    assay(sce)[,1:10] <- 0
    out <- trainSingleR(sce, sce$label)

    for (u in unique(sce$label)) {
        current <- u == sce$label
        expect_identical(out$original.exprs[[u]], assay(sce)[,current])
        expect_s4_class(out$nn.indices[[u]], "BiocNeighborIndex")

        # Zero variance cells are lost.
        expect_identical(nrow(out$nn.indices[[u]]), sum(current[-(1:10)]))
    }
})

test_that("trainSingleR works on raw expression matrices", {
    set.seed(102)
    out <- trainSingleR(training, training$label)

    set.seed(102)
    alt <- trainSingleR(assay(training), training$label)
    expect_identical(out, alt)

    blah <- training
    assays(blah) <- list(stuff=matrix(0, nrow(blah), ncol(blah)), whee=assay(training))

    set.seed(102)
    re.alt <- trainSingleR(blah, blah$label, assay.type="whee")
    expect_identical(out, re.alt)
})

test_that("trainSingleR is invariant to simple transformations", {
    sce <- training
    assay(sce, "shifted") <- assay(sce) + 1
    assay(sce, "scaled") <- assay(sce) * 2

    set.seed(3523)
    out <- trainSingleR(sce, sce$label)
    same.fields <- setdiff(names(out), "original.exprs")

    set.seed(3523)
    alt <- trainSingleR(sce, sce$label, assay.type="shifted")
    expect_identical(out[same.fields], alt[same.fields])

    set.seed(3523)
    alt <- trainSingleR(sce, sce$label, assay.type="scaled")
    expect_identical(out[same.fields], alt[same.fields])

    # DE/SD magnitudes change upon log-transform, so don't test that.
    assay(sce, "logcounts") <- log(assay(sce) + 1)

    set.seed(3523)
    out2 <- trainSingleR(sce, sce$label, genes='all')
    set.seed(3523)
    alt2 <- trainSingleR(sce, sce$label, genes='all', assay.type="logcounts")

    out2$search$extra <- alt2$search$extra <- NULL
    expect_identical(out2[same.fields], alt2[same.fields])
})

test_that("trainSingleR behaves with silly inputs", {
    out <- trainSingleR(training[,0], training$label[0])
    expect_identical(out$common.genes, character(0))
    expect_identical(out$original.exprs, List())
    expect_identical(out$nn.indices, List())

    out <- trainSingleR(training[0,], training$label)
    expect_identical(out$common.genes, character(0))
    expect_identical(names(out$search$extra), unique(training$label))
    expect_identical(names(out$original.exprs), unique(training$label))
    expect_identical(names(out$nn.indices), unique(training$label))

    out <- trainSingleR(training[0,], training$label, genes='sd')
    expect_identical(out$common.genes, character(0))
    expect_identical(colnames(out$search$extra), unique(training$label))
    expect_identical(names(out$original.exprs), unique(training$label))
    expect_identical(names(out$nn.indices), unique(training$label))

    unnamed <- unname(training)
    expect_error(trainSingleR(unnamed, unnamed$label), "must have row names")
})
