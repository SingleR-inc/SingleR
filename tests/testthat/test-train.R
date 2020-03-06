# This tests the trainSingleR function.
# library(testthat); library(SingleR); source("setup.R"); source("test-train.R")

test_that("trainSingleR works correctly for genes='de'", {
    out <- trainSingleR(training, training$label)
    expect_identical(out$search$mode, "de")

    # Checking that the original expression is correctly returned,
    # and that the NN indices are correctly constructed.
    for (u in unique(training$label)) {
        current <- u == training$label
        expect_identical(out$original.exprs[[u]], logcounts(training)[,current]+0)
        expect_s4_class(out$nn.indices[[u]], "BiocNeighborIndex")
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

    # Because the search is equivalent to doing a fixed 'SD' search,
    # there's no DE information for fine-tuning later.
    ref <- trainSingleR(training, training$label, genes='sd')
    expect_identical(ref$search, out$search)
})

test_that("trainSingleR works correctly for a list of lists of genes", {
    set.seed(92) # As NN index construction uses the random seed.
    ref <- trainSingleR(training, training$label, genes='de')

    collected <- SingleR:::.get_genes_by_de(logcounts(training), training$label)
    expect_identical(sort(names(collected)), sort(unique(training$label)))
    in.names <- unique(lapply(collected, names))
    expect_identical(length(in.names), 1L)
    expect_identical(in.names[[1]], names(collected))

    set.seed(92)
    out <- trainSingleR(training, training$label, genes=collected)

    expect_identical(ref, out)
})

test_that("trainSingleR works correctly for a list of genes", {
    collected <- SingleR:::.get_genes_by_de(logcounts(training), training$label)

    set.seed(92)
    ref <- trainSingleR(training, training$label, genes=collected)

    set.seed(92)
    re.collected <- lapply(collected, unlist, use.names=FALSE)
    out <- trainSingleR(training, training$label, genes=re.collected)

    expect_identical(ref$common.genes, out$common.genes)
    expect_identical(names(ref$search$extra), names(out$search$extra))
    expect_identical(lapply(ref$search$extra, names), lapply(out$search$extra, names))
})

test_that("trainSingleR works correctly for a list of lists of genes", {
    set.seed(92) # As NN index construction uses the random seed.
    ref <- trainSingleR(training, training$label, genes='de')

    set.seed(92)
    markers <- SingleR:::.get_genes_by_de(logcounts(training), training$label)
    out <- trainSingleR(training, training$label, genes=markers)
    expect_identical(ref, out)

    # Same results if we get a List of List of genes, which is correctly coerced to ordinary lists.
    set.seed(92)
    markers2 <- List(lapply(markers, List))
    out2 <- trainSingleR(training, training$label, genes=markers)
    expect_identical(ref, out2)

    # Fails when a weird gene set input is provided.
    expect_error(trainSingleR(training, training$label, genes=list(A=list(), B=character(0))), "'genes' must be")
    expect_error(trainSingleR(training, training$label, genes=list(A=list(), B=list())), "for each label")

    empty <- rep(list(list()), length(unique(training$label)))
    names(empty) <- unique(training$label)
    expect_error(trainSingleR(training, training$label, genes=empty), "between each pair of labels")
})

test_that("trainSingleR works correctly for other DE testing methods", {
    # For Wilcox.
    by.t <- scran::pairwiseWilcox(logcounts(training), training$label, direction="up") # do NOT move below set.seed(); some BiocParallel setup changes the seed.
    markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=10)

    set.seed(92) 
    ref <- trainSingleR(training, training$label, genes='de', de.method="wilcox")
    
    set.seed(92) 
    trained <- trainSingleR(training, training$label, genes=markers)

    expect_identical(ref, trained)

    # For t-tests.
    set.seed(102) 
    ref <- trainSingleR(training, training$label, genes='de', de.method="t")

    set.seed(102) 
    by.t <- scran::pairwiseTTests(logcounts(training), training$label, direction="up")
    markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=10)
    trained <- trainSingleR(training, training$label, genes=markers)

    expect_identical(ref, trained)

    # Responds to the requested number of genes.
    set.seed(102) 
    ref <- trainSingleR(training, training$label, genes='de', de.method="t", de.n=20, de.args=list(lfc=1))

    set.seed(102) 
    by.t <- scran::pairwiseTTests(logcounts(training), training$label, direction="up", lfc=1)
    markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=20)
    trained <- trainSingleR(training, training$label, genes=markers)

    expect_identical(ref, trained)
})

test_that("trainSingleR is robust to no-variance cells", {
    sce <- training
    logcounts(sce)[,1:10] <- 0
    out <- trainSingleR(sce, sce$label)

    for (u in unique(sce$label)) {
        current <- u == sce$label
        expect_identical(out$original.exprs[[u]], logcounts(sce)[,current])
        expect_s4_class(out$nn.indices[[u]], "BiocNeighborIndex")
        expect_identical(nrow(out$nn.indices[[u]]), sum(current))
    }
})

test_that("trainSingleR is robust to non-character labels", {
    ids <- sample(1:5, ncol(training), replace=TRUE)

    set.seed(999)
    out <- trainSingleR(training, ids)

    set.seed(999)
    ref <- trainSingleR(training, as.character(ids))

    expect_equal(out, ref)
})

test_that("trainSingleR works on raw expression matrices", {
    set.seed(102)
    out <- trainSingleR(training, training$label)

    set.seed(102)
    alt <- trainSingleR(logcounts(training), training$label)
    expect_identical(out, alt)

    blah <- training
    assays(blah, withDimnames=FALSE) <- list(stuff=matrix(0, nrow(blah), ncol(blah)), whee=logcounts(training))

    set.seed(102)
    re.alt <- trainSingleR(blah, blah$label, assay.type="whee")
    expect_identical(out, re.alt)
})

test_that("trainSingleR is invariant to simple transformations", {
    sce <- training
    assay(sce, "shifted") <- logcounts(sce) + 1
    assay(sce, "scaled") <- logcounts(sce) * 2

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
    assay(sce, "double_log") <- log(logcounts(sce) + 1)

    set.seed(3523)
    out2 <- trainSingleR(sce, sce$label, genes='all')
    set.seed(3523)
    alt2 <- trainSingleR(sce, sce$label, genes='all', assay.type="double_log")

    out2$search$extra <- alt2$search$extra <- NULL
    expect_identical(out2[same.fields], alt2[same.fields])
})

test_that("trainSingleR behaves with NAs", {
    sce <- training
    logcounts(sce)[1,1] <- NA

    set.seed(30101)
    expect_warning(out <- trainSingleR(sce, sce$label), "missing values")
    set.seed(30101)
    ref <- trainSingleR(sce[-1,], sce$label)

    expect_identical(out, ref)
})

test_that("trainSingleR behaves with multiple references", {
    set.seed(1000)
    ref1 <- trainSingleR(training, training$label)
    ref2 <- trainSingleR(training, training$label)
    set.seed(1000)
    out <- trainSingleR(list(training, training), list(training$label, training$label))

    expect_identical(ref1, out[[1]])
    expect_identical(ref2, out[[2]])

    # Checking that the union of common genes are taken correctly 
    # by scrambling the genes and makeing sure that we get the union.
    training1 <- training2 <- training
    training1 <- training1[sample(nrow(training1)),]
    rownames(training1) <- rownames(training)

    set.seed(2000)
    ref1 <- trainSingleR(training1, training1$label)
    ref2 <- trainSingleR(training2, training2$label)
    set.seed(2000)
    out <- trainSingleR(list(training1, training2), list(training1$label, training2$label))

    expect_identical(out[[1]]$search, ref1$search)
    expect_identical(out[[2]]$search, ref2$search)
    expect_false(identical(sort(ref1$common.genes), sort(ref2$common.genes)))
    expect_identical(out[[1]]$common.genes, union(ref1$common.genes, ref2$common.genes))
    expect_identical(out[[1]]$common.genes, out[[2]]$common.genes)

    # Works with pre-specified marker lists.
    markers <- SingleR:::.get_genes_by_de(logcounts(training), training$label)
    markers2 <- SingleR:::.get_genes_by_de(logcounts(training), training$label, de.method="t")
    markers2 <- lapply(markers2, unlist, use.names=FALSE)

    set.seed(2000)
    ref1 <- trainSingleR(training1, training1$label, genes=markers)
    ref2 <- trainSingleR(training2, training2$label, genes=markers2)
    set.seed(2000)
    out <- trainSingleR(list(training1, training2), list(training1$label, training2$label), genes=list(markers, markers2))

    expect_identical(out[[1]]$search, ref1$search)
    expect_identical(out[[2]]$search, ref2$search)
    expect_false(identical(sort(ref1$common.genes), sort(ref2$common.genes)))
    expect_identical(out[[1]]$common.genes, union(ref1$common.genes, ref2$common.genes))
    expect_identical(out[[1]]$common.genes, out[[2]]$common.genes)

    # Throws errors correctly.
    expect_error(trainSingleR(list(training1, training2), list(training1$label)), "same length")
    expect_error(trainSingleR(list(training1, training2), list(training1$label, training2$label), genes=list(training1$label)), "same length")
    expect_error(trainSingleR(list(training1, training2[1:10,]), list(training1$label)), "not identical")
})

test_that("trainSingleR behaves with aggregation turned on", {
    set.seed(10000)
    out <- trainSingleR(training, training$label, aggr.ref=TRUE)
    expect_true(sum(vapply(out$nn.indices, nrow, 0L)) <= ncol(training))
    expect_identical(out$search$mode, "de")

    set.seed(10000)
    out2 <- trainSingleR(ref=list(training, training), label=list(training$label, training$label), aggr.ref=TRUE)
    expect_identical(out2[[1]], out)
    expect_false(identical(out2[[2]], out)) # different k-means initialization.
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

test_that("trainSingleR works when more genes are supplied than in reference", {
    train.sub <- head(training, 90)
    collected <- SingleR:::.get_genes_by_de(logcounts(training), training$label)
    genes <- unique(unlist(collected))
    
    # Make sure more genes than ref
    expect_false(all(genes %in% row.names(train.sub)))
    set.seed(100)
    expect_error(out <- SingleR::trainSingleR(train.sub, training$label, genes = collected), NA)

    # Behaves the same as if those genes were intersected.
    set.seed(100)
    collected2 <- lapply(collected, function(l) lapply(l, intersect, y=rownames(train.sub))) 
    ref <- SingleR::trainSingleR(train.sub, training$label, genes = collected2)
    expect_identical(out, ref)
})
