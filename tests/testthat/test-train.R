# This tests the trainSingleR function.
# library(testthat); library(SingleR); source("setup.R"); source("test-train.R")

test_that("trainSingleR works correctly for genes='de'", {
    out <- trainSingleR(training, training$label)

    expect_identical(out$ref, logcounts(training))
    expect_identical(out$labels$full, training$label)
    expect_identical(out$labels$unique, sort(unique(training$label)))
    expect_identical(sort(out$markers$unique), sort(unique(unlist(out$markers$full))))

    # Checking the structure of the DE gene set.
    expect_identical(names(out$markers$full), sort(unique(training$label)))

    for (u in names(out$markers$full)) {
        expect_identical(names(out$markers$full[[u]]), names(out$markers$full))
        expect_identical(out$markers$full[[u]][[u]], character(0))

        # Genes in opposite directions should not intersect.
        for (j in names(out$markers$full)) {
            combined <- intersect(out$markers$full[[u]][[j]], out$markers$full[[j]][[u]])
            expect_identical(combined, character(0))
        }
    }
})

test_that("trainSingleR works correctly for a list of lists of genes", {
    collected <- SingleR:::.get_genes_by_de(logcounts(training), training$label, de.n=13)
    expect_identical(sort(names(collected)), sort(unique(training$label)))

    in.names <- unique(lapply(collected, names))
    expect_identical(length(in.names), 1L)
    expect_identical(in.names[[1]], names(collected))

    out <- trainSingleR(training, training$label, genes=collected)
    expect_identical(out$markers$full, collected)
})

test_that("trainSingleR works correctly for a list of genes", {
    collected <- SingleR:::.get_genes_by_de(logcounts(training), training$label)
    re.collected <- lapply(collected, unlist, use.names=FALSE)
    out <- trainSingleR(training, training$label, genes=re.collected)

    expect_identical(sort(out$markers$unique), sort(unique(unlist(out$markers$full))))

    for (u in names(out$markers$full)) {
        for (j in names(out$markers$full)) {
            if (u == j) {
                expect_equal(out$markers$full[[u]][[u]], unique(re.collected[[u]]))
            } else {
                expect_identical(out$markers$full[[u]][[j]], character(0))
            }
        }
    }
})

test_that("trainSingleR fails correctly for a list of lists of genes", {
    # Fails when a weird gene set input is provided.
    expect_error(trainSingleR(training, training$label, genes=list(A=list(), B=character(0))), "'genes' must be")
    expect_error(trainSingleR(training, training$label, genes=list(A=list(), B=list())), "for each label")

    empty <- rep(list(list()), length(unique(training$label)))
    names(empty) <- unique(training$label)
    expect_error(trainSingleR(training, training$label, genes=empty), "between each pair of labels")
})

test_that("trainSingleR works correctly for other DE testing methods", {
    # For Wilcox.
    by.t <- scran::pairwiseWilcox(logcounts(training), training$label, direction="up")
    markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=10)

    ref <- trainSingleR(training, training$label, genes='de', de.method="wilcox")
    trained <- trainSingleR(training, training$label, genes=markers)
    expect_identical(ref$markers, trained$markers)

    # For t-tests.
    by.t <- scran::pairwiseTTests(logcounts(training), training$label, direction="up")
    markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=10)

    ref <- trainSingleR(training, training$label, genes='de', de.method="t")
    trained <- trainSingleR(training, training$label, genes=markers)
    expect_identical(ref$markers, trained$markers)

    # Responds to the requested number of genes.
    by.t <- scran::pairwiseTTests(logcounts(training), training$label, direction="up", lfc=1)
    markers <- scran::getTopMarkers(by.t[[1]], by.t[[2]], n=20)

    ref <- trainSingleR(training, training$label, genes='de', de.method="t", de.n=20, de.args=list(lfc=1))
    trained <- trainSingleR(training, training$label, genes=markers)
    expect_identical(ref$markers, trained$markers)
})

test_that("trainSingleR is robust to non-character labels", {
    ids <- sample(1:5, ncol(training), replace=TRUE)
    out <- trainSingleR(training, ids)
    ref <- trainSingleR(training, as.character(ids))
    expect_equal(out$labels, ref$labels)
})

test_that("trainSingleR works on various expression matrices", {
    out <- trainSingleR(training, training$label)
    alt <- trainSingleR(logcounts(training), training$label)
    expect_identical(out$ref, alt$ref)

    # assay.type= works.
    blah <- training
    assays(blah, withDimnames=FALSE) <- list(stuff=matrix(0, nrow(blah), ncol(blah)), whee=logcounts(training))
    re.alt <- trainSingleR(blah, blah$label, assay.type="whee")
    expect_identical(out$ref, re.alt$ref)

    # robust to invariant transformations.
    sce <- training
    assay(sce, "shifted") <- logcounts(sce) + 1
    assay(sce, "scaled") <- logcounts(sce) * 2

    out <- trainSingleR(sce, sce$label)
    alt <- trainSingleR(sce, sce$label, assay.type="shifted")
    expect_identical(out$markers, alt$markers)

    alt <- trainSingleR(sce, sce$label, assay.type="scaled")
    expect_identical(out$markers, alt$markers)
})

test_that("trainSingleR strips out NAs", {
    sce <- training
    logcounts(sce)[1,1] <- NA

    out <- trainSingleR(sce, sce$label)
    ref <- trainSingleR(sce[-1,], sce$label)

    expect_identical(out$ref, ref$ref)
})

test_that("trainSingleR behaves with multiple references, plus recomputation", {
    training1 <- training2 <- training
    training1 <- training1[sample(nrow(training1)),]
    rownames(training1) <- rownames(training)
    
    set.seed(1000)
    ref1 <- trainSingleR(training1, training1$label)
    ref2 <- trainSingleR(training2, training2$label)

    set.seed(1000)
    out <- trainSingleR(list(training1, training2), list(training1$label, training2$label), recompute=TRUE)

    except.built <- setdiff(names(ref1), "built")
    expect_identical(ref1[except.built], out[[1]][except.built])
    expect_identical(ref2[except.built], out[[2]][except.built])
})

test_that("trainSingleR behaves with aggregation turned on", {
    set.seed(10000)
    suppressWarnings(out <- trainSingleR(training, training$label, aggr.ref=TRUE))
    expect_true(ncol(out$ref) <= ncol(training))

    set.seed(10000)
    suppressWarnings(out2 <- trainSingleR(ref=list(training, training), label=list(training$label, training$label), aggr.ref=TRUE))
    expect_identical(out2[[1]]$ref, out$ref)
    expect_false(identical(out2[[2]]$ref, out$ref)) # different k-means initialization.
})

test_that("trainSingleR behaves with silly inputs", {
    expect_error(out <- trainSingleR(training[,0], training$label[0]), "at least one column")

    out <- trainSingleR(training[0,], training$label)
    expect_identical(length(out$markers$unique), 0L)

    unnamed <- unname(training)
    expect_error(trainSingleR(unnamed, unnamed$label), "must have row names")
})

test_that("trainSingleR works when 'genes' contains markers outside of the reference", {
    train.sub <- head(training, 90)
    collected <- SingleR:::.get_genes_by_de(logcounts(training), training$label)
    genes <- unique(unlist(collected))
    
    # Make sure more genes than ref
    expect_false(all(genes %in% row.names(train.sub)))
    expect_error(out <- SingleR::trainSingleR(train.sub, training$label, genes = collected), NA)

    # Behaves the same as if those genes were intersected.
    collected2 <- lapply(collected, function(l) lapply(l, intersect, y=rownames(train.sub))) 
    ref <- SingleR::trainSingleR(train.sub, training$label, genes = collected2)
    expect_identical(out$markers, ref$markers)
})

test_that("trainSingleR works when restricting", {
    keep <- c(letters, head(rownames(training), 90))
    expect_error(out <- SingleR::trainSingleR(training, training$label, restrict=keep), NA)

    # Behaves the same as if those genes were intersected.
    ref <- SingleR::trainSingleR(head(training, 90), training$label)
    expect_identical(out$markers, ref$markers)
})

test_that("trainSingleR auto-eliminates NA labels", {
    populate <- rbinom(length(training$label), 1, 0.2)==1
    training$label[populate] <- NA

    out <- trainSingleR(training, training$label)
    ref <- trainSingleR(training[,!populate], training$label[!populate])
    expect_identical(out$labels, ref$labels)
    expect_identical(out$ref, ref$ref)
})

