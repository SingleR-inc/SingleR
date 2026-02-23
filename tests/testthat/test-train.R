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
    effects <- scrapper::scoreMarkers(logcounts(training), training$label, all.pairwise=TRUE)

    VERIFY <- function(ref.markers, effect.sizes, hard.limit, extra, de.n = 10) {
        sulabels <- sort(unique(training$label))
        expect_identical(sort(names(ref.markers)), sulabels)

        for (n in names(ref.markers)) {
            current.markers <- ref.markers[[n]]
            expect_identical(sort(names(current.markers)), sulabels)

            for (n2 in names(current.markers)) {
                my.markers <- current.markers[[n2]]
                if (n == n2) {
                    expect_identical(length(my.markers), 0L)
                } else {
                    expect_gt(length(my.markers), 0L)

                    my.effects <- effect.sizes[n2, n,]
                    keep <- which(my.effects > hard.limit)
                    o <- order(my.effects[keep], decreasing=TRUE)
                    expect_identical(rownames(training)[keep[head(o, de.n)]], my.markers)

                    is.chosen <- rownames(training) %in% my.markers
                    min.chosen <- min(my.effects[is.chosen])
                    expect_gte(min.chosen, max(my.effects[!is.chosen]))
                    expect_gt(min.chosen, hard.limit)

                    if (!is.null(extra)) {
                        extra(n, n2, my.markers)
                    }
                }
            }
        }
    }

    # For Wilcox.
    ref <- trainSingleR(training, training$label, genes='de', de.method="wilcox")
    VERIFY(ref$markers$full, effects$auc, 0.5, extra=NULL)

    # For t-tests.
    ref <- trainSingleR(training, training$label, genes='de', de.method="t")
    VERIFY(ref$markers$full, effects$cohens.d, 0, extra=function(n, n2, markers) {
        expect_lte(length(markers), 10)
        left <- Matrix::rowMeans(logcounts(training)[markers, training$label == n])
        right <- Matrix::rowMeans(logcounts(training)[markers, training$label == n2])
        expect_true(all(left > right))
    })

    # Behaves with a large number of top genes.
    ref <- trainSingleR(training, training$label, genes='de', de.method="t", de.n=10000)
    VERIFY(
        ref$markers$full,
        effects$cohens.d,
        0,
        extra=function(n, n2, markers) {
            expect_gt(length(markers), 10)
            left <- Matrix::rowMeans(logcounts(training)[markers, training$label == n])
            right <- Matrix::rowMeans(logcounts(training)[markers, training$label == n2])
            expect_true(all(left > right))
        },
        de.n=10000
    )

    # Responds to threshold specification.
    thresh.effects <- scrapper::scoreMarkers(logcounts(training), training$label, threshold=1, all.pairwise=TRUE)
    ref <- trainSingleR(training, training$label, genes='de', de.method="t", de.args=list(threshold=1))
    VERIFY(ref$markers$full, thresh.effects$cohens.d, 0, extra=function(n, n2, markers) {
        left <- Matrix::rowMeans(logcounts(training)[markers, training$label == n])
        right <- Matrix::rowMeans(logcounts(training)[markers, training$label == n2])
        expect_true(all(left > right + 1))
    })
})

test_that("trainSingleR is robust to non-character labels", {
    ids <- sample(1:5, ncol(training), replace=TRUE)
    out <- trainSingleR(training, ids)
    ref <- trainSingleR(training, as.character(ids))
    expect_equal(out$labels, ref$labels)
})

test_that("trainSingleR behaves correctly with gene intersections", {
    random.test <- sample(rownames(test), 500)
    out <- trainSingleR(training, training$label, test.genes=random.test)
    ref <- trainSingleR(training[rownames(training) %in% random.test,], training$label)
    expect_identical(out$labels, ref$labels)
    expect_identical(out$markers$full, ref$markers$full)
    expect_identical(sort(out$markers$unique), sort(ref$markers$unique))
    expect_identical(out$ref[rownames(training) %in% random.test,], ref$ref)
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

    expect_identical(as.matrix(out$ref), ref$ref)
    expect_identical(out$markers, ref$markers)
})

test_that("trainSingleR behaves with multiple references, plus recomputation", {
    training1 <- training2 <- training
    training1 <- training1[sample(nrow(training1)),]
    rownames(training1) <- rownames(training)

    ref1 <- trainSingleR(training1, training1$label)
    ref2 <- trainSingleR(training2, training2$label)
    out <- trainSingleR(list(training1, training2), list(training1$label, training2$label))

    except.built <- setdiff(names(ref1), "built")
    expect_identical(ref1[except.built], out[[1]][except.built])
    expect_identical(ref2[except.built], out[[2]][except.built])

    # Same result with names.
    out <- trainSingleR(list(foo=training1, bar=training2), list(training1$label, training2$label))
    expect_identical(names(out), c("foo", "bar"))
})

test_that("trainSingleR behaves with aggregation turned on", {
    suppressWarnings(out <- trainSingleR(training, training$label, aggr.ref=TRUE))
    expect_true(ncol(out$ref) <= ncol(training))

    suppressWarnings(out2 <- trainSingleR(ref=list(training, training), label=list(training$label, training$label), aggr.ref=TRUE))
    expect_identical(out2[[1]]$ref, out$ref)
    expect_identical(out2[[2]]$ref, out$ref)
})

test_that("trainSingleR behaves with silly inputs", {
    expect_error(out <- trainSingleR(training[,0], training$label[0]), "at least one column")

    # R itself is buggy right now w.r.t. handling dimnames when the extent is zero.
    #out <- trainSingleR(training[0,], training$label)
    #expect_identical(length(out$markers$unique), 0L)

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
    expect_identical(as.matrix(out$ref), ref$ref)
})

test_that("trainSingleR auto-eliminates duplicate row names", {
    expect_identical(anyDuplicated(rownames(training)), 0L) # make sure there weren't any to begin with.

    rownames(training) <- head(rep(rownames(training), each=2), nrow(training))
    out <- trainSingleR(training, training$label)

    keep <- !duplicated(rownames(training))
    ref <- trainSingleR(training[keep,], training$label)
    expect_identical(out$labels, ref$labels)
    expect_identical(as.matrix(out$ref), ref$ref)
})
