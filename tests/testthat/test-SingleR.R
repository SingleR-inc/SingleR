# End-to-end testing of SingleR on a dummy training/test dataset that should be near-perfectly assigned.
# library(testthat); library(SingleR); source("setup.R")

test_that("SingleR works in DE mode", {
    out <- SingleR(test=test, training=training, labels=training$label, genes="de")
    tab <- table(out$labels, test$label)
    expect_true(sum(diag(tab))/sum(tab) > 0.95)
})

test_that("SingleR works in SD mode", {
    out <- SingleR(test=test, training=training, labels=training$label, genes="sd")
    tab <- table(out$labels, test$label)
    expect_true(sum(diag(tab))/sum(tab) > 0.95)
})

test_that("SingleR works with custom gene selection", {
    all.labs <- unique(training$label)
    collected <- rep(list(tail(rownames(training), 100)), length(all.labs))
    names(collected) <- all.labs

    more.collected <- rep(list(collected), length(all.labs))
    names(more.collected) <- all.labs

    out <- SingleR(test=test, training=training, labels=training$label, genes=more.collected)
    tab <- table(out$labels, test$label)
    expect_true(sum(diag(tab))/sum(tab) > 0.95)
    
    # We should get, in this case, the same result with a list of vectors.
    out2 <- SingleR(test=test, training=training, labels=training$label, genes=collected)
    expect_identical(out, out2)
})

