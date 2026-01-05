library(SingleCellExperiment)
set.seed(100)

###########################################
## Mocking up some example training data ##
###########################################

Ngroups <- 5
Ngenes <- 1000
means <- matrix(rnorm(Ngenes*Ngroups), nrow=Ngenes)
means[1:900,] <- 0
colnames(means) <- LETTERS[1:5]

N <- 100
g <- sample(LETTERS[1:5], N, replace=TRUE)
training <- SingleCellExperiment(
    list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
    colData=DataFrame(label=g)
)

rownames(training) <- sprintf("GENE_%s", seq_len(nrow(training)))
training <- scrapper::normalizeRnaCounts.se(training, more.norm.args=list(delayed=FALSE))

##################################################
## Mocking up some test data for classification ##
##################################################

N <- 100
g <- sample(LETTERS[1:5], N, replace=TRUE)
test <- SingleCellExperiment(
    list(counts=matrix(rpois(1000*N, lambda=2^means[,g]), ncol=N)),
    colData=DataFrame(label=g)
)

rownames(test) <- sprintf("GENE_%s", seq_len(nrow(test)))
test <- scrapper::normalizeRnaCounts.se(test, more.norm.args=list(delayed=FALSE))

##################################################
## Setting up flush tests for inadvertent DA BP ##
##################################################

library(BiocParallel)
failgen <- setRefClass("FailParam",
    contains="BiocParallelParam",
    fields=list(),
    methods=list())

FAIL <- failgen()
register(FAIL)

library(DelayedArray)
setAutoBPPARAM(FAIL)
