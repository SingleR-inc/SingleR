\name{SingleRnews}
\title{SingleR News}
\encoding{UTF-8}

\section{Version 2.8.0}{\itemize{
\item Added the \code{test.genes=} argument to \code{trainSingleR()}, to restrict marker detection to only those genes in the test dataset.
This is also checked against \code{rownames(test)} in \code{classifySingleR()} to ensure that the test's feature space is consistent with the space used during training.

\item Restored the \code{BNPARAM=} argument in \code{trainSingleR()}, to enable more fine-grained specification of neighbor search algorithms.
The \code{approximate=} argument is deprecated. 
}}

\section{Version 2.0.0}{\itemize{
\item The format of the output of \code{trainSingleR()} has changed and is no longer back-compatible.
\item \code{recompute=FALSE} in \code{trainSingleR()} does nothing; all integrated analyses are now done with \code{recompute=TRUE}.
To that end, \code{combineCommonResults()} is also deprecated.
\item \code{genes = "sd"} and its associated options in \code{trainSingleR()} are no longer supported.
\item \code{first.labels} is no longer reported in \code{classifySingleR()}.
\item Added another parallelization mechanism via \code{num.threads=} and C++11 threads.
This should be much more memory efficient than using \pkg{BiocParallel}.
\item \code{combineRecomputedScores()} will automatically handle mismatches in the input references by default.
}}

\section{Version 1.6.0}{\itemize{
\item Relaxed the requirements for consistent row names in \code{combineRecomputedResults()}.

\item Support sparse DelayedArray inputs in \code{classifySingleR()}.

\item Parallelize over labels instead of rows in \code{aggregateReference()}, with minor changes in the setting of the seed.
Restrict the PCA to the top 1000 most highly variable genes, for speed.
}}

\section{Version 1.4.0}{\itemize{
\item Migrated all of the dataset getter functions to the \pkg{celldex} package.

\item Streamlined the vignette to point to the book at \url{http://bioconductor.org/books/devel/SingleRBook/}.

\item Added a \code{restrict=} argument to \code{trainSingleR()} and \code{SingleR()} to easily restrict to a subset of features.

\item Deprecated the \code{method=} argument in \code{SingleR()}.

\item Protect against accidental data.frames in \code{ref=} or \code{test=} in all functions.
}}

\section{Version 1.2.0}{\itemize{
\item Added support for consolidating labels from multiple references via \code{combineResults()}.

\item Added mappings to standardized Cell Ontology terms in all \code{*Data()} functions.

\item Changed the name of the \code{labels} input of \code{plotScoreDistribution()} to \code{labels.use} for consistency across functions.

\item Fixed a label from adipocytes to astrocytes in \code{BlueprintEncodeData()}.

\item Removed umlauts from labels (e.g., naive) in \code{NovershternHematopoieticData()} to avoid problems with Windows.

\item Perform PCA before clustering in \code{aggregateReference()} for speed and memory efficiency.

\item Modified \code{genes="all"} behavior in \code{trainSingleR()} to report DE-based markers for fine-tuning only.
}}

\section{Version 1.0.0}{\itemize{
\item New package SingleR for cell type annotation.
}}
