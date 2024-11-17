\name{SingleRnews}
\title{SingleR News}
\encoding{UTF-8}

\section{Version 2.10.0}{\itemize{
\item Switch to \pkg{scrapper} for DE detection when \code{de.method="t"} or \code{de.method="wilcox"} in \code{trainSingleR()}.
This should give similar results to but is faster than the previous \pkg{scran} functions. 

\item Switch to \pkg{scrapper} for variance calculation, PCA and clustering in \code{aggregateReference()}.
This should be faster than the previous \pkg{BiocSingular} and \code{stats::kmeans} functions, and avoids the need to consider the seed.
}}

\section{Version 2.8.0}{\itemize{
\item Added the \code{test.genes=} argument to \code{trainSingleR()}, to restrict marker detection to only those genes in the test dataset.
This is also checked against \code{rownames(test)} in \code{classifySingleR()} to ensure that the test's feature space is consistent with the space used during training.

\item The introduction of \code{test.genes=} means that we no longer need to explicitly subset the rows of the reference dataset (to match the test features) in \code{SingleR()}.
This saves memory by avoiding an unnecessary copy of the reference dataset, but may also slightly alter the marker selection as ties are broken in a different way.
Namely, if the top X genes are used as markers, and the X-th and (X+1)-th gene have the same log-fold change, 
tie breaking will be based on the ordering of the rows in the reference matrix - which is no longer the same as in the previous version of \pkg{SingleR}.
This results in some slight differences in the markers that propagate down to the classification results.

\item Restored the \code{BNPARAM=} argument in \code{trainSingleR()}, to enable more fine-grained specification of neighbor search algorithms.
The \code{approximate=} argument is deprecated. 

\item Soft-deprecated \code{check.missing=} in \code{classifySingleR()} and \code{combineRecomputedResults()}.
This is because any filtering will cause a mismatch between the row names of \code{tests} and the \code{test.genes} in \code{trained}.
Rather, filtering should be done prior to \code{trainSingleR()}, as is done in the main \code{SingleR()} function.

\item \code{combineRecomputedResults()} now supports fine-tuning to resolve closely-related labels from different references.
This is similar to the fine-tuning in \code{classifySingleR()} where the feature space is iterately redefined as the union of markers of labels with near-highest scores.

\item Added the \code{plotMarkerHeatmap()} function to plot a diagnostic heatmap of the most interesting markers for each label.
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

\item Streamlined the vignette to point to the book at \url{https://bioconductor.org/books/devel/SingleRBook/}.

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
