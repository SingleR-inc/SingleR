\name{SingleRnews}
\title{SingleR News}
\encoding{UTF-8}

\section{Version 1.2.0}{\itemize{
\item Added support for consolidating labels from multiple references via \code{combineResults()}.

\item Added mappings to standardized Cell Ontology terms in all \code{*Data()} functions.

\item Changed the name of the \code{labels} input of \code{plotScoreDistribution()} to \code{labels.use} for consistency across functions.

\item Fixed a label from adipocytes to astrocytes in \code{BlueprintEncodeData()}.

\item Removed umlauts from labels (e.g., naive) in \code{NovershternHematopoieticData()} to avoid problems with Windows.

\item Perform PCA before clustering in \code{aggregateReference()} for speed and memory efficiency.
}}

\section{Version 1.0.0}{\itemize{
\item New package SingleR for cell type annotation.
}}
