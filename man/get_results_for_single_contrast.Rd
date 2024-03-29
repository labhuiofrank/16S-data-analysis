% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DESeq2.R
\name{get_results_for_single_contrast}
\alias{get_results_for_single_contrast}
\title{Get result table from DESeq2 object between given variable level vs reference}
\usage{
get_results_for_single_contrast(diagdds, factor_name, ref, level)
}
\arguments{
\item{diagdds}{DESeq2 object}

\item{factor_name}{condition variable to test}

\item{ref}{reference level for variable}

\item{level}{other level for variable}
}
\value{
result dataframe
}
\description{
Get result table from DESeq2 object between given variable level vs reference
}
\examples{
get_results_for_single_contrast(diagdds, "Group", "A", "B")
}
