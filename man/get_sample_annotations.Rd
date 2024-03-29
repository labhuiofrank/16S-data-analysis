% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap-util.R
\name{get_sample_annotations}
\alias{get_sample_annotations}
\title{Format sample metadata from a phyloseq object for the ComplexHeatmap package.}
\usage{
get_sample_annotations(ps, columns = NULL, extra = NULL, colors = NULL)
}
\arguments{
\item{ps}{Phyloseq object}

\item{columns}{Columns from sample metadata to display in the heatmap}

\item{extra}{(list) Extra columns (not in the metadata) to include.
The number of values must match the number of samples in the phyloseq object.}

\item{colors}{(list) Custom colors to use for any of the metadata in the
selected columns or extra metadata. Fields are variable names and values are
named vectors (names: levels, values: colors).}
}
\value{
List object with each variable as factor and an extra "col" field
with custom colors for each of those variable, as provided in the colors
parameter.
}
\description{
Format sample metadata from a phyloseq object for the ComplexHeatmap package.
}
\examples{
get_sample_annotations(ps, c("Site_id", "Season"))
get_sample_annotations(ps, "Site_id", extra=list(pH=rep(7, nsamples(ps))))
get_sample_annotations(
  ps, "Site_id", 
  colors=list(Site_id=c(S1="red", S2="blue", S3="green"))
)
}
