% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading.R
\name{find_pipeline_files}
\alias{find_pipeline_files}
\title{Check if the output folder from the cmaiki pipeline
contains all the necessary files (abund, tax) for a
given OTU identity theshold and returns
the paths in a list}
\usage{
find_pipeline_files(folder, otu_id = 100, meta_file = NULL)
}
\arguments{
\item{folder}{Path to pipeline output folder}

\item{otu_id}{OTU Identity threshold}

\item{meta_file}{Path to metadata file. Look for csv
files in folder if NULL}
}
\value{
list of file paths
}
\description{
Check if the output folder from the cmaiki pipeline
contains all the necessary files (abund, tax) for a
given OTU identity theshold and returns
the paths in a list
}
\examples{
find_pipeline_files("2022-07-01/")
}
