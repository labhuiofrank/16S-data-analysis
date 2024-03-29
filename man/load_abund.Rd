% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading.R
\name{load_abund}
\alias{load_abund}
\title{Load OTU abundance file from mothur
(for either count_table or shared format)}
\usage{
load_abund(file)
}
\arguments{
\item{file}{Path to shared file}
}
\value{
Abundance dataframe (sample x OTU)
}
\description{
Load OTU abundance file from mothur
(for either count_table or shared format)
}
\examples{
load_abundance("abundance.100.shared")
load_abundance("abundance.100.rep.count_table")
}
