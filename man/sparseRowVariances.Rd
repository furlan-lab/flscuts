% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R, R/project_data.R
\name{sparseRowVariances}
\alias{sparseRowVariances}
\title{Sparse Row Variances}
\usage{
sparseRowVariances(m)

sparseRowVariances(m)
}
\arguments{
\item{m}{A sparse matrix of class \code{dgCMatrix}.}
}
\value{
A numeric vector of row variances.
}
\description{
Computes the variances of rows in a sparse matrix efficiently.
}
\references{
Granja, J. M.et al. (2019). Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype
acute leukemia. Nature Biotechnology, 37(12), 1458–1465.
}
\keyword{internal}
