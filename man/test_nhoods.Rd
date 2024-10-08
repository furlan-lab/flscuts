% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/milo.R
\name{test_nhoods}
\alias{test_nhoods}
\title{Test Milo Neighborhoods for Differential Abundance}
\usage{
test_nhoods(sce_milo, formula_str)
}
\arguments{
\item{sce_milo}{A list containing a Milo object and a design matrix, as returned by \code{\link{make_nhoods}}.}

\item{formula_str}{A character string representing the design formula for differential testing (e.g., "~ condition").}
}
\value{
A data frame of differential abundance results for the neighborhoods.
}
\description{
Tests Milo neighborhoods for differential abundance using a specified design formula.
}
\details{
This function performs differential abundance testing on the neighborhoods using the specified design formula.
}
\examples{
# Example usage:
# da_results <- test_nhoods(sce_milo, formula_str = "~ condition")

}
