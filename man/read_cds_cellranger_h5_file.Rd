% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.R
\name{read_cds_cellranger_h5_file}
\alias{read_cds_cellranger_h5_file}
\title{Read Cell Ranger H5 File into Monocle3 cell_data_set}
\usage{
read_cds_cellranger_h5_file(h5_file)
}
\arguments{
\item{h5_file}{Path to the Cell Ranger H5 file.}
}
\value{
A \code{cell_data_set} object containing the expression data and associated metadata.
}
\description{
Reads a Cell Ranger output H5 file and converts it into a \code{cell_data_set} object for use with Monocle3.
}
\details{
This function reads the H5 file produced by Cell Ranger and constructs a \code{cell_data_set} object, including expression matrices, cell metadata, and gene metadata.
}
\examples{
\dontrun{
cds <- read_cds_cellranger_h5_file("path/to/filtered_feature_bc_matrix.h5")
}
}
\keyword{internal}
