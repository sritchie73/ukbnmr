#' Tools for processing Nightingale NMR biomarker data in UK Biobank
#'
#' This package provides utilities for working with the NMR metabolomics data
#' in \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank}.
#'
#' The \code{\link{extract_biomarkers}()} function will take a raw dataset
#' output by \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv},
#' extract the NMR biomarker fields and give them short comprehensible column
#' names. Measurements are also split into multiple rows where a participant has
#' measurements at both baseline and first repeat assessment.
#'
#' The \code{\link{compute_nightingale_ratios}()} function will compute the
#' \href{https://nightingalehealth.com/biomarkers}{81 Nightingale Health biomarker ratios}
#' from the 168 biomarkers available for download from
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank}.
#'
#' This package also provides a \code{data.frame} of biomarker information, loaded
#' as \code{\link{nmr_info}}.
#'
#' @docType package
#' @name ukbnmr
#' @import data.table
#' @import splitstackshape
#' @importFrom stats na.omit
#' @keywords package
NULL
