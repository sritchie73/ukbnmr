#' Tools for processing Nightingale NMR biomarker data in UK Biobank
#'
#'  @description
#' This package provides utilities for working with the NMR metabolomics data
#' in \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank}.
#'
#' @details
#' There are two groups of functions in this package: (1) data extraction functions,
#' and (2) methods for computing biomarker ratios.
#'
#' This package also provides a \code{data.frame} of biomarker information, loaded
#' as \code{\link{nmr_info}}.
#'
#' @section Data Extraction Functions:
#' The \code{\link{extract_biomarkers}()} function will take a raw dataset output
#' by \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv},
#' extract the \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{NMR biomarker fields}
#' and give them short comprehensible column names. Measurements are also split
#' into multiple rows where a participant has measurements at both baseline and first repeat assessment.
#'
#' The \code{\link{extract_biomarker_qc_flags}()} function will take a raw dataset output
#' by \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv},
#' extract the \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221}{Nightingale quality control flags}
#' for each biomarker measurement, returning a single column per biomarker
#' (corresponding to respective columns output by \code{\link{extract_biomarkers}()}).
#'
#' The \code{\link{extract_sample_qc_flags}()}function will take a raw dataset output
#' by \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv},
#' and extract the \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222}{sample quality control tags}
#' for the Nightingale NMR metabolomics data.
#'
#'
#' @section  Methods for computing biomarker ratios:
#' The \code{\link{compute_nightingale_ratios}()} function will compute the
#' \href{https://nightingalehealth.com/biomarkers}{81 Nightingale Health biomarker ratios}
#' from the 168 biomarkers currently available for download from
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank}.
#' A companion function, \code{\link{compute_nightingale_ratio_qc_flags}()} will
#' aggregate the QC flags for the biomarkers underlying each ratio.
#'
#' @docType package
#' @name ukbnmr
#' @import data.table
#' @importFrom stats na.omit
#' @keywords package
NULL
