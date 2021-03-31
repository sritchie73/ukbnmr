#' Tools for processing Nightingale NMR biomarker data in UK Biobank
#'
#' This package provides utilities for working with the NMR metabolomics data
#' in \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank}.
#'
#' The \code{\link{extract_biomarkers}()} function will take a raw dataset
#' output by \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv},
#' extract the NMR biomarker fields and give them short comprehensible column
#' names. Measurements are also split into multiple rows where a participant has
#' measurements at both baseline and first repeat assessment. All functions below
#' will also work with raw data extracted from UK Biobank.
#'
#' The \code{\link{compute_nightingale_ratios}()} function will compute the
#' \href{https://nightingalehealth.com/biomarkers}{81 Nightingale Health biomarker ratios}
#' that were missing prior to the May update of
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank}.
#'
#' The \code{\link{compute_extended_ratios}()} function will compute an extended
#' set of biomarker ratios expanding on the biomarkers available directly from
#' the Nightingale platform.
#'
#' The \code{\link{recompute_derived_biomarkers}()} function will recompute all
#' composite biomarkers and ratios from 107 non-derived biomarkers, which is
#' useful for ensuring data consistency when adjusting for unwanted biological
#' variation.
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
