#' Compute 81 biomarker ratios on the Nightingale platform
#'
#' The Nightingale Health NMR metabolomics biomarker platform quantifies
#' \href{https://research.nightingalehealth.com/biomarkers/}{249 biomarkers}, including 81
#' biomarker ratios. Prior to the August 2021 update, UK Biobank only provided the 168
#' biomarkers that are not ratios for
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{download}.
#' This function will compute the 81 missing ratios.
#'
#' @details
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' All biomarkers in the input data are also returned alongside the ratios computed
#' by this function.
#'
#' @param x \code{data.frame} containing NMR metabolomics data from UK Biobank.
#'   May either be raw field data output by
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   or data with column names corresponding to biomarkers listed in \code{\link{nmr_info}}.
#'
#' @return a \code{data.frame} with the additional computed biomarker ratios.
#'
#' @seealso \code{\link{nmr_info}} for list of computed biomarker ratios,
#'   \code{\link{compute_nightingale_ratio_qc_flags}()} for obtaining an
#'   aggregate of the biomarker QC flags from the biomarkers underlying each
#'   computed ratio, and \code{\link{extract_biomarkers}()} for details on how
#'   raw data from
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   is processed.
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' nmr <- compute_nightingale_ratios(ukb_data)
#'
#' @export
compute_nightingale_ratios <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  # Process data to correct format
  x <- process_data(x, type="biomarkers") # copy of x created if already in right format

  # compute ratios
  x <- nightingale_ratio_compute(x)

  # Return
  returnDT(x)
}

#' Compute extended set of biomarker ratios
#'
#' Computes 76 additional ratios not provided by the Nightingale platform. These
#' include lipid fractions in HDL, LDL, VLDL, and total serum, as well as
#' cholesterol fractions, and omega to polyunsaturated fatty acid ratios. See
#' \code{\link{nmr_info}} for details.
#'
#' @details
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' All biomarkers in the input data are also returned alongside the ratios computed
#' by this function.
#'
#' @param x \code{data.frame} containing NMR metabolomics data from UK Biobank.
#'   May either be raw field data output by
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   or data with column names corresponding to biomarkers listed in \code{\link{nmr_info}}.
#'
#' @return a \code{data.frame} with the additional computed biomarker ratios.
#'
#' @seealso \code{\link{nmr_info}} for list of computed biomarker ratios,
#'   \code{\link{compute_extended_ratio_qc_flags}()} for obtaining an
#'   aggregate of the biomarker QC flags from the biomarkers underlying each
#'   computed ratio, and \code{\link{extract_biomarkers}()} for details on how
#'   raw data from
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   is processed.
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' nmr <- compute_extended_ratios(ukb_data)
#'
#' @export
compute_extended_ratios <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  # Process data to correct format
  x <- process_data(x, type="biomarkers") # copy of x created if already in right format

  # Are the nightingale ratios already computed? If not, also compute these
  Type <- Nightingale <- Biomarker <- NULL # Silence CRAN NOTES for data.table columns
  nightingale_ratios <- ukbnmr::nmr_info[
    (Nightingale) & !(Type %in% c("Non-derived", "Composite")),
    Biomarker]
  has_nightingale_ratios <- (length(intersect(names(x), nightingale_ratios)) > 0)

  # compute ratios
  if (!has_nightingale_ratios)
    x <- nightingale_ratio_compute(x)

  x <- extended_ratios_compute(x)

  # Return
  returnDT(x)
}

#' Recompute composite biomarkers and ratios from the 107 non-derived biomarkers
#'
#' When adjusting biomarkers for unwanted biological covariates, it is desirable
#' to recompute composite biomarkers and ratios to ensure consistency in the
#' adjusted dataset. This function will compute all composite biomarkers and
#' ratios from their parts (see \code{\link{nmr_info}} for biomarker details).
#'
#' @details
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' All biomarkers in the input data are also returned alongside those computed
#' by this function.
#'
#' @param x \code{data.frame} containing NMR metabolomics data from UK Biobank.
#'   May either be raw field data output by
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   or data with column names corresponding to biomarkers listed in \code{\link{nmr_info}}.
#'
#' @return a \code{data.frame} with all composite biomarkers and ratios
#'   (re)computed from the 107 non-derived biomarkers (see \code{\link{nmr_info}}
#'   for details).
#'
#' @seealso \code{\link{nmr_info}} for list of computed biomarker ratios,
#'   \code{\link{recompute_derived_biomarker_qc_flags}()} for obtaining an
#'   aggregate of the biomarker QC flags from the biomarkers underlying each
#'   computed biomarker, and \code{\link{extract_biomarkers}()} for details on how
#'   raw data from
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   is processed.
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' bio_qc <- recompute_derived_biomarkers(ukb_data)
#'
#' @export
recompute_derived_biomarkers <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  # Process data to correct format
  x <- process_data(x, type="biomarkers") # copy of x created if already in right format

  # Composite biomarkers *must* be recomputed before downstream ratios
  x <- nightingale_composite_biomarker_compute(x)
  x <- nightingale_ratio_compute(x)
  x <- extended_ratios_compute(x)

  # Return
  returnDT(x)
}
