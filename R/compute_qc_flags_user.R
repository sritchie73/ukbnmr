#' Aggregate QC Flags when computing Nightingale's 81 biomarker ratios
#'
#' For the 81 biomarker ratios computed by \code{\link{compute_nightingale_ratios}()},
#' aggregates the biomarker QC flags from the biomarkers composing each ratio
#' (see \code{\link{nmr_info}}), which can be useful for determining the reason
#' underlying missing values in the biomarker ratios.
#'
#' @details
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' Biomarker QC Flags in the input data are also returned alongside those
#' aggregated by this function for the computed biomarker ratios.
#'
#' @param x \code{data.frame} containing NMR metabolomics data from UK Biobank.
#'   May either be raw field data output by
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   or data with column names corresponding to biomarkers listed in \code{\link{nmr_info}}.
#'
#' @return a \code{data.frame} with QC flags aggregated for all computed
#'         biomarker ratios.
#'
#' @seealso \code{\link{nmr_info}} for list of computed biomarker ratios and
#'   \code{\link{extract_biomarkers}()} for details on how raw data from
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   is processed.
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' biomarker_qc_flags <- compute_nightingale_ratio_qc_flags(ukb_data)
#'
#' @export
compute_nightingale_ratio_qc_flags <- function(x) {
  # Process data to correct format
  x <- process_data(x, type="biomarker_qc_flags") # copy of x created if already in right format

  # compute ratios
  x <- nightingale_ratio_flags(x)

  # Return
  returnDT(x)
}

#' Aggregate QC Flags when computing the extended set of biomarker ratios
#'
#' For the 76 biomarker ratios computed by \code{\link{compute_extended_ratios}()},
#' aggregates the biomarker QC flags from the biomarkers composing each ratio
#' (see \code{\link{nmr_info}}), which can be useful for determining the reason
#' underlying missing values in the biomarker ratios.
#'
#' @details
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' Biomarker QC Flags in the input data are also returned alongside those
#' aggregated by this function for the computed biomarker ratios.
#'
#' @param x \code{data.frame} containing NMR metabolomics data from UK Biobank.
#'   May either be raw field data output by
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   or data with column names corresponding to biomarkers listed in \code{\link{nmr_info}}.
#'
#' @return a \code{data.frame} with QC flags aggregated for all computed
#'         biomarker ratios.
#'
#' @seealso \code{\link{nmr_info}} for list of computed biomarker ratios and
#'   \code{\link{extract_biomarkers}()} for details on how raw data from
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   is processed.
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' biomarker_qc_flags <- compute_extended_ratio_qc_flags(ukb_data)
#'
#' @export
compute_extended_ratio_qc_flags <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  # Process data to correct format
  x <- process_data(x, type="biomarker_qc_flags") # copy of x created if already in right format

  # Are the nightingale ratios already computed? If not, also compute these
  Type <- Nightingale <- Biomarker <- NULL # Silence CRAN NOTES for data.table columns
  nightingale_ratios <- ukbnmr::nmr_info[
    (Nightingale) & !(Type %in% c("Non-derived", "Composite")),
    Biomarker]
  has_nightingale_ratios <- (length(intersect(names(x), nightingale_ratios)) > 0)

  # compute ratios
  if (!has_nightingale_ratios)
    x <- nightingale_ratio_flags(x)

  x <- extended_ratios_flags(x)

  # Return
  returnDT(x)
}

#' Aggregate QC Flags when recomputing all composite and derived biomarkers
#'
#' For the 61 composite biomarkers, 81 Nightingale biomarker ratios, and 76
#' extended biomarker ratios computed by \code{\link{recompute_derived_biomarkers}()},
#' aggregates the biomarker QC flags from the underlying biomarkers
#' (see \code{\link{nmr_info}}).
#'
#' @details
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' Biomarker QC Flags in the input data are also returned alongside those
#' aggregated by this function for the computed biomarker ratios.
#'
#' @param x \code{data.frame} containing NMR metabolomics data from UK Biobank.
#'   May either be raw field data output by
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   or data with column names corresponding to biomarkers listed in \code{\link{nmr_info}}.
#'
#' @return a \code{data.frame} with QC flags aggregated for all computed
#'         biomarkers and ratios.
#'
#' @seealso \code{\link{nmr_info}} for list of computed biomarker ratios and
#'   \code{\link{extract_biomarkers}()} for details on how raw data from
#'   \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#'   is processed.
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' biomarker_qc_flags <- recompute_derived_biomarker_qc_flags(ukb_data)
#'
#' @export
recompute_derived_biomarker_qc_flags <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  # Process data to correct format
  x <- process_data(x, type="biomarker_qc_flags") # copy of x created if already in right format

  # Composite biomarkers *must* be recomputed before downstream ratios
  x <- nightingale_composite_biomarker_flags(x)
  x <- nightingale_ratio_flags(x)
  x <- extended_ratios_flags(x)

  # Return
  returnDT(x)
}
