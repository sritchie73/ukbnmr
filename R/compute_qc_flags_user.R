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
#' \dontrun{
#' # Subsequent to using the compute_biomarker_ratios function you may wish to
#' # collate into a single field for each biomarker and sample any and all
#' # quality control flags from the biomarkers underlying the ratio for that
#' # sample.
#'
#' # First, if we haven't already done so, we can compute the 81 biomarker
#' # ratios missing prior to the August 2021 update of the UK Biobank NMR
#' # metabolomics biomarker data (see help("compute_nightingale_ratios"))
#' ukb_data <- fread("path/to/extracted.csv")
#' nmr <- compute_nightingale_ratios(ukb_data)
#' fwrite(nmr, file="path/to/nmr_biomarker_data.csv")
#'
#' # Now we can aggregate the QC flags, adding them as additional columns to the
#' # input data
#' biomarker_qc_flags <- compute_nightingale_ratio_qc_flags(nmr)
#' fwrite(biomarker_qc_flags, file="path/to/nmr_biomarker_qc_flags.csv")
#' }
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
#' \dontrun{
#' # Subsequent to using the compute_extended_ratios function you may wish to
#' # collate into a single field for each biomarker and sample any and all
#' # quality control flags from the biomarkers underlying the ratio for that
#' # sample.
#'
#' # First, if we haven't already done so, we can compute the 76 additional
#' # biomarker ratios that are not part of the UK Biobank data release using the
#' # compute_extended_ratios function (see help("compute_extended_ratios"))
#' ukb_data <- fread("path/to/extracted.csv")
#' nmr <- compute_extended_ratios(ukb_data)
#' fwrite(nmr, file="path/to/nmr_biomarker_data.csv")
#'
#' # Now we can aggregate the QC flags, adding them as additional columns to the
#' # input data
#' biomarker_qc_flags <- compute_extended_ratio_qc_flags(nmr)
#' fwrite(biomarker_qc_flags, file="path/to/biomarker_qc_flags.csv")
#' }
#'
#' @export
compute_extended_ratio_qc_flags <- function(x) {
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
#' \dontrun{
#' # Subsequent to using the recompute_derived_biomarkers function you may wish
#' # to collate into a single field for each biomarker and sample any and all
#' # quality control flags from the biomarkers underlying each ratio or derived
#' # biomarker.
#'
#' # If we haven't already done so, we can follow the example workflow outlined
#' # in help("recompute_derived_biomarkers"):
#'
#' # First, if we haven't corrected for unwanted technical variation we do so
#' # using the appropriate function (see help("remove_technical_variation")).
#' ukb_data <- fread("path/to/decoded_ukbiobank_data.csv")
#' processed <- remove_technical_variation(ukb_data)
#' tech_qc <- processed$biomarkers
#' fwrite(tech_qc, file="path/to/nmr_biomarker_data.csv")
#' fwrite(processed$biomarker_qc_flags, file="path/to/nmr_biomarker_qc_flags.csv")
#' fwrite(processed$sample_processing, file="path/to/nmr_sample_qc_flags.csv")
#' fwrite(processed$log_offset, file="path/to/nmr_biomarker_log_offset.csv")
#' fwrite(processed$outlier_plate_detection, file="path/to/outlier_plate_info.csv")
#'
#' # Otherwise assuming we load 'tech_qc' from "path/to/mr_biomarker_data.csv".
#'
#' # We now run code to adjust biomarkers for biological covariates. This code is
#' # not supplied by this package, but for illustrative purposes we assume the user
#' # has written a function to do this:
#' bio_qc <- user_function_to_adjust_biomarkers_for_covariates(tech_qc)
#'
#' # Now we recompute the composite biomarkers and derived ratios after
#' # adjustment for additional biological covariates
#' bio_qc <- recompute_derived_biomarkers(bio_qc)
#' fwrite(bio_qc, file="path/to/nmr_biomarkers_adjusted_for_covariates.csv")
#'
#' # Now we can aggregate and save the quality control flags for each sample
#' # from the biomarkers underlying each derived biomarker or ratio, adding
#' # them as additional columns to the input data.
#' biomarker_qc_flags <- recompute_derived_biomarker_qc_flags(nmr)
#' fwrite(biomarker_qc_flags, file="path/to/biomarker_qc_flags.csv")
#' }
#'
#' @export
recompute_derived_biomarker_qc_flags <- function(x) {
  # Process data to correct format
  x <- process_data(x, type="biomarker_qc_flags") # copy of x created if already in right format

  # Composite biomarkers *must* be recomputed before downstream ratios
  x <- nightingale_composite_biomarker_flags(x)
  x <- nightingale_ratio_flags(x)
  x <- extended_ratios_flags(x)

  # Return
  returnDT(x)
}
