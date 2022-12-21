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
#' \dontrun{
#' # Workflow for computing the 81 biomarker ratios missing prior to the August
#' # 2021 update of the UK Biobank NMR metabolomics biomarker data.
#' ukb_data <- fread("path/to/extracted.csv")
#' nmr <- compute_nightingale_ratios(ukb_data)
#' fwrite(nmr, file="path/to/nmr_biomarker_data.csv")
#'
#' # You may also want to aggregate and save the quality control flags for each
#' # sample from the biomarkers underlying each ratio, adding them as additional
#' # columns to the input data (see help("compute_nightingale_ratio_qc_flags")).
#' biomarker_qc_flags <- compute_nightingale_ratio_qc_flags(nmr)
#' fwrite(biomarker_qc_flags, file="path/to/nmr_biomarker_qc_flags.csv")
#' }
#'
#' @export
compute_nightingale_ratios <- function(x) {
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
#' \dontrun{
#' # Example workflow for computing the extended set of ratios added by this
#' # package to the 249 biomarkers and ratios in the UK Biobank dataset. Note
#' # no correction for unwanted technical variation takes place here, this is
#' # only useful if you want to compare their values before and after removal
#' # of unwanted technical variation.
#' ukb_data <- fread("path/to/decoded_ukbiobank_data.csv")
#' nmr <- compute_extended_ratios(ukb_data)
#' fwrite(nmr, file="path/to/nmr_biomarker_data.csv")
#'
#' # You may also want to aggregate and save the quality control flags for each
#' # sample from the biomarkers underlying each ratio, adding them as additional
#' # columns to the input data (see help("compute_nightingale_ratio_qc_flags")).
#' biomarker_qc_flags <- compute_extended_ratio_qc_flags(nmr)
#' fwrite(biomarker_qc_flags, file="path/to/biomarker_qc_flags.csv")
#' }
#'
#'
#' @export
compute_extended_ratios <- function(x) {
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
#' \dontrun{
#' # Example of a typical workflow for adjusting NMR metabolomics biomarkers for
#' # biological covariates
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
#' # Otherwise assuming we load 'tech_qc' from "path/to/mr_biomarker_data.csv",
#' # We now run code to adjust biomarkers for biological covariates (code not
#' # supplied by this package)
#' bio_qc <- user_function_to_adjust_biomarkers_for_covariates(tech_qc)
#'
#' # Now we recompute the composite biomarkers and derived ratios after
#' # adjustment for additional biological covariates
#' bio_qc <- recompute_derived_biomarkers(bio_qc)
#' fwrite(bio_qc, file="path/to/nmr_biomarkers_adjusted_for_covariates.csv")
#'
#' # You may also want to aggregate and save the quality control flags for each
#' # sample from the biomarkers underlying each derived biomarker or ratio,
#' # adding them as additional columns to the input data (see
#' # help("recompute_derived_biomarker_qc_flags")).
#' biomarker_qc_flags <- recompute_derived_biomarker_qc_flags(nmr)
#' fwrite(biomarker_qc_flags, file="path/to/biomarker_qc_flags.csv")
#' }
#'
#' @export
recompute_derived_biomarkers <- function(x) {
  # Process data to correct format
  x <- process_data(x, type="biomarkers") # copy of x created if already in right format

  # Composite biomarkers *must* be recomputed before downstream ratios
  x <- nightingale_composite_biomarker_compute(x)
  x <- nightingale_ratio_compute(x)
  x <- extended_ratios_compute(x)

  # Return
  returnDT(x)
}
