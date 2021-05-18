#' Extract Nightingale biomarkers from a data.frame of UK Biobank fields
#'
#' Given an input \code{data.frame} loaded from an dataset extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' extracts the \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank fields}
#' corresponding to the
#' \href{https://nightingalehealth.com/biomarkers}{Nightingale Health NMR metabolomics biomarkers}
#' giving them short variable names as listed in the \code{\link{nmr_info}} information data sheet
#' available in this package.
#'
#' @details
#' Data sets extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' have one row per UKB biobank participant whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' have the format "<field_id>-<visit_index>.<repeat_index>", where here <field_id>
#' corresponds to a biomarker of interest, e.g. 23474 for 3-Hydroxybutyrate,
#' <visit_index> corresponds to the assessment time point, e.g. 0 for baseline
#' assessment, 1 for first repeat visit. For the UKB NMR data, the <repeat_index>
#' column is reserved for cases where biomarker measurements have more than
#' one QC Flag (see \code{\link{extract_biomarker_qc_flags}()}).
#'
#' In the returned \code{data.frame} there is single column for each biomarker,
#' with an additional column for the visit index. Rows are uniquely identifiable
#' by the combination of entries in columns "eid" and "visit_index". There are
#' currently no repeat measure data for the NMR biomarker data in UKB,
#' so no repeat_index column is returned.
#'
#' This function will also work with data extracted by the \code{ukbtools} R
#' package.
#'
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' A \code{data.table} will be returned instead of a \code{data.frame} if the
#' the user has loaded the package into their R session.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "23474-0.0", "23474-1.0", \dots, "23467-1.0".
#' @return a \code{data.frame} or \code{data.table} with column names "eid",
#'         and "visit_index", followed by columns for each biomarker
#'        e.g. "bOHbutyrate", \dots, "Valine".
#'
#' @export
extract_biomarkers <- function(x) {
  if (detect_format(x, type="biomarkers") == "processed") {
    return(x) # return as is; process_data() would make a copy first
  } else {
    returnDT(process_data(x, type="biomarkers"))
  }
}

#' Extract QC Flags for Nightingale biomarkers from a data.frame of UK Biobank fields
#'
#' Given an input \code{data.frame} loaded from an dataset extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' extracts the \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221}{UK Biobank fields}
#' corresponding to the
#' \href{https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/nmrm_app4.pdf}{Nightingale Health NMR metabolomics biomarker QC Flags}
#' giving them short variable names as listed in the \code{\link{nmr_info}} information data sheet
#' available in this package. QC Flags are separated by "; " in each column where
#' there are multiple QC Flags for a single measurement.
#'
#' @details
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' have one row per UKB biobank participant whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' have the format "<field_id>-<visit_index>.<repeat_index>", where here <field_id>
#' corresponds to a biomarker of interest, e.g. 23774 for the QC Flags for
#' 3-Hydroxybutyrate, and <visit_index> corresponds to the assessment time point,
#' e.g. 0 for baseline assessment, 1 for first repeat visit.
#'
#' The <repeat_index> field is currently used by UK Biobank to index cases where
#' there are multiple QC Flags for a biomarker for a single participant and
#' visit_index. These are collated into a single entry in the output from this
#' function, with distinct QC Flags separated by commas in the output.
#'
#' In the returned \code{data.frame} there is single column for each biomarker,
#' with an additional column for the visit index. Rows are uniquely identifiable
#' by the combination of entries in columns "eid" and "visit_index".
#'
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "23774-0.0", "23774-1.0", \dots, "23767-1.0".
#' @return a \code{data.frame} or \code{data.table} with column names "eid",
#'        and "vist_index" followed by columns for each biomarker
#'        e.g. "bOHbutyrate", \dots, "Valine".
#'
#' @export
extract_biomarker_qc_flags <- function(x) {
  if (detect_format(x, type="biomarker_qc_flags") == "processed") {
    return(x) # return as is; process_data() would make a copy first
  } else {
    returnDT(process_data(x, type="biomarker_qc_flags"))
  }
}

#' Extract Nightingale sample QC Flags from a data.frame of UK Biobank fields
#'
#' Given an input \code{data.frame} loaded from an dataset extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' extracts the \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222}{UK Biobank fields}
#' corresponding to the
#' \href{https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/nmrm_app4.pdf}{Nightingale Health NMR metabolomics sample QC Flags}
#' giving them short variable names.
#'
#' @details
#' Data sets extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' have one row per UKB biobank participant whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' have the format "<field_id>-<visit_index>.<repeat_index>", where here <field_id>
#' corresponds to a sample QC flag, and <visit_index> corresponds to the assessment
#' time point, e.g. 0 for baseline assessment, 1 for first repeat visit. For the
#' UKB NMR data, the <repeat_index> column is reserved for cases where biomarker
#' measurements have more than one QC Flag (see \code{\link{extract_biomarker_qc_flags}()}).
#'
#' In the returned \code{data.frame} there is single column for each QC Flag,
#' with an additional column for the visit index. Rows are uniquely identifiable
#' by the combination of entries in columns "eid" and "visit_index". There are
#' currently no repeat measure data for the NMR biomarker data in UKB,
#' so no repeat_index column is returned.
#'
#' This function will also work with data extracted by the \code{ukbtools} R
#' package.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "23649-0.0", "23649-1.0", \dots, "23655-1.0".
#' @return a \code{data.frame} or \code{data.table} with column names "eid"
#'        and "visit_index", followed by columns for each sample
#'        QC tag, e.g. "Shipment.Plate", \dots, "Low.Protein".
#'
#' @export
extract_sample_qc_flags <- function(x) {
  if (detect_format(x, type="sample_qc_flags") == "processed") {
    return(x) # return as is; process_data() would make a copy first
  } else {
    returnDT(process_data(x, type="sample_qc_flags"))
  }
}

