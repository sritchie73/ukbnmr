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
#' assessment, 1 for first repeat visit, and <repeat_index> gives a number for
#' repeated measurements at the same time point.
#'
#' In the returned \code{data.frame} there is single column for each biomarker,
#' with additional columns for the instance and array index. Rows are uniquely
#' identifiable by the combination of entries in columns "eid" and "visit_index".
#' Currently there are no repeat measurements in the Nightingale data on
#' UK Biobank.
#'
#' This function will also work with data extracted by the \code{ukbtools} R
#' package.
#'
#' A \code{data.table} will be returned instead of a \code{data.frame} if the
#' the user has loaded the package into their R session.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "23474-0.0", "23474-1.0", \dots, "23467-1.0".
#' @return a \code{data.frame} or \code{data.table} with column names "eid",
#'        "instance", and "array_index", followed by columns for each biomarker
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
#' available in this package.
#'
#' @details
#' Data sets extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' have one row per UKB biobank participant whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' have the format "<field_id>-<visit_index>.<repeat_index>", where here <field_id>
#' corresponds to a QC Flag for the biomarker of interest, e.g. 23774 for the
#' QC Flag for 3-Hydroxybutyrate. <visit_index> corresponds to the assessment
#' time point, e.g. 0 for baseline assessment, 1 for first repeat visit, and
#' <repeat_index> gives a number for repeated measurements at the same time point.
#'
#' In the returned \code{data.frame} there is single column for the QC Flags for
#' each biomarker, with additional columns for the instance and array index.
#' Rows are uniquely identifiable by the combination of entries in columns "eid" and "visit_index".
#' Currently there are no repeat measurements in the Nightingale data on
#' UK Biobank.
#'
#' This function will also work with data extracted by the \code{ukbtools} R
#' package.
#'
#' A \code{data.table} will be returned instead of a \code{data.frame} if the
#' the user has loaded the package into their R session.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "23774-0.0", "23774-1.0", \dots, "23767-1.0".
#' @return a \code{data.frame} or \code{data.table} with column names "eid",
#'        "instance", and "array_index", followed by columns for each biomarker
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
