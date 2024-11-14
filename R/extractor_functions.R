#' Extract NMR metabolomic biomarkers from a data.frame of UK Biobank fields
#'
#' Given an input \code{data.frame} loaded from a dasaset of
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{NMR metabolomics fields}
#' extracted by the Table Exporter tool on the
#' \href{https://ukbiobank.dnanexus.com/landing}{UK Biobank Research Analysis Platform},
#' this function extracts the NMR metabolomics biomarkers giving them short variable
#' names as listed in the \code{\link{nmr_info}} information data sheet
#' available in this package.
#'
#' @details
#' Data sets extracted on the \href{https://ukbiobank.dnanexus.com/landing}{UK Biobank Research Analysis Platform}
#' have one row per UK Biobank participant, whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' follow a \href{https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-data/accessing-phenotypic-data#database-columns}{naming scheme}
#' based on the unique identifier of each field, assessment visit, and (optionally if relevant)
#' repeated measurement of "p<field_id>_i<visit_index>_a<repeat_index>". For example,
#' the measurement of 3-Hydroxybutyrate at baseline assessment has the column name
#' "p23474_i0". For the UKB NMR data, measurements are available at baseline assessment
#' and the first repeat assessment (e.g. "p23474_i1"). For the UKB NMR data, the
#' <repeat_index> is reserved for cases where biomarker measurements have more than
#' one QC Flag (see \code{\link{extract_biomarker_qc_flags}()}).
#'
#' The \code{data.frame} returned by this function gives each field a unique
#' recognizable name, with measurements from baseline and repeat assessment
#' given in separate rows. The "visit_index" column immediately after the "eid"
#' column indicates whether the biomarker measurement was quantified from the
#' blood samples taken at baseline assessment (visit_index == 0) or first repeat
#' assessment (visit_index == 1). Rows are uniquely identifiable
#' by the combination of entries in columns "eid" and "visit_index".
#'
#' This function will also work with data predating the Research Analysis Platform,
#' including data sets extracted by the \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' tool and/or the \code{ukbtools} R package.
#'
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' A \code{data.table} will be returned instead of a \code{data.frame} if the
#' the user has loaded the package into their R session.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "p23474_i0", "p23474_i1", \dots, "p23467_i1".
#' @return a \code{data.frame} or \code{data.table} with column names "eid",
#'         and "visit_index", followed by columns for each biomarker
#'        e.g. "bOHbutyrate", \dots, "Valine".
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' nmr <- extract_biomarkers(ukb_data)
#'
#' @export
extract_biomarkers <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  if (detect_format(x, type="biomarkers") == "processed") {
    return(x) # return as is; process_data() would make a copy first
  } else {
    returnDT(process_data(x, type="biomarkers"))
  }
}

#' Extract NMR biomarker QC flags from a data.frame of UK Biobank fields
#'
#' Given an input \code{data.frame} loaded from a dasaset of
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221}{NMR metabolomics QC indicator fields}
#' extracted by the Table Exporter tool on the
#' \href{https://ukbiobank.dnanexus.com/landing}{UK Biobank Research Analysis Platform},
#' this function extracts the
#' \href{https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/nmrm_app4.pdf}{quality control (QC) flags for the NMR metaolomics biomarkers}
#' giving them short variable names as listed in the \code{\link{nmr_info}} information data sheet
#' available in this package. QC Flags are separated by "; " in each column where
#' there are multiple QC Flags for a single measurement.
#'
#' @details
#' #' Data sets extracted on the \href{https://ukbiobank.dnanexus.com/landing}{UK Biobank Research Analysis Platform}
#' have one row per UK Biobank participant, whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' follow a \href{https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-data/accessing-phenotypic-data#database-columns}{naming scheme}
#' based on the unique identifier of each field, assessment visit, and (optionally if relevant)
#' repeated measurement of "p<field_id>_i<visit_index>_a<repeat_index>". For example,
#' the QC flags for the measurement of 3-Hydroxybutyrate at baseline assessment
#' has the column names "p23774_i0_a0", "p23774_i0_a1", and "p23774_i0_a2"; indicating
#' that the 3-Hydroxybutyrate measurement can have up to three QC flags per
#' sample at baseline assessment. Measurements for blood samples collected at the
#' first repeat assessment have 1 in the visit index, e.g. for 3-Hydroxybutyrate
#' at the first repeat assessment there are three columns "p23774_i1_a0",
#' "p23774_i1_a1", "p23774_i1_a2".
#'
#' The \code{data.frame} returned by this function gives each field a unique
#' recognizable name, with measurements from baseline and repeat assessment
#' given in separate rows. The "visit_index" column immediately after the "eid"
#' column indicates whether the biomarker measurement was quantified from the
#' blood samples taken at baseline assessment (visit_index == 0) or first repeat
#' assessment (visit_index == 1). Where multiple QC flags were present at the
#' same measurement, these are collated into a single entry with the multiple QC
#' flags separated by a "; ". Rows are uniquely identifiable by the combination
#' of entries in columns "eid" and "visit_index".
#'
#' This function will also work with data predating the Research Analysis Platform,
#' including data sets extracted by the \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' tool and/or the \code{ukbtools} R package.
#'
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' A \code{data.table} will be returned instead of a \code{data.frame} if the
#' the user has loaded the package into their R session.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "p23774_i0_a0", \dots "p23774_i1_a0", \dots, "p23767_i1".
#' @return a \code{data.frame} or \code{data.table} with column names "eid",
#'        and "visit_index" followed by columns for each biomarker
#'        e.g. "bOHbutyrate", \dots, "Valine".
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' biomarker_qc_flags <- extract_biomarker_qc_flags(ukb_data)
#'
#' @export
extract_biomarker_qc_flags <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  if (detect_format(x, type="biomarker_qc_flags") == "processed") {
    return(x) # return as is; process_data() would make a copy first
  } else {
    returnDT(process_data(x, type="biomarker_qc_flags"))
  }
}

#' Extract NMR sample QC flags from a data.frame of UK Biobank fields
#'
#' Given an input \code{data.frame} loaded from a dasaset of
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222}{NMR metabolomics processing fields}
#' extracted by the Table Exporter tool on the
#' \href{https://ukbiobank.dnanexus.com/landing}{UK Biobank Research Analysis Platform}, this function
#' extracts the
#' \href{https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/nmrm_app4.pdf}{sample quality control flags for the NMR metabolomics biomarker data} giving them short variable names as listed in the
#' \code{\link{sample_qc_info}} information data sheet available in this package.
#'
#' @details
#' Data sets extracted on the \href{https://ukbiobank.dnanexus.com/landing}{UK Biobank Research Analysis Platform}
#' have one row per UK Biobank participant, whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' follow a \href{https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-data/accessing-phenotypic-data#database-columns}{naming scheme}
#' based on the unique identifier of each field, assessment visit, and (optionally if relevant)
#' repeated measurement of "p<field_id>_i<visit_index>_a<repeat_index>". For example,
#' the Shipment Plate for each sample collected at baseline assessment has the
#' column name "p23649_i0". For the UKB NMR data, measurements are available at
#' baseline assessment and the first repeat assessment (e.g. "p23649_i1"). For the
#' UKB NMR data, the <repeat_index> is reserved for cases where biomarker
#' measurements have more than one QC Flag (see \code{\link{extract_biomarker_qc_flags}()}).
#'
#' The \code{data.frame} returned by this function gives each field a unique
#' recognizable name, with measurements from baseline and repeat assessment
#' given in separate rows. The "visit_index" column immediately after the "eid"
#' column indicates whether the biomarker measurement was quantified from the
#' blood samples taken at baseline assessment (visit_index == 0) or first repeat
#' assessment (visit_index == 1). Rows are uniquely identifiable
#' by the combination of entries in columns "eid" and "visit_index".
#'
#' This function will also work with data predating the Research Analysis Platform,
#' including data sets extracted by the \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' tool and/or the \code{ukbtools} R package.
#'
#' If your UK Biobank project only has access to a subset of biomarkers, then
#' this function will only return the subset of ratios that can be computed from
#' the biomarker data provided.
#'
#' A \code{data.table} will be returned instead of a \code{data.frame} if the
#' the user has loaded the package into their R session.
#'
#' @param x \code{data.frame} with column names "eid" followed by extracted
#'           fields e.g. "p23649_i0", "p23649_i1", \dots, "p23655_i1".
#' @return a \code{data.frame} or \code{data.table} with column names "eid"
#'        and "visit_index", followed by columns for each sample
#'        QC tag, e.g. "Shipment.Plate", \dots, "Low.Protein".
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' sample_qc_flags <- extract_sample_qc_flags(ukb_data)
#'
#' @import bit64
#' @export
extract_sample_qc_flags <- function(x) {
  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  if (detect_format(x, type="sample_qc_flags") == "processed") {
    return(x) # return as is; process_data() would make a copy first
  } else {
    returnDT(process_data(x, type="sample_qc_flags"))
  }
}

