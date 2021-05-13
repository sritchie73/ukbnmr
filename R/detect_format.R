detect_format <- function(x, type) {
  # Silence CRAN NOTES about undefined global variables (columns in data.tables)
  UKB.Field.ID <- QC.Flag.Field.ID <- NULL

  if (type == "biomarkers") {
    field_ids <- ukbnmr::nmr_info[, na.omit(UKB.Field.ID)]
  } else if (type == "biomarker_qc_flags") {
    field_ids <- ukbnmr::nmr_info[, na.omit(QC.Flag.Field.ID)]
  }

  raw_ukb_cols <- c(paste0(field_ids, "-0.0"), paste0(field_ids, "-1.0"))
  ukbtools_cols <- c(paste0(field_ids, "_0_0"), paste0(field_ids, "_1_0"))
  if (!is.data.frame(x)) {
    stop("Input data must be a data.frame or data.table")
  }
  if (length(intersect(names(x), ukbnmr::nmr_info[["Biomarker"]])) > 0) {
    return("processed")
  } else if (length(intersect(names(x), raw_ukb_cols)) > 0) {
    return("raw")
  } else if (length(intersect(gsub(".*_f", "", names(x)), ukbtools_cols)) > 0) {
    return("ukbtools")
  } else {
    if (type == "biomarkers") {
      stop("NMR biomarker fields (#23400-23578) not found in input data")
    } else if (type == "biomarker_qc_flags") {
      stop("NMR biomarker QC Flag fields (#23700-23878) not found in input data")
    } else {
      stop("internal error: 'type' must be one of \"biomarkers\" or \"biomarker_qc_flags\"")
    }
  }
}
