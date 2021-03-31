# Helper function for processing user input data
#
# Columns in the raw UK Biobank data are <field_id>-<instance>.<array_index>,
# where <field_id> corresponds to the biomarker, e.g. 23474 for 3-Hydroxybutyrate
# (given the variable name bOHbutyrate), <instance> corresponds to the timepoint
# of biomarker quantification, 0 for baseline assessment, 1 for first repeat
# assessment, and <array_index> corresponds to a number 0--N for repeated
# measures at the same time point (always present, all 0 if only 1 measure
# taken).
#
# This function splits out these columns into multiple rows, 1 per <instance>
# and <array_index> per participant ("eid"), and maps the field identifiers to
# the short biomarker variable names typically provided by Nightingale Health,
# listed in the Biomarker column in the nmr_info data sheet included with this
# package.
process_data <- function(x) {
  # Silence CRAN NOTES about undefined global variables (columns in data.tables)
  value <- Biomarker <- UKB.Field.ID <- NULL

  # Determine format of data
  data_format <- detect_format(x)

  # Make sure its a data.table
  x <- getDT(x)

  # If already in correct format, no need to do anything extra, return.
  if (data_format == "processed")
    return(x)

  # Give field ids consistent names based on data format
  if (data_format == "raw") {
    setnames(x, gsub("-|\\.", "_", names(x)))
  } else if (data_format == "ukbtools") {
    setnames(x, gsub(".*_f", "", names(x)))
  }

  # Extract field IDs present
  field_ids <- data.table(UKB.Field.ID = names(x))
  field_ids[, UKB.Field.ID := gsub("_.*", "", UKB.Field.ID)]
  field_ids <- unique(field_ids)
  field_ids <- field_ids[UKB.Field.ID %in% na.omit(ukbnmr::nmr_info$UKB.Field.ID)]

  # Map to biomarker variable names
  field_ids[, UKB.Field.ID := as.integer(UKB.Field.ID)]
  field_ids[ukbnmr::nmr_info, on = list(UKB.Field.ID), Biomarker := Biomarker]
  field_ids[, UKB.Field.ID := as.character(UKB.Field.ID)]

  # Split out instance and array index fields
  setkey(x, "eid")
  x <- splitstackshape::merged.stack(x, id.vars="eid", sep="_", keep.all=FALSE,
                                     var.stubs=field_ids$UKB.Field.ID)
  setnames(x, c(".time_1", ".time_2"), c("instance", "array_index"))

  # Rename Field IDs to biomarker variable names
  setnames(x, field_ids$UKB.Field.ID, field_ids$Biomarker)

  # Drop instance and array index combinations with all missing data
  x <- melt(x, id.vars=c("eid", "instance", "array_index"))
  x <- x[!is.na(value)]
  x <- dcast(x, eid + instance + array_index ~ variable, value.var="value")
  setkeyv(x, c("eid", "instance", "array_index"))

  # Finished processing
  return(x)
}
