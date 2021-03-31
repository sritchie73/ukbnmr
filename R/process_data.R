# Helper function for processing user input data
#
# Currently expects only a data.frame loaded from the output of ukbconv, in which
# case it converts the ukb field identifiers to biomarker names, splitting x
# repeat instances and samples per eid.
process_data <- function(x) {
  # The goal here is to convert field codes like "23474-0.0" to "bOHbutyrate"
  # using the nmr_info table. Since columns are of the format
  # <field_id>-<instance>.<array_index> we need to melt to long so we can split
  # the column names into their respective parts. Here, "instance" is the
  # assessment: 0 is baseline assessment, and 1 first repeat assessment.
  # Array_index is used for repeated measures at the same time point.

  # Silence CRAN NOTES about undefined global variables (columns in data.tables)
  variable <- value <- name <- Biomarker <- UKB.Field.ID <- NULL

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

  # Melt wide to long format
  x <- melt(x, id.vars="eid")
  x <- x[!is.na(value)] # drops empty columns

  # Split what were the column names into separate fields
  x[, c("field_id", "instance", "array_index") := lapply(tstrsplit(variable, split="_"), as.integer)]

  # Map field IDs to short biomarker names to use as column names
  x[ukbnmr::nmr_info[!is.na(UKB.Field.ID)],
    on = c("field_id" = "UKB.Field.ID"),
    name := Biomarker]

  # Drop any additional fields not relating to the NMR data or participant ID
  x <- x[!is.na(name)]

  # Make sure values are numeric (may be character if other fields were present)
  if (!is.numeric(x$value))
    x[, value := as.numeric(value)]

  # Convert back to wide format
  x <- dcast(x, eid + instance + array_index ~ name, value.var="value")

  # Finished processing
  return(x)
}
