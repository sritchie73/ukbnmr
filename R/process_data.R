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
process_data <- function(x, type) {
  # Silence CRAN NOTES about undefined global variables (columns in data.tables)
  value <- Biomarker <- UKB.Field.ID <- QC.Flag.Field.ID <- visit_index <-
    repeat_index <- integer_rep <- flag <- variable <- NULL

  # Determine format of data
  data_format <- detect_format(x, type)

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

  # Filter to those we want to extract
  if (type == "biomarkers") {
    field_ids <- field_ids[UKB.Field.ID %in% na.omit(ukbnmr::nmr_info$UKB.Field.ID)]
  } else if (type == "biomarker_qc_flags") {
    field_ids <- field_ids[UKB.Field.ID %in% na.omit(ukbnmr::nmr_info$QC.Flag.Field.ID)]
  } else {
    stop("internal error: 'type' must be one of \"biomarkers\" or \"biomarker_qc_flags\"")
  }

  # Map to biomarker variable names
  if (type == "biomarkers") {
    field_ids[, UKB.Field.ID := as.integer(UKB.Field.ID)]
    field_ids[ukbnmr::nmr_info, on = list(UKB.Field.ID), Biomarker := Biomarker]
    field_ids[, UKB.Field.ID := as.character(UKB.Field.ID)]
  } else if (type == "biomarker_qc_flags") {
    field_ids[, UKB.Field.ID := as.integer(UKB.Field.ID)]
    field_ids[ukbnmr::nmr_info, on = list(UKB.Field.ID=QC.Flag.Field.ID), Biomarker := Biomarker]
    field_ids[, UKB.Field.ID := as.character(UKB.Field.ID)]
  }

  # Split out instance (visit) and array index (repeat measure) fields so they
  # are rows instead of columns
  visit_repeats <- setdiff(unique(gsub("^[0-9]+_", "", names(x))), "eid")
  x <- rbindlist(fill=TRUE, use.names=TRUE, lapply(visit_repeats, function(vr) {
    # Find columns matching this visit repeat pair (e.g. ending in _0_0)
    this_cols <- names(x)[grepl(pattern=paste0(vr, "$"), names(x))]

    # Filter to these columns
    this_x <- x[, .SD, .SDcols=c("eid", this_cols)]

    # Drop repeat visit pair label from column name
    setnames(this_x, this_cols, gsub(paste0("_", vr, "$"), "", this_cols))

    # Get fields present in the visit repeat pair
    this_fields <- intersect(names(this_x), field_ids$UKB.Field.ID)

    # Return empty for this visit repeat pair if no fields of interest present
    if (length(this_fields) == 0) {
      return(NULL)
    }

    # Add columns for visit and repeat index
    this_x[, visit_index := as.integer(gsub("_.*", "", vr))]
    this_x[, repeat_index := as.integer(gsub(".*_", "", vr))]

    # Filter to field IDs of interest
    this_x <- this_x[, .SD, .SDcols=c("eid", "visit_index", "repeat_index", this_fields)]

    # Rename Field IDs to biomarker variable names
    this_names <- field_ids[UKB.Field.ID %in% this_fields]
    setnames(this_x, this_names$UKB.Field.ID, this_names$Biomarker)

    # Drop instance and array index combinations with all missing data
    # eid, visit_index, and array_index always non-missing
    this_x <- this_x[apply(this_x, 1, function(row) { sum(!is.na(row)) > 3L })]

    # Return to rbindlist - which will row-bind all data.tables in the list
    # (one data.table per visit repeat pair)
    return(this_x)
  }))

  # For biomarker QC flags, map integers to flag values:
  # https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=2310
  if (type == "biomarker_qc_flags") {
    x <- melt(x, id.vars=c("eid", "visit_index", "repeat_index"))
    x <- x[!is.na(value)]

    flag_map <- data.table(
      integer_rep = 1L:10L,
      flag = c("Below limit of quantification", "Citrate plasma", "Degraded sample",
               "High ethanol", "Isopropyl alcohol", "Low glutamine or high glutamate",
               "Medium ethanol", "Polysaccharides", "Unknown contamination", "Ethanol")
    )
    x[flag_map, on = .(value=integer_rep), flag := flag]

    x <- dcast(x, eid + visit_index + repeat_index ~ variable, value.var="flag")
  }

  # Set column keys for fast joining by user
  setkeyv(x, c("eid", "visit_index", "repeat_index"))

  # Finished processing
  return(x)
}
