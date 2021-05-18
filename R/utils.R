# For wrapping data.table assign := statements.
#
# Basically runs a command, ignoring any errors produced.
#
# Example:
#   tryAssign(dt[, foo := bar * 2])
#
# Will create column 'foo', or do nothing when 'bar' does not exist as a column
# in 'dt'.
#
# Note, disregards all errors, so will also silently do nothing if 'bar' exists
# but cannot be multiplied (e.g. if 'bar' is class character).
#
tryAssign <- function(...) {
 tryCatch(..., error = function(e) {
   if (!(e %like% "object .* not found"))
     stop(e) # throw error if not related to column not existing
 })
}

# Cast object to data.table, or if already a data.table, make a copy so we don't
# modify by reference
getDT <- function(x) {
  if (is.data.table(x)) {
    copy(x)
  } else {
    as.data.table(x)
  }
}

# If a user is not using the data.table package (i.e. its not loaded by the user
# in their R session) make sure to return a data.frame instead.
returnDT <- function(x) {
  if (is.data.table(x) && !("data.table" %in% .packages())) {
    setDF(x)
  }
  return(x)
}

# Is a vector an integer? It may have been loaded as character / numeric
is.integer <- function(x) {
  if (is.factor(x)) {
    return(FALSE)
  }
  tryCatch({
    as.integer(x)
    return(TRUE)
  }, warning=function(w) {
    return(FALSE)
  })
}

# Paste together QC flags for biomarkers when computing ratios
collate_flags <- function(...) {
  colnames <- sapply(substitute(list(...))[-1], deparse)
  cols <- list(...)

  # Throw appropriate error if column does not exist
  for (ii in seq_along(colnames)) {
    if (is.null(cols[[ii]])) {
      stop(sprintf("object '%s' not found", colnames[ii]))
    }
  }

  tags_with_names <- sapply(seq_along(colnames), function(argIdx) {
    ifelse(is.na(cols[[argIdx]]), NA_character_,
      ifelse(grepl(":", cols[[argIdx]]), gsub(".$", cols[[argIdx]]),
        sprintf("%s: %s", colnames[argIdx], cols[[argIdx]])))
  })
  apply(tags_with_names, 1, function(row) {
    ifelse(all(is.na(row)), NA_character_, paste0(paste(stats::na.omit(row), collapse=". "), "."))
  })
}
