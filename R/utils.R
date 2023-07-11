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
  } else {
    x[,] # Otherwise the first show() call will be empty
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

  # If any of the input biomarkers have QC flags, aggregate, and add the respective
  # biomarker name. If no tags, return NA. If any column is itself an aggregated
  # set of flags (e.g. a ratio of Total_C / Total_L will be aggregated from the
  # sums of all cholesterols and lipids), don't add the biomarker name.
  tags_with_names <- sapply(seq_along(colnames), function(argIdx) {
    ifelse(is.na(cols[[argIdx]]), NA_character_,
      ifelse(grepl(":", cols[[argIdx]]), gsub(".$", "", cols[[argIdx]]),
        sprintf("%s: %s", colnames[argIdx], cols[[argIdx]])))
  })
  # For each row, get a single string of collated flags.
  collated_tags <- apply(tags_with_names, 1, function(row) {
    ifelse(all(is.na(row)), NA_character_, paste0(paste(stats::na.omit(row), collapse=". "), "."))
  })
  # In some cases, biomarkers may appear multiple times where they contribute to
  # e.g. both a numerator and denominator of a ratio. Split out the string and
  # make sure we only return each underlying biomarker's flags once.
  sapply(strsplit(collated_tags, split="\\. ?"), function(row) {
    ifelse(all(is.na(row)), NA_character_, paste0(paste(sort(unique(row)), collapse=". "), "."))
  })
}

# Convert an ITime object to numeric representation of hours since midnight
hours_decimal <- function(time) {
  hour(time) + minute(time)/60 + second(time)/3600
}

# Compute durations between events (hours)
duration_hours <- function(days1, time1, days2, time2) {
  (days2 - days1)*24 + (hours_decimal(time2) - hours_decimal(time1))
}

# Function for converting a set of numbers to a factor, ordering levels by group size
factor_by_size <- function(x) {
  # Suppress CRAN warnings for data.table columns
  N <- NULL

  group_sizes <- as.data.table(table(x))
  group_sizes <- group_sizes[order(-N)]
  factor(x, levels=group_sizes$x)
}

# Split a series of dates into equal size bins, ordered by date
bin_dates <- function(date, version) {
  # Suppress CRAN warnings for data.table columns
  bin <- NULL

  date_as_int <- as.integer(date)
  date_order <- as.integer(factor(date_as_int))
  if (version == 1L) {
    bin_size <- max(date_order)/10
    return(as.integer(ceiling(date_order/bin_size)))
  } else if (version == 2L) {
    n_bins <- floor(length(date)/2000)
    if (n_bins < 2) n_bins <- 2 # test_data not big enough to bin into groups of 2,000 samples, so just bin into 2
    bins <- cut(unique(date_order), n_bins, labels=FALSE)
    bin_map <- data.table(date=unique(date_order), bin=bins)
    return(bin_map[data.table(date=date_order), on=list(date), bin])
  }
}
