#' Remove technical variation from NMR biomarker data in UK Biobank.
#'
#' @param x \code{data.frame} containing a tabular phenotype dataset extracted by the
#' \href{https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/accessing-data/accessing-phenotypic-data#create-a-tsv-or-csv-file-using-table-exporter}{Table Exporter} tool on the
#' \href{https://ukbiobank.dnanexus.com/landing}{UK Biobank Research Analysis Platform}
#' containing the project specific sample id (eid) and all fields (and instances)
#' relating to the NMR metabolomics data (i.e. fields listed in showcase categories
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{220},
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221}{221}, and
#' \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222}{222}.
#' @param remove.outlier.plates logical, when set to \code{FALSE} biomarker
#' concentrations on outlier shipment plates (see details) are not set to
#' missing but simply flagged in the \code{biomarker_qc_flags} \code{data.frame}
#' in the returned \code{list}.
#' @param skip.biomarker.qc.flags logical, when set to \code{TRUE} biomarker QC
#' flags are not processed or returned.
#' @param version version of the QC algorithm to apply. Defaults to 3, the latest
#' version of the algorithm (see details).
#'
#' @details
#' Three versions of the QC algorithm have been developed. Version 1 was designed
#' based on the first phase of data released to the public covering ~120,000
#' UK Biobank participants. Version 2 made several improvements to the algorithm
#' based on the subsequent second public release of data covering an additional
#' ~150,000 participants. Version 3, the default, makes some further minor tweaks
#' primarily so that the algorithm is compatible with the full public data
#' release covering all ~500,000 participants.
#'
#' Details on the impact of these algorithms on technical variation in the latest
#' UK Biobank data are provided in the package vignette and on the github readme.
#'
#' ## Version 1
#' Setting version to 1 applies the algorithm as exactly described in Ritchie S.
#' C. *et al.* _Sci Data_ 2023. In brief, this multi-step procedure applies the
#' following steps in sequence:
#'
#' 1. First biomarker data is filtered to the 107 biomarkers that cannot be
#'    derived from any combination of other biomarkers.
#' 2. Absolute concentrations are log transformed, with a small offset
#'    applied to biomarkers with concentrations of 0.
#' 3. Each biomarker is adjusted for the time between sample preparation
#'    and sample measurement (hours).
#' 4. Each biomarker is adjusted for systematic differences between rows
#'    (A-H) on the 96-well shipment plates.
#' 5. Each biomarker is adjusted for remaining systematic differences
#'    between columns (1-12) on the 96-well shipment plates.
#' 6. Each biomarker is adjusted for drift over time within each of the six
#'    spectrometers. To do so, samples are grouped into 10 bins, within each
#'    spectrometer, by the date the majority of samples on their respective
#'    96-well plates were measured.
#' 7. Regression residuals after the sequential adjustments are
#'    transformed back to absolute concentrations.
#' 8. Samples belonging to shipment plates that are outliers of
#'    non-biological origin are identified and set to missing.
#' 9. The 61 composite biomarkers and 81 biomarker ratios are recomputed
#'    from their adjusted parts.
#' 10. An additional 76 biomarker ratios of potential biological
#'     significance are computed.
#'
#' At each step, adjustment for technical covariates is performed using
#' \link[MASS:rlm]{robust linear regression}. Plate row, plate column, and
#' sample measurement date bin are treated as factors, using the group with the
#' largest sample size as reference in the regression.
#'
#' ## Version 2
#' Version 2 of the algorithm modifies the above procedure to (1)
#' adjust for well and column within each processing batch separately in steps 4
#' and 5, and (2) in step 6 instead of splitting samples into 10 bins per
#' spectrometer uses a fixed bin size of approximately 2,000 samples per bin,
#' ensuring samples measured on the same plate and plates measured on the same
#' day are grouped into the same bin.
#'
#' The first modification was made as applying version 1 of the algorithm to the
#' combined data from the first and second tranche of measurements revealed
#' introduced stratification by well position when examining the corrected
#' concentrations in each data release separately.
#'
#' The second modification was made to ensure consistent bin sizes across data
#' releases when correcting for drift over time. Otherwise, spectrometers used
#' in multiple data releases would have different bin sizes when adjusting
#' different releases. A bin split is also hard coded on spectrometer 5 between
#' plates 0490000006726 and 0490000006714 which correspond to a large change in
#' concentrations akin to a spectrometer recalibration event most strongly
#' observed for alanine concentrations.
#'
#' ## Version 3
#' Version 3 of the algorithm makes two further minor changes:
#'
#' 1. Imputation of missing sample measurement times has been improved. Previously,
#'    any samples missing time of measurement (N=3 in the phase 2 public release)
#'    had their time of measurement set to 00:00. In version 3, the time of measurement
#'    is set to the median time of measurement for that spectrometer on that day,
#'    which is between 12:00-13:00, instead of 00:00.
#' 2. Samples missing sample preparation dates and times (N=182 in the V20 UK
#'    Biobank data release) use their sample centrifugation date and time as a
#'    proxy to allow adjustment for time between sample prep and sample
#'    measurement. Sample centrifugation always takes place a short time after
#'    sample preparation. As sample centrifugation date and time is not made
#'    available via UK Biobank, this data and time is hard coded; all samples
#'    missing sample preparation date and time had a sample centrifugation date
#'    and time of 2022-12-20 06:39:03 in the extended advance access data, so
#'    we use this date and time for any samples missing sample preparation date
#'    and time.
#' 3. Underlying code for adjusting drift over time has been modified to accommodate
#'    the phase 3 public release, which includes one spectrometer with ~2,500 samples.
#'    Version 2 of the algorithm would split this into two bins, whereas version 3
#'    keeps this as a single bin to better match the bin sizes of the rest of the
#'    spectrometers.
#'
#' @return a \code{list} containing three \code{data.frames}: \describe{
#'   \item{biomarkers}{A \code{data.frame} with column names "eid",
#'        and "visit_index", containing project-specific sample identifier and
#'        UK Biobank visit index (0 for baseline assessment, 1 for first repeat
#'        assessment), followed by columns for each biomarker containing their
#'        absolute concentrations (or ratios thereof) adjusted for technical
#'        variation. See \code{\link{nmr_info}} for information on each biomarker.}
#'   \item{biomarker_qc_flags}{A \code{data.frame} with the same format as
#'         \code{biomarkers} with entries corresponding to the
#'         \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221}{quality
#'         control indicators} for each sample. "High plate outlier" and "Low
#'         plate outlier" indicate the value was set to missing due to systematic
#'         abnormalities in the biomarker's concentration on the sample's shipment
#'         plate compared to all other shipment plates (see Details). For
#'         composite and derived biomarkers, quality control flags are aggregates
#'         of any quality control flags for the underlying biomarkers from which
#'         the composite biomarker or ratio is derived.}
#'   \item{sample_processing}{A \code{data.frame} containing the
#'         \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222}{processing
#'         information and quality control indicators} for each sample, including
#'         those derived for removal of unwanted technical variation by this
#'         function. See \code{\link{sample_qc_info}} for details.}
#'   \item{log_offset}{A \code{data.frame} containing diagnostic information on
#'         the offset applied so that biomarkers with concentrations of 0 could
#'         be log transformed, and any right shift applied to prevent negative
#'         concentrations after rescaling adjusted residuals back to absolute
#'         concentrations. Should contain only biomarkers with minimum
#'         concentrations of 0 (in the "Minimum" column). "Minimum.Non.Zero"
#'         gives the smallest non-zero concentration for the biomarker.
#'         "Log.Offset" the small offset added to all samples prior to
#'         log transformation: half the mininum non-zero concentration.
#'         "Right.Shift" gives the small offset added to prevent negative
#'         concentrations that arise after rescaling residuals to log
#'         concentrations: this should be at least one order of magnitude
#'         smaller than the smallest non-zero value (i.e. the offset added
#'         should amount to noise in numeric precision for all samples). See
#'         publication for more details.}
#'   \item{outlier_plate_detection}{A \code{data.frame} containing diagnostic
#'         information and details of outlier plate detection. For each of the
#'         107 non-derived biomarkers, the median concentration on each of the
#'         1,352 plates was calculated, then plates were flagged as outliers if
#'         their median value deviated more than expected from the mean of plate
#'         medians. "Mean.Plate.Medians" gives the mean of the plate medians for
#'         each biomarker. "Lower.Limit" and "Upper.Limit" give the values below
#'         and above which plates are flagged as outliers based on their plate
#'         median. See publication for more details.}
#'   \item{algorithm_version}{Version of the algorithm run, either 1, 2, or 3
#'         (default).}
#' }
#'
#' @examples
#' ukb_data <- ukbnmr::test_data # Toy example dataset for testing package
#' processed <- remove_technical_variation(ukb_data)
#'
#' @importFrom lubridate ymd_hms
#' @importFrom lubridate ymd
#' @importFrom stats coef
#' @importFrom stats sd
#' @importFrom stats qnorm
#' @importFrom stats ppoints
#' @importFrom stats median
#' @importFrom MASS rlm
#'
#' @export
remove_technical_variation <- function(
  x,
  remove.outlier.plates=TRUE,
  skip.biomarker.qc.flags=FALSE,
  version=3L
) {
  # Silence CRAN NOTES about undefined global variables (columns in data.tables)
  Well.Position.Within.Plate <- Well.Row <- Well.Column <- Sample.Measured.Date.and.Time <-
    Sample.Prepared.Date.and.Time <- Sample.Measured.Date <- Sample.Measured.Time <-
    Sample.Prepared.Date <- Sample.Prepared.Time <- Prep.to.Measure.Duration <-
    Spectrometer.Date.Bin <- Spectrometer <- Type <- Biomarker <- eid <- visit_index <-
    Shipment.Plate <- value <- Log.Offset <- Minimum <- Minimum.Non.Zero <-
    log_value <- adj <- Right.Shift <- outlier <- Lower.Limit <- Upper.Limit <-
    Name <- UKB.Field.ID <- N <- Plate.Measured.Date <- i.Sample.Measured.Date <-
    Processing.Batch <- Spectrometer.Group <- max_bin <- i.offset <- i.Processing.Batch <-
    NULL

  # Check for valid algorithm version
  if (version != 1L && version != 2L && version != 3L) stop("'version' must be 1, 2, or 3")

  message("Checking for revelant UKB fields...")

  # Check relevant fields exist
  f1 <- detect_format(x, type="biomarkers")
  f2 <- if (skip.biomarker.qc.flags) { "skipped" } else { detect_format(x, type="biomarker_qc_flags") }
  f3 <- detect_format(x, type="sample_qc_flags")

  # Make sure its not a "processed" dataset returned by extract_biomarkers() and
  # related functions.
  if (f1 == "processed" || f2 == "processed" || f3 == "processed") {
    stop("'x' should be a 'data.frame' containing UK Biobank field ID columns, not one extracted by extract_biomarkers() or related functions'")
  }

  # Check if running on package test_data, and if so, force data.table to be
  # single threaded so that we can avoid a NOTE on CRAN submission due to their
  # misconfigured debian server
  if (isTRUE(all.equal(x, ukbnmr::test_data))) {
    registered_threads <- getDTthreads()
    setDTthreads(1)
    on.exit({ setDTthreads(registered_threads) }) # re-register so no unintended side effects for users
  }

  message("Extracting and pre-processing data...")

  # Extract biomarkers, biomarker QC flags, and sample processing information
  bio <- process_data(x, type="biomarkers")
  if (!skip.biomarker.qc.flags) {
    bio_qc <- process_data(x, type="biomarker_qc_flags")
  }
  sinfo <- process_data(x, type="sample_qc_flags")

  message("Checking for required sample processing fields needed for QC procedure...")

  # Check required sample processing fields exist
  req <- c("Well.Position.Within.Plate", "Sample.Measured.Date.and.Time",
           "Sample.Prepared.Date.and.Time", "Spectrometer", "Shipment.Plate")
  miss_req <- setdiff(req, names(sinfo))
  if (length(miss_req) > 0) {
    err_txt <- ukbnmr::sample_qc_info[Name %in% miss_req]
    err_txt <- err_txt[, sprintf("%s (Field: %s)", Name, UKB.Field.ID)]
    err_txt <- paste(err_txt, collapse=", ")
    stop("Missing required sample processing fields: ", err_txt, ".")
  }

  if (version > 1L && !("Processing.Batch" %in% names(sinfo))) {
    warning("Processing.Batch missing (Field 20282), inferring from Shipment.Plate (Field 23649)")
    sinfo[plate_batch_map, on = list(Shipment.Plate), Processing.Batch := i.Processing.Batch]
    if (any(is.na(sinfo$Processing.Batch))) {
      warning("Inference of Processing.Batch from Shipment.Plate failed, reverting to algorithm version 1")
      version <- 1L
    }
  }

  message("Processing sample processing fields for QC procedure...")

  # Split out row and column information from 96-well plate information
  sinfo[, Well.Row := gsub("[0-9]", "", Well.Position.Within.Plate)]
  sinfo[, Well.Column := as.integer(gsub("[A-Z]", "", Well.Position.Within.Plate))]

  # Split out date and time for sample measurement and sample prep
  sinfo[, c("Sample.Measured.Date", "Sample.Measured.Time") := IDateTime(harmonize_datetime(Sample.Measured.Date.and.Time))]
  sinfo[, c("Sample.Prepared.Date", "Sample.Prepared.Time") := IDateTime(harmonize_datetime(Sample.Prepared.Date.and.Time))]

  # Handle the handful of samples missing a sample measurement time. In the above
  # these default to 00:00:00, instead, impute these to the median time of measurement
  # for that spectrometer on that day
  if (version == 3L) {
    missing_time <- sinfo[, which(is.na(suppressWarnings(ymd_hms(Sample.Measured.Date.and.Time))))]
    for (idx in missing_time) {
      this_date <- sinfo[idx, Sample.Measured.Date]
      this_spectrometer <- sinfo[idx, Spectrometer]
      median_measure_time <- sinfo[Sample.Measured.Date == this_date & Spectrometer == this_spectrometer, median(Sample.Measured.Time)]
      sinfo[idx, Sample.Measured.Time := median_measure_time]
    }
  }

  # In the V20 data release 182 samples are missing both a sample preparation
  # date and sample preparation time, which are needed to compute hours between
  # sample prep and sample measurement. In the advance access data, we also had
  # access to date+time of sample centrifugation, which always takes place a
  # short time after sample prep, so we substituted that date+time in as a proxy
  # for sample preparation date+time to QC the data. These samples all shared a
  # sample centrifugation date and time of 2022-12-20 06:39:03, which makes it
  # possible to hard-code this value in the package to substitute in for these
  # samples, as centrifugation date+time is not released by UK Biobank.
  if (version == 3L) {
    sinfo[is.na(Sample.Prepared.Date), Sample.Prepared.Date := as.IDate("2022-12-20")]
    sinfo[is.na(Sample.Prepared.Time), Sample.Prepared.Time := as.ITime("06:39:03")]
  }

  # Compute hours between sample prep and sample measurement
  sinfo[, Prep.to.Measure.Duration := duration_hours(
    Sample.Prepared.Date, Sample.Prepared.Time,
    Sample.Measured.Date, Sample.Measured.Time
  )]

  # Get majority date of sample measurement per shipment plate
  plate_measure_dates <- sinfo[,.N,by=list(Shipment.Plate, Sample.Measured.Date)]
  plate_measure_dates <- plate_measure_dates[,.SD[which.max(N)],by=Shipment.Plate]
  sinfo[plate_measure_dates, on = list(Shipment.Plate), Plate.Measured.Date := i.Sample.Measured.Date]

  # For version 2 of the algorithm, we need to hard code a bin split. The easiest way
  # to do this is by creating a dummy spectrometer column to split on, by giving the
  # spectrometer before and after the split a different name
  sinfo[, Spectrometer.Group := Spectrometer]
  if (version > 1L && "0490000006726" %in% sinfo[["Shipment.Plate"]]) {
    sinfo <- sinfo[order(Plate.Measured.Date)][order(Spectrometer)]
    sinfo[, Spectrometer.Group := as.character(Spectrometer.Group)]
    spec_to_split <- sinfo[Shipment.Plate == "0490000006726", Spectrometer.Group][1]
    sinfo[, row := .I]
    idx_to_split <- sinfo[Shipment.Plate == "0490000006726", max(row)]
    sinfo[Spectrometer.Group == spec_to_split & row > idx_to_split,
          Spectrometer.Group := paste(Spectrometer.Group, "Group 2")]
    sinfo[, row := NULL]
  }

  # Bin samples into either 10 equal size bins per spectrometer (algorithm
  # version 1) or into bins of roughly 2000 samples each (algorithm version 2)
  sinfo[, Spectrometer.Date.Bin := bin_dates(Plate.Measured.Date, version), by=Spectrometer.Group]

  # Offset date bin by spectrometer number to make bin allocation unique across
  # all spectrometers
  offset <- sinfo[, list(max_bin=max(Spectrometer.Date.Bin)), by=Spectrometer.Group]
  offset <- offset[order(Spectrometer.Group)]
  offset[-1, offset := offset[, cumsum(max_bin)[-.N]]]
  offset[1, offset := 0]
  sinfo[offset, on = list(Spectrometer.Group), Spectrometer.Date.Bin := Spectrometer.Date.Bin + i.offset]

  # Melt to long, filtering to non-derived biomarkers and dropping missing values
  bio <- melt(bio, id.vars=c("eid", "visit_index"), variable.name="Biomarker", na.rm=TRUE,
    measure.vars=intersect(names(bio), ukbnmr::nmr_info[Type == "Non-derived", Biomarker]))

  # add in relevant sample processing information
  if (version == 1L) {
    bio <- sinfo[, list(eid, visit_index, Prep.to.Measure.Duration, Well.Row,
                        Well.Column, Spectrometer.Date.Bin, Spectrometer,
                        Spectrometer.Group, Shipment.Plate)][
                          bio, on = list(eid, visit_index), nomatch=0]
  } else {
    bio <- sinfo[, list(eid, visit_index, Prep.to.Measure.Duration, Well.Row,
                        Well.Column, Spectrometer.Date.Bin, Spectrometer,
                        Spectrometer.Group, Processing.Batch, Shipment.Plate)][
                          bio, on = list(eid, visit_index), nomatch=0]
  }

  message("Determining log offsets for biomarker concentrations...")

  # Determine offset for log transformation (required for variables with measurements == 0):
  # Acetate, Acetoacetate, Albumin, bOHbutyrate, Clinical_LDL_C, Gly, Ile, L_LDL_CE, L_LDL_FC, M_LDL_CE
  log_offset <- bio[!is.na(value), list(Minimum=min(value), Minimum.Non.Zero=min(value[value != 0])),by=Biomarker]
  log_offset[, Log.Offset := ifelse(Minimum == 0, Minimum.Non.Zero / 2, 0)]

  # Get log transformed raw value
  bio[log_offset, on = list(Biomarker), log_value := log(value + Log.Offset)]

  message("Adjusting for time between sample prep and sample measurement...")

  # Adjust for time between sample prep and sample measurement
  bio[, adj := MASS::rlm(log_value ~ log(Prep.to.Measure.Duration))$residuals, by=Biomarker]

  message("Adjusting for within plate structure across 96-well plate rows A-H...")

  # Adjust for within plate structure across 96-well plate rows A-H
  if (version == 1L) {
    bio[, adj := MASS::rlm(adj ~ factor_by_size(Well.Row))$residuals, by=Biomarker]
  } else {
    bio[, adj := MASS::rlm(adj ~ factor_by_size(Well.Row))$residuals, by=list(Processing.Batch, Biomarker)]
  }

  message("Adjusting for within plate structure across 96-well plate columns 1-12...")

  # Adjust for within plate structure across 96-well plate columns 1-12
  if (version == 1L) {
    bio[, adj := MASS::rlm(adj ~ factor_by_size(Well.Column))$residuals, by=Biomarker]
  } else {
    bio[, adj := MASS::rlm(adj ~ factor_by_size(Well.Column))$residuals, by=list(Processing.Batch, Biomarker)]
  }

  message("Adjusting for drift over time within spectrometer...")

  if (version == 3L) {
    # Adjust for drift over time within spectrometer (only applicable to spectrometers with > 1 bin)
    to_adjust <- sinfo[, list(N=length(unique(Spectrometer.Date.Bin))), by=list(Spectrometer)]
    to_adjust <- to_adjust[N > 1, Spectrometer]
    bio[Spectrometer %in% to_adjust, adj := MASS::rlm(adj ~ factor_by_size(Spectrometer.Date.Bin))$residuals, by=list(Biomarker, Spectrometer.Group)]
  } else {
    bio[, adj := MASS::rlm(adj ~ factor_by_size(Spectrometer.Date.Bin))$residuals, by=list(Biomarker, Spectrometer.Group)]
  }

  message("Rescaling adjusted biomarkers to absolute concentrations...")

  # Rescale to absolute units. First, shift the residuals to have the same central
  # parameter estimate (e.g. mean, estimated via robust linear regression) as the
  # log concentrations
  bio[, adj := adj + as.vector(coef(MASS::rlm(log_value ~ 1)))[1], by=Biomarker]

  # Undo the log transform
  bio[, adj := exp(adj)]

  # Remove the log offset
  bio[log_offset, on = list(Biomarker), adj := adj - Log.Offset]

  # Some values that were 0 are now < 0, apply small right shift
  # for these biomarkers (shift is very small, i.e. the impact on
  # the distribution's median is essentially numeric error)
  shift <- bio[, list(Right.Shift=-pmin(0, min(adj))), by=Biomarker]
  log_offset <- log_offset[shift, on = list(Biomarker)]
  bio[log_offset, on = list(Biomarker), adj := adj + Right.Shift]

  message("Identifying outlier plates and setting their concentrations to NA...")

  # Identify and remove outlier plates.
  # Model plate medians as a normal distribution (across all 1,352 plates), then
  # we can flag outlier plates as those > 3.374446 standard deviations from the
  # mean where 3.374446 is derived from the range of the theoretical normal
  # distribution for 1,352 samples
  n_plates <- sinfo[, length(unique(Shipment.Plate))] # 1,352
  sdlim <- max(qnorm(ppoints(n_plates))) # 3.374446
  plate_medians <- bio[, list(value = median(adj)), by=list(Biomarker, Shipment.Plate)]

  outlier_lim <- plate_medians[, list(
    Lower.Limit = mean(value) - sd(value) * sdlim,
    Mean.Plate.Medians = mean(value),
    Upper.Limit = mean(value) + sd(value) * sdlim
  ), by=Biomarker]

  plate_medians[, outlier := "no"]
  plate_medians[outlier_lim, on = list(Biomarker, value < Lower.Limit), outlier := "low"]
  plate_medians[outlier_lim, on = list(Biomarker, value > Upper.Limit), outlier := "high"]

  # Add outlier plate tags to biomarker qc tags
  if (!skip.biomarker.qc.flags) {
    message("Adding outlier plates to measurement QC tags...")

    bio_qc <- melt(bio_qc, id.vars=c("eid", "visit_index"), variable.name="Biomarker", na.rm=TRUE,
                   measure.vars=intersect(names(bio_qc), ukbnmr::nmr_info[Type == "Non-derived", Biomarker]))

    outlier_flags <- bio[plate_medians[outlier != "no"], on = list(Shipment.Plate, Biomarker),
                         list(eid, visit_index, Biomarker, value=ifelse(
                           outlier == "high", "High outlier plate", "Low outlier plate"))]

    bio_qc <- rbind(bio_qc, outlier_flags)
    bio_qc <- bio_qc[, list(value = paste(value, collapse="; ")), by=list(eid, visit_index, Biomarker)]
  }

  # Remove outlier plates from biomarker data
  if (remove.outlier.plates) {
    bio <- bio[!plate_medians[outlier != "no"], on = list(Biomarker, Shipment.Plate)]
  }

  # Convert back to wide format
  bio <- dcast(bio, eid + visit_index ~ Biomarker, value.var="adj")
  if (!skip.biomarker.qc.flags) {
    bio_qc <- dcast(bio_qc, eid + visit_index ~ Biomarker, value.var="value")
  }

  message("Recalculating derived biomarkers...")

  # Compute derived biomarkers and ratios
  bio <- nightingale_composite_biomarker_compute(bio)
  bio <- nightingale_ratio_compute(bio)
  bio <- extended_ratios_compute(bio)

  if (!skip.biomarker.qc.flags) {
    message("Collating measurement QC tags for derived biomarkers...")

    bio_qc <- nightingale_composite_biomarker_flags(bio_qc)
    bio_qc <- nightingale_ratio_flags(bio_qc)
    bio_qc <- extended_ratios_flags(bio_qc)

    # Add back in samples with no QC flags for any biomarkers
    bio_qc <- bio_qc[bio[, list(eid, visit_index)], on = list(eid, visit_index)]
  }

  # Remove dummy Spectrometer.Group column from returned output
  sinfo[, Spectrometer.Group := NULL]

  message("Returning result...")

  # Return list
  if (!skip.biomarker.qc.flags) {
    return(list(
      biomarkers=returnDT(bio),
      biomarker_qc_flags=returnDT(bio_qc),
      sample_processing=returnDT(sinfo),
      log_offset=returnDT(log_offset[Log.Offset != 0 | Right.Shift != 0]),
      outlier_plate_detection=returnDT(outlier_lim),
      algorithm_version=version
    ))
  } else {
    return(list(
      biomarkers=returnDT(bio),
      sample_processing=returnDT(sinfo),
      log_offset=returnDT(log_offset[Log.Offset != 0 | Right.Shift != 0]),
      outlier_plate_detection=returnDT(outlier_lim),
      algorithm_version=version
    ))
  }
}
