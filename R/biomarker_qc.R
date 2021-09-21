#' Remove technical variation from NMR biomarker data in UK Biobank.
#'
#' @param x \code{data.frame} containing a dataset extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' with \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}{UK Biobank fields}
#' containing the \href{https://nightingalehealth.com/biomarkers}{Nightingale Health NMR metabolomics biomarker} data.
#' @param remove.outlier.plates logical, when set to \code{FALSE} biomarker
#' concentrations on outlier shipment plates (see details) are not set to
#' missing but simply flagged in the \code{biomarker_qc_flags} \code{data.frame}
#' in the returned \code{list}.
#'
#' @details
#' A multi-step procedure is applied to the raw biomarker data to remove the
#' effects of technical variation:
#' \itemize{
#'   \item{(1) }{First biomarker data is filtered to the 107 biomarkers that
#'   cannot be derived from any combination of other biomarkers.}
#'   \item{(2) }{Absolute concentrations are log transformed, with a small offset
#'   applied to biomarkers with concentrations of 0.}
#'   \item{(3) }{Each biomarker is adjusted for the time between sample preparation
#'   and sample measurement (hours)}
#'   \item{(4) }{Each biomarker is adjusted for systematic differences between rows
#'   (A-H) on the 96-well shipment plates.}
#'   \item{(5) }{Each biomarker is adjusted for remaining systematic differences
#'   between columns (1-12) on the 96-well shipment plates.}
#'   \item{(6) }{Each biomarker is adjusted for drift over time within each of the
#'   six spectrometers (samples grouped into 10 bins by measurement date in each
#'   spectrometer separately)}
#'   \item{(7) }{Regression residuals after the sequential adjustments are
#'   transformed back to absolute concentrations.}
#'   \item{(8) }{Samples belonging to shipment plates that are outliers of
#'   non-biological origin are identified and set to missing.}
#'   \item{(9) }{The 61 composite biomarkers and 81 biomarker ratios are recomputed
#'   from their adjusted parts.}
#'   \item{(10) }{An additional 76 biomarker ratios of potential biological
#'   significance are computed.}
#' }
#'
#' At each step, adjustment for technical covariates is performed using
#' \link[MASS:rlm]{robust linear regression}. Plate row, plate column, and
#' sample measurement date bin are treated as factors, using the group with the
#' largest sample size as reference in the regression.
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
#'         plate compared to all other shipment plates (see Details).}
#'   \item{sample_processing}{A \code{data.frame} containing the
#'         \href{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222}{processing
#'         information and quality control indicators} for each sample. See
#'         \code{\link{sample_qc_fields}} for details.}
#' }
#'
#' @importFrom stats coef
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom stats qnorm
#' @importFrom stats ppoints
#' @importFrom MASS rlm
#' @export
remove_technical_variation <- function(x, remove.outlier.plates=TRUE) {
  # Silence CRAN NOTES about undefined global variables (columns in data.tables)
  Well.Position.Within.Plate <- Well.Row <- Well.Column <- Sample.Measured.Date.and.Time <-
    Sample.Prepared.Date.and.Time <- Sample.Measured.Date <- Sample.Measured.Time <-
    Sample.Prepared.Date <- Sample.Prepared.Time <- Prep.to.Measure.Duration <-
    Spectrometer.Date.Bin <- Spectrometer <- Type <- Biomarker <- eid <- visit_index <-
    Shipment.Plate <- value <- Log.Offset <- Minimum <- Minimum.Non.Zero <-
    log_value <- adj <- Right.Shift <- outlier <- lower_lim <- upper_lim <-
    Name <- UKB.Field.ID <- NULL

  # Check relevant fields exist
  f1 <- detect_format(x, type="biomarkers")
  f2 <- detect_format(x, type="biomarker_qc_flags")
  f3 <- detect_format(x, type="sample_qc_flags")

  # Make sure its not a "processed" dataset returned by extract_biomarkers() and
  # related functions.
  if (f1 == "processed" || f2 == "processed" || f3 == "processed") {
    stop("'x' should be a 'data.frame' containg UK Biobank field ID columns, not one extracted by extract_biomarkers() or related functions'")
  }

  # Extract biomarkers, biomarker QC flags, and sample processing information
  bio <- process_data(x, type="biomarkers")
  bio_qc <- process_data(x, type="biomarker_qc_flags")
  sinfo <- process_data(x, type="sample_qc_flags")

  # Check required sample processing fields exist
  req <- c("Well.Position.Within.Plate", "Sample.Measured.Date.and.Time",
           "Sample.Prepared.Date.and.Time", "Spectrometer", "Shipment.Plate")
  miss_req <- setdiff(req, names(sinfo))
  if (length(miss_req) > 0) {
    err_txt <- ukbnmr::sample_qc_fields[Name %in% miss_req]
    err_txt <- err_txt[, sprintf("%s (Field: %s)", Name, UKB.Field.ID)]
    err_txt <- paste(err_txt, collapse=", ")
    stop("Missing required sample processing fields: ", err_txt, ".")
  }

  # Split out row and column information from 96-well plate information
  sinfo[, Well.Row := gsub("[0-9]", "", Well.Position.Within.Plate)]
  sinfo[, Well.Column := as.integer(gsub("[A-Z]", "", Well.Position.Within.Plate))]

  # Split out date and time for sample measurement and sample prep
  sinfo[, Sample.Measured.Date := as.IDate(gsub(" .*$", "", Sample.Measured.Date.and.Time))]
  sinfo[, Sample.Measured.Time := as.ITime(gsub("^.* ", "", Sample.Measured.Date.and.Time))]
  sinfo[, Sample.Prepared.Date := as.IDate(gsub(" .*$", "", Sample.Prepared.Date.and.Time))]
  sinfo[, Sample.Prepared.Time := as.ITime(gsub("^.* ", "", Sample.Prepared.Date.and.Time))]

  # Compute hours between sample prep and sample measurement
  sinfo[, Prep.to.Measure.Duration := duration.hours(
    Sample.Prepared.Date, Sample.Prepared.Time,
    Sample.Measured.Date, Sample.Measured.Time
  )]

  # Split sample measurement date into 10 equal size bins per spectrometer
  sinfo[, Spectrometer.Date.Bin := bin_dates(Sample.Measured.Date), by=Spectrometer]

  # Melt to long, filtering to non-derived biomarkers and dropping missing values
  bio <- melt(bio, id.vars=c("eid", "visit_index"), variable.name="Biomarker", na.rm=TRUE,
    measure.vars=intersect(names(bio), ukbnmr::nmr_info[Type == "Non-derived", Biomarker]))

  # add in relevant sample processing information
  bio <- sinfo[, list(eid, visit_index, Prep.to.Measure.Duration, Well.Row,
    Well.Column, Spectrometer.Date.Bin, Spectrometer, Shipment.Plate)][
      bio, on = list(eid, visit_index), nomatch=0]

  # Determine offset for log transformation (required for variables with measurements == 0):
  # Acetate, Acetoacetate, Albumin, bOHbutyrate, Clinical_LDL_C, Gly, Ile, L_LDL_CE, L_LDL_FC, M_LDL_CE
  log_offset <- bio[!is.na(value), list(Minimum=min(value), Minimum.Non.Zero=min(value[value != 0])),by=Biomarker]
  log_offset[, Log.Offset := ifelse(Minimum == 0, Minimum.Non.Zero / 2, 0)]

  # Get log transformed raw value
  bio[log_offset, on = list(Biomarker), log_value := log(value + Log.Offset)]

  # Adjust for time between sample prep and sample measurement
  bio[, adj := MASS::rlm(log_value ~ Prep.to.Measure.Duration)$residuals, by=Biomarker]

  # Adjust for within plate structure across 96-well plate rows A-H
  bio[, adj := MASS::rlm(adj ~ factor_by_size(Well.Row))$residuals, by=Biomarker]

  # Adjust for within plate structure across 96-well plate columns 1-12
  bio[, adj := MASS::rlm(adj ~ factor_by_size(Well.Column))$residuals, by=Biomarker]

  # Adjust for drift over time within spectrometer
  bio[, adj := MASS::rlm(adj ~ factor_by_size(Spectrometer.Date.Bin))$residuals, by=list(Biomarker, Spectrometer)]

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

  # Identify and remove outlier plates.
  # Model plate medians as a normal distribution (across all 1,352 plates), then
  # we can flag outlier plates as those > 3.374446 standard deviations from the
  # mean where 3.374446 is derived from the range of the theoretical normal
  # distribution for 1,352 samples
  n_plates <- sinfo[, length(unique(Shipment.Plate))] # 1,352
  sdlim <- max(qnorm(ppoints(n_plates))) # 3.374446
  plate_medians <- bio[, list(value = median(adj)), by=list(Biomarker, Shipment.Plate)]

  outlier_lim <- plate_medians[, list(
    lower_lim = mean(value) - sd(value) * sdlim,
    upper_lim = mean(value) + sd(value) * sdlim
  ), by=Biomarker]

  plate_medians[, outlier := "no"]
  plate_medians[outlier_lim, on = list(Biomarker, value < lower_lim), outlier := "low"]
  plate_medians[outlier_lim, on = list(Biomarker, value > upper_lim), outlier := "high"]

  # Add outlier plate tags to biomarker qc tags
  bio_qc <- melt(bio_qc, id.vars=c("eid", "visit_index"), variable.name="Biomarker", na.rm=TRUE,
              measure.vars=intersect(names(bio_qc), ukbnmr::nmr_info[Type == "Non-derived", Biomarker]))

  outlier_flags <- bio[plate_medians[outlier != "no"], on = list(Shipment.Plate, Biomarker),
    list(eid, visit_index, Biomarker, value=ifelse(outlier == "high", "High outlier plate", "Low outlier plate"))]

  bio_qc <- rbind(bio_qc, outlier_flags)
  bio_qc <- bio_qc[, list(value = paste(value, collapse="; ")), by=list(eid, visit_index, Biomarker)]

  # Remove outlier plates from biomarker data
  if (remove.outlier.plates) {
    bio <- bio[!plate_medians[outlier != "no"], on = list(Biomarker, Shipment.Plate)]
  }

  # Convert back to wide format
  bio <- dcast(bio, eid + visit_index ~ Biomarker, value.var="adj")
  bio_qc <- dcast(bio_qc, eid + visit_index ~ Biomarker, value.var="value")

  # Compute derived biomarkers and ratios
  bio <- nightingale_composite_biomarker_compute(bio)
  bio <- nightingale_ratio_compute(bio)
  bio <- extended_ratios_compute(bio)

  bio_qc <- nightingale_composite_biomarker_flags(bio_qc)
  bio_qc <- nightingale_ratio_flags(bio_qc)
  bio_qc <- extended_ratios_flags(bio_qc)

  # Return list
  return(list(
    biomarkers=returnDT(bio),
    biomarker_qc_flags=returnDT(bio_qc),
    sample_processing=returnDT(sinfo)
  ))
}
