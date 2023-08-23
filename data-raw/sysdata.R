## Internal data tables used to map QC Flags to their values


# Biomarker QC Flags
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=2310
biomarker_flag_map <- data.table(
  integer_rep = 1L:10L,
  flag = c("Below limit of quantification", "Citrate plasma", "Degraded sample",
           "High ethanol", "Isopropyl alcohol", "Low glutamine or high glutamate",
           "Medium ethanol", "Polysaccharides", "Unknown contamination", "Ethanol")
)

# Measurement.Quality.Flagged field codings
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=2302
measure_quality_map <- data.table(
  integer_rep = c("1", "2"), # char not int so we can join and reassign column
  flag = c("Not enough sample material for measurement", "Solid material")
)

# Hand curated mapping between Shipment.Plate and Processing.Batch in phase 1 + 2 data release
plate_batch_map <- fread("~/Clusters/CSD3/batch_map.csv") # Stored on CSD3

save(biomarker_flag_map, measure_quality_map, plate_batch_map,
     version=2, compress="bzip2", file="R/sysdata.rda")

# Random sampling to create test dataset from real UK data.
setwd("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/NMR_metabolomics")
decoded <- fread("data/raw/ukbiobank/extracted/nmr.csv")
extracted <- fread("output/nmr_techadj.txt")
test_data <- decoded[0][1:50]

# Split data to reflect the presence of data at two timepoints. The split here
# ensures that the simulated samples in the test_data always have data at 1 or
# both time points (depending on the row) rather than a hodgepodge depending on
# the random sampling
ev1 <- extracted[visit_index == 0]
ev2 <- extracted[visit_index == 1]
fields <- unique(gsub("-.*", "", names(decoded)))


for (ff in fields[-1]) {
  if (paste0(ff, "-0.0") %in% names(decoded)) {
    dv1ff <- decoded[eid %in% ev1$eid][[paste0(ff, "-0.0")]]
    test_data[3:50, paste0(ff, "-0.0") := sample(dv1ff, 48, replace=TRUE)]
  }
}

for (ff in fields[-1]) {
  if (paste0(ff, "-1.0") %in% names(decoded)) {
    dv2ff <- decoded[eid %in% ev2$eid][[paste0(ff, "-1.0")]]
    test_data[1:10, paste0(ff, "-1.0") := sample(dv2ff, 10, replace=TRUE)]
  }
}

# Duration between sample prep and measurement comes from two columns, we need
# to sample these together to make sure values are sensible (e.g. measurement
# occurring after sample prep, number of hours between two events is within
# expected range observed in UKB)
dv1sp <- decoded[eid %in% ev1$eid][sample(1:.N, 48, replace=TRUE)]
test_data[3:50, c("23658-0.0") := dv1sp[["23658-0.0"]]]
test_data[3:50, c("23659-0.0") := dv1sp[["23659-0.0"]]]

dv2sp <- decoded[eid %in% ev2$eid][sample(1:.N,10, replace=TRUE)]
test_data[1:10, c("23658-1.0") := dv2sp[["23658-1.0"]]]
test_data[1:10, c("23659-1.0") := dv2sp[["23659-1.0"]]]

# For plate well, sample from a handful of possible values so the groupings are
# large enough for rlm to converge
mkwell <- function(N) {
  sprintf("%s0%s",
          sample(c("A", "C", "F"), N, replace=TRUE),
          sample(c(3, 6, 9), N, replace=TRUE))
}
test_data[3:50, c("23660-0.0") := mkwell(48)]
test_data[1:10, c("23660-1.0") := mkwell(10)]

# For the same reason, set all samples to be from the same spectrometer
test_data[3:50, c("23650-0.0") := 3L]
test_data[1:10, c("23650-1.0") := 3L]

# Set dummy data for Processing Batch
test_data[!is.na(`23649-0.0`), c("20282-0.0") := 1L]
test_data[!is.na(`23649-1.0`), c("20282-1.0") := 1L]

# Create fake sample identifiers and shuffle rows
test_data[, eid := sample(1001:2001, 50)]
test_data <- test_data[sample(1:50, 50, replace=FALSE)]

# Add in some resolved plate swaps
test_data[, c("20283-0.0") := NA_integer_]
test_data[, c("20283-1.0") := NA_integer_]

test_data[!is.na(`23650-0.0`), c("20283-0.0") := sample(c(1L, NA_integer_), .N, replace=TRUE, prob=c(0.2, 0.8))]
test_data[!is.na(`23650-1.0`), c("20283-1.0") := sample(c(1L, NA_integer_), .N, replace=TRUE, prob=c(0.2, 0.8))]

# Add in new biomarkers "Spectrometer corrected alanine" and "Glucose-Lactate"
test_data[!is.na(`23460-0.0`), c("20281-0.0") := sample(`23460-0.0`, .N)]
test_data[!is.na(`23460-1.0`), c("20281-1.0") := sample(`23460-1.0`, .N)]

# For Glucose-lactate we don't have the data available yet to simulate realistic distribution of values,
# but based on the showcase it seems to be a not-quite linear combination of Glucose and lactate (i.e.
# both glucose and lactate have median concentrations of 3.X, whereas Glucose-lactate has median of 5.X)
test_data[!is.na(`23470-0.0`) & !is.na(`23471-0.0`), c("20280-0.0") := (`23470-0.0` + `23471-0.0`)*sample(seq(0.65, 0.75, by=0.01), .N, replace=TRUE)]
test_data[!is.na(`23470-1.0`) & !is.na(`23471-1.0`), c("20280-1.0") := (`23470-1.0` + `23471-1.0`)*sample(seq(0.65, 0.75, by=0.01), .N, replace=TRUE)]

# Order columns
name_order <- data.table(names=names(test_data))
name_order[, field_num := as.integer(gsub("-.*", "", names))]
name_order[, visit_num := as.integer(gsub(".*-", "", gsub("\\..*", "", names)))]
name_order[, instance_num := as.integer(gsub(".*\\.", "", names))]
name_order[names == "eid", c("field_num", "visit_num", "instance_num") := .(0, 0, 0)]
name_order <- name_order[order(instance_num)][order(visit_num)][order(field_num)]
test_data <- test_data[,.SD,.SDcols=name_order$names]

# Extract smaller set of rows that won't cause rlm to to fail to converge to
# make testing code faster
test_data <- test_data[c(25L, 35L, 14L, 13L, 6L, 26L, 9L, 12L, 46L, 43L)]

# Extract smaller subset of columns for illustrative purposes
sinfo_tag_cols <- c(20282, 23658, 23659, 23649, 23650, 23660)
nmr_cols <- c(23465, 23466, 23467, 23464, 23444, 23445, 23446, 23459)
qc_cols <- c(23746, 23751, 23752, 23759, 23764, 23765, 23766, 23767)
fields <- sort(c(sinfo_tag_cols, nmr_cols, qc_cols))
fields <- as.vector(sapply(fields, function(fid) { c(paste0(fid, "-0.0"), paste0(fid, "-1.0")) }))
test_data <- test_data[,.SD,.SDcols=intersect(c("eid", fields), names(test_data))]

# Manually add in some QC flags so they're not all empty
test_data[1, `23765-0.0` := 1L]
test_data[10, `23765-1.0` := 1L]
test_data[c(2,3), `23744-0.0` := 9L]
test_data[7, `23744-0.0` := 8L]
test_data[, `23764-0.0` := as.integer(`23764-0.0`) ]
test_data[, `23767-0.0` := as.integer(`23767-0.0`) ]
test_data[c(1, 4), c("23764-0.0", "23767-0.0") := 4L]

# Manually add in some zeros so log_offset table can be returned
test_data[1, `23465-0.0` := 0]
test_data[10, `23465-1.0` := 0]

# Save
save(test_data, version=2, compress="bzip2", file="data/test_data.rda")
