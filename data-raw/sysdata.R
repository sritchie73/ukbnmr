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

save(biomarker_flag_map, measure_quality_map,
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

# Order columns
name_order <- data.table(names=names(test_data))
name_order[, field_num := as.integer(gsub("-.*", "", names))]
name_order[, visit_num := as.integer(gsub(".*-", "", gsub("\\..*", "", names)))]
name_order[, instance_num := as.integer(gsub(".*\\.", "", names))]
name_order[names == "eid", c("field_num", "visit_num", "instance_num") := .(0, 0, 0)]
name_order <- name_order[order(instance_num)][order(visit_num)][order(field_num)]
test_data <- test_data[,.SD,.SDcols=name_order$names]

# Ensure plate columns are character, not integer64, when saving object
test_data[, c("23649-0.0") := as.character(`23649-0.0`)]
test_data[, c("23649-1.0") := as.character(`23649-1.0`)]

# Save
save(test_data, version=2, compress="bzip2", file="~/test_data.rda")
