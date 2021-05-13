## Internal data tables used to map QC Flags to their values


# Biomarker QC Flags
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=2310
biomarker_flag_map <- data.table(
  integer_rep = 1L:10L,
  flag = c("Below limit of quantification", "Citrate plasma", "Degraded sample",
           "High ethanol", "Isopropyl alcohol", "Low glutamine or high glutamate",
           "Medium ethanol", "Polysaccharides", "Unknown contamination", "Ethanol")
)

# Sample QC Field names
# https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222
sample_qc_fields <- data.table(
  UKB.Field.ID = 23649L:23655L,
  Name = c("Shipment.Plate", "Spectrometer", "Measurement.Quality.Flagged",
           "High.Lactate", "High.Pyruvate", "Low.Glucose", "Low.Protein")
)

# Measurement.Quality.Flagged field codings
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=2302
measure_quality_map <- data.table(
  integer_rep = c(NA_integer_, 1L:2L),
  flag = c(NA_character_, "Not enough sample material for measurement", "Solid material")
)

save(biomarker_flag_map, sample_qc_fields, measure_quality_map,
     version=2, compress="bzip2", file="R/sysdata.rda")
