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
