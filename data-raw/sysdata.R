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
  UKB.Field.ID = as.character(23649L:23655L),
  Name = c("Shipment.Plate", "Spectrometer", "Measurement.Quality.Flagged",
           "High.Lactate", "High.Pyruvate", "Low.Glucose", "Low.Protein")
)

# Measurement.Quality.Flagged field codings
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=2302
measure_quality_map <- data.table(
  integer_rep = c("1", "2"), # char not int so we can join and reassign column
  flag = c("Not enough sample material for measurement", "Solid material")
)

# Biomarker QC Flag fields that are empty in the phase 1 release, and thus are
# not in the dataset downloaded from UKB:
empty_qc_fields <- as.character(c(
  23777L, 23779L, 23740L, 23739L, 23733L, 23732L, 23731L, 23826L,
  23875L, 23812L, 23819L, 23718L, 23827L, 23820L, 23716L, 23781L,
  23730L, 23823L, 23858L, 23830L, 23795L, 23729L, 23865L, 23837L,
  23802L, 23872L, 23844L, 23809L, 23851L, 23788L, 23816L, 23728L,
  23722L, 23828L, 23877L, 23821L, 23720L, 23780L, 23706L, 23765L,
  23771L, 23714L, 23825L, 23832L, 23713L, 23874L, 23846L, 23811L,
  23818L, 23712L, 23702L, 23700L, 23701L, 23727L, 23715L, 23719L,
  23726L, 23824L, 23859L, 23831L, 23725L, 23723L, 23866L, 23838L,
  23803L, 23873L, 23845L, 23810L, 23817L, 23724L, 23711L, 23707L,
  23710L, 23829L, 23836L, 23709L, 23815L, 23822L, 23708L, 23703L
))

save(biomarker_flag_map, sample_qc_fields, measure_quality_map, empty_qc_fields,
     version=2, compress="bzip2", file="R/sysdata.rda")
