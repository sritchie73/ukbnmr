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

# Biomarker QC Flag fields that are empty in the phase 1 release, and thus are
# not in the dataset downloaded from UKB:
empty_qc_fields <- c(
  23700L, 23701L, 23702L, 23703L, 23706L, 23707L, 23708L, 23709L,
  23710L, 23711L, 23712L, 23713L, 23714L, 23715L, 23716L, 23718L,
  23719L, 23720L, 23722L, 23723L, 23724L, 23725L, 23726L, 23727L,
  23728L, 23729L, 23730L, 23731L, 23732L, 23733L, 23739L, 23740L,
  23741L, 23765L, 23771L, 23777L, 23779L, 23780L, 23781L, 23788L,
  23795L, 23802L, 23803L, 23809L, 23810L, 23811L, 23812L, 23815L,
  23816L, 23817L, 23818L, 23819L, 23820L, 23821L, 23822L, 23823L,
  23824L, 23825L, 23826L, 23827L, 23828L, 23829L, 23830L, 23831L,
  23832L, 23836L, 23837L, 23838L, 23839L, 23844L, 23845L, 23846L,
  23851L, 23858L, 23859L, 23865L, 23866L, 23872L, 23873L, 23874L,
  23875L, 23877L, 23899L, 23900L, 23903L, 23904L, 23905L, 23906L,
  23907L, 23908L, 23909L, 23910L, 23911L, 23912L, 23913L, 23914L,
  23918L, 23919L, 23924L, 23944L, 23945L, 23947L
)

save(biomarker_flag_map, measure_quality_map, empty_qc_fields,
     version=2, compress="bzip2", file="R/sysdata.rda")
