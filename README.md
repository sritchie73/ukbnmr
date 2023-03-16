# The ukbnmr R package

This package provides utilities for working with the [UK Biobank NMR metabolomics data](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

There are three groups of functions in this package: (1) data extraction, (2) removal of technical variation, and (3) recomputing derived biomarkers and biomarker ratios.

All functions can be applied directly to raw data extracted from UK Biobank.

This package also provides a `data.frame` of biomarker information, loaded as `nmr_info`, and a `data.frame` of sample processing information, loaded as `sample_qc_info`. See `help("nmr_info")` and `help("sample_qc_info")` for details on column contents.

## Data Extraction Functions

The `extract_biomarkers()` function will take a decoded UK Biobank dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), extract the [NMR metabolomics biomarker fields](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220) and give them short comprehensible column names as described in `nmr_info`. Measurements are also split into multiple rows where a participant has measurements at both baseline and first repeat assessment.

The `extract_biomarker_qc_flags()` function will take a decoded UK Biobank dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), extract the [per-biomarker measurement quality control flags](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221) for each biomarker measurement, returning a single column per biomarker (corresponding to respective columns output by `extract_biomarkers()`).

The `extract_sample_qc_flags()` function will take a decoded UK Biobank dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide) and extract the [sample quality control tags](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222) for the NMR metabolomics data.

An example workflow for extracting these data and saving them for later use:

```R
library(ukbnmr)

decoded <- fread("path/to/decoded_ukbiobank_data.csv") # file save by ukbconv tool

nmr <- extract_biomarkers(decoded)
biomarker_qc_flags <- extract_biomarker_qc_flags(decoded)
sample_qc_flags <- extract_sample_qc_flags(decoded)

fwrite(nmr, file="path/to/nmr_biomarker_data.csv")
fwrite(biomarker_qc_flags, file="path/to/nmr_biomarker_qc_flags.csv")
fwrite(sample_qc_flags, file="path/to/nmr_sample_qc_flags.csv")
```

## Removal of technical variation

The `remove_technical_variation()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), remove the effects of technical variation on biomarker concentrations, and return a list containing the adjusted NMR biomarker data, biomarker QC flags, and sample quality control and processing information.

This applies a multistep process as described in Ritchie *et al.* 2021:

  1. First biomarker data is filtered to the 107 biomarkers that cannot be derived from any combination of other biomarkers.
  2. Absolute concentrations are log transformed, with a small offset applied to biomarkers with concentrations of 0.
  3. Each biomarker is adjusted for the time between sample preparation and sample measurement (hours) on a log scale.
  4. Each biomarker is adjusted for systematic differences between rows (A-H) on the 96-well shipment plates.
  5. Each biomarker is adjusted for remaining systematic differences between columns (1-12) on the 96-well shipment plates.
  6. Each biomarker is adjusted for drift over time within each of the six spectrometers. To do so, samples are grouped into 10
     bins, within each spectrometer, by the date the majority of samples on their respective 96-well plates were measured.
  7. Regression residuals after the sequential adjustments are transformed back to absolute concentrations.
  8. Samples belonging to shipment plates that are outliers of non-biological origin are identified and set to missing.
  9. The 61 composite biomarkers and 81 biomarker ratios are recomputed from their adjusted parts.
  10. An additional 76 biomarker ratios of potential biological significance are computed.

Further details can be found in the preprint Ritchie S. C. *et al.*, Quality control and removal of technical variation of NMR metabolic biomarker data in ~120,000 UK Biobank participants, **medRxiv** (2021). doi: [10.1101/2021.09.24.21264079](https://www.medrxiv.org/content/10.1101/2021.09.24.21264079v1).

This function takes 10-15 minutes to run, and requires at least 16 GB of RAM, so you will want to save the output, rather than incorporate this function into your analysis scripts.

An example workflow for using this function and saving the output for loading into future R sessions or other programs:

```R
library(ukbnmr)
decoded <- fread("path/to/decoded_ukbiobank_data.csv") # file save by ukbconv tool

processed <- remove_technical_variation(decoded)

fwrite(processed$biomarkers, file="path/to/nmr_biomarker_data.csv")
fwrite(processed$biomarker_qc_flags, file="path/to/nmr_biomarker_qc_flags.csv")
fwrite(processed$sample_processing, file="path/to/nmr_sample_qc_flags.csv")
fwrite(processed$log_offset, file="path/to/nmr_biomarker_log_offset.csv")
fwrite(processed$outlier_plate_detection, file="path/to/outlier_plate_info.csv")
```

## Methods for computing derived biomarkers and ratios

Analysts may wish to further adjust data for biological covariates. We provide an additional function, `recompute_derived_biomarkers()` to recompute all composite biomarkers and ratios from 107 non-derived biomarkers, which is useful for ensuring data consistency when adjusting for unwanted biological variation. A companion function, `recompute_derived_biomarker_qc_flags()` will aggregate the QC flags for the biomarkers underlying each composite biomarker and ratio.

Note these functions assume the data has been returned to absolute units after adjusting for technical covariates. For example the ratio of two biomarkers A and B is computed as A/B, which may not be true if the two biomarkers are on different scales (e.g. regression residuals) after adjustment.

If using these functions, please cite Ritchie S. C. *et al.*, Quality control and removal of technical variation of NMR metabolic biomarker data in ~120,000 UK Biobank participants, **Sci Data** *10* 64 (2023). doi: [10.1038/s41597-023-01949-y](https://www.nature.com/articles/s41597-023-01949-y).

An example workflow:

```R
library(ukbnmr)

# First, if we haven't corrected for unwanted technical variation we do so
# using the appropriate function (see help("remove_technical_variation")).
decoded <- fread("path/to/decoded_ukbiobank_data.csv") # file save by ukbconv tool

processed <- remove_technical_variation(decoded)
tech_qc <- processed$biomarkers

fwrite(tech_qc, file="path/to/nmr_biomarker_data.csv")
fwrite(processed$biomarker_qc_flags, file="path/to/nmr_biomarker_qc_flags.csv")
fwrite(processed$sample_processing, file="path/to/nmr_sample_qc_flags.csv")
fwrite(processed$log_offset, file="path/to/nmr_biomarker_log_offset.csv")
fwrite(processed$outlier_plate_detection, file="path/to/outlier_plate_info.csv")

# Otherwise assuming we load 'tech_qc' from "path/to/mr_biomarker_data.csv".

# We now run code to adjust biomarkers for biological covariates. This code is
# not supplied by this package, but for illustrative purposes we assume the user
# has written a function to do this:
bio_qc <- user_function_to_adjust_biomarkers_for_covariates(tech_qc)

# Now we recompute the composite biomarkers and derived ratios after
# adjustment for additional biological covariates
bio_qc <- recompute_derived_biomarkers(bio_qc)
fwrite(bio_qc, file="path/to/nmr_biomarkers_adjusted_for_covariates.csv")

# You may also want to aggregate and save the quality control flags for each
# sample from the biomarkers underlying each derived biomarker or ratio,
# adding them as additional columns to the input data (see
# help("recompute_derived_biomarker_qc_flags")).
biomarker_qc_flags <- recompute_derived_biomarker_qc_flags(nmr)
fwrite(biomarker_qc_flags, file="path/to/biomarker_qc_flags.csv")
```

## Installation

This package can be installed from GitHub with the remotes package:

```R
remotes::install_github("sritchie73/ukbnmr", dependencies = TRUE, build_vignettes = TRUE)
```
