# R package ukbnmr

## Tools for processing Nightingale NMR biomarker data in UK Biobank

This package provides utilities for working with the NMR metabolomics data in [UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

There are three groups of functions in this package: (1) data extraction, (2) removal of technical variation, and (3) recomputing derived biomarkers and biomarker ratios.

All functions can be applied directly to raw data extracted from UK Biobank.

This package also provides a `data.frame` of biomarker information, loaded as `nmr_info`, and a `data.frame` of sample processing information, loaded as `sample_qc_info`.

### Data Extraction Functions

The `extract_biomarkers()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), extract the [NMR biomarker fields](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220) and give them short comprehensible column names as described in `nmr_info`. Measurements are also split into multiple rows where a participant has measurements at both baseline and first repeat assessment.

The `extract_biomarker_qc_flags()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), extract the [Nightingale quality control flags](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221) for each biomarker measurement, returning a single column per biomarker (corresponding to respective columns output by `extract_biomarkers()`).

The `extract_sample_qc_flags()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide) and extract the [sample quality control tags](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222) for the Nightingale NMR metabolomics data.

### Removal of technical variation

The `remove_technical_variation()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), remove the effects of technical variation on biomarker concentrations, and return a list containing the adjusted NMR biomarker data, biomarker QC flags, and sample quality control and processing information.

This applies a multistep process as described in Ritchie *et al.* 2021:

  1. First biomarker data is filtered to the 107 biomarkers that cannot be derived from any combination of other biomarkers.
  2. Absolute concentrations are log transformed, with a small offset applied to biomarkers with concentrations of 0.
  3. Each biomarker is adjusted for the time between sample preparation and sample measurement (hours).
  4. Each biomarker is adjusted for systematic differences between rows (A-H) on the 96-well shipment plates.
  5. Each biomarker is adjusted for remaining systematic differences between columns (1-12) on the 96-well shipment plates.
  6. Each biomarker is adjusted for drift over time within each of the six spectrometers (samples grouped into 10 bins by measurement date in each spectrometer separately.
  7. Regression residuals after the sequential adjustments are transformed back to absolute concentrations.
  8. Samples belonging to shipment plates that are outliers of non-biological origin are identified and set to missing.
  9. The 61 composite biomarkers and 81 biomarker ratios are recomputed from their adjusted parts.
  10. An additional 76 biomarker ratios of potential biological significance are computed.

Further details can be found in the preprint Ritchie S. C. *et al.*, Quality control and removal of technical variation of NMR metabolic biomarker data in ~120,000 UK Biobank participants, **medRxiv** 2021.

This function takes 10-15 minutes to run, and requires at least 16 GB of RAM, so you will want to save the output, rather than incorporate this function into your analysis scripts.

### Methods for computing biomarker ratios.

The `compute_nightingale_ratios()` function will compute the [81 Nightingale Health biomarker ratios](https://nightingalehealth.com/biomarkers) that were missing prior to the August 2021 update of [UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220). A companion function, `compute_nightingale_ratio_qc_flags()` will aggregate the QC flags for the biomarkers underlying each ratio. 

The `compute_extended_ratios()` function will compute an extended set of 76 biomarker ratios expanding on the biomarkers available directly from the Nightingale platform. A companion function, `compute_extended_ratio_qc_flags()`, will aggregate the QC flags for the biomarkers underlying each ratio.

The `recompute_derived_biomarkers()` function will recompute allcomposite biomarkers and ratios from 107 non-derived biomarkers, which is useful for ensuring data consistency when adjusting for unwanted biological variation. A companion function, `recompute_derived_biomarker_qc_flags()` will aggregate the QC flags for the biomarkers underlying each composite biomarker and ratio.

## Installation

This package can be installed from GitHub with the remotes package:

```
remotes::install_github("sritchie73/ukbnmr")
```

The package will be submitted to CRAN subsequent to preprint peer review.
