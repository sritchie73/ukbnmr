# R package ukbnmr

## Tools for processing Nightingale NMR biomarker data in UK Biobank

This package provides utilities for working with the NMR metabolomics data in [UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

There are two groups of functions in this package: (1) data extraction functions,
and (2) methods for computing biomarker ratios.

This package also provides a `data.frame` of biomarker information, loaded
as `nmr_info`.

### Data Extraction Functions

The `extract_biomarkers()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), extract the [NMR biomarker fields](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220) and give them short comprehensible column names. Measurements are also split into multiple rows where a participant has measurements at both baseline and first repeat assessment.

The `extract_biomarker_qc_flags()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), extract the [Nightingale quality control flags](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=221) for each biomarker measurement, returning a single column per biomarker (corresponding to respective columns output by `extract_biomarkers()`).

The `extract_sample_qc_flags()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide) and extract the [sample quality control tags](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222) for the Nightingale NMR metabolomics data.

### Methods for computing biomarker ratios.

The `compute_nightingale_ratios()` function will compute the [81 Nightingale Health biomarker ratios](https://nightingalehealth.com/biomarkers) from the 168 biomarkers currently available for download from [UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

## Installation

This package can be installed from GitHub with the remotes package:

```
remotes::install_github("sritchie73/ukbnmr")
```
