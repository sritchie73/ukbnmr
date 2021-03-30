# R package ukbnmr

## Tools for processing Nightingale NMR biomarker data in UK Biobank

This package provides utilities for working with the NMR metabolomics data in [UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

There are currently two main functions: `extract_biomarkers()` and `compute_nightingale_ratios()`.

The `extract_biomarkers()` function will take a raw dataset output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide), extract the NMR biomarker fields and give them short comprehensible column names. Measurements are also split into multiple rows where a participant has measurements at both baseline and first repeat assessment.

The `compute_nightingale_ratios()` function will compute the [81 Nightingale Health biomarker ratios](https://nightingalehealth.com/biomarkers) from the 168 biomarkers available for download from [UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

This package also provides a `data.frame` of biomarker information, loaded
as `nmr_info`.

## Installation

This package can be installed from GitHub with the remotes package:

```
remotes::install_github("sritchie73/ukbnmr")
```
