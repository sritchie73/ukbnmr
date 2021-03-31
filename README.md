# R package ukbnmr

## Tools for processing Nightingale NMR biomarker data in UK Biobank

This package provides utilities for working with the NMR metabolomics data in [UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

There are currently four main functions:

The `extract_biomarkers()` function will take a raw dataset
output by [ukbconv](https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide),
extract the NMR biomarker fields and give them short comprehensible column
names. Measurements are also split into multiple rows where a participant has
measurements at both baseline and first repeat assessment.

The `compute_nightingale_ratios()` function will compute the
[81 Nightingale Health biomarker ratios](\https://nightingalehealth.com/biomarkers)
that were missing prior to the May update of
[UK Biobank](https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220).

The `compute_extended_ratios()` function will compute an extended
set of biomarker ratios expanding on the biomarkers available directly from
the Nightingale platform.

The `recompute_derived_biomarkers()` function will recompute all
composite biomarkers and ratios from 107 non-derived biomarkers, which is
useful for ensuring data consistency when adjusting for unwanted biological
variation.

All functions can be applied directly to raw data extracted from UK Biobank.

This package also provides a `data.frame` of biomarker information, loaded
as `nmr_info`.

## Installation

This package can be installed from GitHub with the remotes package:

```
remotes::install_github("sritchie73/ukbnmr")
```
