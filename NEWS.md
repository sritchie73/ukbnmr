# NEWS

## Version 1.5

 - Updated citations in documentation to reflect article publication in 
   Scientific Data

## Version 1.4

 - Created example toy dataset for testing package functions and updated 
   documentation
 
 - Removed GitHub README page from package bundle

## Version 1.3
 
 - Fixed URL in DESCRIPTION

## Version 1.2

 - Minor changes to fix DESCRIPTION based on feedback from CRAN maintainers 

## Version 1.1

- Minor changes to fix NOTES and WARNINGS thrown by CRAN: 

  - URLs which have been moved since the initial documentation was written have 
    been fixed.
    
  - Examples have been added to the documentation for each package function.
  
  - A vignette including the most relevant example workflows has been added.

## Version 1.0

- Minor changes to DESCRIPTION, README, and inst/CITATION to prepare for CRAN
  submission following paper acceptance.

## Version 0.7

- `biomarker_qc()` now corrects for sample degradation time on a log scale 
  instead of a linear scale to model exponential decay. This mainly impacts
  histidine concentrations, whose new post-QC values have Pearson correlation
  of 0.974 with those assuming linear decay, while all other biomarkers have
  Pearson correlation > 0.99.

## Version 0.6.3

- `extract_biomarker_qc_flags()` no longer relies on list of hard coded field 
  IDs it expects not to be present due to no data in UK Biobank, making code
  more robust to differences between UK Biobank data releases.

## Version 0.6.2

- `remove_technical_variation()` now gains a `skip.biomarker.qc.flags` argument
  that allows you to skip the collation and curation of biomarker QC flags when
  removing technical variation from the data
  
- `extract_biomarker_qc_flags()` and `remove_technical_variation()` should now
  be more robust to potential changes made by UK Biobank in extracted data 
  formats, particularly for fields that are empty/contain no data in the 
  showcase (or which gain data between showcase updates).
