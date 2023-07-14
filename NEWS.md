# NEWS

## Version 2.0

 - Added updated version of algorithm for removing technical variation, which has
   been modified after exploration of the July 2023 release of the second tranche
   of UK Biobank NMR data covering ~275,000 participants.  
   
 - Updated GitHub README and package vignette to provide details and justification
   for the update algorithm
   
 - Added support for additional sample quality control field 20282 "Processing
   batch from Nightingale Health data", which is used as part of the updated
   algorithm for removing technical variation
   
 - Added support for new biomarker fields 20281 "Spectrometer-corrected alanine"
   and 20280 "Glucose-lactate"

## Version 1.5.3

 - Added dependency to bit64 package to ensure the Shipment.Plate column is always
   correctly interpreted by internal package functions.
 
## Version 1.5.2

 - Added support for new sample quality control field 20283 "Resolved plate swaps" 
 
 - Returned data.tables now behave as expected with respect to printing contents:
   i.e. running a function without storing the result will now show the contents
   of the returned table, and typing the name of the variable storing the result
   and hitting enter will show the contents on first try.
   
 - Fixed bug where columns corresponding to UK Biobank fields not available to
   the user would be filled with NAs rather than missing from the returned 
   results.
 
## Version 1.5.1

 - Fixed error in GitHub README example code, which has now been made consistent 
   with the vignette.
 
 - plate ID and timestamp columns in the package test data have now been set as 
   character class instead of data.table specific representations of the 
   integer64 class (from the bit64 package) and POSIXt to safeguard against 
   intermittent errors arising from incorrect type conversion when running the
   example code without first loading the data.table package.
 
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
