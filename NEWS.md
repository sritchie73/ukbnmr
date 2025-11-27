# NEWS

## Version 3.3.0

 - Updated package to work with UK Biobank V20 data release. The V20 data 
   release included 182 samples missing sample preparation date and time, 
   which are needed during the QC procedure. Fortunately cross-referencing 
   the V20 data release with the advance access data in UK Biobank application 
   30418 provided to us directly by Nightingale Health revealed that all 182
   samples shared the same date+time of sample preparation making it possible
   to hard code this date+time stamp into the QC procedure for missing samples.
   
 - Removed unecessary dependency in DESCRIPTION to roxygen2 package (this 
   package is only required to build the package manuals during package 
   development, and is not required by end users on package installation).
   
 - Updated computation of hours between sample preparation and sample 
   measurement for forward compatibility with future lubridate updates.

## Version 3.2.0

 - Removed broken URLs flagged by CRAN submission
 
 - Removed deprecated functionality that is no longer relevant since the UK 
   Biobank data refresh in July 2023.

## Version 3.1.0

 - Added worked example of adjusting for biological covariates and recomputing
   derived biomarkers (see https://github.com/sritchie73/ukbnmr/issues/7).

## Version 3.0.0

 - Updated README page and vignette with details on technical variation in full
   UK Biobank data release.

 - Fixed bug where empty character strings were incorrectly not interpreted as 
   missing data (see https://github.com/sritchie73/ukbnmr/issues/10).

 - Added dependency to the lubridate package to enable more flexible parsing of
   date-time formats, which differ between data extraction and data reading 
   methods (see https://github.com/sritchie73/ukbnmr/issues/10).

 - Updated package documentation to reflect that the UK Biobank Research Analysis
   Platform data format is what will now be encountered by analysts due to UK 
   Biobank policy changes on data access and downloads.

 - Updated package to work with new data format on the UK Biobank Research 
   Analysis Platform.

## Version 2.2.2

 - Added progress messages to the long-running remove_technical_variation() 
   function.

## Version 2.2.1

 - Fixed unexpected warning about NAs arising due to mislabelling of well 
   position by UK Biobank to use lowercase instead of uppercase letters for
   the row position (see https://github.com/sritchie73/ukbnmr/issues/6)

## Version 2.2

 - Fixed issue where plate numbers only sometimes had a leading zero included in
   their identifier. For consistency with the data released in UK Biobank, plate 
   identifiers now always include a leading 0.

## Version 2.1.3

 - Dramatically reduced test_data size in part due to unsuccessful attempts to 
   address a CRAN NOTE. Test data has now been reduced to 10 rows as before, and
   39 columns covering a subset of 8 biomarkers and 6 required sample processing
   flags - this also makes example output easier to visually inspect.
 
## Version 2.1.2

 - CRAN NOTE resolved for remove_technical_variation() using solution below, but
   similar NOTE was generated for other functions on submission of version 2.1.1
   to CRAN. Work-around code has now been added to all user-facing functions.

## Version 2.1.1

 - Reverted test_data back to 50 row version after cutting down size did not
   resolve the problematic CRAN NOTE (see version 2.1 below)
   
 - R-package-devel suggested CRAN NOTE issue may be due to misconfiguration of
   CRAN's server making data.table use too many threads while running examples.
   As a work around, remove_technical_variation() now explicitly sets the number
   of threads data.table can use to 1 if it is running on ukbnmr::test_data 
   (thanks to Ivan Krylov and Dirk Eddelbuettel if this does solve the issue) 
   
## Version 2.1

 - Processing.Batch now inferred from Shipment.Plate if Processing.Batch missing
   (UK Biobank field #20282) in input data and algorithm version 2 used in 
   remove_technical_variation().

 - Package overview help file is now correctly documented as requested by CRAN
   following breaking changes in Roxygen 7.0.0 that changed the way this help 
   file was internally tagged in the source code.
   
 - Reduced test_data from 50 to 10 rows to (hopefully) get around CRAN NOTE 
   blocking publication on CRAN due to test code exceeding 5s on some very 
   slow CRAN debian servers 

## Version 2.0.1

 - Reduced image filesizes to address CRAN NOTE

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
