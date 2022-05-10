# NEWS

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
