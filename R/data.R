#' Nightingale biomarker information
#'
#' @description {
#' Contains details on the Nightingale biomarkers available in UK Biobank and
#' computed by this package.
#'
#' If you use the estimated fatty acids in lipoprotein measurements please cite:
#' Belkadi A. \emph{et al.} Size-resolved lipoprotein fatty acid content as a
#' novel nuclear magnetic resonance-derived trait specifically associates with
#' genetic variants that control fatty acid metabolism. \emph{J. Proteome Res.}
#'  (2026) doi: \href{https://pubs.acs.org/doi/10.1021/acs.jproteome.6c00107}{10.1021/acs.jproteome.6c00107}
#' }
#'
#' @format A data table with 344 rows and 8 columns:
#' \describe{
#'   \item{Biomarker}{Short column name assigned to the biomarker}
#'   \item{Description}{Biomarker description, matching the description field provided by UK Biobank and Nightingale Health}
#'   \item{Units}{Units of measurement for the biomarker ("mmol/L", "g/L", "nm", "degree", "ratio", or "\%")}
#'   \item{Type}{Biomarker type ("Non-derived", "Composite", "Ratio", "Percentage", or "Estimate")}
#'   \item{Group}{Biomarker group as provided by Nightingale Health}
#'   \item{Sub.Group}{Biomarker sub-group as provided by Nightingale Health}
#'   \item{Nightingale}{Logical. Indicates biomarker is quantified by the Nightingale Health platform}
#'   \item{UKB.Field.ID}{Field ID in UK Biobank, see \url{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}}
#'   \item{QC.Flag.Field.ID}{Field ID in UK Biobank for the biomarker QC Flags, see \url{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=220}}
#'   \item{Full.Formula}{For composite biomarkers and ratios, details formula through which the biomarker can be derived from the 107 non-derived biomarkers}
#'   \item{Simplified.Formula}{Simplified form of the full formula most clearly expressing how each composite biomarker and ratio can be rederived from other biomarkers}
#' }
"nmr_info"

#' Nightingale biomarker sample processing information
#'
#' Contains details on the sample processing and quality control information for
#' the NMR biomarker data in UK Biobank.
#'
#' @format A data table with 18 rows and 3 columns:
#' \describe{
#'   \item{Name}{Column name assigned to the sample information field}
#'   \item{Description}{Brief description of the field contents. Further details
#'   on the Nightingale sample QC columns can be found in the \href{https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=3004}{accompanying resource on the UK Biobank showcase}.}
#'   \item{UKB.Field.ID}{Field ID in UK Biobank, see \url{https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=222}.
#'   Rows missing \code{UKB.Field.ID} entries correspond to additional sample
#'   processing information derived from these fields and returned by \code{\link{remove_technical_variation}}.}
#' }
"sample_qc_info"

#' Data for testing package functions
#'
#' Dataset mimicking structure of decoded UK Biobank dataset of NMR metabolomics
#' biomarker concentrations and associated processing variables for testing
#' package functions.
#'
#' @format A data table with 50 rows and 735 columns with column names "eid"
#' followed by extracted UK Biobank field data of the format "p23649_i0",
#' "p23649_i1", \dots, "p23655_i1".
#'
#' @source Data in each column has been randomly drawn from the distribution
#' present in the UK Biobank dataset. Importantly, random sampling was performed
#' for each column separately, thus no rows represent real observations or
#' participants in UK Biobank.
"test_data"
