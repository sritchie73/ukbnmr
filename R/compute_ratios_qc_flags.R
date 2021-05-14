#' Compute QC Flags for 81 biomarker ratios on the Nightingale platform
#'
#' For the 81 biomarker ratios computed by \code{\link{compute_nightingale_ratios}()},
#' aggregates the biomarker QC flags from the biomarkers composing each ratio
#' (see \code{\link{nmr_info}}), which can be useful for determining the reason
#' underlying missing values in the biomarker ratios.
#'
#' @details
#' Data sets extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' have one row per UKB biobank participant whose project specific sample
#' identifier is given in the first column named "eid". Columns following this
#' have the format "<field_id>-<instance>.<array_index>", "where here <field_id>
#' corresponds to a biomarker of interest, "e.g. 23474 for 3-Hydroxybutyrate,
#' <instance> corresponds to the assessment time point, "e.g. 0 for baseline
#' assessment, "1 for first repeat visit, "and <array_index> gives a number for
#' repeated measurements at the same time point.
#'
#' In the returned \code{data.frame} there is single column for each biomarker,
#' with additional columns for the instance and array index. Rows are uniquely
#' identifiable by the combination of entries in columns "eid", ""instance",
#' and "array index".
#'
#' Input data may be (1) a raw dataset extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' and loaded into R, "(2) a raw dataset extracted by
#' \href{https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide}{ukbconv}
#' and loaded into R using the \code{ukbtools} R package, "or (3) a processed
#' \code{data.frame} of biomarker data obtained from the above using \code{\link{extract_biomarkers}()}.
#'
#' A \code{data.table} will be returned instead of a \code{data.frame} if the
#' the user has loaded the package into their R session.
#'
#' @param x \code{data.frame} containing NMR metabolomics data from UK Biobank.
#'
#' @return a \code{data.frame} or \code{data.table} with column names "eid",
#'        "instance", "and "array_index", "followed by columns for each of the
#'         249 biomarkers and ratios.
#'
#' @export
compute_nightingale_ratio_qc_flags <- function(x) {
  # Silence CRAN NOTES about undefined global variables (columns in the scope of x)
  ApoA1 <- ApoB <- ApoB_by_ApoA1 <- DHA <- DHA_pct <- IDL_C <- IDL_C_pct <-
    IDL_CE <- IDL_CE_pct <- IDL_FC <- IDL_FC_pct <- IDL_L <- IDL_PL <- IDL_PL_pct <-
    IDL_TG <- IDL_TG_pct <- L_HDL_C <- L_HDL_C_pct <- L_HDL_CE <- L_HDL_CE_pct <-
    L_HDL_FC <- L_HDL_FC_pct <- L_HDL_L <- L_HDL_PL <- L_HDL_PL_pct <- L_HDL_TG <-
    L_HDL_TG_pct <- L_LDL_C <- L_LDL_C_pct <- L_LDL_CE <- L_LDL_CE_pct <- L_LDL_FC <-
    L_LDL_FC_pct <- L_LDL_L <- L_LDL_PL <- L_LDL_PL_pct <- L_LDL_TG <- L_LDL_TG_pct <-
    L_VLDL_C <- L_VLDL_C_pct <- L_VLDL_CE <- L_VLDL_CE_pct <- L_VLDL_FC <- L_VLDL_FC_pct <-
    L_VLDL_L <- L_VLDL_PL <- L_VLDL_PL_pct <- L_VLDL_TG <- L_VLDL_TG_pct <- LA <-
    LA_pct <- M_HDL_C <- M_HDL_C_pct <- M_HDL_CE <- M_HDL_CE_pct <- M_HDL_FC <-
    M_HDL_FC_pct <- M_HDL_L <- M_HDL_PL <- M_HDL_PL_pct <- M_HDL_TG <- M_HDL_TG_pct <-
    M_LDL_C <- M_LDL_C_pct <- M_LDL_CE <- M_LDL_CE_pct <- M_LDL_FC <- M_LDL_FC_pct <-
    M_LDL_L <- M_LDL_PL <- M_LDL_PL_pct <- M_LDL_TG <- M_LDL_TG_pct <- M_VLDL_C <-
    M_VLDL_C_pct <- M_VLDL_CE <- M_VLDL_CE_pct <- M_VLDL_FC <- M_VLDL_FC_pct <-
    M_VLDL_L <- M_VLDL_PL <- M_VLDL_PL_pct <- M_VLDL_TG <- M_VLDL_TG_pct <- MUFA <-
    MUFA_pct <- Omega_3 <- Omega_3_pct <- Omega_6 <- Omega_6_by_Omega_3 <- Omega_6_pct <-
    Phosphoglyc <- PUFA <- PUFA_by_MUFA <- PUFA_pct <- S_HDL_C <- S_HDL_C_pct <-
    S_HDL_CE <- S_HDL_CE_pct <- S_HDL_FC <- S_HDL_FC_pct <- S_HDL_L <- S_HDL_PL <-
    S_HDL_PL_pct <- S_HDL_TG <- S_HDL_TG_pct <- S_LDL_C <- S_LDL_C_pct <- S_LDL_CE <-
    S_LDL_CE_pct <- S_LDL_FC <- S_LDL_FC_pct <- S_LDL_L <- S_LDL_PL <- S_LDL_PL_pct <-
    S_LDL_TG <- S_LDL_TG_pct <- S_VLDL_C <- S_VLDL_C_pct <- S_VLDL_CE <- S_VLDL_CE_pct <-
    S_VLDL_FC <- S_VLDL_FC_pct <- S_VLDL_L <- S_VLDL_PL <- S_VLDL_PL_pct <- S_VLDL_TG <-
    S_VLDL_TG_pct <- SFA <- SFA_pct <- TG_by_PG <- Total_FA <- Total_TG <- XL_HDL_C <-
    XL_HDL_C_pct <- XL_HDL_CE <- XL_HDL_CE_pct <- XL_HDL_FC <- XL_HDL_FC_pct <-
    XL_HDL_L <- XL_HDL_PL <- XL_HDL_PL_pct <- XL_HDL_TG <- XL_HDL_TG_pct <- XL_VLDL_C <-
    XL_VLDL_C_pct <- XL_VLDL_CE <- XL_VLDL_CE_pct <- XL_VLDL_FC <- XL_VLDL_FC_pct <-
    XL_VLDL_L <- XL_VLDL_PL <- XL_VLDL_PL_pct <- XL_VLDL_TG <- XL_VLDL_TG_pct <-
    XS_VLDL_C <- XS_VLDL_C_pct <- XS_VLDL_CE <- XS_VLDL_CE_pct <- XS_VLDL_FC <-
    XS_VLDL_FC_pct <- XS_VLDL_L <- XS_VLDL_PL <- XS_VLDL_PL_pct <- XS_VLDL_TG <-
    XS_VLDL_TG_pct <- XXL_VLDL_C <- XXL_VLDL_C_pct <- XXL_VLDL_CE <- XXL_VLDL_CE_pct <-
    XXL_VLDL_FC <- XXL_VLDL_FC_pct <- XXL_VLDL_L <- XXL_VLDL_PL <- XXL_VLDL_PL_pct <-
    XXL_VLDL_TG <- XXL_VLDL_TG_pct <- NULL


  # Process data to correct format
  x <- process_data(x, type="biomarker_qc_flags") # copy of x created if already in right format

  # For each of the 14 lipoprotein subclasses, compute percentage of total
  # lipids composed of cholesteryl esters (CE), free cholesterol (FC), total
  # cholester (C), phospholipids (PL), and triglycerides (TG).
  tryAssign(x[, XXL_VLDL_CE_pct := collate_flags(XXL_VLDL_CE, XXL_VLDL_L)])
  tryAssign(x[, XXL_VLDL_FC_pct := collate_flags(XXL_VLDL_FC, XXL_VLDL_L)])
  tryAssign(x[, XXL_VLDL_C_pct := collate_flags(XXL_VLDL_C, XXL_VLDL_L)])
  tryAssign(x[, XXL_VLDL_PL_pct := collate_flags(XXL_VLDL_PL, XXL_VLDL_L)])
  tryAssign(x[, XXL_VLDL_TG_pct := collate_flags(XXL_VLDL_TG, XXL_VLDL_L)])

  tryAssign(x[, XL_VLDL_CE_pct := collate_flags(XL_VLDL_CE, XL_VLDL_L)])
  tryAssign(x[, XL_VLDL_FC_pct := collate_flags(XL_VLDL_FC, XL_VLDL_L)])
  tryAssign(x[, XL_VLDL_C_pct := collate_flags(XL_VLDL_C, XL_VLDL_L)])
  tryAssign(x[, XL_VLDL_PL_pct := collate_flags(XL_VLDL_PL, XL_VLDL_L)])
  tryAssign(x[, XL_VLDL_TG_pct := collate_flags(XL_VLDL_TG, XL_VLDL_L)])

  tryAssign(x[, L_VLDL_CE_pct := collate_flags(L_VLDL_CE, L_VLDL_L)])
  tryAssign(x[, L_VLDL_FC_pct := collate_flags(L_VLDL_FC, L_VLDL_L)])
  tryAssign(x[, L_VLDL_C_pct := collate_flags(L_VLDL_C, L_VLDL_L)])
  tryAssign(x[, L_VLDL_PL_pct := collate_flags(L_VLDL_PL, L_VLDL_L)])
  tryAssign(x[, L_VLDL_TG_pct := collate_flags(L_VLDL_TG, L_VLDL_L)])

  tryAssign(x[, M_VLDL_CE_pct := collate_flags(M_VLDL_CE, M_VLDL_L)])
  tryAssign(x[, M_VLDL_FC_pct := collate_flags(M_VLDL_FC, M_VLDL_L)])
  tryAssign(x[, M_VLDL_C_pct := collate_flags(M_VLDL_C, M_VLDL_L)])
  tryAssign(x[, M_VLDL_PL_pct := collate_flags(M_VLDL_PL, M_VLDL_L)])
  tryAssign(x[, M_VLDL_TG_pct := collate_flags(M_VLDL_TG, M_VLDL_L)])

  tryAssign(x[, S_VLDL_CE_pct := collate_flags(S_VLDL_CE, S_VLDL_L)])
  tryAssign(x[, S_VLDL_FC_pct := collate_flags(S_VLDL_FC, S_VLDL_L)])
  tryAssign(x[, S_VLDL_C_pct := collate_flags(S_VLDL_C, S_VLDL_L)])
  tryAssign(x[, S_VLDL_PL_pct := collate_flags(S_VLDL_PL, S_VLDL_L)])
  tryAssign(x[, S_VLDL_TG_pct := collate_flags(S_VLDL_TG, S_VLDL_L)])

  tryAssign(x[, XS_VLDL_CE_pct := collate_flags(XS_VLDL_CE, XS_VLDL_L)])
  tryAssign(x[, XS_VLDL_FC_pct := collate_flags(XS_VLDL_FC, XS_VLDL_L)])
  tryAssign(x[, XS_VLDL_C_pct := collate_flags(XS_VLDL_C, XS_VLDL_L)])
  tryAssign(x[, XS_VLDL_PL_pct := collate_flags(XS_VLDL_PL, XS_VLDL_L)])
  tryAssign(x[, XS_VLDL_TG_pct := collate_flags(XS_VLDL_TG, XS_VLDL_L)])

  tryAssign(x[, L_LDL_CE_pct := collate_flags(L_LDL_CE, L_LDL_L)])
  tryAssign(x[, L_LDL_FC_pct := collate_flags(L_LDL_FC, L_LDL_L)])
  tryAssign(x[, L_LDL_C_pct := collate_flags(L_LDL_C, L_LDL_L)])
  tryAssign(x[, L_LDL_PL_pct := collate_flags(L_LDL_PL, L_LDL_L)])
  tryAssign(x[, L_LDL_TG_pct := collate_flags(L_LDL_TG, L_LDL_L)])

  tryAssign(x[, M_LDL_CE_pct := collate_flags(M_LDL_CE, M_LDL_L)])
  tryAssign(x[, M_LDL_FC_pct := collate_flags(M_LDL_FC, M_LDL_L)])
  tryAssign(x[, M_LDL_C_pct := collate_flags(M_LDL_C, M_LDL_L)])
  tryAssign(x[, M_LDL_PL_pct := collate_flags(M_LDL_PL, M_LDL_L)])
  tryAssign(x[, M_LDL_TG_pct := collate_flags(M_LDL_TG, M_LDL_L)])

  tryAssign(x[, S_LDL_CE_pct := collate_flags(S_LDL_CE, S_LDL_L)])
  tryAssign(x[, S_LDL_FC_pct := collate_flags(S_LDL_FC, S_LDL_L)])
  tryAssign(x[, S_LDL_C_pct := collate_flags(S_LDL_C, S_LDL_L)])
  tryAssign(x[, S_LDL_PL_pct := collate_flags(S_LDL_PL, S_LDL_L)])
  tryAssign(x[, S_LDL_TG_pct := collate_flags(S_LDL_TG, S_LDL_L)])

  tryAssign(x[, IDL_CE_pct := collate_flags(IDL_CE, IDL_L)])
  tryAssign(x[, IDL_FC_pct := collate_flags(IDL_FC, IDL_L)])
  tryAssign(x[, IDL_C_pct := collate_flags(IDL_C, IDL_L)])
  tryAssign(x[, IDL_PL_pct := collate_flags(IDL_PL, IDL_L)])
  tryAssign(x[, IDL_TG_pct := collate_flags(IDL_TG, IDL_L)])

  tryAssign(x[, XL_HDL_CE_pct := collate_flags(XL_HDL_CE, XL_HDL_L)])
  tryAssign(x[, XL_HDL_FC_pct := collate_flags(XL_HDL_FC, XL_HDL_L)])
  tryAssign(x[, XL_HDL_C_pct := collate_flags(XL_HDL_C, XL_HDL_L)])
  tryAssign(x[, XL_HDL_PL_pct := collate_flags(XL_HDL_PL, XL_HDL_L)])
  tryAssign(x[, XL_HDL_TG_pct := collate_flags(XL_HDL_TG, XL_HDL_L)])

  tryAssign(x[, L_HDL_CE_pct := collate_flags(L_HDL_CE, L_HDL_L)])
  tryAssign(x[, L_HDL_FC_pct := collate_flags(L_HDL_FC, L_HDL_L)])
  tryAssign(x[, L_HDL_C_pct := collate_flags(L_HDL_C, L_HDL_L)])
  tryAssign(x[, L_HDL_PL_pct := collate_flags(L_HDL_PL, L_HDL_L)])
  tryAssign(x[, L_HDL_TG_pct := collate_flags(L_HDL_TG, L_HDL_L)])

  tryAssign(x[, M_HDL_CE_pct := collate_flags(M_HDL_CE, M_HDL_L)])
  tryAssign(x[, M_HDL_FC_pct := collate_flags(M_HDL_FC, M_HDL_L)])
  tryAssign(x[, M_HDL_C_pct := collate_flags(M_HDL_C, M_HDL_L)])
  tryAssign(x[, M_HDL_PL_pct := collate_flags(M_HDL_PL, M_HDL_L)])
  tryAssign(x[, M_HDL_TG_pct := collate_flags(M_HDL_TG, M_HDL_L)])

  tryAssign(x[, S_HDL_CE_pct := collate_flags(S_HDL_CE, S_HDL_L)])
  tryAssign(x[, S_HDL_FC_pct := collate_flags(S_HDL_FC, S_HDL_L)])
  tryAssign(x[, S_HDL_C_pct := collate_flags(S_HDL_C, S_HDL_L)])
  tryAssign(x[, S_HDL_PL_pct := collate_flags(S_HDL_PL, S_HDL_L)])
  tryAssign(x[, S_HDL_TG_pct := collate_flags(S_HDL_TG, S_HDL_L)])

  # Compute fatty acid percentages
  tryAssign(x[, Omega_3_pct := collate_flags(Omega_3, Total_FA)])
  tryAssign(x[, Omega_6_pct := collate_flags(Omega_6, Total_FA)])
  tryAssign(x[, LA_pct := collate_flags(LA, Total_FA)])
  tryAssign(x[, MUFA_pct := collate_flags(MUFA, Total_FA)])
  tryAssign(x[, PUFA_pct := collate_flags(PUFA, Total_FA)])
  tryAssign(x[, SFA_pct := collate_flags(SFA, Total_FA)])
  tryAssign(x[, DHA_pct := collate_flags(DHA, Total_FA)])

  # Miscellaneous ratios
  tryAssign(x[, Omega_6_by_Omega_3 := collate_flags(Omega_6, Omega_3)])
  tryAssign(x[, ApoB_by_ApoA1 := collate_flags(ApoB, ApoA1)])
  tryAssign(x[, PUFA_by_MUFA := collate_flags(PUFA, MUFA)])
  tryAssign(x[, TG_by_PG := collate_flags(Total_TG, Phosphoglyc)])

  # Return
  returnDT(x)
}
