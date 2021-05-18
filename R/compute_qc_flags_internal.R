# Internal functions for aggregating QC Flags when rederiving biomarkers

nightingale_ratio_flags <- function(x) {
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

  return(x)
}

nightingale_composite_biomarker_flags <- function(x) {
  # Silence CRAN NOTES about undefined global variables (columns in the scope of x)
  XL_VLDL_C <- HDL_C <- HDL_CE <- HDL_FC <- HDL_L <- HDL_P <- HDL_PL <- HDL_TG <-
    IDL_C <- IDL_CE <- IDL_FC <- IDL_L <- IDL_P <- IDL_PL <- IDL_TG <- Ile <-
    L_HDL_C <- L_HDL_CE <- L_HDL_FC <- L_HDL_L <- L_HDL_P <- L_HDL_PL <- L_HDL_TG <-
    L_LDL_C <- L_LDL_CE <- L_LDL_FC <- L_LDL_L <- L_LDL_P <- L_LDL_PL <- L_LDL_TG <-
    L_VLDL_C <- L_VLDL_CE <- L_VLDL_FC <- L_VLDL_L <- L_VLDL_P <- L_VLDL_PL <-
    L_VLDL_TG <- LDL_C <- LDL_CE <- LDL_FC <- LDL_L <- LDL_P <- LDL_PL <- LDL_TG <-
    Leu <- M_HDL_C <- M_HDL_CE <- M_HDL_FC <- M_HDL_L <- M_HDL_P <- M_HDL_PL <-
    M_HDL_TG <- M_LDL_C <- M_LDL_CE <- M_LDL_FC <- M_LDL_L <- M_LDL_P <- M_LDL_PL <-
    M_LDL_TG <- M_VLDL_C <- M_VLDL_CE <- M_VLDL_FC <- M_VLDL_L <- M_VLDL_P <-
    M_VLDL_PL <- M_VLDL_TG <- MUFA <- non_HDL_C <- Omega_3 <- Omega_6 <- PUFA <-
    Remnant_C <- S_HDL_C <- S_HDL_CE <- S_HDL_FC <- S_HDL_L <- S_HDL_P <- S_HDL_PL <-
    S_HDL_TG <- S_LDL_C <- S_LDL_CE <- S_LDL_FC <- S_LDL_L <- S_LDL_P <- S_LDL_PL <-
    S_LDL_TG <- S_VLDL_C <- S_VLDL_CE <- S_VLDL_FC <- S_VLDL_L <- S_VLDL_P <-
    S_VLDL_PL <- S_VLDL_TG <- SFA <- Total_BCAA <- Total_C <- Total_CE <- Total_FA <-
    Total_FC <- Total_L <- Total_P <- Total_PL <- Total_TG <- Val <- VLDL_C <-
    VLDL_CE <- VLDL_FC <- VLDL_L <- VLDL_P <- VLDL_PL <- VLDL_TG <- XL_HDL_C <-
    XL_HDL_CE <- XL_HDL_FC <- XL_HDL_L <- XL_HDL_P <- XL_HDL_PL <- XL_HDL_TG <-
    XL_VLDL_CE <- XL_VLDL_FC <- XL_VLDL_L <- XL_VLDL_P <- XL_VLDL_PL <- XL_VLDL_TG <-
    XS_VLDL_C <- XS_VLDL_CE <- XS_VLDL_FC <- XS_VLDL_L <- XS_VLDL_P <- XS_VLDL_PL <-
    XS_VLDL_TG <- XXL_VLDL_C <- XXL_VLDL_CE <- XXL_VLDL_FC <- XXL_VLDL_L <-
    XXL_VLDL_P <- XXL_VLDL_PL <- XXL_VLDL_TG <- NULL

  # First compute total cholesterol in 14 lipoprotein subclasses: these are the
  # sums of free cholesterol and cholesteryl esters.
  tryAssign(x[, XXL_VLDL_C := collate_flags(XXL_VLDL_CE, XXL_VLDL_FC)])
  tryAssign(x[, XL_VLDL_C := collate_flags(XL_VLDL_CE, XL_VLDL_FC)])
  tryAssign(x[, L_VLDL_C := collate_flags(L_VLDL_CE, L_VLDL_FC)])
  tryAssign(x[, M_VLDL_C := collate_flags(M_VLDL_CE, M_VLDL_FC)])
  tryAssign(x[, S_VLDL_C := collate_flags(S_VLDL_CE, S_VLDL_FC)])
  tryAssign(x[, XS_VLDL_C := collate_flags(XS_VLDL_CE, XS_VLDL_FC)])
  tryAssign(x[, IDL_C := collate_flags(IDL_CE, IDL_FC)])
  tryAssign(x[, L_LDL_C := collate_flags(L_LDL_CE, L_LDL_FC)])
  tryAssign(x[, M_LDL_C := collate_flags(M_LDL_CE, M_LDL_FC)])
  tryAssign(x[, S_LDL_C := collate_flags(S_LDL_CE, S_LDL_FC)])
  tryAssign(x[, XL_HDL_C := collate_flags(XL_HDL_CE, XL_HDL_FC)])
  tryAssign(x[, L_HDL_C := collate_flags(L_HDL_CE, L_HDL_FC)])
  tryAssign(x[, M_HDL_C := collate_flags(M_HDL_CE, M_HDL_FC)])
  tryAssign(x[, S_HDL_C := collate_flags(S_HDL_CE, S_HDL_FC)])

  # Now compute total lipids in lipoprotein subclasses
  # Total Lipids = cholesterol, phospholipids, triglycerides
  tryAssign(x[, XXL_VLDL_L := collate_flags(XXL_VLDL_C, XXL_VLDL_PL, XXL_VLDL_TG)])
  tryAssign(x[, XL_VLDL_L := collate_flags(XL_VLDL_C, XL_VLDL_PL, XL_VLDL_TG)])
  tryAssign(x[, L_VLDL_L := collate_flags(L_VLDL_C, L_VLDL_PL, L_VLDL_TG)])
  tryAssign(x[, M_VLDL_L := collate_flags(M_VLDL_C, M_VLDL_PL, M_VLDL_TG)])
  tryAssign(x[, S_VLDL_L := collate_flags(S_VLDL_C, S_VLDL_PL, S_VLDL_TG)])
  tryAssign(x[, XS_VLDL_L := collate_flags(XS_VLDL_C, XS_VLDL_PL, XS_VLDL_TG)])
  tryAssign(x[, IDL_L := collate_flags(IDL_C, IDL_PL, IDL_TG)])
  tryAssign(x[, L_LDL_L := collate_flags(L_LDL_C, L_LDL_PL, L_LDL_TG)])
  tryAssign(x[, M_LDL_L := collate_flags(M_LDL_C, M_LDL_PL, M_LDL_TG)])
  tryAssign(x[, S_LDL_L := collate_flags(S_LDL_C, S_LDL_PL, S_LDL_TG)])
  tryAssign(x[, XL_HDL_L := collate_flags(XL_HDL_C, XL_HDL_PL, XL_HDL_TG)])
  tryAssign(x[, L_HDL_L := collate_flags(L_HDL_C, L_HDL_PL, L_HDL_TG)])
  tryAssign(x[, M_HDL_L := collate_flags(M_HDL_C, M_HDL_PL, M_HDL_TG)])
  tryAssign(x[, S_HDL_L := collate_flags(S_HDL_C, S_HDL_PL, S_HDL_TG)])

  # Now compute totals for lipoprotein classes.
  # Eg. free cholesterol in LDL (LDL_FC) = sum(free choleseterol in LDL of different sizes).
  tryAssign(x[, VLDL_CE := collate_flags(XXL_VLDL_CE, XL_VLDL_CE, L_VLDL_CE, M_VLDL_CE, S_VLDL_CE, XS_VLDL_CE)])
  tryAssign(x[, VLDL_FC := collate_flags(XXL_VLDL_FC, XL_VLDL_FC, L_VLDL_FC, M_VLDL_FC, S_VLDL_FC, XS_VLDL_FC)])
  tryAssign(x[, VLDL_C := collate_flags(XXL_VLDL_C, XL_VLDL_C, L_VLDL_C, M_VLDL_C, S_VLDL_C, XS_VLDL_C)])
  tryAssign(x[, VLDL_PL := collate_flags(XXL_VLDL_PL, XL_VLDL_PL, L_VLDL_PL, M_VLDL_PL, S_VLDL_PL, XS_VLDL_PL)])
  tryAssign(x[, VLDL_TG := collate_flags(XXL_VLDL_TG, XL_VLDL_TG, L_VLDL_TG, M_VLDL_TG, S_VLDL_TG, XS_VLDL_TG)])
  tryAssign(x[, VLDL_L := collate_flags(XXL_VLDL_L, XL_VLDL_L, L_VLDL_L, M_VLDL_L, S_VLDL_L, XS_VLDL_L)])
  tryAssign(x[, VLDL_P := collate_flags(XXL_VLDL_P, XL_VLDL_P, L_VLDL_P, M_VLDL_P, S_VLDL_P, XS_VLDL_P)])

  tryAssign(x[, LDL_CE := collate_flags(L_LDL_CE, M_LDL_CE, S_LDL_CE)])
  tryAssign(x[, LDL_FC := collate_flags(L_LDL_FC, M_LDL_FC, S_LDL_FC)])
  tryAssign(x[, LDL_C := collate_flags(L_LDL_C, M_LDL_C, S_LDL_C)])
  tryAssign(x[, LDL_PL := collate_flags(L_LDL_PL, M_LDL_PL, S_LDL_PL)])
  tryAssign(x[, LDL_TG := collate_flags(L_LDL_TG, M_LDL_TG, S_LDL_TG)])
  tryAssign(x[, LDL_L := collate_flags(L_LDL_L, M_LDL_L, S_LDL_L)])
  tryAssign(x[, LDL_P := collate_flags(L_LDL_P, M_LDL_P, S_LDL_P)])

  tryAssign(x[, HDL_CE := collate_flags(XL_HDL_CE, L_HDL_CE, M_HDL_CE, S_HDL_CE)])
  tryAssign(x[, HDL_FC := collate_flags(XL_HDL_FC, L_HDL_FC, M_HDL_FC, S_HDL_FC)])
  tryAssign(x[, HDL_C := collate_flags(XL_HDL_C, L_HDL_C, M_HDL_C, S_HDL_C)])
  tryAssign(x[, HDL_PL := collate_flags(XL_HDL_PL, L_HDL_PL, M_HDL_PL, S_HDL_PL)])
  tryAssign(x[, HDL_TG := collate_flags(XL_HDL_TG, L_HDL_TG, M_HDL_TG, S_HDL_TG)])
  tryAssign(x[, HDL_L := collate_flags(XL_HDL_L, L_HDL_L, M_HDL_L, S_HDL_L)])
  tryAssign(x[, HDL_P := collate_flags(XL_HDL_P, L_HDL_P, M_HDL_P, S_HDL_P)])

  # Next compute serum totals for lipids
  tryAssign(x[, Total_CE := collate_flags(VLDL_CE, LDL_CE, IDL_CE, HDL_CE)])
  tryAssign(x[, Total_FC := collate_flags(VLDL_FC, LDL_FC, IDL_FC, HDL_FC)])
  tryAssign(x[, Total_C := collate_flags(VLDL_C, LDL_C, IDL_C, HDL_C)])
  tryAssign(x[, Total_PL := collate_flags(VLDL_PL, LDL_PL, IDL_PL, HDL_PL)])
  tryAssign(x[, Total_TG := collate_flags(VLDL_TG, LDL_TG, IDL_TG, HDL_TG)])
  tryAssign(x[, Total_L := collate_flags(VLDL_L, LDL_L, IDL_L, HDL_L)])
  tryAssign(x[, Total_P := collate_flags(VLDL_P, LDL_P, IDL_P, HDL_P)])

  # Finally miscellaneous composite biomarkers
  tryAssign(x[, PUFA := collate_flags(Omega_3, Omega_6)])
  tryAssign(x[, Total_FA := collate_flags(PUFA, MUFA, SFA)])

  # Miscellaneous composite markers
  tryAssign(x[, Total_BCAA := collate_flags(Leu, Ile, Val)]) # total branched chain amino acids
  tryAssign(x[, non_HDL_C := collate_flags(Total_C, HDL_C)]) # non HDL cholesterol
  tryAssign(x[, Remnant_C := collate_flags(Total_C, HDL_C, LDL_C)]) # remnant cholesterol

  # Finished
  return(x)
}

extended_ratios_flags <- function(x) {
  # Silence CRAN NOTES about undefined global variables (columns in the scope of x)
  HDL_C <- HDL_C_pct <- HDL_CE <- HDL_CE_pct <- HDL_CE_pct_C <- HDL_FC <- HDL_FC_by_CE <-
    HDL_FC_pct <- HDL_FC_pct_C <- HDL_L <- HDL_PL <- HDL_PL_pct <- HDL_TG <-
    HDL_TG_pct <- IDL_C <- IDL_CE <- IDL_CE_pct_C <- IDL_FC <- IDL_FC_by_CE <-
    IDL_FC_pct_C <- L_HDL_C <- L_HDL_CE <- L_HDL_CE_pct_C <- L_HDL_FC <- L_HDL_FC_by_CE <-
    L_HDL_FC_pct_C <- L_LDL_C <- L_LDL_CE <- L_LDL_CE_pct_C <- L_LDL_FC <-
    L_LDL_FC_by_CE <- L_LDL_FC_pct_C <- L_VLDL_C <- L_VLDL_CE <- L_VLDL_CE_pct_C <-
    L_VLDL_FC <- L_VLDL_FC_by_CE <- L_VLDL_FC_pct_C <- LDL_C <- LDL_C_pct <- LDL_CE <-
    LDL_CE_pct <- LDL_CE_pct_C <- LDL_FC <- LDL_FC_by_CE <- LDL_FC_pct <- LDL_FC_pct_C <-
    LDL_L <- LDL_PL <- LDL_PL_pct <- LDL_TG <- LDL_TG_pct <- M_HDL_C <- M_HDL_CE <-
    M_HDL_CE_pct_C <- M_HDL_FC <- M_HDL_FC_by_CE <- M_HDL_FC_pct_C <- M_LDL_C <-
    M_LDL_CE <- M_LDL_CE_pct_C <- M_LDL_FC <- M_LDL_FC_by_CE <- M_LDL_FC_pct_C <-
    M_VLDL_C <- M_VLDL_CE <- M_VLDL_CE_pct_C <- M_VLDL_FC <- M_VLDL_FC_by_CE <-
    M_VLDL_FC_pct_C <- Omega_3 <- Omega_3_pct_PUFA <- Omega_6 <- Omega_6_pct_PUFA <-
    PUFA <- S_HDL_C <- S_HDL_CE <- S_HDL_CE_pct_C <- S_HDL_FC <- S_HDL_FC_by_CE <-
    S_HDL_FC_pct_C <- S_LDL_C <- S_LDL_CE <- S_LDL_CE_pct_C <- S_LDL_FC <-
    S_LDL_FC_by_CE <- S_LDL_FC_pct_C <- S_VLDL_C <- S_VLDL_CE <- S_VLDL_CE_pct_C <-
    S_VLDL_FC <- S_VLDL_FC_by_CE <- S_VLDL_FC_pct_C <- Total_C <- Total_C_pct <-
    Total_CE <- Total_CE_pct <- Total_CE_pct_C <- Total_FC <- Total_FC_by_CE <-
    Total_FC_pct <- Total_FC_pct_C <- Total_L <- Total_PL <- Total_PL_pct <-
    Total_TG <- Total_TG_pct <- VLDL_C <- VLDL_C_pct <- VLDL_CE <- VLDL_CE_pct <-
    VLDL_CE_pct_C <- VLDL_FC <- VLDL_FC_by_CE <- VLDL_FC_pct <- VLDL_FC_pct_C <-
    VLDL_L <- VLDL_PL <- VLDL_PL_pct <- VLDL_TG <- VLDL_TG_pct <- XL_HDL_C <-
    XL_HDL_CE <- XL_HDL_CE_pct_C <- XL_HDL_FC <- XL_HDL_FC_by_CE <- XL_HDL_FC_pct_C <-
    XL_VLDL_C <- XL_VLDL_CE <- XL_VLDL_CE_pct_C <- XL_VLDL_FC <- XL_VLDL_FC_by_CE <-
    XL_VLDL_FC_pct_C <- XS_VLDL_C <- XS_VLDL_CE <- XS_VLDL_CE_pct_C <- XS_VLDL_FC <-
    XS_VLDL_FC_by_CE <- XS_VLDL_FC_pct_C <- XXL_VLDL_C <- XXL_VLDL_CE <-
    XXL_VLDL_CE_pct_C <- XXL_VLDL_FC <- XXL_VLDL_FC_by_CE <-
    XXL_VLDL_FC_pct_C <- NULL

  # Compute percentage of cholesterol from free vs. esterified cholesterol
  tryAssign(x[, XXL_VLDL_FC_pct_C := collate_flags(XXL_VLDL_FC, XXL_VLDL_C)])
  tryAssign(x[, XL_VLDL_FC_pct_C := collate_flags(XL_VLDL_FC, XL_VLDL_C)])
  tryAssign(x[, L_VLDL_FC_pct_C := collate_flags(L_VLDL_FC, L_VLDL_C)])
  tryAssign(x[, M_VLDL_FC_pct_C := collate_flags(M_VLDL_FC, M_VLDL_C)])
  tryAssign(x[, S_VLDL_FC_pct_C := collate_flags(S_VLDL_FC, S_VLDL_C)])
  tryAssign(x[, XS_VLDL_FC_pct_C := collate_flags(XS_VLDL_FC, XS_VLDL_C)])
  tryAssign(x[, L_LDL_FC_pct_C := collate_flags(L_LDL_FC, L_LDL_C)])
  tryAssign(x[, M_LDL_FC_pct_C := collate_flags(M_LDL_FC, M_LDL_C)])
  tryAssign(x[, S_LDL_FC_pct_C := collate_flags(S_LDL_FC, S_LDL_C)])
  tryAssign(x[, IDL_FC_pct_C := collate_flags(IDL_FC, IDL_C)])
  tryAssign(x[, XL_HDL_FC_pct_C := collate_flags(XL_HDL_FC, XL_HDL_C)])
  tryAssign(x[, L_HDL_FC_pct_C := collate_flags(L_HDL_FC, L_HDL_C)])
  tryAssign(x[, M_HDL_FC_pct_C := collate_flags(M_HDL_FC, M_HDL_C)])
  tryAssign(x[, S_HDL_FC_pct_C := collate_flags(S_HDL_FC, S_HDL_C)])

  tryAssign(x[, XXL_VLDL_CE_pct_C := collate_flags(XXL_VLDL_CE, XXL_VLDL_C)])
  tryAssign(x[, XL_VLDL_CE_pct_C := collate_flags(XL_VLDL_CE, XL_VLDL_C)])
  tryAssign(x[, L_VLDL_CE_pct_C := collate_flags(L_VLDL_CE, L_VLDL_C)])
  tryAssign(x[, M_VLDL_CE_pct_C := collate_flags(M_VLDL_CE, M_VLDL_C)])
  tryAssign(x[, S_VLDL_CE_pct_C := collate_flags(S_VLDL_CE, S_VLDL_C)])
  tryAssign(x[, XS_VLDL_CE_pct_C := collate_flags(XS_VLDL_CE, XS_VLDL_C)])
  tryAssign(x[, L_LDL_CE_pct_C := collate_flags(L_LDL_CE, L_LDL_C)])
  tryAssign(x[, M_LDL_CE_pct_C := collate_flags(M_LDL_CE, M_LDL_C)])
  tryAssign(x[, S_LDL_CE_pct_C := collate_flags(S_LDL_CE, S_LDL_C)])
  tryAssign(x[, IDL_CE_pct_C := collate_flags(IDL_CE, IDL_C)])
  tryAssign(x[, XL_HDL_CE_pct_C := collate_flags(XL_HDL_CE, XL_HDL_C)])
  tryAssign(x[, L_HDL_CE_pct_C := collate_flags(L_HDL_CE, L_HDL_C)])
  tryAssign(x[, M_HDL_CE_pct_C := collate_flags(M_HDL_CE, M_HDL_C)])
  tryAssign(x[, S_HDL_CE_pct_C := collate_flags(S_HDL_CE, S_HDL_C)])

  # Compute ratio of free cholesterol to esterified cholesterol
  tryAssign(x[, XXL_VLDL_FC_by_CE := collate_flags(XXL_VLDL_FC, XXL_VLDL_CE)])
  tryAssign(x[, XL_VLDL_FC_by_CE := collate_flags(XL_VLDL_FC, XL_VLDL_CE)])
  tryAssign(x[, L_VLDL_FC_by_CE := collate_flags(L_VLDL_FC, L_VLDL_CE)])
  tryAssign(x[, M_VLDL_FC_by_CE := collate_flags(M_VLDL_FC, M_VLDL_CE)])
  tryAssign(x[, S_VLDL_FC_by_CE := collate_flags(S_VLDL_FC, S_VLDL_CE)])
  tryAssign(x[, XS_VLDL_FC_by_CE := collate_flags(XS_VLDL_FC, XS_VLDL_CE)])
  tryAssign(x[, L_LDL_FC_by_CE := collate_flags(L_LDL_FC, L_LDL_CE)])
  tryAssign(x[, M_LDL_FC_by_CE := collate_flags(M_LDL_FC, M_LDL_CE)])
  tryAssign(x[, S_LDL_FC_by_CE := collate_flags(S_LDL_FC, S_LDL_CE)])
  tryAssign(x[, IDL_FC_by_CE := collate_flags(IDL_FC, IDL_CE)])
  tryAssign(x[, XL_HDL_FC_by_CE := collate_flags(XL_HDL_FC, XL_HDL_CE)])
  tryAssign(x[, L_HDL_FC_by_CE := collate_flags(L_HDL_FC, L_HDL_CE)])
  tryAssign(x[, M_HDL_FC_by_CE := collate_flags(M_HDL_FC, M_HDL_CE)])
  tryAssign(x[, S_HDL_FC_by_CE := collate_flags(S_HDL_FC, S_HDL_CE)])

  # Compute lipid percentages in lipoprotein classes
  tryAssign(x[, VLDL_CE_pct := collate_flags(VLDL_CE, VLDL_L)])
  tryAssign(x[, VLDL_FC_pct := collate_flags(VLDL_FC, VLDL_L)])
  tryAssign(x[, VLDL_C_pct := collate_flags(VLDL_C, VLDL_L)])
  tryAssign(x[, VLDL_PL_pct := collate_flags(VLDL_PL, VLDL_L)])
  tryAssign(x[, VLDL_TG_pct := collate_flags(VLDL_TG, VLDL_L)])

  tryAssign(x[, LDL_CE_pct := collate_flags(LDL_CE, LDL_L)])
  tryAssign(x[, LDL_FC_pct := collate_flags(LDL_FC, LDL_L)])
  tryAssign(x[, LDL_C_pct := collate_flags(LDL_C, LDL_L)])
  tryAssign(x[, LDL_PL_pct := collate_flags(LDL_PL, LDL_L)])
  tryAssign(x[, LDL_TG_pct := collate_flags(LDL_TG, LDL_L)])

  tryAssign(x[, HDL_CE_pct := collate_flags(HDL_CE, HDL_L)])
  tryAssign(x[, HDL_FC_pct := collate_flags(HDL_FC, HDL_L)])
  tryAssign(x[, HDL_C_pct := collate_flags(HDL_C, HDL_L)])
  tryAssign(x[, HDL_PL_pct := collate_flags(HDL_PL, HDL_L)])
  tryAssign(x[, HDL_TG_pct := collate_flags(HDL_TG, HDL_L)])

  # Compute cholesterol percentages
  tryAssign(x[, VLDL_FC_pct_C := collate_flags(VLDL_FC, VLDL_C)])
  tryAssign(x[, LDL_FC_pct_C := collate_flags(LDL_FC, LDL_C)])
  tryAssign(x[, HDL_FC_pct_C := collate_flags(HDL_FC, HDL_C)])

  tryAssign(x[, VLDL_CE_pct_C := collate_flags(VLDL_CE, VLDL_C)])
  tryAssign(x[, LDL_CE_pct_C := collate_flags(LDL_CE, LDL_C)])
  tryAssign(x[, HDL_CE_pct_C := collate_flags(HDL_CE, HDL_C)])

  # Ratios of free cholesterol to esterified cholesterol:
  tryAssign(x[, VLDL_FC_by_CE := collate_flags(VLDL_FC, VLDL_CE)])
  tryAssign(x[, LDL_FC_by_CE := collate_flags(LDL_FC, LDL_CE)])
  tryAssign(x[, HDL_FC_by_CE := collate_flags(HDL_FC, HDL_CE)])

  # Lipid and cholesterol fractions in total serum
  tryAssign(x[, Total_CE_pct := collate_flags(Total_CE, Total_L)])
  tryAssign(x[, Total_FC_pct := collate_flags(Total_FC, Total_L)])
  tryAssign(x[, Total_C_pct := collate_flags(Total_C, Total_L)])
  tryAssign(x[, Total_PL_pct := collate_flags(Total_PL, Total_L)])
  tryAssign(x[, Total_TG_pct := collate_flags(Total_TG, Total_L)])

  tryAssign(x[, Total_FC_pct_C := collate_flags(Total_FC, Total_C)])
  tryAssign(x[, Total_CE_pct_C := collate_flags(Total_CE, Total_C)])
  tryAssign(x[, Total_FC_by_CE := collate_flags(Total_FC, Total_CE)])

  # Omega 3 and Omega 6 percent of PUFA
  tryAssign(x[, Omega_3_pct_PUFA := collate_flags(Omega_3, PUFA)])
  tryAssign(x[, Omega_6_pct_PUFA := collate_flags(Omega_6, PUFA)])

  # Finished
  return(x)
}
