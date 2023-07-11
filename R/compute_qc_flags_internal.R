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
  if (exists("XXL_VLDL_CE", x) && exists("XXL_VLDL_L", x)) {
    x[, XXL_VLDL_CE_pct := collate_flags(XXL_VLDL_CE, XXL_VLDL_L)]
  } else if (!(exists("XXL_VLDL_CE", x) && exists("XXL_VLDL_L", x)) && exists("XXL_VLDL_CE_pct", x)) {
    x[, XXL_VLDL_CE_pct := NULL]
  }

  if (exists("XXL_VLDL_FC", x) && exists("XXL_VLDL_L", x)) {
    x[, XXL_VLDL_FC_pct := collate_flags(XXL_VLDL_FC, XXL_VLDL_L)]
  } else if (!(exists("XXL_VLDL_FC", x) && exists("XXL_VLDL_L", x)) && exists("XXL_VLDL_FC_pct", x)) {
    x[, XXL_VLDL_FC_pct := NULL]
  }

  if (exists("XXL_VLDL_C", x) && exists("XXL_VLDL_L", x)) {
    x[, XXL_VLDL_C_pct := collate_flags(XXL_VLDL_C, XXL_VLDL_L)]
  } else if (!(exists("XXL_VLDL_C", x) && exists("XXL_VLDL_L", x)) && exists("XXL_VLDL_C_pct", x)) {
    x[, XXL_VLDL_C_pct := NULL]
  }

  if (exists("XXL_VLDL_PL", x) && exists("XXL_VLDL_L", x)) {
    x[, XXL_VLDL_PL_pct := collate_flags(XXL_VLDL_PL, XXL_VLDL_L)]
  } else if (!(exists("XXL_VLDL_PL", x) && exists("XXL_VLDL_L", x)) && exists("XXL_VLDL_PL_pct", x)) {
    x[, XXL_VLDL_PL_pct := NULL]
  }

  if (exists("XXL_VLDL_TG", x) && exists("XXL_VLDL_L", x)) {
    x[, XXL_VLDL_TG_pct := collate_flags(XXL_VLDL_TG, XXL_VLDL_L)]
  } else if (!(exists("XXL_VLDL_TG", x) && exists("XXL_VLDL_L", x)) && exists("XXL_VLDL_TG_pct", x)) {
    x[, XXL_VLDL_TG_pct := NULL]
  }


  if (exists("XL_VLDL_CE", x) && exists("XL_VLDL_L", x)) {
    x[, XL_VLDL_CE_pct := collate_flags(XL_VLDL_CE, XL_VLDL_L)]
  } else if (!(exists("XL_VLDL_CE", x) && exists("XL_VLDL_L", x)) && exists("XL_VLDL_CE_pct", x)) {
    x[, XL_VLDL_CE_pct := NULL]
  }

  if (exists("XL_VLDL_FC", x) && exists("XL_VLDL_L", x)) {
    x[, XL_VLDL_FC_pct := collate_flags(XL_VLDL_FC, XL_VLDL_L)]
  } else if (!(exists("XL_VLDL_FC", x) && exists("XL_VLDL_L", x)) && exists("XL_VLDL_FC_pct", x)) {
    x[, XL_VLDL_FC_pct := NULL]
  }

  if (exists("XL_VLDL_C", x) && exists("XL_VLDL_L", x)) {
    x[, XL_VLDL_C_pct := collate_flags(XL_VLDL_C, XL_VLDL_L)]
  } else if (!(exists("XL_VLDL_C", x) && exists("XL_VLDL_L", x)) && exists("XL_VLDL_C_pct", x)) {
    x[, XL_VLDL_C_pct := NULL]
  }

  if (exists("XL_VLDL_PL", x) && exists("XL_VLDL_L", x)) {
    x[, XL_VLDL_PL_pct := collate_flags(XL_VLDL_PL, XL_VLDL_L)]
  } else if (!(exists("XL_VLDL_PL", x) && exists("XL_VLDL_L", x)) && exists("XL_VLDL_PL_pct", x)) {
    x[, XL_VLDL_PL_pct := NULL]
  }

  if (exists("XL_VLDL_TG", x) && exists("XL_VLDL_L", x)) {
    x[, XL_VLDL_TG_pct := collate_flags(XL_VLDL_TG, XL_VLDL_L)]
  } else if (!(exists("XL_VLDL_TG", x) && exists("XL_VLDL_L", x)) && exists("XL_VLDL_TG_pct", x)) {
    x[, XL_VLDL_TG_pct := NULL]
  }


  if (exists("L_VLDL_CE", x) && exists("L_VLDL_L", x)) {
    x[, L_VLDL_CE_pct := collate_flags(L_VLDL_CE, L_VLDL_L)]
  } else if (!(exists("L_VLDL_CE", x) && exists("L_VLDL_L", x)) && exists("L_VLDL_CE_pct", x)) {
    x[, L_VLDL_CE_pct := NULL]
  }

  if (exists("L_VLDL_FC", x) && exists("L_VLDL_L", x)) {
    x[, L_VLDL_FC_pct := collate_flags(L_VLDL_FC, L_VLDL_L)]
  } else if (!(exists("L_VLDL_FC", x) && exists("L_VLDL_L", x)) && exists("L_VLDL_FC_pct", x)) {
    x[, L_VLDL_FC_pct := NULL]
  }

  if (exists("L_VLDL_C", x) && exists("L_VLDL_L", x)) {
    x[, L_VLDL_C_pct := collate_flags(L_VLDL_C, L_VLDL_L)]
  } else if (!(exists("L_VLDL_C", x) && exists("L_VLDL_L", x)) && exists("L_VLDL_C_pct", x)) {
    x[, L_VLDL_C_pct := NULL]
  }

  if (exists("L_VLDL_PL", x) && exists("L_VLDL_L", x)) {
    x[, L_VLDL_PL_pct := collate_flags(L_VLDL_PL, L_VLDL_L)]
  } else if (!(exists("L_VLDL_PL", x) && exists("L_VLDL_L", x)) && exists("L_VLDL_PL_pct", x)) {
    x[, L_VLDL_PL_pct := NULL]
  }

  if (exists("L_VLDL_TG", x) && exists("L_VLDL_L", x)) {
    x[, L_VLDL_TG_pct := collate_flags(L_VLDL_TG, L_VLDL_L)]
  } else if (!(exists("L_VLDL_TG", x) && exists("L_VLDL_L", x)) && exists("L_VLDL_TG_pct", x)) {
    x[, L_VLDL_TG_pct := NULL]
  }


  if (exists("M_VLDL_CE", x) && exists("M_VLDL_L", x)) {
    x[, M_VLDL_CE_pct := collate_flags(M_VLDL_CE, M_VLDL_L)]
  } else if (!(exists("M_VLDL_CE", x) && exists("M_VLDL_L", x)) && exists("M_VLDL_CE_pct", x)) {
    x[, M_VLDL_CE_pct := NULL]
  }

  if (exists("M_VLDL_FC", x) && exists("M_VLDL_L", x)) {
    x[, M_VLDL_FC_pct := collate_flags(M_VLDL_FC, M_VLDL_L)]
  } else if (!(exists("M_VLDL_FC", x) && exists("M_VLDL_L", x)) && exists("M_VLDL_FC_pct", x)) {
    x[, M_VLDL_FC_pct := NULL]
  }

  if (exists("M_VLDL_C", x) && exists("M_VLDL_L", x)) {
    x[, M_VLDL_C_pct := collate_flags(M_VLDL_C, M_VLDL_L)]
  } else if (!(exists("M_VLDL_C", x) && exists("M_VLDL_L", x)) && exists("M_VLDL_C_pct", x)) {
    x[, M_VLDL_C_pct := NULL]
  }

  if (exists("M_VLDL_PL", x) && exists("M_VLDL_L", x)) {
    x[, M_VLDL_PL_pct := collate_flags(M_VLDL_PL, M_VLDL_L)]
  } else if (!(exists("M_VLDL_PL", x) && exists("M_VLDL_L", x)) && exists("M_VLDL_PL_pct", x)) {
    x[, M_VLDL_PL_pct := NULL]
  }

  if (exists("M_VLDL_TG", x) && exists("M_VLDL_L", x)) {
    x[, M_VLDL_TG_pct := collate_flags(M_VLDL_TG, M_VLDL_L)]
  } else if (!(exists("M_VLDL_TG", x) && exists("M_VLDL_L", x)) && exists("M_VLDL_TG_pct", x)) {
    x[, M_VLDL_TG_pct := NULL]
  }


  if (exists("S_VLDL_CE", x) && exists("S_VLDL_L", x)) {
    x[, S_VLDL_CE_pct := collate_flags(S_VLDL_CE, S_VLDL_L)]
  } else if (!(exists("S_VLDL_CE", x) && exists("S_VLDL_L", x)) && exists("S_VLDL_CE_pct", x)) {
    x[, S_VLDL_CE_pct := NULL]
  }

  if (exists("S_VLDL_FC", x) && exists("S_VLDL_L", x)) {
    x[, S_VLDL_FC_pct := collate_flags(S_VLDL_FC, S_VLDL_L)]
  } else if (!(exists("S_VLDL_FC", x) && exists("S_VLDL_L", x)) && exists("S_VLDL_FC_pct", x)) {
    x[, S_VLDL_FC_pct := NULL]
  }

  if (exists("S_VLDL_C", x) && exists("S_VLDL_L", x)) {
    x[, S_VLDL_C_pct := collate_flags(S_VLDL_C, S_VLDL_L)]
  } else if (!(exists("S_VLDL_C", x) && exists("S_VLDL_L", x)) && exists("S_VLDL_C_pct", x)) {
    x[, S_VLDL_C_pct := NULL]
  }

  if (exists("S_VLDL_PL", x) && exists("S_VLDL_L", x)) {
    x[, S_VLDL_PL_pct := collate_flags(S_VLDL_PL, S_VLDL_L)]
  } else if (!(exists("S_VLDL_PL", x) && exists("S_VLDL_L", x)) && exists("S_VLDL_PL_pct", x)) {
    x[, S_VLDL_PL_pct := NULL]
  }

  if (exists("S_VLDL_TG", x) && exists("S_VLDL_L", x)) {
    x[, S_VLDL_TG_pct := collate_flags(S_VLDL_TG, S_VLDL_L)]
  } else if (!(exists("S_VLDL_TG", x) && exists("S_VLDL_L", x)) && exists("S_VLDL_TG_pct", x)) {
    x[, S_VLDL_TG_pct := NULL]
  }


  if (exists("XS_VLDL_CE", x) && exists("XS_VLDL_L", x)) {
    x[, XS_VLDL_CE_pct := collate_flags(XS_VLDL_CE, XS_VLDL_L)]
  } else if (!(exists("XS_VLDL_CE", x) && exists("XS_VLDL_L", x)) && exists("XS_VLDL_CE_pct", x)) {
    x[, XS_VLDL_CE_pct := NULL]
  }

  if (exists("XS_VLDL_FC", x) && exists("XS_VLDL_L", x)) {
    x[, XS_VLDL_FC_pct := collate_flags(XS_VLDL_FC, XS_VLDL_L)]
  } else if (!(exists("XS_VLDL_FC", x) && exists("XS_VLDL_L", x)) && exists("XS_VLDL_FC_pct", x)) {
    x[, XS_VLDL_FC_pct := NULL]
  }

  if (exists("XS_VLDL_C", x) && exists("XS_VLDL_L", x)) {
    x[, XS_VLDL_C_pct := collate_flags(XS_VLDL_C, XS_VLDL_L)]
  } else if (!(exists("XS_VLDL_C", x) && exists("XS_VLDL_L", x)) && exists("XS_VLDL_C_pct", x)) {
    x[, XS_VLDL_C_pct := NULL]
  }

  if (exists("XS_VLDL_PL", x) && exists("XS_VLDL_L", x)) {
    x[, XS_VLDL_PL_pct := collate_flags(XS_VLDL_PL, XS_VLDL_L)]
  } else if (!(exists("XS_VLDL_PL", x) && exists("XS_VLDL_L", x)) && exists("XS_VLDL_PL_pct", x)) {
    x[, XS_VLDL_PL_pct := NULL]
  }

  if (exists("XS_VLDL_TG", x) && exists("XS_VLDL_L", x)) {
    x[, XS_VLDL_TG_pct := collate_flags(XS_VLDL_TG, XS_VLDL_L)]
  } else if (!(exists("XS_VLDL_TG", x) && exists("XS_VLDL_L", x)) && exists("XS_VLDL_TG_pct", x)) {
    x[, XS_VLDL_TG_pct := NULL]
  }


  if (exists("L_LDL_CE", x) && exists("L_LDL_L", x)) {
    x[, L_LDL_CE_pct := collate_flags(L_LDL_CE, L_LDL_L)]
  } else if (!(exists("L_LDL_CE", x) && exists("L_LDL_L", x)) && exists("L_LDL_CE_pct", x)) {
    x[, L_LDL_CE_pct := NULL]
  }

  if (exists("L_LDL_FC", x) && exists("L_LDL_L", x)) {
    x[, L_LDL_FC_pct := collate_flags(L_LDL_FC, L_LDL_L)]
  } else if (!(exists("L_LDL_FC", x) && exists("L_LDL_L", x)) && exists("L_LDL_FC_pct", x)) {
    x[, L_LDL_FC_pct := NULL]
  }

  if (exists("L_LDL_C", x) && exists("L_LDL_L", x)) {
    x[, L_LDL_C_pct := collate_flags(L_LDL_C, L_LDL_L)]
  } else if (!(exists("L_LDL_C", x) && exists("L_LDL_L", x)) && exists("L_LDL_C_pct", x)) {
    x[, L_LDL_C_pct := NULL]
  }

  if (exists("L_LDL_PL", x) && exists("L_LDL_L", x)) {
    x[, L_LDL_PL_pct := collate_flags(L_LDL_PL, L_LDL_L)]
  } else if (!(exists("L_LDL_PL", x) && exists("L_LDL_L", x)) && exists("L_LDL_PL_pct", x)) {
    x[, L_LDL_PL_pct := NULL]
  }

  if (exists("L_LDL_TG", x) && exists("L_LDL_L", x)) {
    x[, L_LDL_TG_pct := collate_flags(L_LDL_TG, L_LDL_L)]
  } else if (!(exists("L_LDL_TG", x) && exists("L_LDL_L", x)) && exists("L_LDL_TG_pct", x)) {
    x[, L_LDL_TG_pct := NULL]
  }


  if (exists("M_LDL_CE", x) && exists("M_LDL_L", x)) {
    x[, M_LDL_CE_pct := collate_flags(M_LDL_CE, M_LDL_L)]
  } else if (!(exists("M_LDL_CE", x) && exists("M_LDL_L", x)) && exists("M_LDL_CE_pct", x)) {
    x[, M_LDL_CE_pct := NULL]
  }

  if (exists("M_LDL_FC", x) && exists("M_LDL_L", x)) {
    x[, M_LDL_FC_pct := collate_flags(M_LDL_FC, M_LDL_L)]
  } else if (!(exists("M_LDL_FC", x) && exists("M_LDL_L", x)) && exists("M_LDL_FC_pct", x)) {
    x[, M_LDL_FC_pct := NULL]
  }

  if (exists("M_LDL_C", x) && exists("M_LDL_L", x)) {
    x[, M_LDL_C_pct := collate_flags(M_LDL_C, M_LDL_L)]
  } else if (!(exists("M_LDL_C", x) && exists("M_LDL_L", x)) && exists("M_LDL_C_pct", x)) {
    x[, M_LDL_C_pct := NULL]
  }

  if (exists("M_LDL_PL", x) && exists("M_LDL_L", x)) {
    x[, M_LDL_PL_pct := collate_flags(M_LDL_PL, M_LDL_L)]
  } else if (!(exists("M_LDL_PL", x) && exists("M_LDL_L", x)) && exists("M_LDL_PL_pct", x)) {
    x[, M_LDL_PL_pct := NULL]
  }

  if (exists("M_LDL_TG", x) && exists("M_LDL_L", x)) {
    x[, M_LDL_TG_pct := collate_flags(M_LDL_TG, M_LDL_L)]
  } else if (!(exists("M_LDL_TG", x) && exists("M_LDL_L", x)) && exists("M_LDL_TG_pct", x)) {
    x[, M_LDL_TG_pct := NULL]
  }


  if (exists("S_LDL_CE", x) && exists("S_LDL_L", x)) {
    x[, S_LDL_CE_pct := collate_flags(S_LDL_CE, S_LDL_L)]
  } else if (!(exists("S_LDL_CE", x) && exists("S_LDL_L", x)) && exists("S_LDL_CE_pct", x)) {
    x[, S_LDL_CE_pct := NULL]
  }

  if (exists("S_LDL_FC", x) && exists("S_LDL_L", x)) {
    x[, S_LDL_FC_pct := collate_flags(S_LDL_FC, S_LDL_L)]
  } else if (!(exists("S_LDL_FC", x) && exists("S_LDL_L", x)) && exists("S_LDL_FC_pct", x)) {
    x[, S_LDL_FC_pct := NULL]
  }

  if (exists("S_LDL_C", x) && exists("S_LDL_L", x)) {
    x[, S_LDL_C_pct := collate_flags(S_LDL_C, S_LDL_L)]
  } else if (!(exists("S_LDL_C", x) && exists("S_LDL_L", x)) && exists("S_LDL_C_pct", x)) {
    x[, S_LDL_C_pct := NULL]
  }

  if (exists("S_LDL_PL", x) && exists("S_LDL_L", x)) {
    x[, S_LDL_PL_pct := collate_flags(S_LDL_PL, S_LDL_L)]
  } else if (!(exists("S_LDL_PL", x) && exists("S_LDL_L", x)) && exists("S_LDL_PL_pct", x)) {
    x[, S_LDL_PL_pct := NULL]
  }

  if (exists("S_LDL_TG", x) && exists("S_LDL_L", x)) {
    x[, S_LDL_TG_pct := collate_flags(S_LDL_TG, S_LDL_L)]
  } else if (!(exists("S_LDL_TG", x) && exists("S_LDL_L", x)) && exists("S_LDL_TG_pct", x)) {
    x[, S_LDL_TG_pct := NULL]
  }


  if (exists("IDL_CE", x) && exists("IDL_L", x)) {
    x[, IDL_CE_pct := collate_flags(IDL_CE, IDL_L)]
  } else if (!(exists("IDL_CE", x) && exists("IDL_L", x)) && exists("IDL_CE_pct", x)) {
    x[, IDL_CE_pct := NULL]
  }

  if (exists("IDL_FC", x) && exists("IDL_L", x)) {
    x[, IDL_FC_pct := collate_flags(IDL_FC, IDL_L)]
  } else if (!(exists("IDL_FC", x) && exists("IDL_L", x)) && exists("IDL_FC_pct", x)) {
    x[, IDL_FC_pct := NULL]
  }

  if (exists("IDL_C", x) && exists("IDL_L", x)) {
    x[, IDL_C_pct := collate_flags(IDL_C, IDL_L)]
  } else if (!(exists("IDL_C", x) && exists("IDL_L", x)) && exists("IDL_C_pct", x)) {
    x[, IDL_C_pct := NULL]
  }

  if (exists("IDL_PL", x) && exists("IDL_L", x)) {
    x[, IDL_PL_pct := collate_flags(IDL_PL, IDL_L)]
  } else if (!(exists("IDL_PL", x) && exists("IDL_L", x)) && exists("IDL_PL_pct", x)) {
    x[, IDL_PL_pct := NULL]
  }

  if (exists("IDL_TG", x) && exists("IDL_L", x)) {
    x[, IDL_TG_pct := collate_flags(IDL_TG, IDL_L)]
  } else if (!(exists("IDL_TG", x) && exists("IDL_L", x)) && exists("IDL_TG_pct", x)) {
    x[, IDL_TG_pct := NULL]
  }


  if (exists("XL_HDL_CE", x) && exists("XL_HDL_L", x)) {
    x[, XL_HDL_CE_pct := collate_flags(XL_HDL_CE, XL_HDL_L)]
  } else if (!(exists("XL_HDL_CE", x) && exists("XL_HDL_L", x)) && exists("XL_HDL_CE_pct", x)) {
    x[, XL_HDL_CE_pct := NULL]
  }

  if (exists("XL_HDL_FC", x) && exists("XL_HDL_L", x)) {
    x[, XL_HDL_FC_pct := collate_flags(XL_HDL_FC, XL_HDL_L)]
  } else if (!(exists("XL_HDL_FC", x) && exists("XL_HDL_L", x)) && exists("XL_HDL_FC_pct", x)) {
    x[, XL_HDL_FC_pct := NULL]
  }

  if (exists("XL_HDL_C", x) && exists("XL_HDL_L", x)) {
    x[, XL_HDL_C_pct := collate_flags(XL_HDL_C, XL_HDL_L)]
  } else if (!(exists("XL_HDL_C", x) && exists("XL_HDL_L", x)) && exists("XL_HDL_C_pct", x)) {
    x[, XL_HDL_C_pct := NULL]
  }

  if (exists("XL_HDL_PL", x) && exists("XL_HDL_L", x)) {
    x[, XL_HDL_PL_pct := collate_flags(XL_HDL_PL, XL_HDL_L)]
  } else if (!(exists("XL_HDL_PL", x) && exists("XL_HDL_L", x)) && exists("XL_HDL_PL_pct", x)) {
    x[, XL_HDL_PL_pct := NULL]
  }

  if (exists("XL_HDL_TG", x) && exists("XL_HDL_L", x)) {
    x[, XL_HDL_TG_pct := collate_flags(XL_HDL_TG, XL_HDL_L)]
  } else if (!(exists("XL_HDL_TG", x) && exists("XL_HDL_L", x)) && exists("XL_HDL_TG_pct", x)) {
    x[, XL_HDL_TG_pct := NULL]
  }


  if (exists("L_HDL_CE", x) && exists("L_HDL_L", x)) {
    x[, L_HDL_CE_pct := collate_flags(L_HDL_CE, L_HDL_L)]
  } else if (!(exists("L_HDL_CE", x) && exists("L_HDL_L", x)) && exists("L_HDL_CE_pct", x)) {
    x[, L_HDL_CE_pct := NULL]
  }

  if (exists("L_HDL_FC", x) && exists("L_HDL_L", x)) {
    x[, L_HDL_FC_pct := collate_flags(L_HDL_FC, L_HDL_L)]
  } else if (!(exists("L_HDL_FC", x) && exists("L_HDL_L", x)) && exists("L_HDL_FC_pct", x)) {
    x[, L_HDL_FC_pct := NULL]
  }

  if (exists("L_HDL_C", x) && exists("L_HDL_L", x)) {
    x[, L_HDL_C_pct := collate_flags(L_HDL_C, L_HDL_L)]
  } else if (!(exists("L_HDL_C", x) && exists("L_HDL_L", x)) && exists("L_HDL_C_pct", x)) {
    x[, L_HDL_C_pct := NULL]
  }

  if (exists("L_HDL_PL", x) && exists("L_HDL_L", x)) {
    x[, L_HDL_PL_pct := collate_flags(L_HDL_PL, L_HDL_L)]
  } else if (!(exists("L_HDL_PL", x) && exists("L_HDL_L", x)) && exists("L_HDL_PL_pct", x)) {
    x[, L_HDL_PL_pct := NULL]
  }

  if (exists("L_HDL_TG", x) && exists("L_HDL_L", x)) {
    x[, L_HDL_TG_pct := collate_flags(L_HDL_TG, L_HDL_L)]
  } else if (!(exists("L_HDL_TG", x) && exists("L_HDL_L", x)) && exists("L_HDL_TG_pct", x)) {
    x[, L_HDL_TG_pct := NULL]
  }


  if (exists("M_HDL_CE", x) && exists("M_HDL_L", x)) {
    x[, M_HDL_CE_pct := collate_flags(M_HDL_CE, M_HDL_L)]
  } else if (!(exists("M_HDL_CE", x) && exists("M_HDL_L", x)) && exists("M_HDL_CE_pct", x)) {
    x[, M_HDL_CE_pct := NULL]
  }

  if (exists("M_HDL_FC", x) && exists("M_HDL_L", x)) {
    x[, M_HDL_FC_pct := collate_flags(M_HDL_FC, M_HDL_L)]
  } else if (!(exists("M_HDL_FC", x) && exists("M_HDL_L", x)) && exists("M_HDL_FC_pct", x)) {
    x[, M_HDL_FC_pct := NULL]
  }

  if (exists("M_HDL_C", x) && exists("M_HDL_L", x)) {
    x[, M_HDL_C_pct := collate_flags(M_HDL_C, M_HDL_L)]
  } else if (!(exists("M_HDL_C", x) && exists("M_HDL_L", x)) && exists("M_HDL_C_pct", x)) {
    x[, M_HDL_C_pct := NULL]
  }

  if (exists("M_HDL_PL", x) && exists("M_HDL_L", x)) {
    x[, M_HDL_PL_pct := collate_flags(M_HDL_PL, M_HDL_L)]
  } else if (!(exists("M_HDL_PL", x) && exists("M_HDL_L", x)) && exists("M_HDL_PL_pct", x)) {
    x[, M_HDL_PL_pct := NULL]
  }

  if (exists("M_HDL_TG", x) && exists("M_HDL_L", x)) {
    x[, M_HDL_TG_pct := collate_flags(M_HDL_TG, M_HDL_L)]
  } else if (!(exists("M_HDL_TG", x) && exists("M_HDL_L", x)) && exists("M_HDL_TG_pct", x)) {
    x[, M_HDL_TG_pct := NULL]
  }


  if (exists("S_HDL_CE", x) && exists("S_HDL_L", x)) {
    x[, S_HDL_CE_pct := collate_flags(S_HDL_CE, S_HDL_L)]
  } else if (!(exists("S_HDL_CE", x) && exists("S_HDL_L", x)) && exists("S_HDL_CE_pct", x)) {
    x[, S_HDL_CE_pct := NULL]
  }

  if (exists("S_HDL_FC", x) && exists("S_HDL_L", x)) {
    x[, S_HDL_FC_pct := collate_flags(S_HDL_FC, S_HDL_L)]
  } else if (!(exists("S_HDL_FC", x) && exists("S_HDL_L", x)) && exists("S_HDL_FC_pct", x)) {
    x[, S_HDL_FC_pct := NULL]
  }

  if (exists("S_HDL_C", x) && exists("S_HDL_L", x)) {
    x[, S_HDL_C_pct := collate_flags(S_HDL_C, S_HDL_L)]
  } else if (!(exists("S_HDL_C", x) && exists("S_HDL_L", x)) && exists("S_HDL_C_pct", x)) {
    x[, S_HDL_C_pct := NULL]
  }

  if (exists("S_HDL_PL", x) && exists("S_HDL_L", x)) {
    x[, S_HDL_PL_pct := collate_flags(S_HDL_PL, S_HDL_L)]
  } else if (!(exists("S_HDL_PL", x) && exists("S_HDL_L", x)) && exists("S_HDL_PL_pct", x)) {
    x[, S_HDL_PL_pct := NULL]
  }

  if (exists("S_HDL_TG", x) && exists("S_HDL_L", x)) {
    x[, S_HDL_TG_pct := collate_flags(S_HDL_TG, S_HDL_L)]
  } else if (!(exists("S_HDL_TG", x) && exists("S_HDL_L", x)) && exists("S_HDL_TG_pct", x)) {
    x[, S_HDL_TG_pct := NULL]
  }


  # Compute fatty acid percentages
  if (exists("Omega_3", x) && exists("Total_FA", x)) {
    x[, Omega_3_pct := collate_flags(Omega_3, Total_FA)]
  } else if (!(exists("Omega_3", x) && exists("Total_FA", x)) && exists("Omega_3_pct", x)) {
    x[, Omega_3_pct := NULL]
  }

  if (exists("Omega_6", x) && exists("Total_FA", x)) {
    x[, Omega_6_pct := collate_flags(Omega_6, Total_FA)]
  } else if (!(exists("Omega_6", x) && exists("Total_FA", x)) && exists("Omega_6_pct", x)) {
    x[, Omega_6_pct := NULL]
  }

  if (exists("LA", x) && exists("Total_FA", x)) {
    x[, LA_pct := collate_flags(LA, Total_FA)]
  } else if (!(exists("LA", x) && exists("Total_FA", x)) && exists("LA_pct", x)) {
    x[, LA_pct := NULL]
  }

  if (exists("MUFA", x) && exists("Total_FA", x)) {
    x[, MUFA_pct := collate_flags(MUFA, Total_FA)]
  } else if (!(exists("MUFA", x) && exists("Total_FA", x)) && exists("MUFA_pct", x)) {
    x[, MUFA_pct := NULL]
  }

  if (exists("PUFA", x) && exists("Total_FA", x)) {
    x[, PUFA_pct := collate_flags(PUFA, Total_FA)]
  } else if (!(exists("PUFA", x) && exists("Total_FA", x)) && exists("PUFA_pct", x)) {
    x[, PUFA_pct := NULL]
  }

  if (exists("SFA", x) && exists("Total_FA", x)) {
    x[, SFA_pct := collate_flags(SFA, Total_FA)]
  } else if (!(exists("SFA", x) && exists("Total_FA", x)) && exists("SFA_pct", x)) {
    x[, SFA_pct := NULL]
  }

  if (exists("DHA", x) && exists("Total_FA", x)) {
    x[, DHA_pct := collate_flags(DHA, Total_FA)]
  } else if (!(exists("DHA", x) && exists("Total_FA", x)) && exists("DHA_pct", x)) {
    x[, DHA_pct := NULL]
  }


  # Miscellaneous ratios
  if (exists("Omega_6", x) && exists("Omega_3", x)) {
    x[, Omega_6_by_Omega_3 := collate_flags(Omega_6, Omega_3)]
  } else if (!(exists("Omega_6", x) && exists("Omega_3", x)) && exists("Omega_6_by_Omega_3", x)) {
    x[, Omega_6_by_Omega_3 := NULL]
  }

  if (exists("ApoB", x) && exists("ApoA1", x)) {
    x[, ApoB_by_ApoA1 := collate_flags(ApoB, ApoA1)]
  } else if (!(exists("ApoB", x) && exists("ApoA1", x)) && exists("ApoB_by_ApoA1", x)) {
    x[, ApoB_by_ApoA1 := NULL]
  }

  if (exists("PUFA", x) && exists("MUFA", x)) {
    x[, PUFA_by_MUFA := collate_flags(PUFA, MUFA)]
  } else if (!(exists("PUFA", x) && exists("MUFA", x)) && exists("PUFA_by_MUFA", x)) {
    x[, PUFA_by_MUFA := NULL]
  }

  if (exists("Total_TG", x) && exists("Phosphoglyc", x)) {
    x[, TG_by_PG := collate_flags(Total_TG, Phosphoglyc)]
  } else if (!(exists("Total_TG", x) && exists("Phosphoglyc", x)) && exists("TG_by_PG", x)) {
    x[, TG_by_PG := NULL]
  }

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
  if (exists("XXL_VLDL_CE", x) && exists("XXL_VLDL_FC", x)) {
    x[, XXL_VLDL_C := collate_flags(XXL_VLDL_CE, XXL_VLDL_FC)]
  } else if (!(exists("XXL_VLDL_CE", x) && exists("XXL_VLDL_FC", x)) && exists("XXL_VLDL_C", x)) {
    x[, XXL_VLDL_C := NULL]
  }

  if (exists("XL_VLDL_CE", x) && exists("XL_VLDL_FC", x)) {
    x[, XL_VLDL_C := collate_flags(XL_VLDL_CE, XL_VLDL_FC)]
  } else if (!(exists("XL_VLDL_CE", x) && exists("XL_VLDL_FC", x)) && exists("XL_VLDL_C", x)) {
    x[, XL_VLDL_C := NULL]
  }

  if (exists("L_VLDL_CE", x) && exists("L_VLDL_FC", x)) {
    x[, L_VLDL_C := collate_flags(L_VLDL_CE, L_VLDL_FC)]
  } else if (!(exists("L_VLDL_CE", x) && exists("L_VLDL_FC", x)) && exists("L_VLDL_C", x)) {
    x[, L_VLDL_C := NULL]
  }

  if (exists("M_VLDL_CE", x) && exists("M_VLDL_FC", x)) {
    x[, M_VLDL_C := collate_flags(M_VLDL_CE, M_VLDL_FC)]
  } else if (!(exists("M_VLDL_CE", x) && exists("M_VLDL_FC", x)) && exists("M_VLDL_C", x)) {
    x[, M_VLDL_C := NULL]
  }

  if (exists("S_VLDL_CE", x) && exists("S_VLDL_FC", x)) {
    x[, S_VLDL_C := collate_flags(S_VLDL_CE, S_VLDL_FC)]
  } else if (!(exists("S_VLDL_CE", x) && exists("S_VLDL_FC", x)) && exists("S_VLDL_C", x)) {
    x[, S_VLDL_C := NULL]
  }

  if (exists("XS_VLDL_CE", x) && exists("XS_VLDL_FC", x)) {
    x[, XS_VLDL_C := collate_flags(XS_VLDL_CE, XS_VLDL_FC)]
  } else if (!(exists("XS_VLDL_CE", x) && exists("XS_VLDL_FC", x)) && exists("XS_VLDL_C", x)) {
    x[, XS_VLDL_C := NULL]
  }

  if (exists("IDL_CE", x) && exists("IDL_FC", x)) {
    x[, IDL_C := collate_flags(IDL_CE, IDL_FC)]
  } else if (!(exists("IDL_CE", x) && exists("IDL_FC", x)) && exists("IDL_C", x)) {
    x[, IDL_C := NULL]
  }

  if (exists("L_LDL_CE", x) && exists("L_LDL_FC", x)) {
    x[, L_LDL_C := collate_flags(L_LDL_CE, L_LDL_FC)]
  } else if (!(exists("L_LDL_CE", x) && exists("L_LDL_FC", x)) && exists("L_LDL_C", x)) {
    x[, L_LDL_C := NULL]
  }

  if (exists("M_LDL_CE", x) && exists("M_LDL_FC", x)) {
    x[, M_LDL_C := collate_flags(M_LDL_CE, M_LDL_FC)]
  } else if (!(exists("M_LDL_CE", x) && exists("M_LDL_FC", x)) && exists("M_LDL_C", x)) {
    x[, M_LDL_C := NULL]
  }

  if (exists("S_LDL_CE", x) && exists("S_LDL_FC", x)) {
    x[, S_LDL_C := collate_flags(S_LDL_CE, S_LDL_FC)]
  } else if (!(exists("S_LDL_CE", x) && exists("S_LDL_FC", x)) && exists("S_LDL_C", x)) {
    x[, S_LDL_C := NULL]
  }

  if (exists("XL_HDL_CE", x) && exists("XL_HDL_FC", x)) {
    x[, XL_HDL_C := collate_flags(XL_HDL_CE, XL_HDL_FC)]
  } else if (!(exists("XL_HDL_CE", x) && exists("XL_HDL_FC", x)) && exists("XL_HDL_C", x)) {
    x[, XL_HDL_C := NULL]
  }

  if (exists("L_HDL_CE", x) && exists("L_HDL_FC", x)) {
    x[, L_HDL_C := collate_flags(L_HDL_CE, L_HDL_FC)]
  } else if (!(exists("L_HDL_CE", x) && exists("L_HDL_FC", x)) && exists("L_HDL_C", x)) {
    x[, L_HDL_C := NULL]
  }

  if (exists("M_HDL_CE", x) && exists("M_HDL_FC", x)) {
    x[, M_HDL_C := collate_flags(M_HDL_CE, M_HDL_FC)]
  } else if (!(exists("M_HDL_CE", x) && exists("M_HDL_FC", x)) && exists("M_HDL_C", x)) {
    x[, M_HDL_C := NULL]
  }

  if (exists("S_HDL_CE", x) && exists("S_HDL_FC", x)) {
    x[, S_HDL_C := collate_flags(S_HDL_CE, S_HDL_FC)]
  } else if (!(exists("S_HDL_CE", x) && exists("S_HDL_FC", x)) && exists("S_HDL_C", x)) {
    x[, S_HDL_C := NULL]
  }


  # Now compute total lipids in lipoprotein subclasses
  # Total Lipids = cholesterol, phospholipids, triglycerides
  if (exists("XXL_VLDL_C", x) && exists("XXL_VLDL_PL", x) && exists("XXL_VLDL_TG", x)) {
    x[, XXL_VLDL_L := collate_flags(XXL_VLDL_C, XXL_VLDL_PL, XXL_VLDL_TG)]
  } else if (!(exists("XXL_VLDL_C", x) && exists("XXL_VLDL_PL", x) && exists("XXL_VLDL_TG", x)) && exists("XXL_VLDL_L", x)) {
    x[, XXL_VLDL_L := NULL]
  }

  if (exists("XL_VLDL_C", x) && exists("XL_VLDL_PL", x) && exists("XL_VLDL_TG", x)) {
    x[, XL_VLDL_L := collate_flags(XL_VLDL_C, XL_VLDL_PL, XL_VLDL_TG)]
  } else if (!(exists("XL_VLDL_C", x) && exists("XL_VLDL_PL", x) && exists("XL_VLDL_TG", x)) && exists("XL_VLDL_L", x)) {
    x[, XL_VLDL_L := NULL]
  }

  if (exists("L_VLDL_C", x) && exists("L_VLDL_PL", x) && exists("L_VLDL_TG", x)) {
    x[, L_VLDL_L := collate_flags(L_VLDL_C, L_VLDL_PL, L_VLDL_TG)]
  } else if (!(exists("L_VLDL_C", x) && exists("L_VLDL_PL", x) && exists("L_VLDL_TG", x)) && exists("L_VLDL_L", x)) {
    x[, L_VLDL_L := NULL]
  }

  if (exists("M_VLDL_C", x) && exists("M_VLDL_PL", x) && exists("M_VLDL_TG", x)) {
    x[, M_VLDL_L := collate_flags(M_VLDL_C, M_VLDL_PL, M_VLDL_TG)]
  } else if (!(exists("M_VLDL_C", x) && exists("M_VLDL_PL", x) && exists("M_VLDL_TG", x)) && exists("M_VLDL_L", x)) {
    x[, M_VLDL_L := NULL]
  }

  if (exists("S_VLDL_C", x) && exists("S_VLDL_PL", x) && exists("S_VLDL_TG", x)) {
    x[, S_VLDL_L := collate_flags(S_VLDL_C, S_VLDL_PL, S_VLDL_TG)]
  } else if (!(exists("S_VLDL_C", x) && exists("S_VLDL_PL", x) && exists("S_VLDL_TG", x)) && exists("S_VLDL_L", x)) {
    x[, S_VLDL_L := NULL]
  }

  if (exists("XS_VLDL_C", x) && exists("XS_VLDL_PL", x) && exists("XS_VLDL_TG", x)) {
    x[, XS_VLDL_L := collate_flags(XS_VLDL_C, XS_VLDL_PL, XS_VLDL_TG)]
  } else if (!(exists("XS_VLDL_C", x) && exists("XS_VLDL_PL", x) && exists("XS_VLDL_TG", x)) && exists("XS_VLDL_L", x)) {
    x[, XS_VLDL_L := NULL]
  }

  if (exists("IDL_C", x) && exists("IDL_PL", x) && exists("IDL_TG", x)) {
    x[, IDL_L := collate_flags(IDL_C, IDL_PL, IDL_TG)]
  } else if (!(exists("IDL_C", x) && exists("IDL_PL", x) && exists("IDL_TG", x)) && exists("IDL_L", x)) {
    x[, IDL_L := NULL]
  }

  if (exists("L_LDL_C", x) && exists("L_LDL_PL", x) && exists("L_LDL_TG", x)) {
    x[, L_LDL_L := collate_flags(L_LDL_C, L_LDL_PL, L_LDL_TG)]
  } else if (!(exists("L_LDL_C", x) && exists("L_LDL_PL", x) && exists("L_LDL_TG", x)) && exists("L_LDL_L", x)) {
    x[, L_LDL_L := NULL]
  }

  if (exists("M_LDL_C", x) && exists("M_LDL_PL", x) && exists("M_LDL_TG", x)) {
    x[, M_LDL_L := collate_flags(M_LDL_C, M_LDL_PL, M_LDL_TG)]
  } else if (!(exists("M_LDL_C", x) && exists("M_LDL_PL", x) && exists("M_LDL_TG", x)) && exists("M_LDL_L", x)) {
    x[, M_LDL_L := NULL]
  }

  if (exists("S_LDL_C", x) && exists("S_LDL_PL", x) && exists("S_LDL_TG", x)) {
    x[, S_LDL_L := collate_flags(S_LDL_C, S_LDL_PL, S_LDL_TG)]
  } else if (!(exists("S_LDL_C", x) && exists("S_LDL_PL", x) && exists("S_LDL_TG", x)) && exists("S_LDL_L", x)) {
    x[, S_LDL_L := NULL]
  }

  if (exists("XL_HDL_C", x) && exists("XL_HDL_PL", x) && exists("XL_HDL_TG", x)) {
    x[, XL_HDL_L := collate_flags(XL_HDL_C, XL_HDL_PL, XL_HDL_TG)]
  } else if (!(exists("XL_HDL_C", x) && exists("XL_HDL_PL", x) && exists("XL_HDL_TG", x)) && exists("XL_HDL_L", x)) {
    x[, XL_HDL_L := NULL]
  }

  if (exists("L_HDL_C", x) && exists("L_HDL_PL", x) && exists("L_HDL_TG", x)) {
    x[, L_HDL_L := collate_flags(L_HDL_C, L_HDL_PL, L_HDL_TG)]
  } else if (!(exists("L_HDL_C", x) && exists("L_HDL_PL", x) && exists("L_HDL_TG", x)) && exists("L_HDL_L", x)) {
    x[, L_HDL_L := NULL]
  }

  if (exists("M_HDL_C", x) && exists("M_HDL_PL", x) && exists("M_HDL_TG", x)) {
    x[, M_HDL_L := collate_flags(M_HDL_C, M_HDL_PL, M_HDL_TG)]
  } else if (!(exists("M_HDL_C", x) && exists("M_HDL_PL", x) && exists("M_HDL_TG", x)) && exists("M_HDL_L", x)) {
    x[, M_HDL_L := NULL]
  }

  if (exists("S_HDL_C", x) && exists("S_HDL_PL", x) && exists("S_HDL_TG", x)) {
    x[, S_HDL_L := collate_flags(S_HDL_C, S_HDL_PL, S_HDL_TG)]
  } else if (!(exists("S_HDL_C", x) && exists("S_HDL_PL", x) && exists("S_HDL_TG", x)) && exists("S_HDL_L", x)) {
    x[, S_HDL_L := NULL]
  }


  # Now compute totals for lipoprotein classes.
  # Eg. free cholesterol in LDL (LDL_FC) = sum(free choleseterol in LDL of different sizes).
  if (exists("XXL_VLDL_CE", x) && exists("XL_VLDL_CE", x) && exists("L_VLDL_CE", x) && exists("M_VLDL_CE", x) && exists("S_VLDL_CE", x) && exists("XS_VLDL_CE", x)) {
    x[, VLDL_CE := collate_flags(XXL_VLDL_CE, XL_VLDL_CE, L_VLDL_CE, M_VLDL_CE, S_VLDL_CE, XS_VLDL_CE)]
  } else if (!(exists("XXL_VLDL_CE", x) && exists("XL_VLDL_CE", x) && exists("L_VLDL_CE", x) && exists("M_VLDL_CE", x) && exists("S_VLDL_CE", x) && exists("XS_VLDL_CE", x)) && exists("VLDL_CE", x)) {
    x[, VLDL_CE := NULL]
  }

  if (exists("XXL_VLDL_FC", x) && exists("XL_VLDL_FC", x) && exists("L_VLDL_FC", x) && exists("M_VLDL_FC", x) && exists("S_VLDL_FC", x) && exists("XS_VLDL_FC", x)) {
    x[, VLDL_FC := collate_flags(XXL_VLDL_FC, XL_VLDL_FC, L_VLDL_FC, M_VLDL_FC, S_VLDL_FC, XS_VLDL_FC)]
  } else if (!(exists("XXL_VLDL_FC", x) && exists("XL_VLDL_FC", x) && exists("L_VLDL_FC", x) && exists("M_VLDL_FC", x) && exists("S_VLDL_FC", x) && exists("XS_VLDL_FC", x)) && exists("VLDL_FC", x)) {
    x[, VLDL_FC := NULL]
  }

  if (exists("XXL_VLDL_C", x) && exists("XL_VLDL_C", x) && exists("L_VLDL_C", x) && exists("M_VLDL_C", x) && exists("S_VLDL_C", x) && exists("XS_VLDL_C", x)) {
    x[, VLDL_C := collate_flags(XXL_VLDL_C, XL_VLDL_C, L_VLDL_C, M_VLDL_C, S_VLDL_C, XS_VLDL_C)]
  } else if (!(exists("XXL_VLDL_C", x) && exists("XL_VLDL_C", x) && exists("L_VLDL_C", x) && exists("M_VLDL_C", x) && exists("S_VLDL_C", x) && exists("XS_VLDL_C", x)) && exists("VLDL_C", x)) {
    x[, VLDL_C := NULL]
  }

  if (exists("XXL_VLDL_PL", x) && exists("XL_VLDL_PL", x) && exists("L_VLDL_PL", x) && exists("M_VLDL_PL", x) && exists("S_VLDL_PL", x) && exists("XS_VLDL_PL", x)) {
    x[, VLDL_PL := collate_flags(XXL_VLDL_PL, XL_VLDL_PL, L_VLDL_PL, M_VLDL_PL, S_VLDL_PL, XS_VLDL_PL)]
  } else if (!(exists("XXL_VLDL_PL", x) && exists("XL_VLDL_PL", x) && exists("L_VLDL_PL", x) && exists("M_VLDL_PL", x) && exists("S_VLDL_PL", x) && exists("XS_VLDL_PL", x)) && exists("VLDL_PL", x)) {
    x[, VLDL_PL := NULL]
  }

  if (exists("XXL_VLDL_TG", x) && exists("XL_VLDL_TG", x) && exists("L_VLDL_TG", x) && exists("M_VLDL_TG", x) && exists("S_VLDL_TG", x) && exists("XS_VLDL_TG", x)) {
    x[, VLDL_TG := collate_flags(XXL_VLDL_TG, XL_VLDL_TG, L_VLDL_TG, M_VLDL_TG, S_VLDL_TG, XS_VLDL_TG)]
  } else if (!(exists("XXL_VLDL_TG", x) && exists("XL_VLDL_TG", x) && exists("L_VLDL_TG", x) && exists("M_VLDL_TG", x) && exists("S_VLDL_TG", x) && exists("XS_VLDL_TG", x)) && exists("VLDL_TG", x)) {
    x[, VLDL_TG := NULL]
  }

  if (exists("XXL_VLDL_L", x) && exists("XL_VLDL_L", x) && exists("L_VLDL_L", x) && exists("M_VLDL_L", x) && exists("S_VLDL_L", x) && exists("XS_VLDL_L", x)) {
    x[, VLDL_L := collate_flags(XXL_VLDL_L, XL_VLDL_L, L_VLDL_L, M_VLDL_L, S_VLDL_L, XS_VLDL_L)]
  } else if (!(exists("XXL_VLDL_L", x) && exists("XL_VLDL_L", x) && exists("L_VLDL_L", x) && exists("M_VLDL_L", x) && exists("S_VLDL_L", x) && exists("XS_VLDL_L", x)) && exists("VLDL_L", x)) {
    x[, VLDL_L := NULL]
  }

  if (exists("XXL_VLDL_P", x) && exists("XL_VLDL_P", x) && exists("L_VLDL_P", x) && exists("M_VLDL_P", x) && exists("S_VLDL_P", x) && exists("XS_VLDL_P", x)) {
    x[, VLDL_P := collate_flags(XXL_VLDL_P, XL_VLDL_P, L_VLDL_P, M_VLDL_P, S_VLDL_P, XS_VLDL_P)]
  } else if (!(exists("XXL_VLDL_P", x) && exists("XL_VLDL_P", x) && exists("L_VLDL_P", x) && exists("M_VLDL_P", x) && exists("S_VLDL_P", x) && exists("XS_VLDL_P", x)) && exists("VLDL_P", x)) {
    x[, VLDL_P := NULL]
  }


  if (exists("L_LDL_CE", x) && exists("M_LDL_CE", x) && exists("S_LDL_CE", x)) {
    x[, LDL_CE := collate_flags(L_LDL_CE, M_LDL_CE, S_LDL_CE)]
  } else if (!(exists("L_LDL_CE", x) && exists("M_LDL_CE", x) && exists("S_LDL_CE", x)) && exists("LDL_CE", x)) {
    x[, LDL_CE := NULL]
  }

  if (exists("L_LDL_FC", x) && exists("M_LDL_FC", x) && exists("S_LDL_FC", x)) {
    x[, LDL_FC := collate_flags(L_LDL_FC, M_LDL_FC, S_LDL_FC)]
  } else if (!(exists("L_LDL_FC", x) && exists("M_LDL_FC", x) && exists("S_LDL_FC", x)) && exists("LDL_FC", x)) {
    x[, LDL_FC := NULL]
  }

  if (exists("L_LDL_C", x) && exists("M_LDL_C", x) && exists("S_LDL_C", x)) {
    x[, LDL_C := collate_flags(L_LDL_C, M_LDL_C, S_LDL_C)]
  } else if (!(exists("L_LDL_C", x) && exists("M_LDL_C", x) && exists("S_LDL_C", x)) && exists("LDL_C", x)) {
    x[, LDL_C := NULL]
  }

  if (exists("L_LDL_PL", x) && exists("M_LDL_PL", x) && exists("S_LDL_PL", x)) {
    x[, LDL_PL := collate_flags(L_LDL_PL, M_LDL_PL, S_LDL_PL)]
  } else if (!(exists("L_LDL_PL", x) && exists("M_LDL_PL", x) && exists("S_LDL_PL", x)) && exists("LDL_PL", x)) {
    x[, LDL_PL := NULL]
  }

  if (exists("L_LDL_TG", x) && exists("M_LDL_TG", x) && exists("S_LDL_TG", x)) {
    x[, LDL_TG := collate_flags(L_LDL_TG, M_LDL_TG, S_LDL_TG)]
  } else if (!(exists("L_LDL_TG", x) && exists("M_LDL_TG", x) && exists("S_LDL_TG", x)) && exists("LDL_TG", x)) {
    x[, LDL_TG := NULL]
  }

  if (exists("L_LDL_L", x) && exists("M_LDL_L", x) && exists("S_LDL_L", x)) {
    x[, LDL_L := collate_flags(L_LDL_L, M_LDL_L, S_LDL_L)]
  } else if (!(exists("L_LDL_L", x) && exists("M_LDL_L", x) && exists("S_LDL_L", x)) && exists("LDL_L", x)) {
    x[, LDL_L := NULL]
  }

  if (exists("L_LDL_P", x) && exists("M_LDL_P", x) && exists("S_LDL_P", x)) {
    x[, LDL_P := collate_flags(L_LDL_P, M_LDL_P, S_LDL_P)]
  } else if (!(exists("L_LDL_P", x) && exists("M_LDL_P", x) && exists("S_LDL_P", x)) && exists("LDL_P", x)) {
    x[, LDL_P := NULL]
  }


  if (exists("XL_HDL_CE", x) && exists("L_HDL_CE", x) && exists("M_HDL_CE", x) && exists("S_HDL_CE", x)) {
    x[, HDL_CE := collate_flags(XL_HDL_CE, L_HDL_CE, M_HDL_CE, S_HDL_CE)]
  } else if (!(exists("XL_HDL_CE", x) && exists("L_HDL_CE", x) && exists("M_HDL_CE", x) && exists("S_HDL_CE", x)) && exists("HDL_CE", x)) {
    x[, HDL_CE := NULL]
  }

  if (exists("XL_HDL_FC", x) && exists("L_HDL_FC", x) && exists("M_HDL_FC", x) && exists("S_HDL_FC", x)) {
    x[, HDL_FC := collate_flags(XL_HDL_FC, L_HDL_FC, M_HDL_FC, S_HDL_FC)]
  } else if (!(exists("XL_HDL_FC", x) && exists("L_HDL_FC", x) && exists("M_HDL_FC", x) && exists("S_HDL_FC", x)) && exists("HDL_FC", x)) {
    x[, HDL_FC := NULL]
  }

  if (exists("XL_HDL_C", x) && exists("L_HDL_C", x) && exists("M_HDL_C", x) && exists("S_HDL_C", x)) {
    x[, HDL_C := collate_flags(XL_HDL_C, L_HDL_C, M_HDL_C, S_HDL_C)]
  } else if (!(exists("XL_HDL_C", x) && exists("L_HDL_C", x) && exists("M_HDL_C", x) && exists("S_HDL_C", x)) && exists("HDL_C", x)) {
    x[, HDL_C := NULL]
  }

  if (exists("XL_HDL_PL", x) && exists("L_HDL_PL", x) && exists("M_HDL_PL", x) && exists("S_HDL_PL", x)) {
    x[, HDL_PL := collate_flags(XL_HDL_PL, L_HDL_PL, M_HDL_PL, S_HDL_PL)]
  } else if (!(exists("XL_HDL_PL", x) && exists("L_HDL_PL", x) && exists("M_HDL_PL", x) && exists("S_HDL_PL", x)) && exists("HDL_PL", x)) {
    x[, HDL_PL := NULL]
  }

  if (exists("XL_HDL_TG", x) && exists("L_HDL_TG", x) && exists("M_HDL_TG", x) && exists("S_HDL_TG", x)) {
    x[, HDL_TG := collate_flags(XL_HDL_TG, L_HDL_TG, M_HDL_TG, S_HDL_TG)]
  } else if (!(exists("XL_HDL_TG", x) && exists("L_HDL_TG", x) && exists("M_HDL_TG", x) && exists("S_HDL_TG", x)) && exists("HDL_TG", x)) {
    x[, HDL_TG := NULL]
  }

  if (exists("XL_HDL_L", x) && exists("L_HDL_L", x) && exists("M_HDL_L", x) && exists("S_HDL_L", x)) {
    x[, HDL_L := collate_flags(XL_HDL_L, L_HDL_L, M_HDL_L, S_HDL_L)]
  } else if (!(exists("XL_HDL_L", x) && exists("L_HDL_L", x) && exists("M_HDL_L", x) && exists("S_HDL_L", x)) && exists("HDL_L", x)) {
    x[, HDL_L := NULL]
  }

  if (exists("XL_HDL_P", x) && exists("L_HDL_P", x) && exists("M_HDL_P", x) && exists("S_HDL_P", x)) {
    x[, HDL_P := collate_flags(XL_HDL_P, L_HDL_P, M_HDL_P, S_HDL_P)]
  } else if (!(exists("XL_HDL_P", x) && exists("L_HDL_P", x) && exists("M_HDL_P", x) && exists("S_HDL_P", x)) && exists("HDL_P", x)) {
    x[, HDL_P := NULL]
  }


  # Next compute serum totals for lipids
  if (exists("VLDL_CE", x) && exists("LDL_CE", x) && exists("IDL_CE", x) && exists("HDL_CE", x)) {
    x[, Total_CE := collate_flags(VLDL_CE, LDL_CE, IDL_CE, HDL_CE)]
  } else if (!(exists("VLDL_CE", x) && exists("LDL_CE", x) && exists("IDL_CE", x) && exists("HDL_CE", x)) && exists("Total_CE", x)) {
    x[, Total_CE := NULL]
  }

  if (exists("VLDL_FC", x) && exists("LDL_FC", x) && exists("IDL_FC", x) && exists("HDL_FC", x)) {
    x[, Total_FC := collate_flags(VLDL_FC, LDL_FC, IDL_FC, HDL_FC)]
  } else if (!(exists("VLDL_FC", x) && exists("LDL_FC", x) && exists("IDL_FC", x) && exists("HDL_FC", x)) && exists("Total_FC", x)) {
    x[, Total_FC := NULL]
  }

  if (exists("VLDL_C", x) && exists("LDL_C", x) && exists("IDL_C", x) && exists("HDL_C", x)) {
    x[, Total_C := collate_flags(VLDL_C, LDL_C, IDL_C, HDL_C)]
  } else if (!(exists("VLDL_C", x) && exists("LDL_C", x) && exists("IDL_C", x) && exists("HDL_C", x)) && exists("Total_C", x)) {
    x[, Total_C := NULL]
  }

  if (exists("VLDL_PL", x) && exists("LDL_PL", x) && exists("IDL_PL", x) && exists("HDL_PL", x)) {
    x[, Total_PL := collate_flags(VLDL_PL, LDL_PL, IDL_PL, HDL_PL)]
  } else if (!(exists("VLDL_PL", x) && exists("LDL_PL", x) && exists("IDL_PL", x) && exists("HDL_PL", x)) && exists("Total_PL", x)) {
    x[, Total_PL := NULL]
  }

  if (exists("VLDL_TG", x) && exists("LDL_TG", x) && exists("IDL_TG", x) && exists("HDL_TG", x)) {
    x[, Total_TG := collate_flags(VLDL_TG, LDL_TG, IDL_TG, HDL_TG)]
  } else if (!(exists("VLDL_TG", x) && exists("LDL_TG", x) && exists("IDL_TG", x) && exists("HDL_TG", x)) && exists("Total_TG", x)) {
    x[, Total_TG := NULL]
  }

  if (exists("VLDL_L", x) && exists("LDL_L", x) && exists("IDL_L", x) && exists("HDL_L", x)) {
    x[, Total_L := collate_flags(VLDL_L, LDL_L, IDL_L, HDL_L)]
  } else if (!(exists("VLDL_L", x) && exists("LDL_L", x) && exists("IDL_L", x) && exists("HDL_L", x)) && exists("Total_L", x)) {
    x[, Total_L := NULL]
  }

  if (exists("VLDL_P", x) && exists("LDL_P", x) && exists("IDL_P", x) && exists("HDL_P", x)) {
    x[, Total_P := collate_flags(VLDL_P, LDL_P, IDL_P, HDL_P)]
  } else if (!(exists("VLDL_P", x) && exists("LDL_P", x) && exists("IDL_P", x) && exists("HDL_P", x)) && exists("Total_P", x)) {
    x[, Total_P := NULL]
  }


  # Finally miscellaneous composite biomarkers
  if (exists("Omega_3", x) && exists("Omega_6", x)) {
    x[, PUFA := collate_flags(Omega_3, Omega_6)]
  } else if (!(exists("Omega_3", x) && exists("Omega_6", x)) && exists("PUFA", x)) {
    x[, PUFA := NULL]
  }

  if (exists("PUFA", x) && exists("MUFA", x) && exists("SFA", x)) {
    x[, Total_FA := collate_flags(PUFA, MUFA, SFA)]
  } else if (!(exists("PUFA", x) && exists("MUFA", x) && exists("SFA", x)) && exists("Total_FA", x)) {
    x[, Total_FA := NULL]
  }


  # Miscellaneous composite markers
  if (exists("Leu", x) && exists("Ile", x) && exists("Val", x)) {
    x[, Total_BCAA := collate_flags(Leu, Ile, Val)] # total branched chain amino acids
  } else if (!(exists("Leu", x) && exists("Ile", x) && exists("Val", x)) && exists("Total_BCAA", x)) {
    x[, Total_BCAA := NULL]
  }

  if (exists("Total_C", x) && exists("HDL_C", x)) {
    x[, non_HDL_C := collate_flags(Total_C, HDL_C)] # non HDL cholesterol
  } else if (!(exists("Total_C", x) && exists("HDL_C", x)) && exists("non_HDL_C", x)) {
    x[, non_HDL_C := NULL]
  }

  if (exists("Total_C", x) && exists("HDL_C", x) && exists("LDL_C", x)) {
    x[, Remnant_C := collate_flags(Total_C, HDL_C, LDL_C)] # remnant cholesterol
  } else if (!(exists("Total_C", x) && exists("HDL_C", x) && exists("LDL_C", x)) && exists("Remnant_C", x)) {
    x[, Remnant_C := NULL]
  }

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
  if (exists("XXL_VLDL_FC", x) && exists("XXL_VLDL_C", x)) {
    x[, XXL_VLDL_FC_pct_C := collate_flags(XXL_VLDL_FC, XXL_VLDL_C)]
  } else if (!(exists("XXL_VLDL_FC", x) && exists("XXL_VLDL_C", x)) && exists("XXL_VLDL_FC_pct_C", x)) {
    x[, XXL_VLDL_FC_pct_C := NULL]
  }

  if (exists("XL_VLDL_FC", x) && exists("XL_VLDL_C", x)) {
    x[, XL_VLDL_FC_pct_C := collate_flags(XL_VLDL_FC, XL_VLDL_C)]
  } else if (!(exists("XL_VLDL_FC", x) && exists("XL_VLDL_C", x)) && exists("XL_VLDL_FC_pct_C", x)) {
    x[, XL_VLDL_FC_pct_C := NULL]
  }

  if (exists("L_VLDL_FC", x) && exists("L_VLDL_C", x)) {
    x[, L_VLDL_FC_pct_C := collate_flags(L_VLDL_FC, L_VLDL_C)]
  } else if (!(exists("L_VLDL_FC", x) && exists("L_VLDL_C", x)) && exists("L_VLDL_FC_pct_C", x)) {
    x[, L_VLDL_FC_pct_C := NULL]
  }

  if (exists("M_VLDL_FC", x) && exists("M_VLDL_C", x)) {
    x[, M_VLDL_FC_pct_C := collate_flags(M_VLDL_FC, M_VLDL_C)]
  } else if (!(exists("M_VLDL_FC", x) && exists("M_VLDL_C", x)) && exists("M_VLDL_FC_pct_C", x)) {
    x[, M_VLDL_FC_pct_C := NULL]
  }

  if (exists("S_VLDL_FC", x) && exists("S_VLDL_C", x)) {
    x[, S_VLDL_FC_pct_C := collate_flags(S_VLDL_FC, S_VLDL_C)]
  } else if (!(exists("S_VLDL_FC", x) && exists("S_VLDL_C", x)) && exists("S_VLDL_FC_pct_C", x)) {
    x[, S_VLDL_FC_pct_C := NULL]
  }

  if (exists("XS_VLDL_FC", x) && exists("XS_VLDL_C", x)) {
    x[, XS_VLDL_FC_pct_C := collate_flags(XS_VLDL_FC, XS_VLDL_C)]
  } else if (!(exists("XS_VLDL_FC", x) && exists("XS_VLDL_C", x)) && exists("XS_VLDL_FC_pct_C", x)) {
    x[, XS_VLDL_FC_pct_C := NULL]
  }

  if (exists("L_LDL_FC", x) && exists("L_LDL_C", x)) {
    x[, L_LDL_FC_pct_C := collate_flags(L_LDL_FC, L_LDL_C)]
  } else if (!(exists("L_LDL_FC", x) && exists("L_LDL_C", x)) && exists("L_LDL_FC_pct_C", x)) {
    x[, L_LDL_FC_pct_C := NULL]
  }

  if (exists("M_LDL_FC", x) && exists("M_LDL_C", x)) {
    x[, M_LDL_FC_pct_C := collate_flags(M_LDL_FC, M_LDL_C)]
  } else if (!(exists("M_LDL_FC", x) && exists("M_LDL_C", x)) && exists("M_LDL_FC_pct_C", x)) {
    x[, M_LDL_FC_pct_C := NULL]
  }

  if (exists("S_LDL_FC", x) && exists("S_LDL_C", x)) {
    x[, S_LDL_FC_pct_C := collate_flags(S_LDL_FC, S_LDL_C)]
  } else if (!(exists("S_LDL_FC", x) && exists("S_LDL_C", x)) && exists("S_LDL_FC_pct_C", x)) {
    x[, S_LDL_FC_pct_C := NULL]
  }

  if (exists("IDL_FC", x) && exists("IDL_C", x)) {
    x[, IDL_FC_pct_C := collate_flags(IDL_FC, IDL_C)]
  } else if (!(exists("IDL_FC", x) && exists("IDL_C", x)) && exists("IDL_FC_pct_C", x)) {
    x[, IDL_FC_pct_C := NULL]
  }

  if (exists("XL_HDL_FC", x) && exists("XL_HDL_C", x)) {
    x[, XL_HDL_FC_pct_C := collate_flags(XL_HDL_FC, XL_HDL_C)]
  } else if (!(exists("XL_HDL_FC", x) && exists("XL_HDL_C", x)) && exists("XL_HDL_FC_pct_C", x)) {
    x[, XL_HDL_FC_pct_C := NULL]
  }

  if (exists("L_HDL_FC", x) && exists("L_HDL_C", x)) {
    x[, L_HDL_FC_pct_C := collate_flags(L_HDL_FC, L_HDL_C)]
  } else if (!(exists("L_HDL_FC", x) && exists("L_HDL_C", x)) && exists("L_HDL_FC_pct_C", x)) {
    x[, L_HDL_FC_pct_C := NULL]
  }

  if (exists("M_HDL_FC", x) && exists("M_HDL_C", x)) {
    x[, M_HDL_FC_pct_C := collate_flags(M_HDL_FC, M_HDL_C)]
  } else if (!(exists("M_HDL_FC", x) && exists("M_HDL_C", x)) && exists("M_HDL_FC_pct_C", x)) {
    x[, M_HDL_FC_pct_C := NULL]
  }

  if (exists("S_HDL_FC", x) && exists("S_HDL_C", x)) {
    x[, S_HDL_FC_pct_C := collate_flags(S_HDL_FC, S_HDL_C)]
  } else if (!(exists("S_HDL_FC", x) && exists("S_HDL_C", x)) && exists("S_HDL_FC_pct_C", x)) {
    x[, S_HDL_FC_pct_C := NULL]
  }


  if (exists("XXL_VLDL_CE", x) && exists("XXL_VLDL_C", x)) {
    x[, XXL_VLDL_CE_pct_C := collate_flags(XXL_VLDL_CE, XXL_VLDL_C)]
  } else if (!(exists("XXL_VLDL_CE", x) && exists("XXL_VLDL_C", x)) && exists("XXL_VLDL_CE_pct_C", x)) {
    x[, XXL_VLDL_CE_pct_C := NULL]
  }

  if (exists("XL_VLDL_CE", x) && exists("XL_VLDL_C", x)) {
    x[, XL_VLDL_CE_pct_C := collate_flags(XL_VLDL_CE, XL_VLDL_C)]
  } else if (!(exists("XL_VLDL_CE", x) && exists("XL_VLDL_C", x)) && exists("XL_VLDL_CE_pct_C", x)) {
    x[, XL_VLDL_CE_pct_C := NULL]
  }

  if (exists("L_VLDL_CE", x) && exists("L_VLDL_C", x)) {
    x[, L_VLDL_CE_pct_C := collate_flags(L_VLDL_CE, L_VLDL_C)]
  } else if (!(exists("L_VLDL_CE", x) && exists("L_VLDL_C", x)) && exists("L_VLDL_CE_pct_C", x)) {
    x[, L_VLDL_CE_pct_C := NULL]
  }

  if (exists("M_VLDL_CE", x) && exists("M_VLDL_C", x)) {
    x[, M_VLDL_CE_pct_C := collate_flags(M_VLDL_CE, M_VLDL_C)]
  } else if (!(exists("M_VLDL_CE", x) && exists("M_VLDL_C", x)) && exists("M_VLDL_CE_pct_C", x)) {
    x[, M_VLDL_CE_pct_C := NULL]
  }

  if (exists("S_VLDL_CE", x) && exists("S_VLDL_C", x)) {
    x[, S_VLDL_CE_pct_C := collate_flags(S_VLDL_CE, S_VLDL_C)]
  } else if (!(exists("S_VLDL_CE", x) && exists("S_VLDL_C", x)) && exists("S_VLDL_CE_pct_C", x)) {
    x[, S_VLDL_CE_pct_C := NULL]
  }

  if (exists("XS_VLDL_CE", x) && exists("XS_VLDL_C", x)) {
    x[, XS_VLDL_CE_pct_C := collate_flags(XS_VLDL_CE, XS_VLDL_C)]
  } else if (!(exists("XS_VLDL_CE", x) && exists("XS_VLDL_C", x)) && exists("XS_VLDL_CE_pct_C", x)) {
    x[, XS_VLDL_CE_pct_C := NULL]
  }

  if (exists("L_LDL_CE", x) && exists("L_LDL_C", x)) {
    x[, L_LDL_CE_pct_C := collate_flags(L_LDL_CE, L_LDL_C)]
  } else if (!(exists("L_LDL_CE", x) && exists("L_LDL_C", x)) && exists("L_LDL_CE_pct_C", x)) {
    x[, L_LDL_CE_pct_C := NULL]
  }

  if (exists("M_LDL_CE", x) && exists("M_LDL_C", x)) {
    x[, M_LDL_CE_pct_C := collate_flags(M_LDL_CE, M_LDL_C)]
  } else if (!(exists("M_LDL_CE", x) && exists("M_LDL_C", x)) && exists("M_LDL_CE_pct_C", x)) {
    x[, M_LDL_CE_pct_C := NULL]
  }

  if (exists("S_LDL_CE", x) && exists("S_LDL_C", x)) {
    x[, S_LDL_CE_pct_C := collate_flags(S_LDL_CE, S_LDL_C)]
  } else if (!(exists("S_LDL_CE", x) && exists("S_LDL_C", x)) && exists("S_LDL_CE_pct_C", x)) {
    x[, S_LDL_CE_pct_C := NULL]
  }

  if (exists("IDL_CE", x) && exists("IDL_C", x)) {
    x[, IDL_CE_pct_C := collate_flags(IDL_CE, IDL_C)]
  } else if (!(exists("IDL_CE", x) && exists("IDL_C", x)) && exists("IDL_CE_pct_C", x)) {
    x[, IDL_CE_pct_C := NULL]
  }

  if (exists("XL_HDL_CE", x) && exists("XL_HDL_C", x)) {
    x[, XL_HDL_CE_pct_C := collate_flags(XL_HDL_CE, XL_HDL_C)]
  } else if (!(exists("XL_HDL_CE", x) && exists("XL_HDL_C", x)) && exists("XL_HDL_CE_pct_C", x)) {
    x[, XL_HDL_CE_pct_C := NULL]
  }

  if (exists("L_HDL_CE", x) && exists("L_HDL_C", x)) {
    x[, L_HDL_CE_pct_C := collate_flags(L_HDL_CE, L_HDL_C)]
  } else if (!(exists("L_HDL_CE", x) && exists("L_HDL_C", x)) && exists("L_HDL_CE_pct_C", x)) {
    x[, L_HDL_CE_pct_C := NULL]
  }

  if (exists("M_HDL_CE", x) && exists("M_HDL_C", x)) {
    x[, M_HDL_CE_pct_C := collate_flags(M_HDL_CE, M_HDL_C)]
  } else if (!(exists("M_HDL_CE", x) && exists("M_HDL_C", x)) && exists("M_HDL_CE_pct_C", x)) {
    x[, M_HDL_CE_pct_C := NULL]
  }

  if (exists("S_HDL_CE", x) && exists("S_HDL_C", x)) {
    x[, S_HDL_CE_pct_C := collate_flags(S_HDL_CE, S_HDL_C)]
  } else if (!(exists("S_HDL_CE", x) && exists("S_HDL_C", x)) && exists("S_HDL_CE_pct_C", x)) {
    x[, S_HDL_CE_pct_C := NULL]
  }


  # Compute ratio of free cholesterol to esterified cholesterol
  if (exists("XXL_VLDL_FC", x) && exists("XXL_VLDL_CE", x)) {
    x[, XXL_VLDL_FC_by_CE := collate_flags(XXL_VLDL_FC, XXL_VLDL_CE)]
  } else if (!(exists("XXL_VLDL_FC", x) && exists("XXL_VLDL_CE", x)) && exists("XXL_VLDL_FC_by_CE", x)) {
    x[, XXL_VLDL_FC_by_CE := NULL]
  }

  if (exists("XL_VLDL_FC", x) && exists("XL_VLDL_CE", x)) {
    x[, XL_VLDL_FC_by_CE := collate_flags(XL_VLDL_FC, XL_VLDL_CE)]
  } else if (!(exists("XL_VLDL_FC", x) && exists("XL_VLDL_CE", x)) && exists("XL_VLDL_FC_by_CE", x)) {
    x[, XL_VLDL_FC_by_CE := NULL]
  }

  if (exists("L_VLDL_FC", x) && exists("L_VLDL_CE", x)) {
    x[, L_VLDL_FC_by_CE := collate_flags(L_VLDL_FC, L_VLDL_CE)]
  } else if (!(exists("L_VLDL_FC", x) && exists("L_VLDL_CE", x)) && exists("L_VLDL_FC_by_CE", x)) {
    x[, L_VLDL_FC_by_CE := NULL]
  }

  if (exists("M_VLDL_FC", x) && exists("M_VLDL_CE", x)) {
    x[, M_VLDL_FC_by_CE := collate_flags(M_VLDL_FC, M_VLDL_CE)]
  } else if (!(exists("M_VLDL_FC", x) && exists("M_VLDL_CE", x)) && exists("M_VLDL_FC_by_CE", x)) {
    x[, M_VLDL_FC_by_CE := NULL]
  }

  if (exists("S_VLDL_FC", x) && exists("S_VLDL_CE", x)) {
    x[, S_VLDL_FC_by_CE := collate_flags(S_VLDL_FC, S_VLDL_CE)]
  } else if (!(exists("S_VLDL_FC", x) && exists("S_VLDL_CE", x)) && exists("S_VLDL_FC_by_CE", x)) {
    x[, S_VLDL_FC_by_CE := NULL]
  }

  if (exists("XS_VLDL_FC", x) && exists("XS_VLDL_CE", x)) {
    x[, XS_VLDL_FC_by_CE := collate_flags(XS_VLDL_FC, XS_VLDL_CE)]
  } else if (!(exists("XS_VLDL_FC", x) && exists("XS_VLDL_CE", x)) && exists("XS_VLDL_FC_by_CE", x)) {
    x[, XS_VLDL_FC_by_CE := NULL]
  }

  if (exists("L_LDL_FC", x) && exists("L_LDL_CE", x)) {
    x[, L_LDL_FC_by_CE := collate_flags(L_LDL_FC, L_LDL_CE)]
  } else if (!(exists("L_LDL_FC", x) && exists("L_LDL_CE", x)) && exists("L_LDL_FC_by_CE", x)) {
    x[, L_LDL_FC_by_CE := NULL]
  }

  if (exists("M_LDL_FC", x) && exists("M_LDL_CE", x)) {
    x[, M_LDL_FC_by_CE := collate_flags(M_LDL_FC, M_LDL_CE)]
  } else if (!(exists("M_LDL_FC", x) && exists("M_LDL_CE", x)) && exists("M_LDL_FC_by_CE", x)) {
    x[, M_LDL_FC_by_CE := NULL]
  }

  if (exists("S_LDL_FC", x) && exists("S_LDL_CE", x)) {
    x[, S_LDL_FC_by_CE := collate_flags(S_LDL_FC, S_LDL_CE)]
  } else if (!(exists("S_LDL_FC", x) && exists("S_LDL_CE", x)) && exists("S_LDL_FC_by_CE", x)) {
    x[, S_LDL_FC_by_CE := NULL]
  }

  if (exists("IDL_FC", x) && exists("IDL_CE", x)) {
    x[, IDL_FC_by_CE := collate_flags(IDL_FC, IDL_CE)]
  } else if (!(exists("IDL_FC", x) && exists("IDL_CE", x)) && exists("IDL_FC_by_CE", x)) {
    x[, IDL_FC_by_CE := NULL]
  }

  if (exists("XL_HDL_FC", x) && exists("XL_HDL_CE", x)) {
    x[, XL_HDL_FC_by_CE := collate_flags(XL_HDL_FC, XL_HDL_CE)]
  } else if (!(exists("XL_HDL_FC", x) && exists("XL_HDL_CE", x)) && exists("XL_HDL_FC_by_CE", x)) {
    x[, XL_HDL_FC_by_CE := NULL]
  }

  if (exists("L_HDL_FC", x) && exists("L_HDL_CE", x)) {
    x[, L_HDL_FC_by_CE := collate_flags(L_HDL_FC, L_HDL_CE)]
  } else if (!(exists("L_HDL_FC", x) && exists("L_HDL_CE", x)) && exists("L_HDL_FC_by_CE", x)) {
    x[, L_HDL_FC_by_CE := NULL]
  }

  if (exists("M_HDL_FC", x) && exists("M_HDL_CE", x)) {
    x[, M_HDL_FC_by_CE := collate_flags(M_HDL_FC, M_HDL_CE)]
  } else if (!(exists("M_HDL_FC", x) && exists("M_HDL_CE", x)) && exists("M_HDL_FC_by_CE", x)) {
    x[, M_HDL_FC_by_CE := NULL]
  }

  if (exists("S_HDL_FC", x) && exists("S_HDL_CE", x)) {
    x[, S_HDL_FC_by_CE := collate_flags(S_HDL_FC, S_HDL_CE)]
  } else if (!(exists("S_HDL_FC", x) && exists("S_HDL_CE", x)) && exists("S_HDL_FC_by_CE", x)) {
    x[, S_HDL_FC_by_CE := NULL]
  }


  # Compute lipid percentages in lipoprotein classes
  if (exists("VLDL_CE", x) && exists("VLDL_L", x)) {
    x[, VLDL_CE_pct := collate_flags(VLDL_CE, VLDL_L)]
  } else if (!(exists("VLDL_CE", x) && exists("VLDL_L", x)) && exists("VLDL_CE_pct", x)) {
    x[, VLDL_CE_pct := NULL]
  }

  if (exists("VLDL_FC", x) && exists("VLDL_L", x)) {
    x[, VLDL_FC_pct := collate_flags(VLDL_FC, VLDL_L)]
  } else if (!(exists("VLDL_FC", x) && exists("VLDL_L", x)) && exists("VLDL_FC_pct", x)) {
    x[, VLDL_FC_pct := NULL]
  }

  if (exists("VLDL_C", x) && exists("VLDL_L", x)) {
    x[, VLDL_C_pct := collate_flags(VLDL_C, VLDL_L)]
  } else if (!(exists("VLDL_C", x) && exists("VLDL_L", x)) && exists("VLDL_C_pct", x)) {
    x[, VLDL_C_pct := NULL]
  }

  if (exists("VLDL_PL", x) && exists("VLDL_L", x)) {
    x[, VLDL_PL_pct := collate_flags(VLDL_PL, VLDL_L)]
  } else if (!(exists("VLDL_PL", x) && exists("VLDL_L", x)) && exists("VLDL_PL_pct", x)) {
    x[, VLDL_PL_pct := NULL]
  }

  if (exists("VLDL_TG", x) && exists("VLDL_L", x)) {
    x[, VLDL_TG_pct := collate_flags(VLDL_TG, VLDL_L)]
  } else if (!(exists("VLDL_TG", x) && exists("VLDL_L", x)) && exists("VLDL_TG_pct", x)) {
    x[, VLDL_TG_pct := NULL]
  }


  if (exists("LDL_CE", x) && exists("LDL_L", x)) {
    x[, LDL_CE_pct := collate_flags(LDL_CE, LDL_L)]
  } else if (!(exists("LDL_CE", x) && exists("LDL_L", x)) && exists("LDL_CE_pct", x)) {
    x[, LDL_CE_pct := NULL]
  }

  if (exists("LDL_FC", x) && exists("LDL_L", x)) {
    x[, LDL_FC_pct := collate_flags(LDL_FC, LDL_L)]
  } else if (!(exists("LDL_FC", x) && exists("LDL_L", x)) && exists("LDL_FC_pct", x)) {
    x[, LDL_FC_pct := NULL]
  }

  if (exists("LDL_C", x) && exists("LDL_L", x)) {
    x[, LDL_C_pct := collate_flags(LDL_C, LDL_L)]
  } else if (!(exists("LDL_C", x) && exists("LDL_L", x)) && exists("LDL_C_pct", x)) {
    x[, LDL_C_pct := NULL]
  }

  if (exists("LDL_PL", x) && exists("LDL_L", x)) {
    x[, LDL_PL_pct := collate_flags(LDL_PL, LDL_L)]
  } else if (!(exists("LDL_PL", x) && exists("LDL_L", x)) && exists("LDL_PL_pct", x)) {
    x[, LDL_PL_pct := NULL]
  }

  if (exists("LDL_TG", x) && exists("LDL_L", x)) {
    x[, LDL_TG_pct := collate_flags(LDL_TG, LDL_L)]
  } else if (!(exists("LDL_TG", x) && exists("LDL_L", x)) && exists("LDL_TG_pct", x)) {
    x[, LDL_TG_pct := NULL]
  }


  if (exists("HDL_CE", x) && exists("HDL_L", x)) {
    x[, HDL_CE_pct := collate_flags(HDL_CE, HDL_L)]
  } else if (!(exists("HDL_CE", x) && exists("HDL_L", x)) && exists("HDL_CE_pct", x)) {
    x[, HDL_CE_pct := NULL]
  }

  if (exists("HDL_FC", x) && exists("HDL_L", x)) {
    x[, HDL_FC_pct := collate_flags(HDL_FC, HDL_L)]
  } else if (!(exists("HDL_FC", x) && exists("HDL_L", x)) && exists("HDL_FC_pct", x)) {
    x[, HDL_FC_pct := NULL]
  }

  if (exists("HDL_C", x) && exists("HDL_L", x)) {
    x[, HDL_C_pct := collate_flags(HDL_C, HDL_L)]
  } else if (!(exists("HDL_C", x) && exists("HDL_L", x)) && exists("HDL_C_pct", x)) {
    x[, HDL_C_pct := NULL]
  }

  if (exists("HDL_PL", x) && exists("HDL_L", x)) {
    x[, HDL_PL_pct := collate_flags(HDL_PL, HDL_L)]
  } else if (!(exists("HDL_PL", x) && exists("HDL_L", x)) && exists("HDL_PL_pct", x)) {
    x[, HDL_PL_pct := NULL]
  }

  if (exists("HDL_TG", x) && exists("HDL_L", x)) {
    x[, HDL_TG_pct := collate_flags(HDL_TG, HDL_L)]
  } else if (!(exists("HDL_TG", x) && exists("HDL_L", x)) && exists("HDL_TG_pct", x)) {
    x[, HDL_TG_pct := NULL]
  }


  # Compute cholesterol percentages
  if (exists("VLDL_FC", x) && exists("VLDL_C", x)) {
    x[, VLDL_FC_pct_C := collate_flags(VLDL_FC, VLDL_C)]
  } else if (!(exists("VLDL_FC", x) && exists("VLDL_C", x)) && exists("VLDL_FC_pct_C", x)) {
    x[, VLDL_FC_pct_C := NULL]
  }

  if (exists("LDL_FC", x) && exists("LDL_C", x)) {
    x[, LDL_FC_pct_C := collate_flags(LDL_FC, LDL_C)]
  } else if (!(exists("LDL_FC", x) && exists("LDL_C", x)) && exists("LDL_FC_pct_C", x)) {
    x[, LDL_FC_pct_C := NULL]
  }

  if (exists("HDL_FC", x) && exists("HDL_C", x)) {
    x[, HDL_FC_pct_C := collate_flags(HDL_FC, HDL_C)]
  } else if (!(exists("HDL_FC", x) && exists("HDL_C", x)) && exists("HDL_FC_pct_C", x)) {
    x[, HDL_FC_pct_C := NULL]
  }


  if (exists("VLDL_CE", x) && exists("VLDL_C", x)) {
    x[, VLDL_CE_pct_C := collate_flags(VLDL_CE, VLDL_C)]
  } else if (!(exists("VLDL_CE", x) && exists("VLDL_C", x)) && exists("VLDL_CE_pct_C", x)) {
    x[, VLDL_CE_pct_C := NULL]
  }

  if (exists("LDL_CE", x) && exists("LDL_C", x)) {
    x[, LDL_CE_pct_C := collate_flags(LDL_CE, LDL_C)]
  } else if (!(exists("LDL_CE", x) && exists("LDL_C", x)) && exists("LDL_CE_pct_C", x)) {
    x[, LDL_CE_pct_C := NULL]
  }

  if (exists("HDL_CE", x) && exists("HDL_C", x)) {
    x[, HDL_CE_pct_C := collate_flags(HDL_CE, HDL_C)]
  } else if (!(exists("HDL_CE", x) && exists("HDL_C", x)) && exists("HDL_CE_pct_C", x)) {
    x[, HDL_CE_pct_C := NULL]
  }


  # Ratios of free cholesterol to esterified cholesterol:
  if (exists("VLDL_FC", x) && exists("VLDL_CE", x)) {
    x[, VLDL_FC_by_CE := collate_flags(VLDL_FC, VLDL_CE)]
  } else if (!(exists("VLDL_FC", x) && exists("VLDL_CE", x)) && exists("VLDL_FC_by_CE", x)) {
    x[, VLDL_FC_by_CE := NULL]
  }

  if (exists("LDL_FC", x) && exists("LDL_CE", x)) {
    x[, LDL_FC_by_CE := collate_flags(LDL_FC, LDL_CE)]
  } else if (!(exists("LDL_FC", x) && exists("LDL_CE", x)) && exists("LDL_FC_by_CE", x)) {
    x[, LDL_FC_by_CE := NULL]
  }

  if (exists("HDL_FC", x) && exists("HDL_CE", x)) {
    x[, HDL_FC_by_CE := collate_flags(HDL_FC, HDL_CE)]
  } else if (!(exists("HDL_FC", x) && exists("HDL_CE", x)) && exists("HDL_FC_by_CE", x)) {
    x[, HDL_FC_by_CE := NULL]
  }


  # Lipid and cholesterol fractions in total serum
  if (exists("Total_CE", x) && exists("Total_L", x)) {
    x[, Total_CE_pct := collate_flags(Total_CE, Total_L)]
  } else if (!(exists("Total_CE", x) && exists("Total_L", x)) && exists("Total_CE_pct", x)) {
    x[, Total_CE_pct := NULL]
  }

  if (exists("Total_FC", x) && exists("Total_L", x)) {
    x[, Total_FC_pct := collate_flags(Total_FC, Total_L)]
  } else if (!(exists("Total_FC", x) && exists("Total_L", x)) && exists("Total_FC_pct", x)) {
    x[, Total_FC_pct := NULL]
  }

  if (exists("Total_C", x) && exists("Total_L", x)) {
    x[, Total_C_pct := collate_flags(Total_C, Total_L)]
  } else if (!(exists("Total_C", x) && exists("Total_L", x)) && exists("Total_C_pct", x)) {
    x[, Total_C_pct := NULL]
  }

  if (exists("Total_PL", x) && exists("Total_L", x)) {
    x[, Total_PL_pct := collate_flags(Total_PL, Total_L)]
  } else if (!(exists("Total_PL", x) && exists("Total_L", x)) && exists("Total_PL_pct", x)) {
    x[, Total_PL_pct := NULL]
  }

  if (exists("Total_TG", x) && exists("Total_L", x)) {
    x[, Total_TG_pct := collate_flags(Total_TG, Total_L)]
  } else if (!(exists("Total_TG", x) && exists("Total_L", x)) && exists("Total_TG_pct", x)) {
    x[, Total_TG_pct := NULL]
  }


  if (exists("Total_FC", x) && exists("Total_C", x)) {
    x[, Total_FC_pct_C := collate_flags(Total_FC, Total_C)]
  } else if (!(exists("Total_FC", x) && exists("Total_C", x)) && exists("Total_FC_pct_C", x)) {
    x[, Total_FC_pct_C := NULL]
  }

  if (exists("Total_CE", x) && exists("Total_C", x)) {
    x[, Total_CE_pct_C := collate_flags(Total_CE, Total_C)]
  } else if (!(exists("Total_CE", x) && exists("Total_C", x)) && exists("Total_CE_pct_C", x)) {
    x[, Total_CE_pct_C := NULL]
  }

  if (exists("Total_FC", x) && exists("Total_CE", x)) {
    x[, Total_FC_by_CE := collate_flags(Total_FC, Total_CE)]
  } else if (!(exists("Total_FC", x) && exists("Total_CE", x)) && exists("Total_FC_by_CE", x)) {
    x[, Total_FC_by_CE := NULL]
  }


  # Omega 3 and Omega 6 percent of PUFA
  if (exists("Omega_3", x) && exists("PUFA", x)) {
    x[, Omega_3_pct_PUFA := collate_flags(Omega_3, PUFA)]
  } else if (!(exists("Omega_3", x) && exists("PUFA", x)) && exists("Omega_3_pct_PUFA", x)) {
    x[, Omega_3_pct_PUFA := NULL]
  }

  if (exists("Omega_6", x) && exists("PUFA", x)) {
    x[, Omega_6_pct_PUFA := collate_flags(Omega_6, PUFA)]
  } else if (!(exists("Omega_6", x) && exists("PUFA", x)) && exists("Omega_6_pct_PUFA", x)) {
    x[, Omega_6_pct_PUFA := NULL]
  }

  # Finished
  return(x)
}
