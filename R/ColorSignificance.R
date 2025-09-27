# @title ColorSignificance internal function
# @description This function colors statistically significant values.
#' @import hwriter
#' @import kableExtra
#' @import DT
# @return returns ColoredTable, sig_matrix.

ColorSignificance <- function(tabletoColor, Klimek = FALSE) {
  asc <- function(z) as.character(z)

  KSimNames <- paste0("KSim", c("I", "E", "alpha", "turnout", "winprop", "sigma", "stdAtt", "theta"))
  KLikNames <- paste0("KLik", c("I", "E", "alpha", "turnout", "winprop", "sigma", "stdAtt", "theta", "LL"))
  MethNames <- c("_2BL", "LastC", "P05s", "C05s", "Skew", "Kurt", "DipT", "Sobyanin", "Obs")

  tab <- tabletoColor

  testfunc <- function(z, i, j) {
    TF <- FALSE
    if (z %in% c("_2BL", "LastC", "P05s", "C05s", "Skew", "Kurt")) {
      z <- gsub("_", "X", z)
      cmd <- paste("c", asc(tab[i + 1, j]), sep = "")
      if (cmd != "c--" & cmd != "c(0, 0)") {
        CI <- eval(parse(text = cmd))
        m <- switch(z,
                    X2BL = 4.187, LastC = 4.5, P05s = 0.2, C05s = 0.2,
                    Skew = 0, Kurt = 3, Correlation = 0.4
        )
        TF <- CI[2] <= m | CI[1] >= m
      }
    } else if (z == "DipT") {
      cmd <- asc(tab[i, j])
      if (cmd != "--") {
        est <- eval(parse(text = cmd))
        TF <- est < 1 - pnorm(qnorm(0.95))
      }
    } else if (z == "Sobyanin") {
      cmd <- asc(tab[i, j])
      if (cmd != "--") {
        est <- eval(parse(text = cmd))
        TF <- abs(est) > qnorm(0.95)
      }
    } else if (z == "Correlation") {
      cmd <- asc(tab[i, j])
      if (cmd != "--") {
        est <- eval(parse(text = cmd))
        TF <- abs(est) > 0.4
      }
    } else if (z %in% paste0("KSim", c("I", "E"))) {
      cmd <- asc(tab[i, j])
      cmd2 <- asc(tab[i + 1, j])
      if (cmd != "--" & cmd2 != "--") {
        est <- eval(parse(text = cmd))
        sdev <- eval(parse(text = cmd2))
        if (sdev != 0) {
          TF <- 1 - pnorm(abs(est / sdev)) < 1 - pnorm(qnorm(0.975))
        }
      }
    } else if (z %in% paste0("KLik", c("I", "E", "alpha", "turnout", "winprop", "sigma", "stdAtt", "theta"))) {
      TF <- NA
      if (!(tab[i, 2] %in% c("Turnout", ""))) {
        ztol <- 1e-9
        Icol <- which(colnames(tab) == "KLikI")
        Ecol <- which(colnames(tab) == "KLikE")
        Iidx <- eval(parse(text = asc(tab[i, Icol]))) > ztol
        Eidx <- eval(parse(text = asc(tab[i, Ecol]))) > ztol
        nparms <- 4 - ifelse(Iidx & Eidx, 0,
                             ifelse(Iidx, 1, ifelse(Eidx, 2, 4)))
        LLcol <- which(colnames(tab) == "KLikLL")
        if (tab[i, LLcol] != "--") {
          LL <- 2 * eval(parse(text = asc(tab[i, LLcol])))
          TF <- ifelse(nparms == 0 | LL < 0, FALSE,
                       1 - pchisq(LL, nparms) < 1 - pnorm(qnorm(0.975)))
        }
      }
    }
    return(TF)
  }

  ColorTabSignificance <- function(data_table, sig_matrixF) {
    # Convert to data frame and remove any completely empty rows
    df <- as.data.frame(data_table, stringsAsFactors = FALSE)
    df <- df[rowSums(df == "" | is.na(df)) != ncol(df), ]  # Remove rows that are all empty/NA

    # Create HTML version with marked cells
    marked_df <- df
    for (i in 1:min(nrow(sig_matrixF), ceiling(nrow(df)/2))) {  # Ensure we don't exceed matrix bounds
      for (j in 1:min(ncol(sig_matrixF), ncol(df))) {        # Ensure we don't exceed matrix bounds
        if (!is.na(sig_matrixF[i,j]) && sig_matrixF[i,j] == 1) {
          # Mark coefficient row (odd rows)
          row_idx <- i*2 - 1
          if (row_idx <= nrow(df)) {  # Additional safety check
            marked_df[row_idx, j] <- sprintf(
              '<span style="color:red; background-color:#FFCCCC">%s</span>',
              df[row_idx, j]
            )
          }
        }
      }
    }

    # Create datatable
    datatable(
      marked_df,
      escape = FALSE,
      options = list(
        paging = FALSE,
        searching = TRUE,
        info = FALSE,
        # Additional option to prevent empty row creation
        drawCallback = DT::JS("function(settings) {
                var api = this.api();
                $(api.table().node()).removeClass('no-footer');
            }")
      )
    )
  }

  # --- Generate significance matrix ---
  estrows <- seq(1, nrow(tab), 2)
  sig_matrix <- matrix(NA, nrow = length(estrows), ncol = ncol(tab))
  cnames <- colnames(tab)

  for (i in seq_along(estrows)) {
    for (j in seq_len(ncol(tab) - 1)) {
      if (!Klimek) {
        sig_matrix[i, j] <- ifelse(testfunc(cnames[j], estrows[i], j), 1, 0)
      } else {
        sig_matrix[i, j] <- FALSE
      }
    }
  }

  sig_matrixF <- sig_matrix[rep(seq_len(nrow(sig_matrix)), each = 2), ]
  sig_matrixF[estrows + 1, ] <- 0

  ColoredTable <- ColorTabSignificance(tab, sig_matrixF)
  resColoredTable<-list(ColoredTable, sig_matrixF)

  return(resColoredTable)}
