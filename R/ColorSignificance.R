# @title ColorSignificance internal function
# @description This function colors statistically significant values.
#' @import hwriter
#' @import kableExtra
# @return returns ColoredTable, sig_matrix.

ColorSignificance <- function(tabletoColor=tabletoColor, Klimek=FALSE) {
  asc <- function(z) { as.character(z) }
  KSimNames <-
    paste("KSim", c("I","E","alpha","turnout","winprop","sigma","stdAtt","theta"),sep="");
  KLikNames <-
    paste("KLik", c("I","E","alpha","turnout","winprop","sigma","stdAtt","theta","LL"),sep="");
  MethNames<-c("_2BL", "LastC", "P05s", "C05s", "Skew", "Kurt", "DipT", "Sobyanin", "Obs");

  tab <- tabletoColor;
  #testfunc(cnames[j],estrows[i],j)
  #z=cnames[j]; i=estrows[i];j
  testfunc <- function(z,i,j) {
    TF <- FALSE;
    if (z %in% c("_2BL", "LastC", "P05s", "C05s", "Skew", "Kurt")) {
      z <- gsub("_","X",z);
      cmd <- paste("c",asc(tab[i+1,j]),sep="");
      if (cmd != "c--" & cmd != "c(0, 0)") {
        CI <- eval(parse(text=cmd));
        m <- switch(z,
                    X2BL = 4.187,
                    LastC = 4.5,
                    P05s = .2,
                    C05s = .2,
                    Skew = 0,
                    Kurt = 3,
                    Correlation = 0.4);
        TF <- CI[2] <= m | CI[1] >= m;
      }
    } else if (z == "DipT") {
      cmd <- asc(tab[i,j]);
      if (cmd != "--") {
        est <- eval(parse(text=cmd));
        TF <- est < 1-pnorm(qnorm(.95));
      }
    } else if (z == "Sobyanin") {
      cmd <- asc(tab[i,j]);
      if (cmd != "--") {
        est <- eval(parse(text=cmd));
        TF <- abs(est) > qnorm(.95);
      }
    } else if (z == "Correlation") {
      cmd <- asc(tab[i,j]);
      if (cmd != "--") {
        est <- eval(parse(text=cmd));
        TF <- abs(est) > 0.4;
      }
    } else if (z %in% paste("KSim", c("I","E"), sep="")) {
      cmd <- asc(tab[i,j]);
      cmd2 <- asc(tab[i+1,j]);
      if (cmd != "--" & cmd2 != "--") {
        est <- eval(parse(text=cmd));
        sdev <- eval(parse(text=cmd2));
        if (sdev != 0) {
          TF <- 1-pnorm(abs(est/sdev)) < 1-pnorm(qnorm(.975));
        }
      }
    } else if (z %in% paste("KLik",
                            c("I","E","alpha","turnout","winprop","sigma","stdAtt","theta"), sep="")) {
      TF <- NA;
      if (!(tab[i,2] %in% c("Turnout",""))) {
        ztol <- 1e-9;
        Icol <- (1:dim(tabletoColor)[2])[dimnames(tabletoColor)[[2]] == "KLikI"];
        Ecol <- (1:dim(tabletoColor)[2])[dimnames(tabletoColor)[[2]] == "KLikE"];
        Iidx <- eval(parse(text=asc(tab[i,Icol]))) > ztol;
        Eidx <- eval(parse(text=asc(tab[i,Ecol]))) > ztol;
        nparms <- 4 - ifelse(Iidx&Eidx, 0,
                             ifelse(Iidx, 1,
                                    ifelse(Eidx, 2, 4)));
        LLcol <- (1:dim(tabletoColor)[2])[dimnames(tabletoColor)[[2]] == "KLikLL"];
        if (tab[i,LLcol]!="--") {
          LL <- 2*eval(parse(text=asc(tab[i,LLcol])));
          TF <- ifelse(nparms==0 | LL<0, FALSE, 1-pchisq(LL,nparms) < 1-pnorm(qnorm(.975)));
        }
      }
    }
    TF;
  };

  point_mrs <- tabletoColor[seq(1,dim(tabletoColor)[1],2),]
  sig_matrix <- matrix(NA,dim(point_mrs)[1],dim(point_mrs)[2])

  cnames <- dimnames(tabletoColor)[[2]];
  estrows <- seq(1,dim(tabletoColor)[1],2);
  for (i in 1:dim(sig_matrix)[1]) {
    for (j in 1:(dim(tabletoColor)[2]-1)) {
      if(Klimek==FALSE){sig_matrix[i,j] <- ifelse(testfunc(cnames[j],estrows[i],j),1,0);
      }else{sig_matrix[i,j]=FALSE}
    }
  }

  sig_matrixF <- sig_matrix[rep(1:nrow(sig_matrix),each=2),]
  sig_matrixF[estrows+1,] <- 0;
  sig_matrixF[sig_matrixF==1] <- "color:red"

  ColoredTable <-hwrite(tabletoColor, row.style=list('font-weight:bold'),
                        style=sig_matrixF,table.style='border-collapse:collapse',
                        table.cellpadding='15px')
  ColoredTable<-kable(ColoredTable, escape = F, caption="Election Forensics Results")
  ColoredTable<-gsub("> x <", "> <", ColoredTable)
  ColoredTable<-kable_styling(ColoredTable, "striped", full_width = F)
  ColoredTable
  resColoredTable<-list(ColoredTable, sig_matrixF)
  return(resColoredTable)
}
