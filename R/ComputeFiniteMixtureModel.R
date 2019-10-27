#' @title Finite Mixture Model
#' @description This function implements the Finite Mixture Model developed and coded by Walter Mebane.
#' @usage ComputeFiniteMixtureModel(dat, MainCandidate="Votes", TotalReg="NVoters",
#'                                  TotalVotes="NValid", cores=2, itstartmax=1)
#' @param dat electoral data.
#' @param MainCandidate variable name for major/incumbent candidate
#' @param TotalReg variable name for the total number of eligible voters
#' @param TotalVotes variable name for the total number of ballots cast
#' @param cores number of cores for parallel computing (2 cores by default)
#' @param itstartmax maximum number of iterations
#' @export
#' @import doParallel
#' @import rgenoud
#' @return A list of FMM estimates.
#' @examples
#' library(EFToolkit)
#'
#' #Russia examples
#' dat<-read.csv(system.file("ruspres2018.csv", package="EFToolkit"))
#' dat<-subset(dat, select=c("region", "NVoters", "NValid", "Votes"))
#' datc<-dat[dat$region=="Volgogradskaya Oblast`",]
#' #Do not run: takes long time to compute
#' #res<-ComputeFiniteMixtureModel(datc, MainCandidate="Votes", TotalReg="NVoters", TotalVotes="NValid")

###################################
###2015/7/18
##Klimek-like Methods
##  Walter Mebane and Naoki Egami
###################################

ComputeFiniteMixtureModel <- function(dat, MainCandidate="Votes", TotalReg="NVoters", TotalVotes="NValid", cores=2, itstartmax=1){
  iterations <<- 5

  cleandata<-function(dat){
    dat$Votes <- dat[, names(dat)%in%MainCandidate]
    dat$Voters <- dat[, names(dat)%in%TotalReg]
    dat$NValid <- dat[, names(dat)%in%TotalVotes]

    kidx <- !is.na(dat$NVoters) & (dat$NVoters >= dat$NValid) &
      (dat$NValid > 0) & (dat$NValid >= dat$Votes);
    if (any(is.na(kidx))) kidx[is.na(kidx)] <- FALSE;

    if(length(kidx)!=0){dat <- dat[kidx,]}
    return(dat)}


  #FMM functions
  IterationN <- function(iterations,data1,cores,cluster,pop=5000,itprint=0,itstartmax=10) {

    # absolute value of Normal density
    danorm <- function(x, m, s) {
      s2 <- s^2;
      1/(s*sqrt(2*pi))* (exp(-(-x-m)^2/(2*s2)) + exp(-(x-m)^2/(2*s2)));
    }

    # absolute value of Normal density with zero mean (half-Normal)
    danorm0 <- function(x, s) {
      2/(s*sqrt(2*pi))* exp(-x^2/(2*s^2));
    }

    erf <- function(x) {
      2*pnorm(sqrt(2)*x) - 1
    }

    # absolute value of Normal cdf with zero mean (half-Normal)
    Fanorm0 <- function(x, s) {
      #    erf(x/sqrt(2*s));
      2*pnorm(sqrt(2)*(x/sqrt(2*s)))-1;
    }

    # Normal by half-normal density
    # for x >= 0
    # Epstein (1948) equation (2);  Curtiss (1941)
    dconnHNI <- function(x, NV, m1, s1, m2, s2){
      F0 <- pnorm(1,mean=m2,sd=s2)-pnorm(0,mean=m2,sd=s2);
      F1 <- Fanorm0(1, s1);
      f <- function(y, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2) {
        ifelse(y==1,0,dnorm(x/(1-y) ,mean=NV*(1-m2),sd=NV*s2) * danorm(y,m1,s1)/(1-y)/(F0*F1));
      }
      intf <- integrate(f, 0, 1, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    # Normal by 1-half-normal density
    # for x >= 0
    # Epstein (1948) equation (2);  Curtiss (1941)
    dconnHNE <- function(x, NV, m1, s1, m2, s2){
      F0 <- pnorm(1,mean=m2,sd=s2)-pnorm(0,mean=m2,sd=s2);
      F1 <- Fanorm0(1, s1);
      f <- function(y, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2) {
        ifelse(y==0,0,dnorm(x/y ,mean=NV*(1-m2),sd=NV*s2) * danorm(y,m1,s1)/y/(F0*F1));
      }
      intf <- integrate(f, 0, 1, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    # no "fraud" joint density
    dconv <- function(xW, xA, NV, m21, s21, m22, s22) {
      F0 <- pnorm(1,mean=m21,sd=s21)-pnorm(0,mean=m21,sd=s21);
      F2 <- pnorm(1,mean=m22,sd=s22)-pnorm(0,mean=m22,sd=s22);
      F02 <- F0*F2^2;
      a <- 1-xA/NV;
      v <-
        dnorm(xW/a, NV*m21, NV*s21) *
        dnorm(a, m22, s22) *
        dnorm(xA, NV*(1-m22), NV*s22)/a/F02;
      ifelse(is.finite(v),v,0);
    }

    # Inc "fraud" joint density
    dconvI <- function(xW, NV, m1, s1, m21, s21, m22, s22, e) {
      F0 <- pnorm(1,mean=m21,sd=s21)-pnorm(0,mean=m21,sd=s21);
      F1 <- Fanorm0(1, s1);
      F2 <- pnorm(1,mean=m22,sd=s22)-pnorm(0,mean=m22,sd=s22);
      F02 <- F0*F2;
      g <- function(z,y, xW, NV, m21, s21, m22, s22, ye, ye1, dy){
        # X = L*A*(1 - U^q) + U*(1-A) + U^q*A
        zye <- z*ye;
        v <- ifelse(z==0|y==1,0,
                    dnorm(xW/zye, NV*(m21 + y*(1-z)/zye + ye1), NV*s21) *
                      dy * dnorm(z, m22, s22)/zye/F02);
        ifelse(is.finite(v),v,0);
      }
      f <- function(y,xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e){
        ye <- 1-y^e;
        ye1 <- y^(e)/(1-y^(e));
        dy <- (1-y)^(-1+1/e)*danorm0(y, s1)/e/F1;
        intg <- integrate(g, 0, 1,
                          y=y, xW=xW,NV=NV,m21=m21,s21=s21,m22=m22,s22=s22,ye=ye,ye1=ye1,dy=dy,
                          stop.on.error=FALSE);
        ret <- ifelse(intg$message=="OK",intg$value, 0);
        if (length(ret)==1) ret <- rep(ret,length(y));
        ret
      }
      intf <- integrate(f, 0, 1,
                        xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    # Ext "fraud" joint density
    dconvE <- function(xW, NV, m1, s1, m21, s21, m22, s22, e) {
      F0 <- pnorm(1,mean=m21,sd=s21)-pnorm(0,mean=m21,sd=s21);
      F1 <- Fanorm0(1, s1);
      F2 <- pnorm(1,mean=m22,sd=s22)-pnorm(0,mean=m22,sd=s22);
      F02 <- F0*F2;
      g <- function(z,y, xW, NV, m21, s21, m22, s22, ye, ye1, dy){
        # X = L*A*(1 - U^q) + U*(1-A) + U^q*A
        zye <- z*ye;
        v <- ifelse(z==0|y==0,0,
                    dnorm(xW/zye, NV*(m21 + (1-y)*(1-z)/zye + ye1), NV*s21) *
                      dy * dnorm(z, m22, s22)/zye/F02);
        ifelse(is.finite(v),v,0);
      }
      f <- function(y,xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e){
        ye <- 1-(1-y)^e;
        ye1 <- (1-y)^(e)/(1-(1-y)^(e));
        dy <- (1-y)^(-1+1/e)*danorm0(y, s1)/e/F1;
        intg <- integrate(g, 0, 1,
                          y=y, xW=xW,NV=NV,m21=m21,s21=s21,m22=m22,s22=s22,ye=ye,ye1=ye1,dy=dy,
                          stop.on.error=FALSE);
        ret <- ifelse(intg$message=="OK",intg$value, 0);
        if (length(ret)==1) ret <- rep(ret,length(y));
        ret
      }
      intf <- integrate(f, 0, 1,
                        xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    EvalFit <- function(par,mdat) {
      alpha        <- par[1];
      p_Att        <- par[2];
      lambda_fraud <- par[3];
      sigma        <- par[4];
      stdAtt       <- par[5];
      theta        <- par[6];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      fitI <-
        sapply(loc,function(i) {
          dconvI(mdat[i,1], NV[i], 0, theta, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNI(mdat[i,3], NV[i], 0, theta, p_Att, stdAtt)
        });

      fitE <-
        sapply(loc,function(i) {
          dconvE(mdat[i,1], NV[i], 0, 0.075, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNE(mdat[i,3], NV[i], 0, 0.075, p_Att, stdAtt)
        });

      return(list(fit0=fit0,fitI=fitI,fitE=fitE));
    }

    EvalFitI <- function(par,mdat) {
      alpha        <- par[1];
      p_Att        <- par[2];
      lambda_fraud <- par[3];
      sigma        <- par[4];
      stdAtt       <- par[5];
      theta        <- par[6];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      fitI <-
        sapply(loc,function(i) {
          dconvI(mdat[i,1], NV[i], 0, theta, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNI(mdat[i,3], NV[i], 0, theta, p_Att, stdAtt)
        });

      return(list(fit0=fit0,fitI=fitI));
    }

    EvalFitE <- function(par,mdat) {
      alpha        <- par[1];
      p_Att        <- par[2];
      lambda_fraud <- par[3];
      sigma        <- par[4];
      stdAtt       <- par[5];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      fitE <-
        sapply(loc,function(i) {
          dconvE(mdat[i,1], NV[i], 0, 0.075, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNE(mdat[i,3], NV[i], 0, 0.075, p_Att, stdAtt)
        });

      return(list(fit0=fit0,fitE=fitE));
    }

    EvalFit0 <- function(par,mdat) {
      p_Att        <- par[1];
      lambda_fraud <- par[2];
      sigma        <- par[3];
      stdAtt       <- par[4];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      return(list(fit0=fit0));
    }

    GetfStarts <- function(fitlist, itprint) {
      #    while (TRUE) {
      #      fi <- runif(1,.1,.95);
      #      fe <- runif(1,.1,.3);
      #      if (fi+fe<1) break;
      #    }
      fi <- fe <- 0;
      if (itprint>1) cat("0starts:\n",fi,fe,"\n");
      fit <- (1-fi-fe)*fitlist$fit0 + fi*fitlist$fitI + fe*fitlist$fitE;
      h0 <- (1-fi-fe)*fitlist$fit0/fit;
      hI <- fi*fitlist$fitI/fit;
      hE <- fe*fitlist$fitE/fit;
      #    fi <- mean(hI);
      #    fe <- mean(hE);
      if (itprint>1) cat("1starts:\n",fi,fe,"\n");
      return(list(fi=fi,fe=fe,h0=h0,hI=hI,hE=hE));
    }

    GetStartsIteration <- function(itstartmax,bmin,bmax,data1,mdat, itprint) {
      fit <- fitmax <- (-1e+50);
      bstart0 <- (bmin+bmax)/2;
      fstartlist0 <- list();
      itmax <- 10^2;
      for (itstart in 1:itstartmax) {
        for (i in 1:itmax) {
          for (j in 1:itmax) {
            alpha <- runif(1,bmin[1],bmax[1]);
            #          p_Att <- runif(1,bmin[2],bmax[2]);
            #          lambda_fraud <- runif(1,bmin[3],bmax[3]);
            t <- data1$NValid/data1$NVoters;
            p_Att <- .99*min(mean(t),median(t));
            t <- data1$Votes/data1$NValid;
            lambda_fraud <- .99*min(mean(t),median(t));
            lambda_fraud <- c(lambda_fraud, 1-lambda_fraud);
            sigma <- runif(1,bmin[4],bmax[4]);
            stdAtt <- runif(1,bmin[5],bmax[5]);
            theta <- runif(1,bmin[6],bmax[6]);
            bstart <- c(alpha,p_Att,lambda_fraud[1],sigma,stdAtt,theta);

            fstartlist <- GetfStarts(EvalFit(bstart,mdat), itprint=itprint);
            h0 <- fstartlist$h0;
            hI <- fstartlist$hI;
            hE <- fstartlist$hE;
            #          fi <- fstartlist$fi;
            #          fe <- fstartlist$fe;
            #          if (is.finite(fi) & is.finite(fe)) break;
            if (all(is.finite(h0))) break;
          }
          fitlist <- EvalFit(bstart,mdat);
          #        fit <- sum(log(h0*fitlist$fit0 + hI*fitlist$fitI + hE*fitlist$fitE));
          fit <- sum(log(fitlist$fit0));
          if (is.finite(fit)) { break; }
        }
        if (is.finite(fit) && is.finite(fitmax) && fit > fitmax) {
          if (itprint>0) print(c(itstart,fit));
          fitmax <- fit;
          bstart0 <- bstart;
          fstartlist0 <- fstartlist;
        }
      }
      return(list(bstart=bstart0, fit=fitmax, fstartlist=fstartlist0));
    }

    GetStarts <- function(cores,bmin,bmax,data1,mdat, itprint, itstartmax=10) {
      Icores <- cores;  # number of processors to use (processes to spawn)
      if (itstartmax > Icores) {
        Ibase <- trunc(itstartmax/Icores);
        Ires  <- itstartmax %% Icores;
        Ivec  <- rep(Ibase,Icores);
        Ivec[1:Ires] <- Ivec[1:Ires] + 1;
      } else {
        Ivec <- rep(1,itstartmax);
      }
      while (TRUE) {
        # prepare to use foreach()
        registerDoParallel(cores=cores);
        StartsList <- foreach(itstartMAX=Ivec, .export='GetStartsIteration') %dopar%{
          GetStartsIteration(itstartMAX,bmin,bmax,data1,mdat, itprint)};
        stopImplicitCluster();
        if (itprint>0) { print("starts report"); }
        jmax <- 1;
        for (j in 1:itstartmax) {
          if (itprint>0) { print(c(j,StartsList[[j]]$fit)); }
          if (StartsList[[j]]$fit > StartsList[[jmax]]$fit) jmax <- j;
        }
        fitmax <- StartsList[[jmax]]$fit;
        if (fitmax > (-1e+50)) break;
      }
      bstart0 <- StartsList[[jmax]]$bstart;
      fstartlist0 <- StartsList[[jmax]]$fstartlist;
      cat("starting fit:  ",fitmax,"\n");
      return(list(bstart=bstart0, fit=fitmax, fstartlist=fstartlist0));
    }

    ###Set FraudFits
    FraudFits <- matrix(NA,nrow=8+2,ncol=0);

    tol <- 1e-3;
    ztol <- 1e-9;
    if (itprint>0) cat("zero tolerance: ",ztol,"\n");
    mdat <- cbind(data1$Votes,data1$NValid-data1$Votes,data1$NVoters-data1$NValid);
    AttendanceQ <- quantile(Attendance <- data1$NValid/data1$NVoters);
    WinpropQ <- quantile(Winprop <- data1$Votes/data1$NValid);
    sdwin <- 2*sqrt(var(Winprop[Winprop<=WinpropQ[3]]));
    sdbal <- 2*sqrt(var(Attendance[Attendance<=AttendanceQ[3]]));
    bmin <- c(0.01,AttendanceQ[1],WinpropQ[1],0.01*sdwin,0.01*sdbal,1e-5);
    bmax <- c(1.5,AttendanceQ[3],WinpropQ[3],sdwin,sdbal,sqrt(2/pi));
    df <- length(data1$NVoters)-length(bmin);

    startlist <- GetStarts(cores,bmin,bmax,data1,mdat, itprint=itprint, itstartmax=itstartmax);
    bstart <- startlist$bstart;
    h0 <- startlist$fstartlist$h0;
    hI <- startlist$fstartlist$hI;
    hE <- startlist$fstartlist$hE;
    fi <- startlist$fstartlist$fi;
    fe <- startlist$fstartlist$fe;
    b1 <- bstart;
    b0 <- 0*b1;
    val1 <- startlist$fit;
    val0 <- (-1e+31);
    wrongway <- 0;

    etype <- 0;
    for (it in 1:iterations) {
      if (itprint>0) {
        cat("iteration ",it," starts:\n",fi,fe,b1,"\n");
        cat("bmax: ",bmax,"\n");
      }
      fife0 <- c(fi,fe);
      if ((val0-val1)/abs(val0+1e-16) > 1e-4) {
        wrongway <- wrongway + 1;
        if (itprint>0) {
          cat("wrong way: ",val0," > ",val1,"\n");
          #        if (fe < ztol & b1[1]>2) {
          #          print("setting fi to zero");
          #          fi <- 0;
          #	}
          #        wrongway <- 0;
        }
      } else {
        wrongway <- 0;
      }
      if (fi > ztol) {
        par0 <- c(fife0,b0);
        par1 <- c(fi,fe,b1);
      } else if (fe > ztol) {  # drop convergence test theta because unindentifiable with fi=0
        par0 <- c(fife0,b0[-6]);
        par1 <- c(fi,fe,b1[-6]);
      } else {  # drop conv. test theta and alpha because unindentifiable with fi=fe=0
        par0 <- c(fife0,b0[-c(1,6)]);
        par1 <- c(fi,fe,b1[-c(1,6)]);
      }
      if (itprint>0) cat("max par change: ",
                         max((abs(par0-par1)/(abs(par0)+1e-16))[abs(par1)>ztol]),"\n");
      if (it>1 & (all(abs(par0-par1)/(abs(par0)+1e-16)<=tol | abs(par1)<ztol))) {
        break;
      }
      fife0 <- c(fi,fe);

      if (fi > ztol & fe > ztol) {
        etype <- 0;
        fSV <- function(b) {
          fitlist <- EvalFit(b,mdat);
          fit <- sum(log(h0*fitlist$fit0 + hI*fitlist$fitI + hE*fitlist$fitE));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      } else if (fi > ztol & fe <= ztol) {
        etype <- 1;
        fSV <- function(b) {
          fitlist <- EvalFitI(b,mdat);
          fit <- sum(log(h0*fitlist$fit0 + hI*fitlist$fitI));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      } else if (fi <= ztol & fe > ztol) {
        etype <- 2;
        fSV <- function(b) {
          fitlist <- EvalFitE(b,mdat);
          fit <- sum(log(h0*fitlist$fit0 + hE*fitlist$fitE));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      } else {
        etype <- 3;
        fSV <- function(b) {
          fitlist <- EvalFit0(b,mdat);
          fit <- sum(log(fitlist$fit0));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      }

      #Fraud Model
      bstart <- b1;
      clarg <- cluster;
      if (length(cluster)==1) clarg <- FALSE;
      if (etype==0 | etype==1) {
        oSV <- genoud(fSV, 6, max=TRUE, Domains=cbind(bmin,bmax),
                      pop.size=pop,
                      starting.values=bstart, boundary.enforcement=2,
                      BFGS=FALSE, gradient.check=FALSE, print.level=0, cluster=clarg);
      } else if (etype==2) {
        oSV <- genoud(fSV, 5, max=TRUE, Domains=cbind(bmin,bmax)[-6,],
                      pop.size=pop,
                      starting.values=bstart[-6], boundary.enforcement=2,
                      BFGS=FALSE, gradient.check=FALSE, print.level=0, cluster=clarg);
      } else {
        oSV <- genoud(fSV, 4, max=TRUE, Domains=cbind(bmin,bmax)[-c(1,6),],
                      pop.size=pop,
                      starting.values=bstart[-c(1,6)], boundary.enforcement=2,
                      BFGS=FALSE, gradient.check=FALSE, print.level=0, cluster=clarg);
      }

      if (itprint>=2) {
        print(oSV$value, digits=12);
        print(oSV$par);
      }

      b0 <- b1;
      if (etype==0) {
        b1 <- oSV$par;
        ##Compute Statistics
        fitlist <- EvalFit(b1,mdat);
        fit <- (1-fi-fe)*fitlist$fit0 + fi*fitlist$fitI + fe*fitlist$fitE;
        h0 <- (1-fi-fe)*fitlist$fit0/fit;
        hI <- fi*fitlist$fitI/fit;
        hE <- fe*fitlist$fitE/fit;
        fi <- mean(hI);
        fe <- mean(hE);

        #      alpha <- oSV$par[1];
        #      p_Att <- oSV$par[2];
        #      lambda_fraud <- c(oSV$par[3],1-oSV$par[3]);
        #      sigma <- oSV$par[4];
        #      stdAtt <- oSV$par[5];
        #      theta <- oSV$par[6];
      } else if (etype==1) {
        b1 <- oSV$par;
        ##Compute Statistics
        fitlist <- EvalFitI(b1,mdat);
        fit <- (1-fi)*fitlist$fit0 + fi*fitlist$fitI;
        h0 <- (1-fi)*fitlist$fit0/fit;
        hI <- fi*fitlist$fitI/fit;
        fi <- mean(hI);
        fe <- 0;
      } else if (etype==2) {
        b1 <- c(oSV$par,0);
        ##Compute Statistics
        fitlist <- EvalFitE(b1[-6],mdat);
        fit <- (1-fe)*fitlist$fit0 + fe*fitlist$fitE;
        h0 <- (1-fe)*fitlist$fit0/fit;
        hE <- fe*fitlist$fitE/fit;
        fi <- 0;
        fe <- mean(hE);
      } else {
        b1 <- c(0,oSV$par,0);
        ##Compute Statistics
        fitlist <- EvalFit0(b1[-c(1,6)],mdat);
        fit <- fitlist$fit0;
        fi <- 0;
        fe <- 0;
      }
      bmax[1] <- ifelse (oSV$par[1]>bmax[1]*.95,oSV$par[1]*1.5,bmax[1]);
      val1 <- oSV$value;
      val0 <- max(val0,val1);

      if (itprint>=1) {
        print(paste("Iteration",it));
        print(oSV$value,digits=12);
        print(oSV$par);
      }
      if (itprint>0) cat("fi, fe, value: ",fi,fe,val1,"\n");
      Ffits <- c(fi,fe,b1,oSV$value,df);

      FraudFits <- cbind(FraudFits,Ffits);
    }

    return(list(FraudFits=FraudFits,h0=h0,hI=hI,hE=hE))
  }
  #
  #
  #
  #
  #

  ###Function 4
  ###Iteration Function
  Iteration <- function(iterations,data1,cores,cluster,
                        pop=5000,itprint=0,itstartmax=10, pstarts=NULL) {

    # absolute value of Normal density
    danorm <- function(x, m, s) {
      s2 <- s^2;
      1/(s*sqrt(2*pi))* (exp(-(-x-m)^2/(2*s2)) + exp(-(x-m)^2/(2*s2)));
    }

    # absolute value of Normal density with zero mean (half-Normal)
    danorm0 <- function(x, s) {
      2/(s*sqrt(2*pi))* exp(-x^2/(2*s^2));
    }

    erf <- function(x) {
      2*pnorm(sqrt(2)*x) - 1
    }

    # absolute value of Normal cdf with zero mean (half-Normal)
    Fanorm0 <- function(x, s) {
      #    erf(x/sqrt(2*s));
      2*pnorm(sqrt(2)*(x/sqrt(2*s)))-1;
    }

    # Normal by half-normal density
    # for x >= 0
    # Epstein (1948) equation (2);  Curtiss (1941)
    dconnHNI <- function(x, NV, m1, s1, m2, s2){
      F0 <- pnorm(1,mean=m2,sd=s2)-pnorm(0,mean=m2,sd=s2);
      F1 <- Fanorm0(1, s1);
      f <- function(y, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2) {
        ifelse(y==1,0,dnorm(x/(1-y) ,mean=NV*(1-m2),sd=NV*s2) * danorm(y,m1,s1)/(1-y)/(F0*F1));
      }
      intf <- integrate(f, 0, 1, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    # Normal by 1-half-normal density
    # for x >= 0
    # Epstein (1948) equation (2);  Curtiss (1941)
    dconnHNE <- function(x, NV, m1, s1, m2, s2){
      F0 <- pnorm(1,mean=m2,sd=s2)-pnorm(0,mean=m2,sd=s2);
      F1 <- Fanorm0(1, s1);
      f <- function(y, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2) {
        ifelse(y==0,0,dnorm(x/y ,mean=NV*(1-m2),sd=NV*s2) * danorm(y,m1,s1)/y/(F0*F1));
      }
      intf <- integrate(f, 0, 1, x = x, NV=NV, m1 = m1, s1 = s1, m2 = m2, s2 = s2,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    # no "fraud" joint density
    dconv <- function(xW, xA, NV, m21, s21, m22, s22) {
      F0 <- pnorm(1,mean=m21,sd=s21)-pnorm(0,mean=m21,sd=s21);
      F2 <- pnorm(1,mean=m22,sd=s22)-pnorm(0,mean=m22,sd=s22);
      F02 <- F0*F2^2;
      a <- 1-xA/NV;
      v <-
        dnorm(xW/a, NV*m21, NV*s21) *
        dnorm(a, m22, s22) *
        dnorm(xA, NV*(1-m22), NV*s22)/a/F02;
      ifelse(is.finite(v),v,0);
    }

    # Inc "fraud" joint density
    dconvI <- function(xW, NV, m1, s1, m21, s21, m22, s22, e) {
      F0 <- pnorm(1,mean=m21,sd=s21)-pnorm(0,mean=m21,sd=s21);
      F1 <- Fanorm0(1, s1);
      F2 <- pnorm(1,mean=m22,sd=s22)-pnorm(0,mean=m22,sd=s22);
      F02 <- F0*F2;
      g <- function(z,y, xW, NV, m21, s21, m22, s22, ye, ye1, dy){
        # X = L*A*(1 - U^q) + U*(1-A) + U^q*A
        zye <- z*ye;
        v <- ifelse(z==0|y==1,0,
                    dnorm(xW/zye, NV*(m21 + y*(1-z)/zye + ye1), NV*s21) *
                      dy * dnorm(z, m22, s22)/zye/F02);
        ifelse(is.finite(v),v,0);
      }
      f <- function(y,xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e){
        ye <- 1-y^e;
        ye1 <- y^(e)/(1-y^(e));
        dy <- (1-y)^(-1+1/e)*danorm0(y, s1)/e/F1;
        intg <- integrate(g, 0, 1,
                          y=y, xW=xW,NV=NV,m21=m21,s21=s21,m22=m22,s22=s22,ye=ye,ye1=ye1,dy=dy,
                          stop.on.error=FALSE);
        ret <- ifelse(intg$message=="OK",intg$value, 0);
        if (length(ret)==1) ret <- rep(ret,length(y));
        ret
      }
      intf <- integrate(f, 0, 1,
                        xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    # Ext "fraud" joint density
    dconvE <- function(xW, NV, m1, s1, m21, s21, m22, s22, e) {
      F0 <- pnorm(1,mean=m21,sd=s21)-pnorm(0,mean=m21,sd=s21);
      F1 <- Fanorm0(1, s1);
      F2 <- pnorm(1,mean=m22,sd=s22)-pnorm(0,mean=m22,sd=s22);
      F02 <- F0*F2;
      g <- function(z,y, xW, NV, m21, s21, m22, s22, ye, ye1, dy){
        # X = L*A*(1 - U^q) + U*(1-A) + U^q*A
        zye <- z*ye;
        v <- ifelse(z==0|y==0,0,
                    dnorm(xW/zye, NV*(m21 + (1-y)*(1-z)/zye + ye1), NV*s21) *
                      dy * dnorm(z, m22, s22)/zye/F02);
        ifelse(is.finite(v),v,0);
      }
      f <- function(y,xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e){
        ye <- 1-(1-y)^e;
        ye1 <- (1-y)^(e)/(1-(1-y)^(e));
        dy <- (1-y)^(-1+1/e)*danorm0(y, s1)/e/F1;
        intg <- integrate(g, 0, 1,
                          y=y, xW=xW,NV=NV,m21=m21,s21=s21,m22=m22,s22=s22,ye=ye,ye1=ye1,dy=dy,
                          stop.on.error=FALSE);
        ret <- ifelse(intg$message=="OK",intg$value, 0);
        if (length(ret)==1) ret <- rep(ret,length(y));
        ret
      }
      intf <- integrate(f, 0, 1,
                        xW=xW,NV=NV,m1=m1,s1=s1,m21=m21,s21=s21,m22=m22,s22=s22,e=e,
                        stop.on.error=FALSE);
      ifelse(intf$message=="OK",intf$value, 0);
    }

    EvalFit <- function(par,mdat) {
      alpha        <- par[1];
      p_Att        <- par[2];
      lambda_fraud <- par[3];
      sigma        <- par[4];
      stdAtt       <- par[5];
      theta        <- par[6];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      fitI <-
        sapply(loc,function(i) {
          dconvI(mdat[i,1], NV[i], 0, theta, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNI(mdat[i,3], NV[i], 0, theta, p_Att, stdAtt)
        });

      fitE <-
        sapply(loc,function(i) {
          dconvE(mdat[i,1], NV[i], 0, 0.075, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNE(mdat[i,3], NV[i], 0, 0.075, p_Att, stdAtt)
        });

      return(list(fit0=fit0,fitI=fitI,fitE=fitE));
    }

    EvalFitI <- function(par,mdat) {
      alpha        <- par[1];
      p_Att        <- par[2];
      lambda_fraud <- par[3];
      sigma        <- par[4];
      stdAtt       <- par[5];
      theta        <- par[6];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      fitI <-
        sapply(loc,function(i) {
          dconvI(mdat[i,1], NV[i], 0, theta, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNI(mdat[i,3], NV[i], 0, theta, p_Att, stdAtt)
        });

      return(list(fit0=fit0,fitI=fitI));
    }

    EvalFitE <- function(par,mdat) {
      alpha        <- par[1];
      p_Att        <- par[2];
      lambda_fraud <- par[3];
      sigma        <- par[4];
      stdAtt       <- par[5];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      fitE <-
        sapply(loc,function(i) {
          dconvE(mdat[i,1], NV[i], 0, 0.075, lambda_fraud, sigma, p_Att, stdAtt,
                 10*alpha)*
            dconnHNE(mdat[i,3], NV[i], 0, 0.075, p_Att, stdAtt)
        });

      return(list(fit0=fit0,fitE=fitE));
    }

    EvalFit0 <- function(par,mdat) {
      p_Att        <- par[1];
      lambda_fraud <- par[2];
      sigma        <- par[3];
      stdAtt       <- par[4];

      loc <- 1:dim(mdat)[1];
      NV <- apply(mdat,1,sum);

      fit0 <-
        sapply(loc,function(i) {
          dconv(mdat[i,1], mdat[i,3], NV[i], lambda_fraud, sigma, p_Att, stdAtt)
        });

      return(list(fit0=fit0));
    }

    GetfStarts <- function(fitlist, itprint) {
      while (TRUE) {
        fi <- runif(1,.1,.95);
        fe <- runif(1,.1,.3);
        if (fi+fe<1) break;
      }
      if (itprint>1) cat("0starts:\n",fi,fe,"\n");
      fit <- (1-fi-fe)*fitlist$fit0 + fi*fitlist$fitI + fe*fitlist$fitE;
      h0 <- (1-fi-fe)*fitlist$fit0/fit;
      hI <- fi*fitlist$fitI/fit;
      hE <- fe*fitlist$fitE/fit;
      fi <- mean(hI);
      fe <- mean(hE);
      if (itprint>1) cat("1starts:\n",fi,fe,"\n");
      return(list(fi=fi,fe=fe,h0=h0,hI=hI,hE=hE));
    }

    GetStartsGrid <- function(cores,pstarts,bmin,bmax,data1,mdat, itprint) {
      fit <- fitmax <- (-1e+50);
      bstart0 <- pstarts;
      fstartlist0 <- list();
      OKlist <- list();
      itmax <- 10;
      agrid <- seq(bmin[1],bmax[1],length.out=itmax);
      tgrid <- seq(bmin[6],bmax[6],length.out=itmax);
      p_Att <- pstarts[2]
      lambda_fraud <- pstarts[3];
      lambda_fraud <- c(lambda_fraud, 1-lambda_fraud);
      sigma <- pstarts[4];
      stdAtt <- pstarts[5];
      # prepare to use foreach()
      registerDoParallel(cores=min(itmax,cores));
      for (i in 1:itmax) {
        alpha <- agrid[i];
        sublist <- foreach(j=1:itmax, .export=c('GetfStarts', 'EvalFit')) %dopar% {
          theta <- tgrid[j];
          bstart <- c(alpha,p_Att,lambda_fraud[1],sigma,stdAtt,theta);
          fstartlist <- GetfStarts(EvalFit(bstart,mdat), itprint=itprint);
          fi <- fstartlist$fi;
          fe <- fstartlist$fe;
          slist <- NULL;
          if (is.finite(fi) & is.finite(fe)) slist <- list(i,j,fstartlist);
          slist;
        }
        OKlist[paste(i, 1:itmax)] <- sublist;
      }
      stopImplicitCluster();
      OKfits <- list();
      for (j in names(OKlist)) {
        if (is.null(OKlist[[j]])) next;
        bstart <- pstarts;
        bstart[c(1,6)] <- c(agrid[OKlist[[j]][[1]]],tgrid[OKlist[[j]][[2]]]);
        h0 <- OKlist[[j]][[3]]$h0;
        hI <- OKlist[[j]][[3]]$hI;
        hE <- OKlist[[j]][[3]]$hE;
        fitlist <- EvalFit(bstart,mdat);
        fit <- sum(log(h0*fitlist$fit0 + hI*fitlist$fitI + hE*fitlist$fitE));
        if (is.finite(fit)) {
          OKfits[[j]] <- list(bstart=bstart, fit=fit, fstartlist=OKlist[[j]][[3]]);
        }
      }
      if (length(OKfits)>0) {
        for (j in names(OKfits)) {
          if (OKfits[[j]]$fit > fitmax) {
            fitmax <- OKfits[[j]]$fit;
            bstart0 <- OKfits[[j]]$bstart;
            fstartlist0 <- OKfits[[j]]$fstartlist;
          }
        }
      }
      if (itprint>0) cat("grid fit start:",fitmax,"\n");
      return(list(bstart=bstart0, fit=fitmax, fstartlist=fstartlist0));
    }

    GetStartsIteration <- function(itstartmax,bmin,bmax,data1,mdat, itprint) {
      fit <- fitmax <- (-1e+50);
      bstart0 <- (bmin+bmax)/2;
      fstartlist0 <- list();
      itmax <- 10^2;
      for (itstart in 1:itstartmax) {
        for (i in 1:itmax) {
          for (j in 1:itmax) {
            alpha <- runif(1,bmin[1],bmax[1]);
            p_Att <- runif(1,bmin[2],bmax[2]);
            lambda_fraud <- runif(1,bmin[3],bmax[3]);
            lambda_fraud <- c(lambda_fraud, 1-lambda_fraud);
            sigma <- runif(1,bmin[4],bmax[4]);
            stdAtt <- runif(1,bmin[5],bmax[5]);
            theta <- runif(1,bmin[6],bmax[6]);
            bstart <- c(alpha,p_Att,lambda_fraud[1],sigma,stdAtt,theta);

            fstartlist <- GetfStarts(EvalFit(bstart,mdat), itprint=itprint);
            h0 <- fstartlist$h0;
            hI <- fstartlist$hI;
            hE <- fstartlist$hE;
            fi <- fstartlist$fi;
            fe <- fstartlist$fe;
            if (is.finite(fi) & is.finite(fe)) break;
          }
          fitlist <- EvalFit(bstart,mdat);
          fit <- sum(log(h0*fitlist$fit0 + hI*fitlist$fitI + hE*fitlist$fitE));
          if (is.finite(fit)) { break; }
        }
        if (is.finite(fit) && is.finite(fitmax) && fit > fitmax) {
          if (itprint>0) print(c(itstart,fit));
          fitmax <- fit;
          bstart0 <- bstart;
          fstartlist0 <- fstartlist;
        }
      }
      return(list(bstart=bstart0, fit=fitmax, fstartlist=fstartlist0));
    }
    #cores,bmin,bmax,data1,mdat, itprint=itprint, itstartmax=itstartmax
    GetStarts <- function(cores,bmin,bmax,data1,mdat, itprint, itstartmax=10) {
      Icores <- cores;  # number of processors to use (processes to spawn)
      if (itstartmax > Icores) {
        Ibase <- trunc(itstartmax/Icores);
        Ires  <- itstartmax %% Icores;
        Ivec  <- rep(Ibase,Icores);
        Ivec[1:Ires] <- Ivec[1:Ires] + 1;
      } else {
        Ivec <- rep(1,itstartmax);
      }
      while (TRUE) {
        # prepare to use foreach()
        registerDoParallel(cores=cores);
        StartsList <- foreach(itstartMAX=Ivec) %dopar%
          GetStartsIteration(itstartMAX,bmin,bmax,data1,mdat, itprint);
        stopImplicitCluster();
        if (itprint>0) { print("starts report"); }
        jmax <- 1;
        for (j in 1:itstartmax) {
          if (itprint>0) { print(c(j,StartsList[[j]]$fit)); }
          if (StartsList[[j]]$fit > StartsList[[jmax]]$fit) jmax <- j;
        }
        fitmax <- StartsList[[jmax]]$fit;
        if (fitmax > (-1e+50)) break;
      }
      bstart0 <- StartsList[[jmax]]$bstart;
      fstartlist0 <- StartsList[[jmax]]$fstartlist;
      cat("starting fit:  ",fitmax,"\n");
      return(list(bstart=bstart0, fit=fitmax, fstartlist=fstartlist0));
    }

    ###Set FraudFits
    FraudFits <- matrix(NA,nrow=8+2,ncol=0);

    tol <- 1e-3;
    ztol <- 1e-9;
    if (itprint>0) cat("zero tolerance: ",ztol,"\n");
    mdat <- cbind(data1$Votes,data1$NValid-data1$Votes,data1$NVoters-data1$NValid);
    AttendanceQ <- quantile(Attendance <- data1$NValid/data1$NVoters);
    WinpropQ <- quantile(Winprop <- data1$Votes/data1$NValid);
    sdwin <- 2*sqrt(var(Winprop[Winprop<=WinpropQ[3]]));
    sdbal <- 2*sqrt(var(Attendance[Attendance<=AttendanceQ[3]]));
    bmin <- c(0.01,AttendanceQ[1],WinpropQ[1],0.01*sdwin,0.01*sdbal,1e-5);
    bmax <- c(1.5,AttendanceQ[3],WinpropQ[3],sdwin,sdbal,sqrt(2/pi));
    df <- length(data1$NVoters)-length(bmin);

    if (!is.null(pstarts)) {
      startlist <- GetStartsGrid(cores,pstarts,bmin,bmax,data1,mdat, itprint=itprint);
      if (startlist$fit > (-1e+50)) {
        cat("starting fit:  ",startlist$fit,"\n");
      } else {
        pstarts <- NULL;
      }
    }
    if (is.null(pstarts)) {
      startlist <- GetStarts(cores,bmin,bmax,data1,mdat, itprint=itprint, itstartmax=itstartmax);
    }
    bstart <- startlist$bstart;
    h0 <- startlist$fstartlist$h0;
    hI <- startlist$fstartlist$hI;
    hE <- startlist$fstartlist$hE;
    fi <- startlist$fstartlist$fi;
    fe <- startlist$fstartlist$fe;
    b1 <- bstart;
    b0 <- 0*b1;
    val1 <- startlist$fit;
    val0 <- (-1e+31);
    wrongway <- 0;

    etype <- 0;
    for (it in 1:iterations) {
      if (itprint>0) {
        cat("iteration ",it," starts:\n",fi,fe,b1,"\n");
        cat("bmax: ",bmax,"\n");
      }
      fife0 <- c(fi,fe);
      if ((val0-val1)/abs(val0+1e-16) > 1e-4) {
        wrongway <- wrongway + 1;
        if (itprint>0) {
          cat("wrong way: ",val0," > ",val1,"\n");
          #        if (fe < ztol & b1[1]>2) {
          #          print("setting fi to zero");
          #          fi <- 0;
          #	}
          #        wrongway <- 0;
        }
      } else {
        wrongway <- 0;
      }
      if (fi > ztol) {
        par0 <- c(fife0,b0);
        par1 <- c(fi,fe,b1);
      } else if (fe > ztol) {  # drop convergence test theta because unindentifiable with fi=0
        par0 <- c(fife0,b0[-6]);
        par1 <- c(fi,fe,b1[-6]);
      } else {  # drop conv. test theta and alpha because unindentifiable with fi=fe=0
        par0 <- c(fife0,b0[-c(1,6)]);
        par1 <- c(fi,fe,b1[-c(1,6)]);
      }
      if (itprint>0) cat("max par change: ",
                         max((abs(par0-par1)/(abs(par0)+1e-16))[abs(par1)>ztol]),"\n");
      if (it>1 & (all(abs(par0-par1)/(abs(par0)+1e-16)<=tol | abs(par1)<ztol))) {
        break;
      }
      fife0 <- c(fi,fe);

      if (fi > ztol & fe > ztol) {
        etype <- 0;
        fSV <- function(b) {
          fitlist <- EvalFit(b,mdat);
          fit <- sum(log(h0*fitlist$fit0 + hI*fitlist$fitI + hE*fitlist$fitE));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      } else if (fi > ztol & fe <= ztol) {
        etype <- 1;
        fSV <- function(b) {
          fitlist <- EvalFitI(b,mdat);
          fit <- sum(log(h0*fitlist$fit0 + hI*fitlist$fitI));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      } else if (fi <= ztol & fe > ztol) {
        etype <- 2;
        fSV <- function(b) {
          fitlist <- EvalFitE(b,mdat);
          fit <- sum(log(h0*fitlist$fit0 + hE*fitlist$fitE));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      } else {
        etype <- 3;
        fSV <- function(b) {
          fitlist <- EvalFit0(b,mdat);
          fit <- sum(log(fitlist$fit0));
          if (is.na(fit)) fit <- (-1e+10);
          return(fit);
        }
      }

      #Fraud Model
      bstart <- b1;
      clarg <- cluster;
      if (length(cluster)==1) clarg <- FALSE;
      if (etype==0 | etype==1) {
        oSV <- genoud(fSV, 6, max=TRUE, Domains=cbind(bmin,bmax),
                      pop.size=pop,
                      starting.values=bstart, boundary.enforcement=2,
                      BFGS=FALSE, gradient.check=FALSE, print.level=0, cluster=clarg);
      } else if (etype==2) {
        oSV <- genoud(fSV, 5, max=TRUE, Domains=cbind(bmin,bmax)[-6,],
                      pop.size=pop,
                      starting.values=bstart[-6], boundary.enforcement=2,
                      BFGS=FALSE, gradient.check=FALSE, print.level=0, cluster=clarg);
      } else {
        oSV <- genoud(fSV, 4, max=TRUE, Domains=cbind(bmin,bmax)[-c(1,6),],
                      pop.size=pop,
                      starting.values=bstart[-c(1,6)], boundary.enforcement=2,
                      BFGS=FALSE, gradient.check=FALSE, print.level=0, cluster=clarg);
      }

      if (itprint>=2) {
        print(oSV$value, digits=12);
        print(oSV$par);
      }

      b0 <- b1;
      if (etype==0) {
        b1 <- oSV$par;
        ##Compute Statistics
        fitlist <- EvalFit(b1,mdat);
        fit <- (1-fi-fe)*fitlist$fit0 + fi*fitlist$fitI + fe*fitlist$fitE;
        h0 <- (1-fi-fe)*fitlist$fit0/fit;
        hI <- fi*fitlist$fitI/fit;
        hE <- fe*fitlist$fitE/fit;
        fi <- mean(hI);
        fe <- mean(hE);

        #      alpha <- oSV$par[1];
        #      p_Att <- oSV$par[2];
        #      lambda_fraud <- c(oSV$par[3],1-oSV$par[3]);
        #      sigma <- oSV$par[4];
        #      stdAtt <- oSV$par[5];
        #      theta <- oSV$par[6];
      } else if (etype==1) {
        b1 <- oSV$par;
        ##Compute Statistics
        fitlist <- EvalFitI(b1,mdat);
        fit <- (1-fi)*fitlist$fit0 + fi*fitlist$fitI;
        h0 <- (1-fi)*fitlist$fit0/fit;
        hI <- fi*fitlist$fitI/fit;
        fi <- mean(hI);
        fe <- 0;
      } else if (etype==2) {
        b1 <- c(oSV$par,0);
        ##Compute Statistics
        fitlist <- EvalFitE(b1[-6],mdat);
        fit <- (1-fe)*fitlist$fit0 + fe*fitlist$fitE;
        h0 <- (1-fe)*fitlist$fit0/fit;
        hE <- fe*fitlist$fitE/fit;
        fi <- 0;
        fe <- mean(hE);
      } else {
        b1 <- c(0,oSV$par,0);
        ##Compute Statistics
        fitlist <- EvalFit0(b1[-c(1,6)],mdat);
        fit <- fitlist$fit0;
        fi <- 0;
        fe <- 0;
      }
      bmax[1] <- ifelse (oSV$par[1]>bmax[1]*.95,oSV$par[1]*1.5,bmax[1]);
      val1 <- oSV$value;
      val0 <- max(val0,val1);

      if (itprint>=1) {
        print(paste("Iteration",it));
        print(oSV$value,digits=12);
        print(oSV$par);
      }
      if (itprint>0) cat("fi, fe, value: ",fi,fe,val1,"\n");
      Ffits <- c(fi,fe,b1,oSV$value,df);

      FraudFits <- cbind(FraudFits,Ffits);
    }

    return(list(FraudFits=FraudFits,h0=h0,hI=hI,hE=hE))
  }

  ElectionFitter <- function(data1, rations=2, itprint=0, pop=5000, cores=1, itstartmax=1) {
    ##Initiitealize some parameters

    ####################
    cluster <- rep("localhost",cores);
    FFlist <- IterationN(iterations,data1,cores,cluster,
                         pop=pop,itprint=itprint,itstartmax=itstartmax);
    FraudFits <- FFlist$FraudFits;

    FFnames <- c("incremental","extreme","alpha","turnout","winprop","sigma","stdAtt","theta",
                 "loglik","df");
    if (length(dim(FraudFits))>1) {
      FFuse <- 1:(dim(FraudFits)[2]);
      if (length(FFuse)>2) FFuse <- FFuse[(dim(FraudFits)[2])-0:1];
      FF <- cbind(FraudFits[,dim(FraudFits)[2]],
                  sqrt(apply(FraudFits[,FFuse,drop=FALSE],1,var)));
      dimnames(FF) <- list(FFnames,c("est","sdev"));
    } else  {
      FF <- FraudFits;
      names(FF) <- FFnames;
    }

    results<-list()
    results$FF_null <- FF;
    results$FFlist_null <- FFlist;

    print("null model results:");
    print(results$FF_null);

    pstarts <- results$FF_null[3:8,1];

    FFlist <- Iteration(iterations,data1,cores,cluster,
                        pop=pop,itprint=itprint,itstartmax=itstartmax, pstarts=pstarts);
    FraudFits <- FFlist$FraudFits;

    FFnames <- c("incremental","extreme","alpha","turnout","winprop","sigma","stdAtt","theta",
                 "loglik","df");
    if (length(dim(FraudFits))>1) {
      FFuse <- 1:(dim(FraudFits)[2]);
      if (length(FFuse)>2) FFuse <- FFuse[(dim(FraudFits)[2])-0:1];
      FF <- cbind(FraudFits[,dim(FraudFits)[2]],
                  sqrt(apply(FraudFits[,FFuse,drop=FALSE],1,var)));
      dimnames(FF) <- list(FFnames,c("est","sdev"));
    }  else  {
      FF <- FraudFits;
      names(FF) <- FFnames;
    }

    results$FF <- FF;
    results$FFlist <- FFlist;
    return(results)
  }

  datC <- cleandata(dat)

  resultFMM <- ElectionFitter(datC, cores, itstartmax)
return(resultFMM)}
