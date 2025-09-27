#' @title Basic Election Forensics function
#' @description This function helps to run basic election forensics methods: Benford's 2nd Digit, Last Digit (Counts),
#' Counts (05s), Multimodality Test, Sobyanin-Sukhovolsky test and Correlation test.
#' @usage BasicElectionForensics(data, Candidates, Level,
#'                              TotalReg, TotalVotes, Methods, R=1000, cores=2)
#' @param data data frame
#' @param Candidates variable name referring to vote counts for candidates/parties
#' @param Level variable name depicting the level of analysis ("National", i.e. whole dataset by default)
#' @param TotalReg variable name for the total number of eligible voters
#' @param TotalVotes variable name for the total number of ballots cast
#' @param Methods
#' \itemize{
#'   \item "_2BL" - second-digit mean test (to be compared to the mean value Benford's Law)
#'   \item "LastC" - last-digit mean test (to be compared to the mean value implied by uniformly distributed last digits)
#'   \item "P05s" - percentage last-digit 0/5 indicator mean (to be compared to the mean value implied by uniformly distributed percentages)
#'   \item "C05s" - count last-digit 0/5 indicator mean count (last-digit 0/5 indicator mean)
#'   \item "Sobyanin" - Sobynin-Sukhovolsky measure
#'   \item "Skew" - skewness (the extent to which a variable departs from a normal distribution by being asymmetric)
#'   \item "Kurt" - Kurt (the extent to which a variable departs from a normal distribution by being spread out too much or not enough)
#'   \item "DipT" - Unimodality test (tests whether the distribution of a variable departs from unimodality)
#'   \item "Correlation" - correlation coefficient between turnout and vote share
#' }
#' @param R number of simulations (1000  simulations by default)
#' @param cores number of cores for parallel computing (2 cores by default)
#' @export
#' @importFrom xtable xtable
#' @import boot
#' @import stats
#' @importFrom moments skewness kurtosis
#' @import diptest
#' @import doParallel
#' @import hwriter
#' @import kableExtra
#' @return A list with the following parameters:
#' \itemize{
#'   \item table - dataframe with results
#'   \item html -  html table with results
#'   \item tex - tex table with results
#'   \item sigMatrix - significance matrix
#' }
#'
#' @examples
#' library(EFToolkit)
#'
#' dat<-read.csv(system.file("extdata/Albania2013.csv", package="EFToolkit"))
#'
#' #NB! R=100 to speed up computations for this example.
#' eldata<-BasicElectionForensics(dat,
#'                           Candidates=c("C035", "C050"),
#'                           Level="Prefectures", TotalReg="Registered",
#'                           TotalVotes="Ballots",
#'                           Methods=c("P05s", "C05s", "_2BL", "Sobyanin",
#'                          "DipT", "Skew", "Kurt", "Correlation"), cores=2, R=100)
#' eldata
#'


############################################################
##               Election Forensics Toolkit               ##
##  25sep2015, 24oct2019, 23sep2025                       ##
##  Kirill Kalinin and Walter R. Mebane, Jr               ##
############################################################

BasicElectionForensics<-function(data, Candidates, Level="National", TotalReg, TotalVotes, Methods, R=1000, cores=2){

  closeAllConnections()

  MethNames<-c("_2BL", "LastC", "P05s", "C05s", "Sobyanin", "Skew", "Kurt", "DipT", "Correlation", "Obs");
  MethNames2Boot <- c("_2BL", "LastC", "P05s", "C05s", "Skew", "Kurt")
  nullvar="-"

  #Functions
  SecondDigitGetter <- function(num,digit) {
    s <- as.character(num);
    idx <- sapply(s, function(x){ nchar(x) >= digit; });
    dout <- ifelse(idx, NA, "");
    if (any(idx)) {
      dout[idx] <- sapply(s[idx], function(x){ substr(x,digit,digit) })}
    return(dout);
  }

  LastDigitGetter <- function(num) {
    s <- as.character(num);
    dout <-  as.numeric(substr(num, nchar(num),nchar(num)));
    return(dout)
  }

  sobyanin_sukhovolsky<-function(turnout, voteshare){
    est_reg<-lm(voteshare~turnout)
    mean_voteshare<-mean(voteshare, na.rm=TRUE)
    est_avs<-(coef(est_reg)[2]-mean_voteshare) #large positive difference = fraud
    #tstat<- est_avs/coefficients(summary(est_reg))[1,2]
    return(est_avs)}

  correlation_coef<-function(turnout, voteshare){
    cor_result<-cor(voteshare, turnout, use="pairwise.complete.obs")
    return(cor_result)}

  NPAnal <- function(methods,pct=NULL,vot=NULL, mnames=MethNames) {
    outvec <- rep(NA,length(methods));
    names(outvec) <- methods;
    if (!is.null(pct)) {
      pct100 <- pct*100;
      for (j in methods) {
        if (j == "_2BL") outvec["_2BL"] <-
            mean(as.numeric(SecondDigitGetter(vot,2)), na.rm=TRUE);
        if (j == "LastC") outvec["LastC"] <-
            mean(LastDigitGetter(vot), na.rm=TRUE);
        if (j == "P05s") outvec["P05s"] <-
            sum(LastDigitGetter(round(pct100,0))==5 |
                  LastDigitGetter(round(pct100,0))==0,
                na.rm=TRUE)/(length(pct100)-sum(is.na(pct100)));
        if (j == "C05s") outvec["C05s"] <-
            mean(LastDigitGetter(vot)==5|LastDigitGetter(vot)==0,
                 na.rm=TRUE);
        if (j == "Skew") {
          outvec["Skew"] <- skewness(pct,na.rm=TRUE);
        }
        if (j == "Kurt") {
          outvec["Kurt"] <- kurtosis(pct,na.rm=TRUE);
        }
      }
    } else {
      for (j in methods) {
        if (j == "_2BL") outvec["_2BL"] <-
            mean(as.numeric(SecondDigitGetter(vot,2)), na.rm=TRUE);
        if (j == "LastC") outvec["LastC"] <-
            mean(LastDigitGetter(vot), na.rm=TRUE);
        if (j == "P05s") outvec["P05s"] <- NA;
        if (j == "C05s") outvec["C05s"] <-
            mean(LastDigitGetter(vot)==5|LastDigitGetter(vot)==0,
                 na.rm=TRUE);
        if (j == "Skew")  outvec["Skew"] <- NA;
        if (j == "Kurt")  outvec["Kurt"] <- NA;
      }
    }
    outvec;
  }

  bootAnalT <- function(dat,idx, methods="_2BL") {
    dat <- dat[idx,,drop=FALSE];
    eligible<-dat$totalReg
    totalVoters<-dat$totalVotes
    turnoutPercent<-totalVoters/eligible
    turnoutPercent<-ifelse(turnoutPercent>=1,1,turnoutPercent)

    NPAnal(methods,pct=turnoutPercent,vot=totalVoters);
  }

  bootAnalV <- function(dat,idx, methods="_2BL", cand=NULL) {
    dat <- dat[idx,,drop=FALSE];
    voteCandidate<-dat[,names(dat) %in% cand]
    if (!is.null(dat$totalVotes)) {
      totalVoters<-dat$totalVotes
      percentCandidate<-voteCandidate/totalVoters
      #Controlling for over 1
      percentCandidate<-ifelse(percentCandidate>=1,1,percentCandidate)
      NPAnal(methods,pct=percentCandidate, vot=voteCandidate);
    }
    else {
      NPAnal(methods,vot=voteCandidate);
    }
  }

  statsArray2Table <- function(arr) {
    nLevel <- dim(arr)[1];
    nCands <- dim(arr)[2];
    nParms <- dim(arr)[3];
    nRow <- nLevel*nCands;
    nCol <- nParms + 3;
    Parms <- dimnames(arr)[3];
    tab <- data.frame(matrix(NA, nRow*2, nCol));
    names(tab) <- c("Level","Candidate",dimnames(arr)[[3]],"Obs");
    tab$Level[2*(1:nRow)-1] <- dimnames(arr)[[1]];
    for (j in 1:nCands) {
      tab$Candidate[(j-1)*2*nLevel + 2*(1:nLevel)-1] <- dimnames(arr)[[2]][j];
      for (k in 1:nParms) {
        tab[(j-1)*2*nLevel + 2*(1:nLevel)-1,2+k] <- round(arr[,j,k,1],3);
        for (L in 1:nLevel) {
          if (all(!is.na(arr[L,j,k,2:3]))) {
            tab[(j-1)*2*nLevel + 2*L,2+k] <-
              paste("(",paste(round(arr[L,j,k,2:3],3),collapse=", "),
                    ")",sep="",collapse="");
          }
          else if (!is.na(arr[L,j,k,2])) {
            tab[(j-1)*2*nLevel + 2*L,2+k] <-
              paste("(",paste(round(arr[L,j,k,2],3),collapse=", "),
                    ")",sep="",collapse="");
          }
        }
        if (regexpr("KSim|KLik",Parms[k])<0 & any(!is.na(arr[,j,k,4]))) {
          tab$Obs[(j-1)*2*nLevel + 2*(1:nLevel)-1] <- arr[,j,k,4];  # "Obs"
        }
      }
    }
    tab;
  }

  ###Main
  if (TotalReg != nullvar & TotalVotes != nullvar & TotalReg != TotalVotes) {
    totalReg<-data[,which(names(data)%in%TotalReg)]
    totalVotes<-data[,which(names(data)%in%TotalVotes)]
    turnoutPercent<-totalVotes/totalReg
    data$Turnout<-turnoutPercent
    Candidates<-c("Turnout",Candidates)
  }else if (TotalVotes != nullvar) {
    totalVotes <- data[,which(names(data)%in%TotalVotes)];
    totalReg <- turnoutPercent <- data$Turnout <- NULL;
  }else {
    totalReg <- totalVotes <- turnoutPercent <- data$Turnout <- NULL;}

  Candidates[grepl("All Leaders",Candidates)]<-"AllLeaders"

  #Obtain Compact dataset
  data$AllLeaders<-1
  candidates<-data[, names(data)%in%Candidates, drop=FALSE];
  num_candidates<-length(Candidates)
  if (!is.null(totalReg) & !is.null(totalVotes)) {
    dataC<-data.frame(totalReg,totalVotes, candidates)
    colnames(dataC)<-c("totalReg","totalVotes", colnames(candidates))
  }else if (!is.null(totalReg)) {
    dataC<-data.frame(totalReg, candidates)
    colnames(dataC)<-c("totalReg", colnames(candidates))
  }else if (!is.null(totalVotes)) {
    dataC<-data.frame(totalVotes, candidates)
    colnames(dataC)<-c("totalVotes", colnames(candidates))
  }else {
    dataC <- data.frame(candidates);
    colnames(dataC) <- c(colnames(candidates));
  }

  levels<-data[, names(data)%in%Level];
  sampleSize <- dim(data)[1];

  if(length(levels)==0) {
    levels<-rep(1,dim(dataC)[1]);
    data$National<-levels
  }
  dataCC<-split(dataC, levels)
  levellength<-length(names(dataCC))
  if(levellength=='1') {
    names(dataCC)<-"National"
  }
  levelnames<-names(dataCC)

  whichCandidates<-regexpr("[CP]\\d",names(data))>0;
  matrixCandidates<- split(data[,whichCandidates,drop=FALSE],levels)

  AllLeaders<-rep(NA, dim(data)[1])
  if(("All Leaders" %in% Candidates) | ("AllLeaders" %in% Candidates)) {
    resElections <- lapply(matrixCandidates, function(x)
      names(data[,whichCandidates,drop=FALSE])[which.max(colSums(x, na.rm=TRUE))])
    for(ii in 1:length(names(resElections))) {
      AllLeaders[levels==names(resElections)[ii]] <-
        data[levels==names(resElections)[ii],names(data)%in%resElections[ii]]
    }
    dataC$AllLeaders<-AllLeaders
  }

  dataCC<-split(dataC, levels)
  levellength<-length(names(dataCC))
  if(levellength=='1') {
    names(dataCC)<-"National"
  }
  levelnames<-names(dataCC)

  #Create holder for results
  resultsLength<-num_candidates*levellength*2
  resultsWidth <- length(Methods[Methods%in%MethNames]);
  resultsMatrix <- matrix(NA,resultsLength,resultsWidth);

  #Loops
  resultsArray <- array(NA,dim=c(levellength,num_candidates,resultsWidth,4));
  dimnames(resultsArray) <- list(NULL, NULL, NULL, NULL);
  dimnames(resultsArray)[[1]] <- levelnames;
  dimnames(resultsArray)[[2]] <- Candidates;
  dimnames(resultsArray)[[3]] <- MethNames[MethNames%in%Methods]
  dimnames(resultsArray)[[4]] <- c("point","lo","up","obs");
  sourceArray <- resultsArray;
  sourceArray[,,,] <- TRUE;

  varssArray <- array(NA,dim=c(levellength,num_candidates));

  varReg <- ifelse("totalReg" %in% names(dataCC[[1]]), "totalReg", "");
  varVotes <- ifelse("totalVotes" %in% names(dataCC[[1]]),"totalVotes","")


  if (any(MethNames %in% Methods)) {
    print('Toolkit is working.')
    for(candidID in 1:num_candidates){  #Loop over Candidates
      for(levelID in 1:levellength) {     #Loop over Levels
        print(paste("Analyzing ", Candidates[candidID], ".", sep=""),quote = FALSE)
        dataSubset<-dataCC[[levelID]]
        numberUnits<-dim(dataSubset)[1]
        if (dim(dataSubset)[1]==0) next;

        voteCandidateG <- dataSubset[,names(dataSubset)%in%Candidates[candidID]]
        if (!is.null(dataSubset$totalReg)) {
          eligibleG <- dataSubset$totalReg;  #var1
        }else{
          eligibleG <- NULL;
        }
        if (!is.null(dataSubset$totalVotes)) {
          totalVotersG <- dataSubset$totalVotes;  #var2
          percentCandidateG <- voteCandidateG/totalVotersG
        }else{
          totalVotersG <- percentCandidateG <- NULL;
        }
        if (!is.null(eligibleG) & !is.null(totalVotersG)) {
          #turnoutPercentG <- na.omit(totalVotersG/eligibleG);  #var4
          turnoutPercentG <- totalVotersG/eligibleG
        }else{
          turnoutPercentG <- NULL;
        }

        #Nonparametric bootstrap
        if (Candidates[candidID] == "Turnout" &&
            !is.null(dataSubset$totalReg) && !is.null(dataSubset$totalVotes)) {
          if (any(methodsidx <- Methods %in% MethNames2Boot)) {
            bootT <- boot(dataSubset, bootAnalT, R, parallel="multicore", ncpus=cores, methods=Methods);
            resultsArray[levelID,candidID,names(bootT$t0),1] <- bootT$t0;
            resultsArray[levelID,candidID,names(bootT$t0),4] <- dim(bootT$data)[1];
            for (j in 1:length(bootT$t0)) {
              if (!is.na(bootT$t0[j])){
                compute_ci<-tryCatch(boot.ci(bootT, index=j, type="basic")$basic[4:5], error = function(e) e)
                if(inherits(compute_ci,  "error")|is.null(compute_ci)) compute_ci<-c(NA,NA)
                resultsArray[levelID,candidID,names(bootT$t0)[j],2:3] <-  compute_ci
              }
            }
          }
          if ("DipT" %in% Methods) {
            if (length(turnoutPercentG)==0) next;
            resultsArray[levelID,candidID,"DipT",1] <- dip.test(turnoutPercentG)$p.value
            resultsArray[levelID,candidID,"DipT",4] <- length(turnoutPercentG);
          }
        } else {
          if (any(methodsidx <- Methods %in% MethNames2Boot)) {
            bootV <- boot(dataSubset, bootAnalV, R, parallel="multicore", ncpus=cores, cand=Candidates[candidID], methods=Methods);
            resultsArray[levelID,candidID,names(bootV$t0),1] <- bootV$t0;
            resultsArray[levelID,candidID,names(bootV$t0),4] <- dim(bootV$data)[1];
            for (j in 1:length(bootV$t0)) {
              if (!is.na(bootV$t0[j])){
                compute_ci<-tryCatch(boot.ci(bootV, index=j, type="basic")$basic[4:5], error = function(e) e)
                if(inherits(compute_ci,  "error")|is.null(compute_ci)) compute_ci<-c(NA,NA)
                resultsArray[levelID,candidID,names(bootV$t0)[j],2:3] <-  compute_ci
              }
            }
          }
          if ("DipT" %in% Methods) {
            resultsArray[levelID,candidID,"DipT",1] <- dip.test(percentCandidateG)$p.value
            resultsArray[levelID,candidID,"DipT",4] <- length(percentCandidateG);
          }

          if ("Sobyanin" %in% Methods) {
            resultsArray[levelID,candidID,"Sobyanin",1] <- sobyanin_sukhovolsky(turnoutPercentG, percentCandidateG)
            resultsArray[levelID,candidID,"Sobyanin",4] <- length(percentCandidateG);
          }


          if ("Correlation" %in% Methods) {
            resultsArray[levelID,candidID,"Correlation",1] <- correlation_coef(turnoutPercentG, percentCandidateG)
            resultsArray[levelID,candidID,"Correlation",4] <- length(percentCandidateG);
          }
        }
      }
    } #loop over candidates
  }  #with progress

  CandidatesNames<-rep(Candidates,each=levellength*2)
  LevelNames<-rep(levelnames,length(Candidates), each=2)
  SelectedMethods <- statsArray2Table(resultsArray);
  SelectedMethods[,dim(SelectedMethods)[2]][seq(2,dim(SelectedMethods)[1],2)]<-""
  SelectedMethods[,1][seq(2,length(SelectedMethods[,1]),2)]<-""
  SelectedMethods[,2][seq(2,length(SelectedMethods[,2]),2)]<-""
  SelectedMethods[is.na(SelectedMethods)] <- "--";
  colnames(SelectedMethods) <- c("Level", "Candidate's Name", MethNames[MethNames %in% Methods], "Obs");

  coloredResults<-ColorSignificance(SelectedMethods)
  coloredResults[[1]]
  SigMatrix <- coloredResults[[2]]
  texResults <- xtable(SelectedMethods)
  coloredResults <- coloredResults[[1]]
  Results<-list(SelectedMethods,coloredResults,texResults, SigMatrix)

  names(Results) <- c("table", "html", "tex", "sigMatrix")
  return(Results)
}
