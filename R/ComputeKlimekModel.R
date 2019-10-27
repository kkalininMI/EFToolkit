#' @title Compute Klimek Model function
#' @description This function helps to estimate the Klimek et al.(2012) model.
#' @usage ComputeKlimekModel(data, Candidates, Level, TotalReg, TotalVotes, R=1000, cores=2)
#' @param data data frame
#' @param Candidates variable name referring to vote counts for candidates/parties
#' @param Level variable name depicting the level of analysis (whole dataset by default)
#' @param TotalReg variable name for the total number of eligible voters
#' @param TotalVotes  variable name for the total number of ballots cast
#' @param R number of simulations (1000  simulations by default)
#' @param cores number of cores for parallel computing (2 cores by default)
#' @export
#' @import doParallel
#' @importFrom graphics hist
#' @importFrom xtable xtable
#' @return A list with the following parameters:
#' \itemize{
#'   \item table - data frame with results
#'   \item html -  html table with results
#'   \item tex - tex table with results
#'   \item sigMatrix - significance matrix
#' }
#' @examples
#' library(EFToolkit)
#'
#' dat<-read.csv(system.file("Albania2013.csv", package="EFToolkit"))
#'
#' #NB! R=100 to speed up computations for this example.
#' klimek<-ComputeKlimekModel(dat, Candidates="C050", Level="National",
#'                         TotalReg="Registered", TotalVotes="Ballots", cores=1, R=100)
#'


############################################################
##               Election Forensics Toolkit               ##
##  2014/2/5, 08mar2014, 26jul2015                        ##
##  Walter Mebane and Naoki Egami                         ##
##  Klimek Methods                                        ##
##  Pure Conversion from Matlab ElectionFritter.m to R.   ##
##  Make Iteration part as function                       ##
##  Make this as a big function                           ##
############################################################
ComputeKlimekModel<-function(data, Candidates,
                             Level="National",
                             TotalReg, TotalVotes, R=1000, cores=2){

  nullvar="-"

  ###############
  ###Function 1
  Estimate<-function(bins,thres,Attendance,v){
    ##Estimate v.
    n<-hist(v[,1],breaks=bins,plot=F)$counts
    #  x<-seq(from=0,to=1,by=0.01*(5/thres))###change the code!!
    x <- bins;
    sl<-diff(n[10:length(n)])
    v_max<-x[which(n==max(n))][1]  ###Probably I need to put [1]
    lfthres <- thres;
    while (lfthres>1) {
      lambda_fraud<-c(x[which(sl< -lfthres)[1]+9], 1-x[which(sl<  -lfthres)[1]+9])
      if (any(is.na(lambda_fraud))) {
        lfthres <- lfthres - 0.5;
      }
      else {
        break;
      }
    }
    ##Estimate a.
    n1<-hist(Attendance[is.finite(Attendance)],breaks=bins,plot=F)$counts
    #  x1<-seq(from=0,to=1,by=0.01*(5/thres))
    x1 <- bins;
    sl1<-diff(n1[15:length(n1)])
    p_Att<-x1[which(sl1<  -lfthres)[1]+14]

    ###Compute the variance
    s1<-vector();
    s2<-vector();
    s3<-vector();

    for(i in 1:dim(v)[1]){
      if (v[i,1]<lambda_fraud[1] & Attendance[i]<p_Att){
        s1 <- c(s1,Attendance[i]-p_Att);     ### I slightly change the code for R.
      }
      if(v[i,1]<lambda_fraud[1]){
        s2 <- c(s2,v[i,1]-lambda_fraud[1])
      }
      if(v[i,1]>v_max){
        s3 <- c(s3,v[i,1]-v_max);
      }
    }
    stdAtt<-sqrt(sum(s1^2)/length(s1))
    theta<-sqrt(sqrt(sum(s3^2)/length(s3)))      ##OK
    sigma<-c(sqrt(2*sum(s2^2)/length(s2)),sqrt(2*sum(s2^2)/length(s2))) ##OK
    sigma_fraud<-sigma    ##OK

    First_Estimated<-c(lambda_fraud,p_Att,stdAtt,theta,sigma_fraud)
    names(First_Estimated) <-
      c("lambda_fraud","lambda_fraud","p_Att","stdAtt","theta","sigma_fraud","sigma_fraud")
    return(First_Estimated)
  }
  ###############
  ### Function 2
  ##Simulate FraudVotes for each set of fi, fe, fz.
  Sim_Vote<-function(p_Att,stdAtt,lambda_fraud,sigma_fraud,NVoters,N,f2,f1,theta,alpha){
    FraudVotes<-matrix(0,nrow=N,ncol=2)

    ##Model the turnout a and vote rate parametre l for the fair election.
    ##a
    a <- p_Att+rnorm(N,sd=stdAtt);
    while (any(atest <- a<0 | a>1)) a[atest] <- p_Att+rnorm(sum(atest),sd=stdAtt);

    l <-
      cbind(rnorm(N,mean=lambda_fraud[1],sd=sigma_fraud),rnorm(N,mean=lambda_fraud[2],sd=sigma_fraud));
    l <- l/apply(l,1,sum);
    while (any(ltest <- apply(l<0,1,sum) | apply(l>1,1,sum))) {
      l[ltest,] <-
        cbind(rnorm(sum(ltest),mean=lambda_fraud[1],sd=sigma_fraud),
              rnorm(sum(ltest),mean=lambda_fraud[2],sd=sigma_fraud));
      l[ltest,] <- l[ltest,]/apply(l[ltest,,drop=FALSE],1,sum);
    }

    FraudVotes <- round(NVoters*a*l)

    FraudFlag <- runif(N) < f2;
    subflag <- runif(N) < f1;
    IncreFlag <- FraudFlag & subflag;
    ExtreFlag <- FraudFlag & !subflag;

    us <- abs(rnorm(sum(IncreFlag),sd=theta));    ##xi
    while (any(utest <- us>=1 | us<=0)) us[utest] <- abs(rnorm(sum(utest),sd=theta));
    use <- 1-abs(rnorm(sum(ExtreFlag),sd=0.075)); ##yi
    while (any(utest <- use>=1 | use<=0)) use[utest] <- abs(rnorm(sum(utest),sd=0.075));
    us0 <- rep(0.0, N);
    if (any(IncreFlag)) us0[IncreFlag] <- us;
    if (any(ExtreFlag)) us0[ExtreFlag] <- use;
    cv0 <- us0^alpha;
    FraudVotes[,1] <-
      FraudVotes[,1] + floor(us0*(NVoters-(FraudVotes[,1]+FraudVotes[,2]))) +  ##From non-attended.
      round(cv0*FraudVotes[,2]);                                               ##From opponents.
    FraudVotes[,2] <- FraudVotes[,2] - round(cv0*FraudVotes[,2])

    return(FraudVotes)
  }
  ##############
  ##Function 3
  ##Histogram-type-Vote Count  Function
  Sim.Histo<-function(FraudVotes,bins,N){
    vf<-matrix(0,nrow=N,ncol=2)
    for(i in 1:N){
      vf[i,1:2]<-FraudVotes[i,1:2]/sum(FraudVotes[i,])
    }
    Sim.H.Vote<-hist(vf[,1],breaks=bins,plot=F)$counts
    Sim.result<-list()
    Sim.result$Sim.H.Vote<-Sim.H.Vote
    Sim.result$vf<-vf
    return(Sim.result)
  }
  ###############
  ###Function 4
  ###Iteration Function
  Iteration_sim<-function(iterations,f1range,f2range,arange,p_Att,stdAtt,lambda_fraud,sigma_fraud,
                          data1,N,theta,bins,Obs.H.Vote,v, cores){
    ###Set FraudFits
    FraudFits<-matrix(NA,nrow=8,ncol=0)   ##I changed the code.
    ###Set Iterations-Loop

    registerDoParallel(cores=cores);

    for(it in 1:iterations){
      Ffits<-matrix(0,nrow=8,ncol=length(f1range)*length(f2range)*length(arange))
      #Fraud Model
      ##Create Comb, which contains every possible combination of fi, fe, and alpha.
      Comb<-expand.grid(f1range=f1range,f2range=f2range,arange=arange)
      #Set fi,fe, and alpha
      fi<-c()
      fe<-c()
      alpha<-c()
      f2<-c()
      f1<-c()
      FitsList <- foreach(z=1:(length(f1range)*length(f2range)*length(arange)), .export=c('Sim_Vote', 'Sim.Histo')) %dopar% {
        #      for(z in 1:(length(f1range)*length(f2range)*length(arange))){
        fi<-Comb[z,1]
        fe<-Comb[z,2]
        alpha<-Comb[z,3]
        Alpha<-Comb[z,3]
        f2<-fi+fe
        f1<-fi/(fi+fe)
        FraudVotes<-Sim_Vote(p_Att,stdAtt,lambda_fraud,sigma_fraud,data1$NVoters,N,f2,f1,theta,Alpha)
        ###we have simulated data for each district/ This is for each (fi,fe,alpha)set
        Sim.R<-Sim.Histo(FraudVotes,bins,N)
        ##Compute Statistics
        c(fi,fe,alpha,
          sum(((Sim.R$Sim.H.Vote-Obs.H.Vote)/(Obs.H.Vote+1))^2),
          sum(((Obs.H.Vote-Sim.R$Sim.H.Vote)^2)/(Sim.R$Sim.H.Vote+1)),
          sum(((data1$Votes-FraudVotes[,1])^2)/(FraudVotes[,1]+0.5)),
          sum(((data1$Votes-FraudVotes[,1])^2)/
                (FraudVotes[,1]+0.5)) +
            sum(((data1$NValid-data1$Votes-FraudVotes[,2])^2)/
                  (FraudVotes[,2]+0.5)) +
            sum(((-data1$NValid + apply(FraudVotes,1,sum))^2)/
                  (data1$NVoters-apply(FraudVotes,1,sum)+0.5)),
          length(v[,1])-9)
      }
      # stopImplicitCluster();

      for (j in 1:length(FitsList)) {
        Ffits[1:8,j] <- FitsList[[j]];
      }

      ind<-which(Ffits[4,]==min(Ffits[4,]))    ##I changed the code
      if (length(ind)>1) {
        FraudFits<-cbind(FraudFits,apply(Ffits[,ind],1,mean));
      }
      else {
        FraudFits<-cbind(FraudFits,Ffits[,ind])       ##I changed the code
      }
    }
    return(FraudFits)
  }

  ###Main Function
  ElectionFitter_sim <- function(data1,
                                 thres=5,
                                 f1range=seq(from=0,to=1,by=0.25),
                                 f2range=seq(from=0,to=1,by=0.25),
                                 arange=seq(from=0.5,to=2,by=0.5),
                                 iterations=2, cores=2) {
    ##Initialize some parameters
    N<-length(data1$Votes)
    Attendance<-data1$NValid/data1$NVoters

    ##Compute votes for the party and others.
    v<-matrix(NA,nrow=N,ncol=2)
    for(i in 1:N){
      v[i,1]<-data1$Votes[i]/data1$NValid[i]
      v[i,2]<-(data1$NValid[i]-data1$Votes[i])/data1$NValid[i]
    }
    ##Detect the local Maximum
    ##For Histgram
    while (TRUE) {
      x0<- -0.005*(5/thres);
      h<- 0.01*(5/thres);
      if (0.8*N > 1/h) break;
      if (thres>0.5) {
        thres <- thres - 0.5;
      }
      else {
        break;
      }
    }
    bins<-seq(from=x0,to=1-x0,by=0.01*(5/thres))
    #########

    First_Estimated<-Estimate(bins,thres,Attendance,v)
    ###Just Naming; I will fix this later
    lambda_fraud<-First_Estimated[1:2]
    p_Att<-First_Estimated[3]
    stdAtt<-First_Estimated[4]
    theta<-First_Estimated[5]
    sigma_fraud<-First_Estimated[6:7]

    #######Here, Changing the vector version successful
    ###Empirical vote distribution is stored in Obs.Vote
    Obs.H.Vote<-hist(v[,1],breaks=bins,plot=F)$counts

    FraudFits <- Iteration_sim(iterations,f1range,f2range,arange,p_Att,stdAtt,lambda_fraud,sigma_fraud,
                               data1,N,theta,bins,Obs.H.Vote,v, cores)

    FFnames <- c("incremental","extreme","alpha",
                 "Winner.HFit.Klimek*","Winner.HFit.chi2","Winner.Fit.chi2","chi2","df");
    if (length(dim(FraudFits))>1) {
      FF <- cbind(apply(FraudFits,1,mean), sqrt(apply(FraudFits,1,var)));
      dimnames(FF) <- list(FFnames,c("mean","sdev"));
    }
    else  {
      FF <- FraudFits;
      names(FF) <- FFnames;
    }

    results<-list()
    results$FF <- FF;
    results$FraudFits<-round(FraudFits,3)
    results$First_Estimated<-First_Estimated
    return(results)
  }

  Methods=7
  closeAllConnections()
  KSimNames <- paste("KSim", c("I","E","alpha","turnout","winprop","sigma","stdAtt","theta"),sep="");
  nullvar="-"
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
    totalReg <- totalVotes <- turnoutPercent <- data$Turnout <- NULL;
  }
  Candidates[grepl("All Leaders",Candidates)]<-"AllLeaders"

  #Obtain Compact dataset    #FOR KLIMEK
  data$AllLeaders<-1
  candidates<-data[, names(data)%in%Candidates, drop=FALSE];
  num_candidates<-length(Candidates)                      #FOR KLIMEK
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

  kidx <- !is.na(dataC$totalReg) & (dataC$totalReg >= dataC$totalVotes) &
    (dataC$totalVotes > 0) & (dataC$totalVotes >= dataC$Votes);
  if (any(is.na(kidx))) kidx[is.na(kidx)] <- FALSE;

  if(length(kidx)!=0) dataC <- dataC[kidx,];

  statsArray2Table <- function(arr) {
    #  resultsArray[1:levellength,1:num_candidates,resultsWidth,4]
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


  #FOR KLIMEK
  levels<-data[, names(data)%in%Level];
  methods<-as.numeric(Methods)
  sampleSize <- dim(data)[1];
  #Divide by level  1 for national
  if(length(levels)==0) {
    levels<-rep(1,dim(dataC)[1]);
    dataC$National<-levels
  }
  dataCC<-split(dataC, levels)
  levellength<-length(names(dataCC))
  if(levellength=='1') {
    names(dataCC)<-"National"
  }
  levelnames<-names(dataCC)

  dataCC<-split(dataC, levels)
  levellength<-length(names(dataCC))
  if(levellength=='1') {
    names(dataCC)<-"National"
  }
  levelnames<-names(dataCC)

  #Create holder for results   FOR KLIMEK
  resultsLength<-num_candidates*levellength*2
  resultsWidth <- length(KSimNames);  #I CHANGED LENGTH
  resultsMatrix <- matrix(NA,resultsLength,resultsWidth); #+2 to account for Level and Candiates columns

  #Loops
  resultsArray <-
    array(NA,dim=c(levellength,num_candidates,resultsWidth,4));

  dimnames(resultsArray) <- list(NULL, NULL, NULL, NULL);
  dimnames(resultsArray)[[1]] <- levelnames;
  dimnames(resultsArray)[[2]] <- Candidates;
  dimnames(resultsArray)[[3]] <- KSimNames   #CHANGED WORDS
  dimnames(resultsArray)[[4]] <- c("point","lo","up","obs");
  sourceArray <- resultsArray;
  sourceArray[,,,] <- TRUE;

  varssArray <- array(NA,dim=c(levellength,num_candidates));

  varReg <- ifelse("totalReg" %in% names(dataCC[[1]]),c("totalReg"),"");
  varVotes <- ifelse("totalVotes" %in% names(dataCC[[1]]),"totalVotes","")

  print('Toolkit is working.')


  ###Klimek et al. simulation approach
  f1r <- seq(from=0.0,to=1.0,by=0.2);       # min=0.0, max=1.0
  f2r <- seq(from=0.0,to=0.3,by=0.1);     # min=0.0, max=1.0
  ar <- seq(from=0.5,to=1,by=0.5);        # min=0.0
  iter <- 10 #!!!

  CandidatesKlimek <- Candidates[Candidates!="Turnout"]
  KlimekMatrix<- matrix(NA,resultsLength,2);
  thresSet <- c(5+c(0:9)/10);
  print("Computing Klimek et al. simulation.", quote=FALSE)
  for (candidID in 1:length(CandidatesKlimek)) {  #Loop over Candidates
    for (levelID in 1:levellength) {     #Loop over Levels
      print(paste("Analyzing ", CandidatesKlimek[candidID], ".", sep=""), quote=FALSE)
      dataSubset<-dataCC[[levelID]]
      numberUnits<-dim(dataSubset)[1]
      NVoters <- dataSubset$totalReg
      NValid  <- dataSubset$totalVotes
      Votes   <- dataSubset[,names(dataSubset) %in% CandidatesKlimek[candidID]]
      KlData<-data.frame(NVoters,NValid,Votes)
      KlData<-na.omit(KlData)
      KlData<-KlData[!c(KlData$NVoters==0|KlData$NValid==0|KlData$Votes==0),]
      for (k in 1:length(thresSet)) {
        KlimekResults <-
          try(ElectionFitter_sim(data=KlData, thres=thresSet[k],
                                 f1range=f1r, f2range=f2r, arange=ar, iterations=iter,
                                 cores=cores));
        if (!is.character(KlimekResults)) break;
        if (k==length(thresSet) && is.character(KlimekResults)) KlimekResults <- NULL;
      }
      if (is.null(KlimekResults)) next;
      resultsArray[levelnames[levelID],CandidatesKlimek[candidID],KSimNames[1:3],1] <-
        KlimekResults$FF[1:3,1];
      resultsArray[levelnames[levelID],CandidatesKlimek[candidID],KSimNames[1:3],2] <-
        KlimekResults$FF[1:3,2];
      resultsArray[levelnames[levelID],CandidatesKlimek[candidID],KSimNames[4:8],1] <-
        c(KlimekResults$First_Estimated[c(1,3,6,4,5)]);
      resultsArray[levelnames[levelID],CandidatesKlimek[candidID],KSimNames,4] <-
        KlimekResults$FF[8,1]+9;
    }
  }
  # mK <- unlist(MethNamesList[methods[methods == 7]]);
  # Cidx <- Candidates!="Turnout";
  #  sourceArray[,!Cidx,mK,] <- FALSE;
  #  if (any(sourceArray[,Cidx,mK,,drop=FALSE])) {
  #    statsArray2XML(dataID, dataNAME, dataDIM,
  #      resultsArray[,Cidx,mK,,drop=FALSE], sourceArray[,Cidx,mK,,drop=FALSE],
  #      regvar=varReg, totalvar=varVotes, leveltype=Level);
  #  }

  CandidatesNames<-rep(Candidates,each=levellength*2)
  LevelNames<-rep(levelnames,length(Candidates), each=2)
  SelectedMethods <- statsArray2Table(resultsArray);
  SelectedMethods[,dim(SelectedMethods)[2]][seq(2,dim(SelectedMethods)[1],2)]<-""
  SelectedMethods[,1][seq(2,length(SelectedMethods[,1]),2)]<-""
  SelectedMethods[,2][seq(2,length(SelectedMethods[,2]),2)]<-""
  SelectedMethods[is.na(SelectedMethods)] <- "--";
  #colnames(SelectedMethods) <- c("Level", "Candidate's Name",
  #  MethNames[MethodsMarkers %in% c(methods,0)]);


  coloredResults<-ColorSignificance(SelectedMethods, Klimek=TRUE)
  coloredResults
  texResults=xtable(SelectedMethods)
  Results<-list(SelectedMethods,coloredResults,texResults)
  names(Results) <- c("table", "html", "tex")
  return(Results)
}

