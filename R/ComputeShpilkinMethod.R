#' @title Shpilkin's Method
#' @description This function implements the revised version of Shpilkin's method (NB! It doesn't replicate original method fully).
#' @usage ComputeShpilkinMethod(data, Candidates, CandidatesText=NULL, MainCandidate, TotalReg, Level=NULL,
#'                              Methodmax="M1", MaxtThreshold=0.8, WindowSize=5, FigureName)
#' @param data electoral data.
#' @param Candidates vector of variable names for all candidates/parties participated in the election
#' @param CandidatesText vector of candidates/parties' names participated in the election used to draw the figure
#' @param MainCandidate variable name for main/incumbent candidate
#' @param TotalReg  variable name for the total number of eligible voters
#' @param Level variable name depicting the level of analysis ("National", i.e. whole dataset by default)
#' @param Methodmax clean peak search
#' \itemize{
#'   \item M1 - absolute clean peak search on the left handside from official turnout
#'   \item M2 - relative search on the left handside from official turnout within a range defined by WindowSize
#' }
#' @param MaxtThreshold  anomalous turnout threshold (by default 0.8)
#' @param WindowSize define WindowSize for M1.  Algorithm searches for max change within prespecified WindowSize (by default WindowSize=5\%)
#' @param FigureName  figure's name
#' @export
#' @import ggplot2
#' @return list containing results of analysis
#' \itemize{
#'   \item list_graphs - list of graphs
#'   \item stats_table -  table with results
#'   \item Level -  external parameter
#'   \item creationdate - date/time of analysis
#' }
#' @examples
#' library(EFToolkit)
#'
#' dat<-read.csv(system.file("ruspres2018.csv", package="EFToolkit"))
#' res<-ComputeShpilkinMethod(dat, Candidates=c("P1", "P2",  "P3",  "P4", "P5", "P6", "P7", "P8" ),
#'                     CandidatesText=c("Baburin", "Grudinin",  "Zhirinovsky",  "Putin", "Sobchak",
#'                              "Suraikin", "Titov", "Yavlinsky"),
#'                     MainCandidate="P4",
#'                     TotalReg="NVoters",
#'                     Methodmax="M1",
#'                     FigureName="Russian Presidential Elections, 2018",
#'                     Level="region",
#'                     MaxtThreshold=0.85)



############################################################
##               Election Forensics Toolkit               ##
##  24oct2019 by Kirill Kalinin                           ##
##  Kirill Kalinin and Walter R. Mebane, Jr               ##
############################################################


ComputeShpilkinMethod<-function(data, Candidates, CandidatesText = NULL, MainCandidate, TotalReg, Level = NULL,
                         Methodmax = "M1", MaxtThreshold = 0.8, WindowSize = 5, FigureName){


  estimate_fraud<-function(data,Candidates,CandidatesText, MainCandidate, TotalReg, Methodmax,
                           FigureName, WindowSize, MaxtThreshold, whole.real.turnout = NULL){

    MainCandidate.v <- data[,names(data) %in% MainCandidate];
    acandidates.v <- data[,names(data) %in% Candidates];
    mcandidateM <- Candidates %in% MainCandidate

    if(is.null(CandidatesText)) CandidatesText<- Candidates
    ccandidate <- CandidatesText[mcandidateM]

    valid.v <- apply(acandidates.v, 1, function(x) sum(x, na.rm = TRUE))
    TotalReg.v <- data[,names(data) %in% TotalReg];
    mdat <- data.frame(valid.v, MainCandidate.v,   TotalReg.v, acandidates.v)
    turnout = mdat$valid.v/mdat$TotalReg.v
    official.turnout = sum(mdat$valid.v, na.rm = TRUE)/sum(mdat$TotalReg.v,na.rm = TRUE)
    incumbent.support = sum(mdat$MainCandidate.v , na.rm = TRUE)/sum(mdat$valid.v, na.rm = TRUE)

    i = 0
    ruler<-seq(0,1.01,.01)
    res_matrix<-matrix(NA,103,length(Candidates))
    while(i <= 102){
      res_matrix[i+1,]<-sapply(1:length(Candidates), function(j) {
        sum(mdat[,3+j][turnout >= ruler[i]&turnout<ruler[i+1]], na.rm = TRUE)})
      i = i + 1
    }
    #adjustment
    res_matrix<-res_matrix[-c(1, 103),]; ruler = ruler[-101]

    percentUR <- res_matrix[,mcandidateM]/apply(res_matrix,1,function(x) sum(x, na.rm = TRUE))
    peakyL <- ruler[which(sign(diff(percentUR)) == -1)]
    peakyR <- ruler[which(sign(diff(percentUR)) == 1)]
    max.peaky <- which.max(apply(res_matrix,1,sum))/100
    peaky <- c(peakyL[peakyL <= max.peaky], peakyR[peakyR >= max.peaky])

    stepp <- c(0,diff(peaky))
    max_height <- res_matrix[which(ruler %in% peaky),1]

    peaks_matrix <- data.frame(peaky, stepp, max_height)
    peaks_matrix_less <- peaks_matrix[peaky < official.turnout,]


    if(Methodmax == "M1"){
      max.peak.turnout <- peaks_matrix_less$peaky[which.max(peaks_matrix_less$max_height)]
    }

    if(Methodmax == "M2"){
      i=0
      l=dim(peaks_matrix_less)[1]

      if(length(peaks_matrix_less$peaky) < WindowSize){
        WindowSize = length(peaks_matrix_less$peaky);
        warning("adjusting for a Window size")}

      k <- which.max(peaks_matrix_less$max_height)

      while((l - i - WindowSize + 1) > 0){
        selected = seq(l - i, l - i - WindowSize + 1, -1)
        if(any(peaks_matrix_less$stepp[selected] > 0.01)){
          o <- which.max(peaks_matrix_less$stepp[selected])
          p <- o[length(o)]
          k <- selected[p]
          break
        }
        i = i + 1
      }
      max.peak.turnout = peaks_matrix_less$peaky[k]
    }

    if(official.turnout > MaxtThreshold & !is.null(whole.real.turnout)){
                                            max.peak.turnout <- whole.real.turnout}

    mat100 <- length(res_matrix[,1])
    clean.votes <- res_matrix[which(ruler == max.peak.turnout),mcandidateM]
    clean.valid.votes <- sum(res_matrix[which(ruler == max.peak.turnout),], na.rm=TRUE)
    clean.all.butUR <- sum(res_matrix[which(ruler == max.peak.turnout), !mcandidateM], na.rm=TRUE)
    clean.ur.prop <- clean.votes/clean.valid.votes
    clean.all.butUR.prop <- 1 - clean.ur.prop
    inflation.factor <- clean.votes / clean.all.butUR
    clean.votes.vector <- apply(res_matrix[, !mcandidateM],1,function(x) sum(x, na.rm=TRUE)) * inflation.factor
    clean.votes.total <- sum(clean.votes.vector, na.rm=TRUE)

    if (max.peak.turnout > 1){k1 = 1}else{k1 = 100}

    magnitude.election.fraud<-
      sum(res_matrix[(round(max.peak.turnout * k1 + 1, 0)):101, mcandidateM], na.rm=TRUE) -
      sum(clean.votes.vector[(round(max.peak.turnout * k1 + 1, 0)):101], na.rm=TRUE)

    realURsupport <- round(percentUR[round(max.peak.turnout * k1, 0) + 1]*100, 0)

    ballot_stuffing <- round((official.turnout - max.peak.turnout) * sum(mdat$TotalReg.v,na.rm = TRUE), 0)
    if(ballot_stuffing > magnitude.election.fraud){ballot_stuffing <- round(magnitude.election.fraud, 0)}
    ballot_switching <- round((magnitude.election.fraud - ballot_stuffing) / 2, 0)

    #Setting up the table
    colors_vector = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    colors <- sample(colors_vector,  dim(res_matrix)[2])
    colors[mcandidateM] <- "chartreuse4"
    colors <- c("grey", colors)
    colV <- as.vector(matrix(rep(colors,each = mat100),nrow = mat100))

    #party vector
    if(length(clean.votes.vector) == 0){clean.votes.vector = 0}
    ggparties <- cbind(clean.votes.vector, res_matrix)
    votV <- as.vector(as.matrix(ggparties))

    rulV <- rep(ruler*100, times = dim(ggparties)[2])

    namV <- as.vector(
                    matrix(rep(c("Clean votes", paste("P/C: ", CandidatesText, sep = "")),
                               each = mat100),nrow = mat100))
    #CI
    u <- ifelse(namV == "Clean votes", votV, NA)
    l <- ifelse(namV == paste("P/C: ", ccandidate, sep = ""), votV, NA);
    lna <- rep(NA,length(namV)); lna[1:101] <- l[!is.na(l)]
    l <- lna

    ggdata <- data.frame(namV, rulV, colV, votV, u, l)
    o.turnout <- round(official.turnout*100,digits = 1)
    r.turnout <- round(max.peak.turnout*100,digits = 1)
    turnout.info <- paste("Official turnout: ",
                          round(o.turnout,digits = 1), "%",
                          "   Real turnout: ",
                          round(r.turnout,digits = 1), "%", sep = "")
    support.info <- paste("Official support: ",
                          round(incumbent.support * 100, 1), "%",
                          "   Real support: ",
                          round(realURsupport, digits = 1), "%", sep = "")
    m.fraud.info <- paste("Fraud: ", round(magnitude.election.fraud, 0), sep = "")

    drawfig <- ggplot(ggdata, aes(x = rulV, y = votV, color = namV, group = namV)) +
      scale_x_discrete(limits = seq(0,100,10)) +
      theme_bw() +
      geom_vline(aes(xintercept = o.turnout), color="blue", linetype="dashed", size = 1, show.legend = F) +
      geom_vline(aes(xintercept = r.turnout), color="grey", linetype="dashed", size = 1, show.legend = F) +
      geom_line(size = 2) +
      geom_ribbon(aes(ymin = l,ymax = u), fill = "grey", alpha = 0.5, show.legend = F, color = NA)+
      theme(legend.title = element_blank(),
            legend.text = element_text(size = 10),
            legend.key = element_blank(),
            legend.background = element_blank()) +
      labs(title = FigureName,
           y = "Number of Votes", x = "Turnout") +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(caption = paste("\n", turnout.info,"\n",
                           support.info,"\n",
                           m.fraud.info, sep = ""))

    fraud_measure <- round(magnitude.election.fraud, 0)
    votes_incumbent <- sum(mdat$MainCandidate.v)

    results <- list(official_turnout = o.turnout, real_turnout = r.turnout,
                    official_support = round(incumbent.support*100, 1), real_support = round(realURsupport, 1),
                    ballot_stuffing = ballot_stuffing, ballot_switching = ballot_switching,
                    total_fraud = fraud_measure, drawfig = drawfig)
    return(results)}

  gresults <- list()

  if(!is.null(Level)){
    gresults[['Whole dataset']] <- estimate_fraud(data,Candidates,CandidatesText,
                                                  MainCandidate, TotalReg, Methodmax,
                                                  FigureName, WindowSize, MaxtThreshold)

    whole.real.turnout <- gresults[['Whole dataset']]$real_turnout/100

    if(Level!="National"){
      splby <- data[,names(data)%in%Level]
      mdat <- split(data, splby)

      for(i in 1:length(names(mdat))){
        mitem <- mdat[[i]]
        FigureName <- names(mdat)[i]
        gresults[[FigureName]] <- estimate_fraud(mitem,Candidates,CandidatesText,
                                                 MainCandidate, TotalReg, Methodmax,
                                                 FigureName, WindowSize = 5,
                                                 whole.real.turnout = whole.real.turnout,
                                                 MaxtThreshold)
      }
      }

    }else{

      gresults[['Whole dataset']] <- estimate_fraud(data,Candidates,CandidatesText,
                                                    MainCandidate, TotalReg, Methodmax,
                                                    FigureName, WindowSize, MaxtThreshold)
      }

  lst_lngth <- length(gresults[[1]])
  list_graphs  <-  lapply(gresults, function(x) x[[lst_lngth]])
  stats_table  <-  do.call(rbind.data.frame,
                         lapply(gresults, function(x) x[1:(lst_lngth-1)]))

  list_results <- list(list_graphs = list_graphs, stats_table = stats_table,
                       Level = Level, creationdate = Sys.time())
}
