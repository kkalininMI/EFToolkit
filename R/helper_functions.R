#' Internal functions for election fraud detection
#'
#' Various functions for detecting election anomalies using statistical methods
#'
#' @name internal-functions
#' @keywords internal
#' @export
#' @import survey
#' @import ggplot2
#' @import stringi
#' @import plotly
#' @import reshape2
#' @import mixtools
#' @import cluster
#' @importFrom boot boot
#' @importFrom pracma findpeaks
#' @importFrom robustbase covMcd
#' @importFrom dbscan dbscan
#' @importFrom dbscan kNNdist
#' @return Returns various data outputs.


var_addition <- function(data,
                         Level,
                         MainCandidate,
                         Candidates,
                         TotalVotes = NULL,
                         CandidatesText = NULL,
                         TotalReg) {

  if(Level!="National"){
    region <- as.character(data[,names(data) %in% Level])
    data$id <- paste(region, ave(seq_len(nrow(data)), region, FUN = seq_along), sep = ":")}

  if(Level=="National"){
    data$id <- ave(seq_len(nrow(data)), FUN = seq_along)}

  mcandidate_idx <- which(Candidates == MainCandidate)
  if (length(mcandidate_idx) == 0) stop("MainCandidate not found in Candidates list")

  # Compute Valid votes
  if (!is.null(TotalVotes)) {
    NValid <- if (length(TotalVotes) > 1) {
      rowSums(data[, TotalVotes], na.rm = TRUE)
    } else {
      data[[TotalVotes]]
    }
  } else {
    NValid <- rowSums(data[,names(data) %in% Candidates], na.rm = TRUE)
  }

  Votes <- data[, MainCandidate, drop = TRUE]
  NVoters <- data[, TotalReg, drop = TRUE]

  # Compute per-row and global metrics
  turnout <- NValid / NVoters
  incumbent <- Votes / NValid  #incumbent.support.v

  # Compute official metrics
  official.turnout <- sum(NValid, na.rm = TRUE) / sum(NVoters, na.rm = TRUE)
  incumbent.support <- sum(Votes, na.rm = TRUE) / sum(NValid, na.rm = TRUE)


  data$Votes <- Votes
  data$NVoters <- NVoters
  data$NValid <- NValid
  data$turnout <- turnout
  data$incumbent <- incumbent

  # Attach summary attributes
  attr(data, "official.turnout") <- official.turnout
  attr(data, "incumbent.support") <- incumbent.support

  return(data)
}


estimate_fraud<-function(data, Candidates, CandidatesText, MainCandidate,
                         FigureName, MaxThreshold,
                         colors.v, precinct.level, mode_search, man_turnout){

  cand_vo <- data[, names(data) %in% Candidates];
  mcand_id <- Candidates %in% MainCandidate
  turnout <- data$turnout
  incumbent<- data$incumbent
  incumbent.support = sum(data$Votes, na.rm = TRUE) / sum(data$NValid, na.rm = TRUE) #attr(data, "incumbent.support")
  official.turnout = sum(data$NValid, na.rm = TRUE) / sum(data$NVoters, na.rm = TRUE) #attr(data, "official.turnout")

  completeDat = TRUE
  columns_to_check <- c(MainCandidate, Candidates)

  for (col in columns_to_check) {
    if (all(is.na(data[[col]]))) {
      print(paste("Incomplete data found in column:", col))
      completeDat = FALSE
    }
  }

  if (completeDat==FALSE){
    results <- list(official_turnout = NA, real_turnout = NA, official_support = NA,
                    real_support = NA, ballot_stuffing = NA, ballot_switching = NA,
                    total_fraud = NA, prop_fraud = NA, drawfig = NA,
                    precinct_fraud = list(id = NA, clean.votes = NA, fraud.votes = NA))

    estimate_fraud_res <- list(
      official_turnout   = NA,
      real_turnout       = NA,
      official_support   = NA,
      real_support       = NA,
      ballot_stuffing    = NA,
      ballot_switching   = NA,
      total_fraud        = NA,
      prop_fraud         = NA,
      drawfig            = NA,
      precinct_fraud     = data.frame(
        id = NA,
        clean.votes = NA,
        fraud.votes = NA
      )
    )
    return(results)}


  if(is.null(CandidatesText)) CandidatesText<- Candidates
  ccandidate <- CandidatesText[mcand_id]

  mdat <- data.frame(
    subset(data, select=c("NValid", "Votes", "NVoters")), cand_vo)

  candidate_matrix <- matrix(NA,103,length(Candidates))
  valid100 <- registered100 <- mcandidate100 <- rep(NA,103)

  i = 0
  ruler <- seq(0,1.01,.01)
  while(i <= 102){
    candidate_matrix[i+1,]<-sapply(1:length(Candidates), function(j) {
      sum(mdat[,3+j][turnout > ruler[i] & turnout<=ruler[i+1]], na.rm = TRUE)})
    valid100[i+1] <- sum(mdat[,1][turnout > ruler[i] & turnout <= ruler[i+1]], na.rm = TRUE)
    mcandidate100[i+1] <- sum(mdat[,2][turnout > ruler[i]&turnout <= ruler[i+1]], na.rm = TRUE)
    registered100[i+1] <- sum(mdat[,3][turnout > ruler[i]&turnout <= ruler[i+1]], na.rm = TRUE)
    i = i + 1
  }

  candidate_matrix <- candidate_matrix[-c(1, 103),]; ruler = ruler[-101]
  valid100 <- valid100[-c(1, 103)];
  registered100 <- registered100[-c(1, 103)];
  mcandidate100 <- mcandidate100[-c(1, 103)];

  percentMainCandidateSupport <-  mcandidate100 / valid100

  peaky <- ruler[!is.na(percentMainCandidateSupport)]
  max_height <-   mcandidate100[which(ruler %in% peaky)]
  peaks_matrix <- data.frame(peaky, max_height)
  #peaks_matrix_less <- peaks_matrix[peaky < official.turnout,]
  #peaks_matrix_less <- peaks_matrix

  if (!is.null(mode_search) & is.null(man_turnout)){

    if (mode_search$pick_by %in% c("height", "area", "cluster")){
      max.peak.turnout = peak_search_basic(peaks_matrix, mode_search, official.turnout, peaky)
    } else if (mode_search$pick_by %in% "quantile"){
      max.peak.turnout = peak_search_quantile(turnout)
    }

    if (mode_search$pick_by %in% c("height", "area", "cluster", "quantile")) {
      peak_row_idx <- which(as.character(ruler) == as.character(max.peak.turnout))
      clean.votes <- candidate_matrix[peak_row_idx,mcand_id]

      clean.valid.votes <- sum(candidate_matrix[peak_row_idx, ], na.rm = TRUE)
      clean.all.butMainCandidate <- valid100[peak_row_idx] - clean.votes

      clean.ur.prop <- clean.votes / clean.valid.votes
      clean.all.butMainCandidate.prop <- 1 - clean.ur.prop

      inflation.factor <- clean.votes / clean.all.butMainCandidate

      other_candidate_votes <- valid100 - candidate_matrix[,mcand_id]
      clean.votes.vector <- other_candidate_votes * inflation.factor
      clean.votes.total <- sum(clean.votes.vector, na.rm = TRUE)

      incumbent.share = mdat$Votes/mdat$NValid
      votes_peak<-incumbent.share[as.character(round(turnout, 2))==as.character(max.peak.turnout)]
      clean_SD <- sd(votes_peak, na.rm=TRUE)
    }

    if (mode_search$pick_by %in% "elipse") {
      elipseMCD <- computeMCDforKmeans(data$Votes, data$NValid, data$NVoters)
      draw_elipse <- data.frame(elipseMCD$elipsedat)
      max.peak.turnout <- elipseMCD$clean_mean[1]/100

      clean.votes <- elipseMCD$clean_rawmean[1]
      clean.valid.votes <- elipseMCD$clean_rawmean[2]
      clean.all.butMainCandidate <- clean.valid.votes - clean.votes
      clean.ur.prop <- clean.votes/clean.valid.votes
      clean.all.butMainCandidate.prop <- 1 - clean.ur.prop

      inflation.factor <- clean.votes / clean.all.butMainCandidate

      clean.votes.vector <- (valid100 - candidate_matrix[, mcand_id]) * inflation.factor
      clean.votes.total <- sum(clean.votes.vector, na.rm=TRUE)
      clean_SD <- sqrt(diag(elipseMCD$covmat))}

    }else if (!is.null(man_turnout)){
      if (man_turnout > 1) man_turnout = man_turnout/100
      max.peak.turnout <- man_turnout

      peak_row_idx <- which(as.character(ruler) == as.character(max.peak.turnout))
      clean.votes <- candidate_matrix[peak_row_idx,mcand_id]

      clean.valid.votes <- sum(candidate_matrix[peak_row_idx, ], na.rm = TRUE)
      clean.all.butMainCandidate <- valid100[peak_row_idx] - clean.votes

      clean.ur.prop <- clean.votes / clean.valid.votes
      clean.all.butMainCandidate.prop <- 1 - clean.ur.prop

      inflation.factor <- clean.votes / clean.all.butMainCandidate

      other_candidate_votes <- valid100 - candidate_matrix[,mcand_id]
      clean.votes.vector <- other_candidate_votes * inflation.factor
      clean.votes.total <- sum(clean.votes.vector, na.rm = TRUE)

      incumbent.share = mdat$Votes/mdat$NValid
      votes_peak<-incumbent.share[as.character(round(turnout, 2))==as.character(max.peak.turnout)]
      clean_SD <- sd(votes_peak, na.rm=TRUE)
    }

  if (is.na(max.peak.turnout)) {max.peak.turnout = official.turnout}
  if (max.peak.turnout > 1){k1 = 1}else{k1 = 100}

  magnitude.election.fraud <-
    sum(candidate_matrix[(round(max.peak.turnout * k1 + 1, 0)):101, mcand_id], na.rm=TRUE) -
    sum(clean.votes.vector[(round(max.peak.turnout * k1 + 1, 0)):101], na.rm=TRUE)

  realMainCandidateSupport <-
    round(percentMainCandidateSupport[round(max.peak.turnout * k1, 0) + 1] * 100, 0)

  ballot_stuffing <- round((official.turnout - max.peak.turnout) * sum(mdat$NVoters,na.rm = TRUE), 0)

  if(ballot_stuffing > magnitude.election.fraud){ballot_stuffing <- round(magnitude.election.fraud, 0)}

  ballot_switching <- magnitude.election.fraud - ballot_stuffing #round((magnitude.election.fraud - ballot_stuffing) / 2, 0)

  mat100 <- length(candidate_matrix[,1])
  colV <- as.vector(matrix(rep(colors.v,each = mat100),nrow = mat100))

  if(length(clean.votes.vector) == 0){clean.votes.vector = 0}
  ggparties <- cbind(clean.votes.vector, candidate_matrix)
  votV <- as.vector(as.matrix(ggparties))

  rulV <- rep(ruler*100, times = dim(ggparties)[2])

  namV <- as.vector(
    matrix(rep(c("Clean votes", paste("P/C: ", CandidatesText, sep = "")),
               each = mat100), nrow = mat100))
  #CI
  u <- ifelse(namV == "Clean votes", votV, NA)
  l <- ifelse(namV == paste("P/C: ", ccandidate, sep = ""), votV, NA);
  lna <- rep(NA,length(namV)); if(!all(is.na(l))) lna[1:101] <- l[!is.na(l)]
  l <- lna

  ggdata <- data.frame(namV, rulV, colV, votV, u, l)
  o.turnout <- round(official.turnout * 100, digits = 1)
  r.turnout <- round(max.peak.turnout * 100, digits = 1)
  turnout.info <- paste("Official turnout: ",
                        round(o.turnout,digits = 1), "%",
                        "   Real turnout: ",
                        round(r.turnout,digits = 1), "%", sep = "")
  support.info <- paste("Official support: ",
                        round(incumbent.support * 100, 1), "%",
                        "   Real support: ",
                        round(realMainCandidateSupport, digits = 1), "%", sep = "")
  m.fraud.info <- paste("Anomalies: ", round(magnitude.election.fraud, 0), sep = "")

  names(colors.v) <-unique(namV)
  colScale <- scale_colour_manual(name = "namV", values = colors.v)

  drawfig <- ggplot2::ggplot(ggdata, aes(x = rulV, y = votV, color = namV, group = namV)) +
    scale_x_continuous(breaks = seq(0,100,10)) +
    theme_bw() +
    geom_vline(aes(xintercept = o.turnout), color="blue", linetype="dashed", size = 1, show.legend = F) +
    geom_vline(aes(xintercept = r.turnout), color="grey", linetype="dashed", size = 1, show.legend = F) +
    geom_line(size = 2, na.rm=TRUE) +
    colScale +
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
  votes_incumbent <- sum(mdat$Votes)

  if(mode_search$pick_by %in% "elipse"){

    turnoutP <- data$NValid/data$NVoters * 100
    voteshareMainCandidate <- data$Votes/data$NValid * 100
    voteshareCandidates <- apply(data[,names(data) %in% Candidates], 2, function(x) x/data$NValid)
    votV <- as.vector(as.matrix(voteshareCandidates))*100
    rulV <- rep(turnoutP, times = ncol(voteshareCandidates))
    namV2 <- as.vector(rep(paste("P/C: ", CandidatesText, sep = ""), each = nrow(data)))
    colV <- rep(colors.v[-1],each = nrow(data))

    ggdata <- data.frame(namV2, rulV, colV, votV)

    colors.v <- colors.v[-1]
    names(colors.v) <- unique(namV2)
    stats.info <- paste("Clean cluster      Turnout ", round(elipseMCD$clean_mean[1], 1),
                        "  Vote percentage ", round(elipseMCD$clean_mean[2], 1), sep = "")

    colScale <- scale_colour_manual(name = "namV2", values = colors.v)

    drawfig2 <- ggplot2::ggplot(NULL) +
      scale_x_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0,100,10), expand = c(0, 0)) +
      theme_bw() +
      ggplot2::geom_point(data = ggdata, aes(x = rulV, y = votV, color = namV2, group = namV2), size=1, shape=".") +
      theme(legend.title = element_blank(),
            legend.text = element_blank(),
            legend.key = element_blank(),
            legend.background = element_blank()) +
      colScale +
      ggplot2::geom_point(data = draw_elipse, aes(x = X1, y = X2), color ="black", size = 4, shape=".") +
      labs(caption = paste("\n", stats.info,"\n", sep = "")) +
      labs(title = FigureName,
           y = "Vote percentage", x = "Turnout") +
      theme(plot.title = element_text(hjust = 0.5))

    ggdata2 <- data.frame(turnoutP, voteshareMainCandidate)/100

    drawfig3 <- ggplot2::ggplot(data = ggdata2, aes(x = turnoutP, y = voteshareMainCandidate)) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
      scale_fill_distiller(palette=4, direction=-1) +
      scale_x_continuous(breaks = seq(0,1,0.1), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq(0,1,0.1), expand = c(0, 0)) +
      theme(legend.position='none') +
      labs(title = FigureName, y = "Vote share", x = "Turnout")

    drawfig = list(drawfig, drawfig2, drawfig3)

  }else{
    drawfig = list(drawfig)
  }
  if(precinct.level){

    precinct.fraud <- compute_precinctf(data = data, clean.votes.vector = clean.votes.vector,
                                        magnitude.election.fraud = magnitude.election.fraud, postStratify = TRUE)


  }else{
    precinct.fraud <- NULL
  }

  prop.fraud <- fraud_measure/sum(data$Votes, na.rm=TRUE)

  results <- list(official_turnout = o.turnout, real_turnout = r.turnout,
                  official_support = round(incumbent.support*100, 1), real_support = round(realMainCandidateSupport, 1),
                  ballot_stuffing = ballot_stuffing, ballot_switching = ballot_switching,
                  total_fraud = fraud_measure, prop_fraud = prop.fraud, drawfig = drawfig, precinct_fraud = precinct.fraud)

  return(results)}



area_estimator<-function(x_mind, elmat){
  elmat2<-elmat
  elmat2$ruler <- elmat2$ruler * 100
  y_mind <- elmat2$mcandidate100[which(as.character(elmat2$ruler)%in%as.character(x_mind))]

  reg_m <- elmat2[1:which(as.character(elmat2$ruler)%in%as.character(x_mind)),]
  reg_m <-reg_m[order(reg_m$ruler, decreasing = TRUE),]

  s<-2; if(nrow(reg_m)==1) s<-1

  for (i in s:nrow(reg_m)){
    num_peak_start <- reg_m$ruler[i]
    if (reg_m$mcandidate100[reg_m$ruler%in%num_peak_start] >= y_mind) {
      x_boundary <- num_peak_start
      y_boundary <- reg_m$mcandidate100[reg_m$ruler%in%num_peak_start]
      full_rect_sq <- (x_mind - x_boundary + 1) * y_mind
      part_rect <- sum(reg_m$mcandidate100[c(which(as.character(reg_m$ruler)%in%as.character(x_mind)):
                                               (which(reg_m$ruler%in%x_boundary)-1),
                                             which(as.character(reg_m$ruler)%in%as.character(x_mind)))], na.rm=TRUE)
      resL <- full_rect_sq - part_rect
      break
    }else{
      resL <- 0
    }
  }

  reg_m <- elmat2[which(as.character(elmat2$ruler)%in%as.character(x_mind)):length(elmat2$ruler),]
  reg_m <-reg_m[order(reg_m$ruler, decreasing = FALSE),]


  s<-2; if(nrow(reg_m)==1) s<-1

  for (i in s:nrow(reg_m)){
    num_peak_start <- reg_m$ruler[i]
    if (reg_m$mcandidate100[reg_m$ruler%in%num_peak_start] >= y_mind) {
      x_boundary <- num_peak_start
      y_boundary <- reg_m$mcandidate100[reg_m$ruler%in%num_peak_start]
      full_rect_sq <- (x_boundary - x_mind + 1) * y_mind
      part_rect <- sum(reg_m$mcandidate100[c(which(as.character(reg_m$ruler)%in%as.character(x_mind)):
                                               (which(reg_m$ruler%in%x_boundary)-1),
                                             which(as.character(reg_m$ruler)%in%as.character(x_mind)))], na.rm=TRUE)
      resR <- full_rect_sq - part_rect
      break;}else{
        resR <- 0
      }
  }
  res <- list(resL = resL, resR = resR)
  return(res)}


computeMCDforKmeans <- function(Votes, Valid, NVoters){

  turnoutP <- Valid/NVoters * 100
  voteshareMainCandidate <- Votes/Valid * 100

  dat <- data.frame(Votes, Valid, NVoters, turnoutP, voteshareMainCandidate)
  filter_dat <- apply(is.na(dat), 1, any) #na.omit(dat)
  dat<-dat[!filter_dat,]


  kmres <- stats::kmeans(subset(dat, select=c("turnoutP", "voteshareMainCandidate")),
                  centers=3, iter.max = 10, nstart = 1, trace=FALSE)

  clean_id <- which.min(kmres$centers[,"turnoutP"])
  clean_dat <- dat[kmres$cluster==clean_id,]
  clean_center <- kmres$centers[clean_id,]

  extrem_id <- which.max(kmres$centers[,"turnoutP"])
  extrem_dat <- dat[kmres$cluster==extrem_id,]

  increm_id <- rownames(kmres$centers)[!rownames(kmres$centers)%in%as.character(c(clean_id, extrem_id))]
  increm_dat <- dat[kmres$cluster==increm_id,]

  cmcd <- covMcd(subset(dat, select=c("turnoutP", "voteshareMainCandidate")), nrow(dat), alpha=.50)

  covmat <- cmcd$cov
  evals <- eigen(covmat)$values
  evecs <- eigen(covmat)$vectors

  a <- seq(0, 2*pi, len=1000)
  c2 <- qchisq(0.5, 2)
  c <- sqrt(c2)

  xT <- c * sqrt(evals[1]) * cos(a)
  yT <- c * sqrt(evals[2]) * sin(a)
  M <- cbind(xT, yT)

  transM <- evecs %*% t(M)
  transM <- t(transM)

  elipsedat <- transM + rep(as.vector(clean_center), each = nrow(transM))

  datc <- dat[cmcd$best,]
  meanTV <- apply(datc, 2, function(x) sum(x, na.rm=TRUE))
  turnoutCL <- meanTV[2]/meanTV[3]*100
  voteperCL <- meanTV[1]/meanTV[2]*100

  result <- list(clean_mean = c(turnoutCL, voteperCL),
                 clean_rawmean = meanTV,
                 best = cmcd$best,
                 raw.weights = cmcd$raw.weights,
                 clean_center = clean_center,
                 covmat=covmat,
                 elipsedat=elipsedat)

  return (result)}


ComputeSimulationParametric <- function(data, Candidates, CandidatesText,
                                        MainCandidate, FigureName, MaxThreshold,
                                        colors.v,
                                        precinct.level, sims,
                                        mode_search,
                                        man_turnout,
                                        grid_type = "2D"){

  data_filter <- is.na(data$turnout)|is.na(data$incumbent)|
                 apply(data[Candidates]/data$NValid, 1, function(row) any(is.na(row)))
  if(any(data_filter==TRUE)){data <- data[!data_filter,]}
  candidates <- data[Candidates] / data$NValid
  id <- data$id
  n <- nrow(data)
  NVoters <- data$NVoters

  res_parameters <- matrix(NA, nrow = 4,  ncol = sims)
  res_precincts <- matrix(NA, nrow = n, ncol = sims)

  for (i in 1:sims) {
    turnoutV <- rbinom(n, size = data$NVoters, prob = data$turnout)
    candidates_sim <- apply(candidates, 2, function(x)  rbinom(n, size = turnoutV, prob = x))

    NValid_sim <- apply(candidates_sim, 1, function(x) sum(x, na.rm=TRUE))
    Votes_sim <- candidates_sim[,as.character(MainCandidate)]

    turnout_sim <- NValid_sim / data$NVoters
    incumbent_sim <- Votes_sim / NValid_sim  #incumbent.support.v

    official.turnout_sim <- sum(NValid_sim, na.rm = TRUE) / sum(data$NVoters, na.rm = TRUE)
    incumbent.support_sim <- sum(Votes_sim, na.rm = TRUE) / sum(NValid_sim, na.rm = TRUE)
    dat_sim <-data.frame(id, candidates_sim, NValid_sim, Votes_sim, turnout_sim, incumbent_sim, NVoters)

    colnames(dat_sim)<-gsub("_sim", "", colnames(dat_sim))
    attr(dat_sim, "official.turnout") <- official.turnout_sim
    attr(dat_sim, "incumbent.support") <- incumbent.support_sim

    if (grid_type=="1D"){
      estimate_fraud_res <- estimate_fraud(dat_sim, Candidates, CandidatesText, MainCandidate,
                                           FigureName, MaxThreshold, colors.v,
                                           precinct.level,
                                           mode_search = mode_search,
                                           man_turnout = man_turnout)
      }else{
        estimate_fraud_res <- estimate_fraud_2d(dat_sim, Candidates, CandidatesText, MainCandidate,
                                           FigureName, MaxThreshold, colors.v,
                                           precinct.level,
                                           mode_search = mode_search,
                                           man_turnout = man_turnout,
                                           graphsoff = TRUE)}

    res_parameters[,i]<- c(estimate_fraud_res$real_turnout, estimate_fraud_res$real_support,
                           estimate_fraud_res$total_fraud, estimate_fraud_res$prop_fraud)
    res_precincts[,i] <- estimate_fraud_res$precinct_fraud$fraud.votes
  }

  param_mean<-apply(res_parameters,1, function(x) mean(x, na.rm=TRUE))
  param_sd <- apply(res_parameters,1, function(x) sd(x, na.rm=TRUE))
  param_stats <- c(param_mean[1], param_sd[1], param_mean[2], param_sd[2], param_mean[3], param_sd[3], param_mean[4], param_sd[4])
  names(param_stats) <- c("real_turnoutM", "real_turnoutSD", "real_supportM","real_supportSD",
                          "total_fraudM", "total_fraudSD", "prop_fraudM", "prop_fraudSD")
  precinct_mean <- apply(res_precincts,1, function(x) mean(x, na.rm=TRUE))
  precinct_sd <- apply(res_precincts,1, function(x) sd(x, na.rm=TRUE))
  precinct_stats <-  data.frame(id=id, precinct_mean, precinct_sd)

  res <- list(param_stats=param_stats, precinct_stats=precinct_stats)
  return(res)}


gen_colors <- function(Candidates, MainCandidate){
  colors.v = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  colors.v <- sample(colors.v,  length(Candidates))
  colors.v[which(Candidates%in%MainCandidate)] <- "chartreuse4"
  colors.v <- c("grey", colors.v)
  return(colors.v)}

peak_search_quantile <- function(turnout) {
  turnout <- turnout[!is.na(turnout)]

  #3-component mixture model
  res_ml <- tryCatch(
    mixtools::normalmixEM(turnout,
                mu = quantile(turnout, c(0.25, 0.5, 0.75)),
                sigma = 0.01),
    error = function(e) e
  )

  #try 2-component mixture model
  if (inherits(res_ml, "error")) {
    res_ml <- tryCatch(
      mixtools::normalmixEM(turnout,
                  mu = quantile(turnout, c(0.25, 0.75)),
                  sigma = 0.01),
      error = function(e) e
    )
  }

  # Find component with minimum mean and return rounded value
  cind <- which.min(as.numeric(res_ml$mu))
  max.peak.turnout <- round(as.numeric(res_ml$mu)[cind], 2)

  return(max.peak.turnout)
}


ComputeSimulationNonparametric <- function(data, Candidates,
                                           CandidatesText, MainCandidate,
                                           FigureName, MaxThreshold, colors.v,
                                           precinct.level = precinctLevel, sims=sims,
                                           mode_search = mode_search,
                                           man_turnout = man_turnout,
                                           grid_type = "2D"){

  n_precincts <- nrow(data)

  bootDat <- function(dat, idx, level, grid_type) {
    dat <- dat[idx, ]

    if (grid_type == "1D") {
      estimate_fraud_res <- estimate_fraud(dat, Candidates, CandidatesText, MainCandidate,
                                           FigureName, MaxThreshold, colors.v,
                                           precinct.level = precinct.level,
                                           mode_search = mode_search,
                                           man_turnout = man_turnout)
    } else {
      estimate_fraud_res <- estimate_fraud_2d(dat, Candidates, CandidatesText, MainCandidate,
                                              FigureName, MaxThreshold, colors.v,
                                              precinct.level = precinct.level,
                                              mode_search = mode_search,
                                              man_turnout = man_turnout,
                                              graphsoff = TRUE)
    }

    main_params <- c(estimate_fraud_res$total_fraud,
                     estimate_fraud_res$real_turnout,
                     estimate_fraud_res$real_support)

    if (precinct.level && !is.null(estimate_fraud_res$precinct_fraud) &&
        !is.null(estimate_fraud_res$precinct_fraud$fraud.votes)) {

      precinct_results <- estimate_fraud_res$precinct_fraud$fraud.votes

      # ensure consistent length (pad with NA if necessary)
      if (length(precinct_results) < n_precincts) {
        precinct_results <- c(precinct_results, rep(NA, n_precincts - length(precinct_results)))
      } else if (length(precinct_results) > n_precincts) {
        precinct_results <- precinct_results[1:n_precincts]
      }
    } else {
      # if no precinct analysis, fill with zeros or NAs
      precinct_results <- rep(0, n_precincts)
    }

    # combine main parameters with precinct results
    res_parameters <- c(main_params, precinct_results)

    return(res_parameters)
  }

  tryCatch({
    results <- boot::boot(data=data, statistic=bootDat, R=sims,
                    level=precinct.level, grid_type=grid_type)

    boot.sd <- apply(results$t, 2, sd, na.rm=TRUE)
    boot.mean <- results$t0

    # extract parameter statistics (first 3 elements)
    param_stats <- c(boot.mean[2], boot.sd[2],  # real_turnout
                     boot.mean[3], boot.sd[3],   # real_support
                     boot.mean[1], boot.sd[1])   # total_fraud
    names(param_stats) <- c("real_turnoutM", "real_turnoutSD",
                            "real_supportM", "real_supportSD",
                            "total_fraudM", "total_fraudSD")

    # extract precinct statistics (elements 4 onwards)
    if (precinct.level && ncol(results$t) > 3) {
      precinct_indices <- 4:ncol(results$t)
      precinct_sd <- apply(results$t[, precinct_indices, drop=FALSE], 2,
                           function(x) sd(x, na.rm=TRUE))
      precinct_mean <- results$t0[precinct_indices]

      # ensure we have the right number of precincts
      if (length(precinct_mean) >= n_precincts) {
        precinct_stats <- data.frame(
          id = data$id[1:n_precincts],
          precinct_mean = precinct_mean[1:n_precincts],
          precinct_sd = precinct_sd[1:n_precincts]
        )
      } else {
        precinct_stats <- data.frame(
          id = data$id,
          precinct_mean = c(precinct_mean, rep(NA, n_precincts - length(precinct_mean))),
          precinct_sd = c(precinct_sd, rep(NA, n_precincts - length(precinct_sd)))
        )
      }
    } else {
      # no precinct analysis
      precinct_stats <- data.frame(
        id = data$id,
        precinct_mean = rep(0, n_precincts),
        precinct_sd = rep(0, n_precincts)
      )
    }

    res <- list(param_stats=param_stats, precinct_stats=precinct_stats)

  }, error = function(e) {
    warning("Bootstrap failed: ", e$message)

    # Return default results in case of error
    param_stats <- c(0, 0, 0, 0, 0, 0)
    names(param_stats) <- c("real_turnoutM", "real_turnoutSD",
                            "real_supportM", "real_supportSD",
                            "total_fraudM", "total_fraudSD")

    precinct_stats <- data.frame(
      id = data$id,
      precinct_mean = rep(0, n_precincts),
      precinct_sd = rep(0, n_precincts)
    )

    res <- list(param_stats=param_stats, precinct_stats=precinct_stats)
  })

  return(res)
}


peak_search_basic <- function(peaks_matrix,
                              mode_search, official.turnout = 0,
                              peaky) {

  peaks_matrix_less <- peaks_matrix[peaky < official.turnout,]
  if (nrow(peaks_matrix_less) == 0) {

    peaks_matrix_less <- peaks_matrix
  }

  # extract parameters with defaults
  pick_by <- ifelse(is.null(mode_search$pick_by), "height", mode_search$pick_by)

  peaks_raw <- tryCatch({
    pracma::findpeaks(peaks_matrix_less$max_height,
              npeaks = mode_search$npeaks,
              sortstr = mode_search$sortstr,
              minpeakdistance = mode_search$minpeakdistance)
    }, error = function(e) {
      peaks_raw = NULL
      })

  # handle degenerate case
  if (is.null(peaks_raw)) {
    # fallback: choose the max from c$max_height as a dummy peak
    max_index <- which.max(peaks_matrix_less$max_height)
    peak_height <- peaks_matrix_less$max_height[max_index]
    # construct dummy peak row [peak height, index, left base, right base]
    peaks_raw <- matrix(c(peak_height, max_index, max_index, max_index), nrow = 1, byrow = TRUE)
  }
  if (is.null(peaks_raw)) peaks_raw <- matrix(c(0, 0, 0, 0), nrow = 1, byrow = TRUE)
  #print(peaks_raw)

  # Fix overlapping base indices
  tryCatch({
    if (!(is.null(peaks_raw)) | nrow(peaks_raw) > 1) {
      used_indices <- c()
      for (i in 1:nrow(peaks_raw)) {
        for (col in 3:4) {  # columns 3 and 4 are Left_Base and Right_Base
          while (peaks_raw[i, col] %in% used_indices) {
            peaks_raw[i, col] <- peaks_raw[i, col] + 1
          }
          used_indices <- c(used_indices, peaks_raw[i, col])
        }
      }
    }
  }, error = function(e) {
    message("Error caught in peak adjustment block:")
    message(e$message)
  })

  # calculate peak areas
  peak_areas <- apply(peaks_raw, 1, function(row) {
    sum(peaks_matrix_less$max_height[row[3]:row[4]], na.rm = TRUE)
  })

  # create final peak matrix
  peaks_with_area <- cbind(peaks_raw, peak_areas)
  colnames(peaks_with_area) <- c("Peak_Height", "Peak_Index", "Left_Base", "Right_Base", "Peak_Area")

  # get turnout values for each peak
  peak_turnouts <- peaks_matrix_less$peaky[peaks_with_area[, "Peak_Index"]]

  # filter peaks below turnout threshold
  valid_mask <- peak_turnouts < official.turnout

  if (!any(valid_mask)) {
    warning("No clean peak found with turnout below threshold.")
    return(NA)
  }

  # Select best peak based on criteria
  valid_peaks <- peaks_with_area[valid_mask, , drop = FALSE]
  valid_turnouts <- peak_turnouts[valid_mask]

  if (pick_by == "cluster") {
    if (nrow(valid_peaks) <= 3) {
      # few peaks: use highest peak
      best_idx <- which.max(valid_peaks[, "Peak_Height"])
      peak_index <- valid_peaks[best_idx, "Peak_Index"]
      return(peaks_matrix_less$peaky[peak_index])
    } else {
      # many peaks: use clustering
      num_clusters <- 3
      cl.fit <- stats::kmeans(valid_turnouts, centers = num_clusters)

      # get peaks from median cluster
      cluster_vals <- cl.fit$cluster
      median_cluster <- median(cluster_vals)
      cluster_peaks <- valid_turnouts[cluster_vals == median_cluster]
      return(cluster_peaks[1])
    }
  } else {
    # original area/height selection
    best_idx <- if (pick_by == "area") {
      which.max(valid_peaks[, "Peak_Area"])
    } else {
      which.max(valid_peaks[, "Peak_Height"])
    }
    # return turnout value
    peak_index <- valid_peaks[best_idx, "Peak_Index"]
    return(peaks_matrix_less$peaky[peak_index])
  }
}


compute_precinctf <- function(data, clean.votes.vector,
         magnitude.election.fraud = NULL,
         postStratify = TRUE) {

  Votes <- data$Votes
  NValid <- data$NValid
  NVoters <- data$NVoters

  turn.scaler <- 0:100
  clean.votes.vector.c <- clean.votes.vector
  clean.dat <- data.frame(turn.scaler, clean.votes.vector.c, stringsAsFactors = FALSE)
  colnames(clean.dat) <- c("turn", "Freq")
  clean.dat[] <- lapply(clean.dat, function(x) as.numeric(as.character(x)))

  weight.incumbent <- Votes
  turn.scalerDat <- round(NValid / NVoters, 2) * 100
  real.dat <- data.frame(turn.scalerDat, Votes, weight.incumbent)
  colnames(real.dat) <- c("turn", "Freq", "w")

  real.datC <- real.dat[!(real.dat$turn == 0 | real.dat$Freq == 0), ]
  clean.dat <- clean.dat[!(clean.dat$turn == 0 | clean.dat$Freq == 0), ]
  diff.turn <- unique(c(setdiff(unique(clean.dat$turn), real.datC$turn),
                        setdiff(real.datC$turn, unique(clean.dat$turn))))

  real.datR <- real.datC[!real.datC$turn %in% diff.turn, ]
  clean.datR <- clean.dat[!clean.dat$turn %in% diff.turn, ]

  design <- survey::svydesign(ids = ~1, data = real.datR, weights = real.datR$w)
  clean.datR$turn <- factor(clean.datR$turn)
  real.datR$turn <- factor(real.datR$turn)

  if (postStratify) {
    ps.strat <- tryCatch(postStratify(design, ~turn, clean.datR), error = function(e) e)
    if (inherits(ps.strat, "error")) print("Error in postStratify")
    weights <- attr(ps.strat$postStrata[[1]], 'weights')
  } else {
    ps.strat <- rake(design, list(~turn), list(clean.datR))
    weights <- attr(ps.strat$postStrata[[1]][[1]], 'weights')
  }

  new.dat <- data.frame(real.datR, weights, stringsAsFactors = FALSE)
  new.dat[] <- lapply(new.dat, function(x) as.numeric(as.character(x)))
  new.dat$weights <- ifelse(new.dat$Freq < new.dat$weights, new.dat$Freq, new.dat$weights)

  clean.votes <- real.dat$Freq
  merg.dat <- data.frame(real.dat, clean.votes)
  merg.dat[rownames(merg.dat) %in% rownames(new.dat), "clean.votes"] <- new.dat$weights

  # compute initial fraud votes
  fraud.votes <- merg.dat$Freq - merg.dat$clean.votes

  # normalize total fraud to match magnitude.election.fraud (if provided)
  if (!is.null(magnitude.election.fraud)) {
    fraud.total <- sum(fraud.votes, na.rm=TRUE)
    if (fraud.total > 0) {
      scale.factor <- magnitude.election.fraud / fraud.total
      fraud.votes <- fraud.votes * scale.factor
      # adjust clean.votes accordingly
      merg.dat$clean.votes <- merg.dat$Freq - fraud.votes
    } else {
      warning("Total estimated fraud is zero or negative; skipping normalization.")
    }
  }

  res <- data.frame(id = data$id, clean.votes = merg.dat$clean.votes, fraud.votes = fraud.votes)
  return(res)
}


estimate_fraud_2d <- function(data, Candidates, CandidatesText, MainCandidate,
                              FigureName, MaxtThreshold, colors.v,
                              precinct.level = FALSE, mode_search = NULL,
                              man_turnout = NULL, official.turnout = NULL,
                              incumbent.support = NULL, graphsoff = FALSE) {

  # Load required libraries (with error handling)
  if(!require(ggplot2, quietly = TRUE)) stop("ggplot2 package required")
  if(!require(plotly, quietly = TRUE)) stop("plotly package required")
  if(!require(reshape2, quietly = TRUE)) stop("reshape2 package required")

  # Optional packages
  mixtools_available <- require(mixtools, quietly = TRUE)
  dbscan_available <- require(dbscan, quietly = TRUE)
  cluster_available <- require(cluster, quietly = TRUE)

  # Prepare data
  Votes <- data$Votes
  NValid <- data$NValid
  NVoters <- data$NVoters
  cand_vo <- data[, names(data) %in% Candidates]
  mcand_id <- which(Candidates == MainCandidate)
  turnout <- data$turnout * 100  # Turnout in percentage
  support <- data$incumbent * 100    # Support in percentage
  incumbent.support = sum(data$Votes, na.rm = TRUE) / sum(data$NValid, na.rm = TRUE) #attr(data, "incumbent.support")
  official.turnout = sum(data$NValid, na.rm = TRUE) / sum(data$NVoters, na.rm = TRUE) #attr(data, "official.turnout")

  # Set default values if not provided
  if(is.null(official.turnout)) official.turnout <- median(turnout, na.rm = TRUE)/100
  if(is.null(incumbent.support)) incumbent.support <- median(support, na.rm = TRUE)/100

  # Create 2D grid (1% bins)
  xbins <- seq(0, 100, by = 1)
  ybins <- seq(0, 100, by = 1)

  # Initialize vote count matrices for each candidate
  vote_grids <- lapply(1:length(Candidates), function(i) {
    matrix(0, nrow = length(xbins)-1, ncol = length(ybins)-1)
  })

  # Populate the grids
  for (i in 1:(length(xbins)-1)) {
    for (j in 1:(length(ybins)-1)) {
      in_bin <- turnout >= xbins[i] & turnout < xbins[i+1] &
        support >= ybins[j] & support < ybins[j+1]
      for (k in 1:length(Candidates)) {
        vote_grids[[k]][i,j] <- sum(cand_vo[in_bin, k], na.rm = TRUE)
      }
    }
  }

  for (i in 1:(length(xbins)-1)) {
    for (j in 1:(length(ybins)-1)) {
      in_xbin <- if (i == length(xbins) - 1) {
        turnout >= xbins[i] & turnout <= xbins[i+1]
      } else {
        turnout >= xbins[i] & turnout < xbins[i+1]
      }

      in_ybin <- if (j == length(ybins) - 1) {
        support >= ybins[j] & support <= ybins[j+1]
      } else {
        support >= ybins[j] & support < ybins[j+1]
      }

      in_bin <- in_xbin & in_ybin
      for (k in 1:length(Candidates)) {
        vote_grids[[k]][i, j] <- sum(cand_vo[in_bin, k], na.rm = TRUE)
      }
    }
  }

  # Get main candidate grid
  main_cand_grid <- vote_grids[[mcand_id]]
  rownames(main_cand_grid) <- xbins[-length(xbins)]
  colnames(main_cand_grid) <- ybins[-length(ybins)]

  # Melt for ggplot
  melted_grid <- reshape2::melt(main_cand_grid)
  names(melted_grid) <- c("Turnout", "Support", "Votes")


  # IMPROVED CLUSTERING USING GRID DATA
  cluster_results <- find_clusters_from_grid(main_cand_grid, target_clusters = 3)

  # Extract cluster information
  clean_turnout <- cluster_results$clean_cluster_turnout
  clean_support <- cluster_results$clean_cluster_support
  cluster_centers <- cluster_results$cluster_centers
  cluster_assignments <- cluster_results$cluster_assignments

  # Create cluster visualization data
  cluster_df <- cluster_results$cluster_df

  if (!graphsoff) {
    cluster_plot <- ggplot2::ggplot(cluster_df, aes(x = turnout, y = support,
                                           color = cluster, size = votes)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_point(data = data.frame(turnout = cluster_centers[,1],
                                   support = cluster_centers[,2],
                                   cluster = factor(1:nrow(cluster_centers))),
                 aes(x = turnout, y = support, color = cluster),
                 size = 5, shape = 4, stroke = 2) +
      ggplot2::geom_vline(xintercept = official.turnout*100, linetype = "dashed", color = "blue") +
      ggplot2::geom_vline(xintercept = clean_turnout, linetype = "dashed", color = "red") +
      ggplot2::geom_hline(yintercept = clean_support, linetype = "dashed", color = "red") +
      scale_size_continuous(range = c(0.5, 3)) +
      ggtitle(paste("Cluster Analysis (", FigureName,")", sep="")) +
      labs(subtitle = paste("Clean cluster turnout:", round(clean_turnout, 1), "%")) +
      theme_minimal()
  } else {
    cluster_plot <- NULL
  }

  # Melt for ggplot
  melted_grid <- reshape2::melt(main_cand_grid)
  names(melted_grid) <- c("Turnout", "Support", "Votes")

  # Create complete grid for 3D visualization
  complete_grid <- expand.grid(
    Turnout = seq(0, 99, by = 1),
    Support = seq(0, 99, by = 1)
  )
  merged_grid <- merge(complete_grid, melted_grid, by = c("Turnout", "Support"), all.x = TRUE)
  merged_grid$Votes[is.na(merged_grid$Votes)] <- 0

  # Create z_matrix for precinct-level analysis
  z_matrix <- reshape2::acast(merged_grid, Turnout ~ Support, value.var = "Votes")

  if (!graphsoff) {
    heatmap_plot <- ggplot2::ggplot(merged_grid, aes(x = Turnout, y = Support, fill = Votes)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      geom_vline(xintercept = official.turnout*100, linetype = "dashed", color = "blue") +
      geom_vline(xintercept = clean_turnout, linetype = "dashed", color = "red") +
      geom_hline(yintercept = clean_support, linetype = "dashed", color = "red") +
      ggtitle(paste("2D Vote Distribution for Incumbent (", FigureName,")", sep="")) +
      theme_minimal()
  } else {
    heatmap_plot <- NULL
  }

  if (!graphsoff) {
    surface_plot <- tryCatch({
      plotly::plot_ly(
        x = as.numeric(rownames(z_matrix)),
        y = as.numeric(colnames(z_matrix)),
        z = z_matrix,
        type = "surface",
        colorscale = list(c(0, 1), c("white", "red")),
        hoverinfo = "x+y+z",
        showscale = TRUE
      ) %>%
        layout(
          scene = list(
            xaxis = list(title = "Turnout (%)"),
            yaxis = list(title = "Support (%)"),
            zaxis = list(title = "Votes"),
            camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
          ),
          title = paste("2D Vote Distribution for Incumbent (", FigureName,")", sep="")
        )
    }, error = function(e) NULL)
  } else {
    surface_plot <- NULL
  }

  # Calculate fraud estimates using clean cluster parameters
  high_turnout <- turnout > clean_turnout
  if(sum(high_turnout, na.rm = TRUE) > 0) {
    expected_votes <- sum(NValid[high_turnout], na.rm = TRUE) * (clean_support/100)
    actual_votes <- sum(Votes[high_turnout], na.rm = TRUE)
    fraud_estimate <- max(0, actual_votes - expected_votes, na.rm = TRUE)
    prop_fraud <- ifelse(sum(Votes, na.rm = TRUE) > 0,
                         fraud_estimate / sum(Votes, na.rm = TRUE), 0)
  } else {
    fraud_estimate <- 0
    prop_fraud <- 0
  }
  ##########################
  #elipseMCD <- computeMCDforKmeans(Votes, data$Valid, data$NVoters)
  #draw_elipse <- data.frame(elipseMCD$elipsedat)
  #max.peak.turnout <- elipseMCD$clean_mean[1]/100
  #total votes
  # total_votes_grid contains the sum across all candidates
  #rownames(total_votes_grid) <- xbins[-length(xbins)]
  #colnames(total_votes_grid) <- ybins[-length(ybins)]

  #main candidate grid
  main_cand_grid <- vote_grids[[mcand_id]]
  melted_grid_ur <- reshape2::melt(main_cand_grid)
  names(melted_grid_ur) <- c("Turnout", "Support", "Votes")

  #valid votes
  total_votes_grid <- Reduce(`+`, vote_grids)
  melted_valid_grid <- reshape2::melt(total_votes_grid)
  names(melted_valid_grid) <- c("Turnout", "Support", "Votes")

  #Valid minus main
  other_candidate_votes <- total_votes_grid - main_cand_grid
  clean.votes.v <- melted_grid_ur$Votes[melted_grid_ur$Turnout==clean_turnout]
  clean.valid.votes.v <- melted_valid_grid$Votes[melted_valid_grid$Turnout==clean_turnout]
  clean.all.butMainCandidate.v <- clean.valid.votes.v  - clean.votes.v

  clean.ur.prop.v <- clean.votes.v / clean.valid.votes.v
  clean.all.butMainCandidate.prop.v <- 1 - clean.ur.prop.v
  inflation.factor.v <- clean.votes.v / (clean.all.butMainCandidate.v)

  inflation.factor.v[is.nan(inflation.factor.v)] <- 1

  #148*(inflation.factor.v[inflation.factor.v>=100]-100)/100<-100

  # Replace NaN with 0 (or another value if preferred)
  #inflation.factor.v[is.nan(inflation.factor.v)] <- 1

  # Scale only if max exceeds 100
  max_val <- max(inflation.factor.v, na.rm = TRUE)
  if (max_val > 100) {
    scale_factor <- 100 / max_val
    inflation.factor.v <- inflation.factor.v * scale_factor
  }

  #inflation.factor.v[is.nan(inflation.factor.v)] <- 1
  inflation.factor.v[inflation.factor.v<1] <- 1
  #inflation.factor.v[inflation.factor.v>=100] <- 100

  clean.votes.mat <- other_candidate_votes * inflation.factor.v

  clean.mat <- reshape2::melt(clean.votes.mat)
  names(clean.mat) <- c("Turnout", "Support", "Votes")

  clean.votes.total <- sum(clean.votes.mat, na.rm=TRUE)

  magnitude.election.fraud <- sum(main_cand_grid, na.rm=TRUE) - clean.votes.total
  magnitude.election.fraud_mat <- main_cand_grid - clean.votes.mat

  ballot_stuffing <- round((official.turnout*100 - clean_turnout)/100 * sum(NVoters, na.rm = TRUE), 0)
  if(ballot_stuffing > magnitude.election.fraud) ballot_stuffing <- round(magnitude.election.fraud, 0)
  ballot_switching <- round((magnitude.election.fraud - ballot_stuffing) / 2, 0)

  # Calculate ballot stuffing and switching
  ballot_stuffing <- round((official.turnout*100 - clean_turnout)/100 * sum(NVoters, na.rm = TRUE), 0)
  if(ballot_stuffing > fraud_estimate) ballot_stuffing <- round(fraud_estimate, 0)
  ballot_switching <- round((fraud_estimate - ballot_stuffing) / 2, 0)

  # Precinct-level fraud calculation using 2D approach
  if(precinct.level) {
    precinct.fraud <- tryCatch({
      compute_precinctf_2d(
        data = data,
        clean.votes.vector = clean.mat,
        magnitude.election.fraud = magnitude.election.fraud,
        postStratify = TRUE
      )
    }, error = function(e) {
      warning("Precinct-level fraud calculation failed: ", e$message)
      NULL
    })
  } else {
    precinct.fraud <- NULL
  }

  # 3D Plot for Clean Votes Matrix
  if (!graphsoff) {
    clean_votes_surface <- tryCatch({
      # Prepare clean votes matrix
      clean_melted <- reshape2::melt(clean.votes.mat)
      names(clean_melted) <- c("Turnout", "Support", "Votes")

      # Create complete grid
      clean_complete_grid <- expand.grid(
        Turnout = seq(0, 99, by = 1),
        Support = seq(0, 99, by = 1)
      )
      clean_merged <- merge(clean_complete_grid, clean_melted,
                            by = c("Turnout", "Support"), all.x = TRUE)
      clean_merged$Votes[is.na(clean_merged$Votes)] <- 0

      # Create z-matrix
      clean_z_matrix <- reshape2::acast(clean_merged, Turnout ~ Support, value.var = "Votes")

      # Generate plot
      plotly::plot_ly(
        x = as.numeric(rownames(clean_z_matrix)),
        y = as.numeric(colnames(clean_z_matrix)),
        z = clean_z_matrix,
        type = "surface",
        colorscale = list(c(0, 1), c("white", "green")),
        hoverinfo = "x+y+z",
        showscale = TRUE
      ) %>%
        layout(
          scene = list(
            xaxis = list(title = "Turnout (%)"),
            yaxis = list(title = "Support (%)"),
            zaxis = list(title = "Clean Votes"),
            camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
          ),
          title = paste("Estimated Clean Votes for", MainCandidate, "(", FigureName, ")")
        )
    }, error = function(e) NULL)
  } else {
    clean_votes_surface <- NULL
  }

  # 3D Plot for Fraud Magnitude Matrix
  if (!graphsoff) {
    fraud_surface <- tryCatch({
      # Prepare fraud magnitude matrix
      fraud_melted <- reshape2::melt(magnitude.election.fraud_mat)
      names(fraud_melted) <- c("Turnout", "Support", "Votes")

      # Create complete grid
      fraud_complete_grid <- expand.grid(
        Turnout = seq(0, 99, by = 1),
        Support = seq(0, 99, by = 1)
      )
      fraud_merged <- merge(fraud_complete_grid, fraud_melted,
                            by = c("Turnout", "Support"), all.x = TRUE)
      fraud_merged$Votes[is.na(fraud_merged$Votes)] <- 0

      # Create z-matrix
      fraud_z_matrix <- reshape2::acast(fraud_merged, Turnout ~ Support, value.var = "Votes")

      # Generate plot (using divergent colors)
      plotly::plot_ly(
        x = as.numeric(rownames(fraud_z_matrix)),
        y = as.numeric(colnames(fraud_z_matrix)),
        z = fraud_z_matrix,
        type = "surface",
        colorscale = list(
          list(0, "blue"),
          list(0.5, "white"),
          list(1, "red")
        ),
        hoverinfo = "x+y+z",
        showscale = TRUE
      ) %>%
        layout(
          scene = list(
            xaxis = list(title = "Turnout (%)"),
            yaxis = list(title = "Support (%)"),
            zaxis = list(title = "Fraud Magnitude"),
            camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
          ),
          title = paste("Estimated Fraud Magnitude for", MainCandidate, "(", FigureName, ")")
        )
    }, error = function(e) NULL)
  } else {
    fraud_surface <- NULL
  }

  # Create drawfig list for plots (matching expected format)
  #drawfig <- list(
  #  cluster_plot = cluster_plot,
  #  heatmap_plot = heatmap_plot,
  #  surface_plot = surface_plot
  #)

  drawfig <- list(
    cluster_plot = cluster_plot,
    heatmap_plot = heatmap_plot,
    surface_plot = surface_plot,
    clean_votes_surface = clean_votes_surface,
    fraud_surface = fraud_surface
  )


  # Set variable names to match the target format
  o.turnout <- official.turnout * 100
  r.turnout <- clean_turnout
  realMainCandidateSupport <- clean_support
  fraud_measure <- round(fraud_estimate, 0)
  prop.fraud <- prop_fraud

  # Create results list matching the exact target format
  results <- list(official_turnout = o.turnout, real_turnout = r.turnout,
                  official_support = round(incumbent.support*100, 1), real_support = round(realMainCandidateSupport, 1),
                  ballot_stuffing = ballot_stuffing, ballot_switching = ballot_switching,
                  total_fraud = fraud_measure, prop_fraud = prop.fraud, drawfig = drawfig, precinct_fraud = precinct.fraud)

  return(results)
}


find_clusters_from_grid <- function(grid_matrix, target_clusters = 3, min_votes = 5) {

  # Check package availability
  mixtools_available <- requireNamespace("mixtools", quietly = TRUE)
  dbscan_available <- requireNamespace("dbscan", quietly = TRUE)
  cluster_available <- requireNamespace("cluster", quietly = TRUE)

  # helper function for multivariate normal density (for GMM)
  dmvnorm_simple <- function(x, mu, sigma) {
    tryCatch({
      k <- length(mu)
      det_sigma <- det(sigma)
      if(det_sigma <= 0) return(1e-10)

      inv_sigma <- solve(sigma)
      diff <- x - mu
      exp_part <- -0.5 * t(diff) %*% inv_sigma %*% diff

      (2 * pi)^(-k/2) * det_sigma^(-0.5) * exp(exp_part)
    }, error = function(e) {
      1e-10  # Return small positive value on error
    })
  }

  # convert grid to data points weighted by vote counts
  grid_data <- expand.grid(
    turnout = as.numeric(rownames(grid_matrix)),
    support = as.numeric(colnames(grid_matrix))
  )
  grid_data$votes <- as.vector(grid_matrix)

  # filter out cells with very few votes to reduce noise
  significant_cells <- grid_data[grid_data$votes >= min_votes, ]

  if(nrow(significant_cells) < target_clusters) {
    # fallback to all non-zero cells
    significant_cells <- grid_data[grid_data$votes > 0, ]
  }

  if(nrow(significant_cells) < target_clusters) {
    # ultimate fallback - use raw grid positions
    significant_cells <- grid_data
    significant_cells$votes[significant_cells$votes == 0] <- 1
  }

  # try multiple clustering methods
  clustering_results <- list()

  # method 1: Weighted K-means (votes as weights)
  tryCatch({
    # create weighted dataset by replicating points proportional to votes
    weighted_data <- significant_cells[rep(1:nrow(significant_cells),
                                           pmax(1, round(significant_cells$votes/max(significant_cells$votes) * 10))), ]

    kmeans_result <- stats::kmeans(weighted_data[, c("turnout", "support")],
                            centers = target_clusters, nstart = 25)

    # assign clusters back to original significant cells
    cluster_assignments <- rep(NA, nrow(significant_cells))
    for(i in 1:nrow(significant_cells)) {
      distances <- sqrt(rowSums((t(kmeans_result$centers) - c(significant_cells$turnout[i], significant_cells$support[i]))^2))
      cluster_assignments[i] <- which.min(distances)
    }

    clustering_results[["weighted_kmeans"]] <- list(
      clusters = cluster_assignments,
      centers = kmeans_result$centers,
      method = "Weighted K-means"
    )
  }, error = function(e) {
    message("Weighted K-means failed: ", e$message)
  })

  # method 2: Gaussian Mixture Model with vote weighting
  tryCatch({
    if(mixtools_available && nrow(significant_cells) >= target_clusters * 3) {
      # create weighted dataset by replicating points based on vote counts
      vote_weights <- pmax(1, round(significant_cells$votes / max(significant_cells$votes) * 3))
      weighted_coords <- significant_cells[rep(1:nrow(significant_cells), vote_weights),
                                           c("turnout", "support")]

      # only proceed if we have enough points
      if(nrow(weighted_coords) >= target_clusters * 5) {
        # fit GMM with very conservative settings
        gmm_result <- mixtools::normalmixEM(as.matrix(weighted_coords), k = target_clusters,
                                  verb = FALSE, maxit = 25, epsilon = 1e-2)

        # check if GMM converged properly
        if(!is.null(gmm_result$mu) && length(gmm_result$mu) == target_clusters) {
          # extract centers - mu is a list of vectors
          centers <- matrix(NA, nrow = target_clusters, ncol = 2)
          for(i in 1:target_clusters) {
            centers[i, ] <- gmm_result$mu[[i]]
          }

          # assign original points to closest center
          original_coords <- as.matrix(significant_cells[, c("turnout", "support")])
          cluster_assignments <- apply(original_coords, 1, function(point) {
            distances <- sqrt(rowSums((t(centers) - point)^2))
            which.min(distances)
          })

          colnames(centers) <- c("turnout", "support")

          clustering_results[["gmm"]] <- list(
            clusters = cluster_assignments,
            centers = centers,
            method = "Gaussian Mixture Model"
          )
        }
      }
    }
  }, error = function(e) {
    message("GMM failed: ", e$message)
  })

  # method 3: DBSCAN for density-based clustering (if available)
  tryCatch({
    if(dbscan_available) {
      # use votes as density weights by replicating points
      vote_weights <- pmax(1, round(significant_cells$votes / max(significant_cells$votes) * 5))
      weighted_coords <- significant_cells[rep(1:nrow(significant_cells), vote_weights),
                                           c("turnout", "support")]

      # estimate eps parameter
      if(nrow(weighted_coords) > 10) {
        knn_dist <- kNNdist(as.matrix(weighted_coords), k = 4)
        eps <- quantile(knn_dist, 0.8)

        dbscan_result <- dbscan::dbscan(as.matrix(weighted_coords), eps = eps,
                                minPts = max(3, nrow(weighted_coords)/30))

        # map back to original points
        cluster_assignments <- rep(NA, nrow(significant_cells))
        weighted_clusters <- dbscan_result$cluster

        idx <- 1
        for(i in 1:nrow(significant_cells)) {
          # take the mode of cluster assignments for this original point
          point_clusters <- weighted_clusters[idx:(idx + vote_weights[i] - 1)]
          cluster_assignments[i] <- as.numeric(names(sort(table(point_clusters), decreasing = TRUE)[1]))
          idx <- idx + vote_weights[i]
        }

        if(max(cluster_assignments, na.rm = TRUE) >= 2) {
          # calculate cluster centers
          centers <- matrix(NA, nrow = max(cluster_assignments, na.rm = TRUE), ncol = 2)
          for(i in 1:max(cluster_assignments, na.rm = TRUE)) {
            if(sum(cluster_assignments == i, na.rm = TRUE) > 0) {
              cluster_points <- as.matrix(significant_cells[cluster_assignments == i, c("turnout", "support")])
              cluster_weights <- significant_cells$votes[cluster_assignments == i]
              centers[i, ] <- colSums(cluster_points * cluster_weights) / sum(cluster_weights)
            }
          }
          colnames(centers) <- c("turnout", "support")

          clustering_results[["dbscan"]] <- list(
            clusters = cluster_assignments,
            centers = centers,
            method = "DBSCAN"
          )
        }
      }
    }
  }, error = function(e) {
    message("DBSCAN failed: ", e$message)
  })

  # method 4: hierarchical clustering with vote weighting
  tryCatch({
    coords <- as.matrix(significant_cells[, c("turnout", "support")])

    # create distance matrix
    dist_matrix <- dist(coords)
    hclust_result <- stats::hclust(dist_matrix, method = "ward.D2")
    cluster_assignments <- cutree(hclust_result, k = target_clusters)

    # calculate weighted centers
    centers <- matrix(NA, nrow = target_clusters, ncol = 2)
    for(i in 1:target_clusters) {
      cluster_mask <- cluster_assignments == i
      if(sum(cluster_mask) > 0) {
        cluster_points <- coords[cluster_mask, , drop = FALSE]
        cluster_weights <- significant_cells$votes[cluster_mask]
        if(nrow(cluster_points) == 1) {
          centers[i, ] <- cluster_points[1, ]
        } else {
          centers[i, ] <- colSums(cluster_points * cluster_weights) / sum(cluster_weights)
        }
      }
    }
    colnames(centers) <- c("turnout", "support")

    clustering_results[["hclust"]] <- list(
      clusters = cluster_assignments,
      centers = centers,
      method = "Hierarchical Clustering"
    )
  }, error = function(e) {
    message("Hierarchical clustering failed: ", e$message)
  })

  # select best clustering result
  if(length(clustering_results) == 0) {
    # Ultimate fallback - simple k-means
    kmeans_result <- stats::kmeans(significant_cells[, c("turnout", "support")],
                            centers = target_clusters, nstart = 10)
    best_result <- list(
      clusters = kmeans_result$cluster,
      centers = kmeans_result$centers,
      method = "Simple K-means (fallback)"
    )
  } else {
    # choose the method that best separates turnout ranges
    best_score <- -Inf
    best_result <- NULL

    for(method_name in names(clustering_results)) {
      result <- clustering_results[[method_name]]

      # calculate separation score based on turnout ranges
      turnout_ranges <- sapply(1:max(result$clusters, na.rm = TRUE), function(i) {
        cluster_turnouts <- significant_cells$turnout[result$clusters == i]
        if(length(cluster_turnouts) > 0) {
          max(cluster_turnouts) - min(cluster_turnouts)
        } else {
          Inf
        }
      })

      # score: prefer methods with distinct turnout clusters
      score <- -mean(turnout_ranges, na.rm = TRUE)

      if(score > best_score) {
        best_score <- score
        best_result <- result
      }
    }
  }

  # find the clean cluster (with lowest mean turnout)
  cluster_turnout_means <- sapply(1:max(best_result$clusters, na.rm = TRUE), function(i) {
    cluster_turnouts <- significant_cells$turnout[best_result$clusters == i]
    cluster_votes <- significant_cells$votes[best_result$clusters == i]
    if(length(cluster_turnouts) > 0) {
      # Weighted mean
      sum(cluster_turnouts * cluster_votes) / sum(cluster_votes)
    } else {
      Inf
    }
  })

  clean_cluster_id <- which.min(cluster_turnout_means)

  # get the highest turnout in the clean cluster
  clean_cluster_points <- significant_cells[best_result$clusters == clean_cluster_id, ]
  clean_turnout <- max(clean_cluster_points$turnout)

  # get weighted support for clean cluster
  clean_support <- sum(clean_cluster_points$support * clean_cluster_points$votes) /
    sum(clean_cluster_points$votes)

  # create visualization dataframe
  cluster_df <- data.frame(
    turnout = significant_cells$turnout,
    support = significant_cells$support,
    votes = significant_cells$votes,
    cluster = factor(best_result$clusters)
  )

  return(list(
    clean_cluster_turnout = clean_turnout,
    clean_cluster_support = clean_support,
    cluster_centers = best_result$centers,
    cluster_assignments = best_result$clusters,
    cluster_df = cluster_df,
    method_used = best_result$method,
    all_methods_tried = names(clustering_results)
  ))
}

compute_precinctf_2d <- function(data, clean.votes.vector,
                                 magnitude.election.fraud = NULL,
                                 postStratify = TRUE) {

  Votes <- data$Votes
  NValid <- data$NValid
  NVoters <- data$NVoters
  id <- data$id

  # prepare clean data - ensure proper aggregation
  clean.dat <- as.data.frame(clean.votes.vector)
  colnames(clean.dat) <- c("turnout", "support", "Freq")
  clean.dat$turnout <- round(clean.dat$turnout)
  clean.dat$support <- round(clean.dat$support)

  # aggregation - sum frequencies by unique turnout-support pairs
  clean.dat <- aggregate(Freq ~ turnout + support, data = clean.dat, sum)

  turnout.obs <- round(NValid / NVoters * 100)
  support.obs <- round(Votes / NValid * 100)
  real.dat <- data.frame(
    id = id,
    turnout = turnout.obs,
    support = support.obs,
    Freq = Votes,
    w = Votes  # weights
  )

  # remove zero/NA frequencies
  clean.dat <- clean.dat[clean.dat$Freq > 0 & !is.na(clean.dat$Freq), ]
  real.dat <- real.dat[real.dat$Freq > 0 & !is.na(real.dat$Freq), ]

  # find common (turnout,support) pairs - using UNIQUE combinations
  common_pairs <- merge(
    unique(clean.dat[, c("turnout", "support")]),
    unique(real.dat[, c("turnout", "support")]),
    by = c("turnout", "support")
  )

  if (nrow(common_pairs) == 0) {
    stop("No overlapping (turnout, support) levels between observed and expected data")
  }

  # filter to common levels - PROPERLY JOIN CLEAN DATA
  clean.datR <- merge(common_pairs, clean.dat, by = c("turnout", "support"))
  real.datR <- merge(common_pairs, real.dat, by = c("turnout", "support"))

  # verify clean data has unique combinations
  if (any(duplicated(clean.datR[, c("turnout", "support")]))) {
    stop("Duplicate (turnout,support) pairs in clean data after merging")
  }

  # convert to factors
  clean.datR$turnout <- factor(clean.datR$turnout)
  clean.datR$support <- factor(clean.datR$support)
  real.datR$turnout <- factor(real.datR$turnout)
  real.datR$support <- factor(real.datR$support)

  # create survey design
  design <- survey::svydesign(ids = ~id, data = real.datR, weights = real.datR$w)

  # calculate weights
  if (postStratify) {
    ps.strat <- tryCatch(
      postStratify(design, ~turnout + support, clean.datR),
      error = function(e) {
        warning("Post-stratification failed, using raking: ", e$message)
        rake(design, list(~turnout, ~support), list(clean.datR, clean.datR))
      }
    )
  } else {
    ps.strat <- rake(design, list(~turnout, ~support), list(clean.datR, clean.datR))
  }

  # extract weights
  weights <- weights(ps.strat)

  result <- data.frame(
    id = real.datR$id,
    turnout = real.datR$turnout,
    support = real.datR$support,
    observed_votes = real.datR$Freq,
    clean_votes = pmin(weights, real.datR$Freq), # Can't exceed observed
    fraud.votes = pmax(real.datR$Freq - weights, 0) # Can't be negative
  )

  # if magnitude specified, scale fraud votes accordingly
  if (!is.null(magnitude.election.fraud)) {
    current_fraud <- sum(result$fraud.votes)
    if (current_fraud > 0) {
      scale_factor <- magnitude.election.fraud / current_fraud
      result$fraud.votes <- result$fraud.votes * scale_factor
      result$clean_votes <- result$observed_votes - result$fraud.votes
    } else {
      warning("No fraud detected - cannot scale to requested magnitude")
    }
  }

  final_result <- merge(
    data.frame(id = id),
    result,
    by = "id",
    all.x = TRUE
  )

  # fill NA values (unmatched precincts)
  final_result$clean_votes[is.na(final_result$clean_votes)] <-
    final_result$observed_votes[is.na(final_result$clean_votes)]
  final_result$fraud.votes[is.na(final_result$fraud.votes)] <- 0

  return(final_result)
}
