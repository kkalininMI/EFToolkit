#' @title Nonparametric Method
#' @description This function implements the revised version of Shpilkin's method (NB! It doesn't replicate original method fully).
#' @usage NonparamElectionForensics(data, Candidates, CandidatesText,
#'                                MainCandidate, TotalReg, TotalVotes,
#'                                Level, MaxThreshold,
#'                                FigureName, setcolors,
#'                                precinctLevel, computeSD, sims,
#'                                mode_search = list(npeaks, sortstr,
#'                                              minpeakdistance, pick_by),
#'                                man_turnout,
#'                                grid_type)
#' @param data A data.frame containing precinct-level election results
#' @param Candidates Character vector of candidate variable names in the dataset
#' @param CandidatesText Character vector of text labels for candidates (optional)
#' @param MainCandidate Character name of the main candidate (incumbent) to analyze
#' @param TotalReg Character name of the variable for total registered voters
#' @param TotalVotes Character name of the variable for total valid votes (optional)
#' @param Level Character analysis level ("National" or regional variable name)
#' @param MaxThreshold Numeric maximum turnout threshold for analysis (default: 0.8)
#' @param FigureName Character name for output figures
#' @param setcolors Character vector of custom colors for candidates (optional)
#' @param precinctLevel Logical indicating whether to compute precinct-level fraud estimates (default: TRUE)
#' @param computeSD Character method for uncertainty estimation ("parametric", "nonparametric", or NULL)
#' @param sims Integer number of bootstrap simulations (default: 10)
#' @param mode_search List of parameters for peak detection algorithm (default: list(npeaks = 5, sortstr = TRUE, minpeakdistance = 1, pick_by = "height"))
#' @param man_turnout Numeric manual turnout value to override automatic detection (optional)
#' @param grid_type Character analysis grid type ("1D" for turnout only or "2D" for turnout-support) (default: "1D")
#'
#' @return A list containing:
#' \itemize{
#'   \item list_graphs - Collection of generated plots
#'   \item base_stats - Basic fraud statistics for the whole dataset
#'   \item sim_all_stats - Simulation statistics for the whole dataset (if computeSD specified)
#'   \item sim_hetero_stats_base - Base statistics for regional analyses
#'   \item sim_hetero_stats_sims - Simulation statistics for regional analyses
#'   \item fraud_precinct_data - Precinct-level fraud estimates with uncertainty
#'   \item data - Original input data with computed variables
#'   \item Level - Analysis level used
#'   \item creationdate - Timestamp of creation
#' }
#'
#' @export
#' @import ggplot2
#' @import robustbase
#' @importFrom stats aggregate na.omit sd
#' @importFrom utils head tail
#' @return list containing results of analysis
#' \itemize{
#'   \item list_graphs - list of graphs
#'   \item stats_summary -  table with results for analysis of the whole dataset (for Level!=NULL: data are also
#'                          summed across the units, aggregation error is computed)
#'   \item Level -  external parameter
#'   \item creationdate - date/time of analysis
#'   list_graphs = list_graphs, stats_summary = stats_table, stats_level=stats_level
#' }
#' @examples
#' library(EFToolkit)
#'
#' dat<-read.csv(system.file("extdata/ruspres2000.csv", package="EFToolkit"))
#'
#'# Estimate fraud using the whole dataset with 1D grid analysis
#'res1 <- NonparamElectionForensics(
#'        dat, Candidates = paste("P", 1:12, sep=""),
#'        CandidatesText=c("Stanislav Govorukhin", "Umar Dzhabrailov",
#'        "Vladimir Zhirinovsky", "Gennady Zuganov",
#'        "Ella Pamfilova", "Alexei Podberezkin",
#'        "Vladimir Putin", "Yuri Skuratov",
#'        "Konstantin Titov", "Aman Tuleev",
#'        "Grigorii Yavlinsky", "Against All"),
#'        MainCandidate = "P7",
#'        TotalReg = "NVoters",
#'        TotalVotes = "NValid",
#'        Level = "National",
#'        MaxThreshold = 0.8,
#'        mode_search = list(npeaks = 5, sortstr = TRUE,
#'                minpeakdistance = 1, pick_by = "height"),
#'        FigureName = "Russian Presidential Elections, 2000",
#'        setcolors = c("royalblue2", "springgreen1","blue",
#'                      "red", "green","brown2",
#'                      "darkgreen", "yellow", "lawngreen",
#'                      "purple","chartreuse1", "orange"),
#'        precinctLevel = TRUE, computeSD="nonparametric",
#'        sims=2, grid_type = "1D")
#'
#'# Sum precinct-level fraud estimates
#'sum(res1$fraud_precinct_data$base.fraud.votes, na.rm=TRUE)
#'# 2685246
#'
#'# Sum precinct-level fraud (only statistically significant values)
#'sum(res1$fraud_precinct_data$sim.precinct_mean[res1$fraud_precinct_data$sim.sig_all==TRUE], na.rm=TRUE)
#'# 811500.5
#'
#'# Estimate fraud using region-level analysis with 1D grid
#'res2 <- NonparamElectionForensics(
#'        dat, Candidates = paste("P", 1:12, sep=""),
#'        CandidatesText=c("Stanislav Govorukhin", "Umar Dzhabrailov",
#'        "Vladimir Zhirinovsky", "Gennady Zuganov",
#'        "Ella Pamfilova", "Alexei Podberezkin",
#'        "Vladimir Putin", "Yuri Skuratov",
#'        "Konstantin Titov", "Aman Tuleev",
#'        "Grigorii Yavlinsky", "Against All"),
#'        MainCandidate = "P7",
#'        TotalReg = "NVoters",
#'        TotalVotes = "NValid",
#'        Level = "regname",
#'        MaxThreshold = 0.8,
#'        mode_search = list(npeaks = 5, sortstr = TRUE,
#'                minpeakdistance = 1, pick_by = "height"),
#'        FigureName = "Russian Presidential Elections, 2000",
#'        setcolors = c("royalblue2", "springgreen1","blue",
#'                      "red", "green","brown2",
#'                      "darkgreen", "yellow", "lawngreen",
#'                      "purple","chartreuse1", "orange"),
#'        precinctLevel = TRUE, computeSD="nonparametric",
#'        sims=2, grid_type = "1D")
#'
#'# Sum precinct-level fraud estimates
#'sum(res2$fraud_precinct_data$precinct_mean_hetero, na.rm=TRUE)
#'# 1693742
#'
#'# Sum precinct-level fraud (only statistically significant values)
#'sum(res2$fraud_precinct_data$precinct_mean_hetero[res2$fraud_precinct_data$sim.sig_all==TRUE], na.rm=TRUE)
#'# 1125151
#'
#'# Estimate fraud using the whole dataset with 2D grid analysis
#'res3 <- NonparamElectionForensics(
#'        dat, Candidates = paste("P", 1:12, sep=""),
#'        CandidatesText=c("Stanislav Govorukhin", "Umar Dzhabrailov",
#'        "Vladimir Zhirinovsky", "Gennady Zuganov",
#'        "Ella Pamfilova", "Alexei Podberezkin",
#'        "Vladimir Putin", "Yuri Skuratov",
#'        "Konstantin Titov", "Aman Tuleev",
#'        "Grigorii Yavlinsky", "Against All"),
#'        MainCandidate = "P7",
#'        TotalReg = "NVoters",
#'        TotalVotes = "NValid",
#'        Level = "regname",
#'        MaxThreshold = 0.8,
#'        mode_search = list(npeaks = 5, sortstr = TRUE,
#'                minpeakdistance = 1, pick_by = "height"),
#'        FigureName = "Russian Presidential Elections, 2000",
#'        setcolors = c("royalblue2", "springgreen1","blue",
#'                      "red", "green","brown2",
#'                      "darkgreen", "yellow", "lawngreen",
#'                      "purple","chartreuse1", "orange"),
#'        precinctLevel = TRUE, computeSD="parametric",
#'        sims=2, grid_type = "2D")
#'
#'# Sum precinct-level fraud estimates
#'sum(res3$fraud_precinct_data$sim.precinct_mean, na.rm=TRUE)
#'# -15724742
#'
#'# Sum precinct-level fraud (only statistically significant values)
#' sum(
#' res3$fraud_precinct_data$sim.precinct_mean[
#' res3$fraud_precinct_data$sim.sig_all == TRUE
#' ],
#' na.rm = TRUE)
#'# 0
#'
# Display regional statistics table
#'print(round(res3$sim_hetero_stats_base, 3))
#'
#'# Estimate fraud using region-level analysis with 2D grid
#'#NB! Takes long time to compute
#'\dontrun{
#'res4 <- NonparamElectionForensics(
#'        dat, Candidates = paste("P", 1:12, sep=""),
#'        CandidatesText=c("Stanislav Govorukhin", "Umar Dzhabrailov",
#'        "Vladimir Zhirinovsky", "Gennady Zuganov",
#'        "Ella Pamfilova", "Alexei Podberezkin",
#'        "Vladimir Putin", "Yuri Skuratov",
#'        "Konstantin Titov", "Aman Tuleev",
#'        "Grigorii Yavlinsky", "Against All"),
#'        MainCandidate = "P7",
#'        TotalReg = "NVoters",
#'        TotalVotes = "NValid",
#'        Level = "regname",
#'        MaxThreshold = 0.8,
#'        mode_search = list(npeaks = 5, sortstr = TRUE,
#'                minpeakdistance = 1, pick_by = "height"),
#'        FigureName = "Russian Presidential Elections, 2000",
#'        setcolors = c("royalblue2", "springgreen1","blue",
#'                      "red", "green","brown2",
#'                      "darkgreen", "yellow", "lawngreen",
#'                      "purple","chartreuse1", "orange"),
#'        precinctLevel = TRUE, computeSD="parametric",
#'        sims=2, grid_type = "2D")
#'
#'# Sum precinct-level fraud estimates
#'sum(res4$fraud_precinct_data$precinct_mean_hetero, na.rm=TRUE)
#'# -8503879
#'
#'# Sum precinct-level fraud (only statistically significant values)
#'sum(
#'  res4$sim_hetero_stats_sims[
#'  res4$fraud_precinct_data$precinct_sd_hetero_sig==TRUE
#'  ], na.rm=TRUE)
#'# 0
#'}

NonparamElectionForensics<-function(data,
                                Candidates,
                                CandidatesText = NULL,
                                MainCandidate,
                                TotalReg,
                                TotalVotes=NULL,
                                Level = NULL,
                                MaxThreshold = 0.8,
                                FigureName,
                                setcolors = NULL,
                                precinctLevel = TRUE,
                                computeSD = NULL,
                                sims = 10,
                                mode_search = list(npeaks = 5,
                                                   sortstr = TRUE,
                                                   minpeakdistance = 1,
                                                   pick_by = "height"),
                                man_turnout = NULL,
                                grid_type = "1D"){

  if(is.null(Level)) Level <-"National"
  graph_results <- NULL
  base_precinct <- base_sim_precinct <- sim_all_precinct <- sim_hetero_precinct <- NULL
  base_stats <- sim_all_stats <- sim_hetero_stats1 <- sim_hetero_stats2 <- NULL
  sim_hetero_stats_base <- sim_hetero_stats_sims <- NULL

  stats_table <- NULL

  m1 <- m2 <- m3 <- d1 <- d2 <- d3 <- d4 <- NULL

  mat_result <- NULL

  completeDat = TRUE

  data <- var_addition(data = data, Level=Level,
                       MainCandidate = MainCandidate, Candidates = Candidates,
                       TotalVotes = TotalVotes, TotalReg = TotalReg)

  columns_to_check <- c(MainCandidate, Candidates, TotalVotes, TotalReg)

  for (col in columns_to_check) {
    if (all(is.na(data[[col]]))) {
      stop(paste("Incomplete data found in column:", col, "- exiting the script."))
    }
  }

  if(is.null(setcolors)){
    setcolors.v <- gen_colors(Candidates=Candidates, MainCandidate=MainCandidate)
    }else{
    if(length(setcolors)!=length(Candidates)){
    setcolors.v <- gen_colors(Candidates=Candidates, MainCandidate=MainCandidate)
    }else{
    setcolors.v <- c("grey", setcolors)
    }
      }

   if(Level=="National"){

       if(grid_type == "1D"){
         estimate_fraud_res <- estimate_fraud(data, Candidates, CandidatesText, MainCandidate,
                                              FigureName, MaxThreshold, colors.v = setcolors.v,
                                              precinct.level = precinctLevel,
                                              mode_search = mode_search,
                                              man_turnout = man_turnout)
        }else if(grid_type == "2D"){
          estimate_fraud_res <- estimate_fraud_2d(data, Candidates, CandidatesText, MainCandidate,
                                                  FigureName, MaxThreshold, colors.v = setcolors.v,
                                                  precinct.level = precinctLevel,
                                                  mode_search = mode_search,
                                                  man_turnout = man_turnout,
                                                  graphsoff = FALSE)
          }

       precinct_fraud <- estimate_fraud_res$precinct_fraud$fraud.votes
       base_precinct[['Whole dataset']] <- subset(estimate_fraud_res$precinct_fraud, select=c("id", "fraud.votes"))
       graph_results[['Whole dataset']] <- estimate_fraud_res$drawfig
       base_stats[['Whole dataset']] <- unlist(estimate_fraud_res[1:8])

    if(computeSD=="parametric"){
        uncertainty_res <-  ComputeSimulationParametric(data, Candidates,
                                                        CandidatesText, MainCandidate,
                                                        FigureName, MaxThreshold, colors.v = setcolors.v,
                                                        precinct.level = precinctLevel, sims=sims,
                                                        mode_search = mode_search,
                                                        man_turnout = man_turnout,
                                                        grid_type = grid_type)
      sim_all_stats [['Whole dataset']] <- uncertainty_res$param_stats
      uncertainty_res$precinct_stats$sig_all <- c(abs(uncertainty_res$precinct_stats['precinct_mean']/uncertainty_res$precinct_stats['precinct_sd'])>=1.96)
      sim_all_precinct[['Whole dataset']] <- uncertainty_res$precinct_stats
    }

    if(computeSD=="nonparametric"){
        uncertainty_res <- ComputeSimulationNonparametric(data, Candidates,
                                                          CandidatesText, MainCandidate,
                                                          FigureName, MaxThreshold, colors.v = setcolors.v,
                                                          precinct.level = precinctLevel, sims=sims,
                                                          mode_search = mode_search,
                                                          man_turnout = man_turnout,
                                                          grid_type = grid_type)
      sim_all_stats [['Whole dataset']] <- uncertainty_res$param_stats
      uncertainty_res$precinct_stats$sig_all <- c(abs(uncertainty_res$precinct_stats['precinct_mean']/uncertainty_res$precinct_stats['precinct_sd'])>=1.96)
      sim_all_precinct[['Whole dataset']] <- uncertainty_res$precinct_stats
    }
       }

  if(Level!="National"){

      if(grid_type == "1D"){
        estimate_fraud_res <- estimate_fraud(data, Candidates, CandidatesText, MainCandidate,
                                             FigureName, MaxThreshold, colors.v = setcolors.v,
                                             precinct.level = precinctLevel,
                                             mode_search = mode_search,
                                             man_turnout = man_turnout)
      }else if(grid_type == "2D"){
        estimate_fraud_res <- estimate_fraud_2d(data, Candidates, CandidatesText, MainCandidate,
                                                FigureName, MaxThreshold, colors.v = setcolors.v,
                                                precinct.level = precinctLevel,
                                                mode_search = mode_search,
                                                man_turnout = man_turnout,
                                                graphsoff = FALSE)
      }

      base_stats[['Whole dataset']] <- unlist(estimate_fraud_res[1:8])
      base_precinct[['Whole dataset']] <- subset(estimate_fraud_res$precinct_fraud, select=c("id", "fraud.votes"))
      graph_results[['Whole dataset']] <- estimate_fraud_res

      if(computeSD=="parametric"){
        uncertainty_res <- ComputeSimulationParametric(data, Candidates,
                                                        CandidatesText, MainCandidate,
                                                        FigureName, MaxThreshold, colors.v = setcolors.v,
                                                        precinct.level = precinctLevel, sims=sims,
                                                        mode_search = mode_search,
                                                        man_turnout = man_turnout,
                                                        grid_type = grid_type)

        sim_all_stats [['Whole dataset']] <- uncertainty_res$param_stats
        uncertainty_res$precinct_stats$sig_all <- c(abs(uncertainty_res$precinct_stats['precinct_mean']/uncertainty_res$precinct_stats['precinct_sd'])>=1.96)
        sim_all_precinct[['Whole dataset']] <- uncertainty_res$precinct_stats
        }

      if(computeSD=="nonparametric"){
        uncertainty_res <- ComputeSimulationNonparametric(data, Candidates,
                                                          CandidatesText, MainCandidate,
                                                          FigureName, MaxThreshold, colors.v = setcolors.v,
                                                          precinct.level = precinctLevel, sims=sims,
                                                          mode_search = mode_search,
                                                          man_turnout = man_turnout,
                                                          grid_type = grid_type)
        sim_all_stats [['Whole dataset']] <- uncertainty_res$param_stats
        uncertainty_res$precinct_stats$sig_all <- c(abs(uncertainty_res$precinct_stats['precinct_mean']/uncertainty_res$precinct_stats['precinct_sd'])>=1.96)
        sim_all_precinct[['Whole dataset']] <- uncertainty_res$precinct_stats
      }

      splby <- data[,names(data) %in% Level]
      mdat <- split(data, splby)

      for(i in 1:length(names(mdat))){
        completeDat = TRUE
        mitem <- mdat[[i]]
        FigureName <- names(mdat)[i]
        print(FigureName)

        # check if any of the specified columns are fully NA
        for (col in columns_to_check) {
          if (all(is.na(mitem[[col]]))) {
            print(paste("Incomplete data found in column:", col))
            completeDat = FALSE
          }
        }

        if (completeDat){
          if(grid_type == "1D"){
            estimate_fraud_res <- tryCatch(estimate_fraud(mitem, Candidates, CandidatesText, MainCandidate,
                                                 FigureName, MaxThreshold, colors.v = setcolors.v,
                                                 precinct.level = precinctLevel,
                                                 mode_search = mode_search,
                                                 man_turnout = man_turnout),
                                           error = function(e) e)

          }else if(grid_type == "2D"){
            estimate_fraud_res <- estimate_fraud_2d(mitem, Candidates, CandidatesText, MainCandidate,
                                                    FigureName, MaxThreshold, colors.v = setcolors.v,
                                                    precinct.level = precinctLevel,
                                                    mode_search = mode_search,
                                                    man_turnout = man_turnout,
                                                    graphsoff = FALSE)
            }

          precinct_fraud <- estimate_fraud_res$precinct_fraud$fraud.votes
          base_sim_precinct[[FigureName]] <- subset(estimate_fraud_res$precinct_fraud,
                                                    select=c("id", "fraud.votes"))
          graph_results[[FigureName]] <- estimate_fraud_res
          sim_hetero_stats1[[FigureName]] <- unlist(estimate_fraud_res[1:8])

          if(computeSD=="parametric"){
            uncertainty_res <- tryCatch(ComputeSimulationParametric(mitem, Candidates,
                                                           CandidatesText, MainCandidate,
                                                           FigureName, MaxThreshold, colors.v = setcolors.v,
                                                           precinct.level = precinctLevel, sims=sims,
                                                           mode_search = mode_search,
                                                           man_turnout = man_turnout,
                                                           grid_type = grid_type),
                                        error = function(e) e)

            uncertainty_res$precinct_stats$sig_all <- c(abs(uncertainty_res$precinct_stats['precinct_mean']/uncertainty_res$precinct_stats['precinct_sd'])>=1.96)
            sim_hetero_precinct[[FigureName]] <- uncertainty_res$precinct_stats
            sim_hetero_stats2[[FigureName]] <- uncertainty_res$param_stats
          }

          if(computeSD=="nonparametric"){
            uncertainty_res <- ComputeSimulationNonparametric(mitem, Candidates,
                                                              CandidatesText, MainCandidate,
                                                              FigureName, MaxThreshold, colors.v = setcolors.v,
                                                              precinct.level = precinctLevel, sims=sims,
                                                              mode_search = mode_search,
                                                              man_turnout = man_turnout,
                                                              grid_type = grid_type)

            uncertainty_res$precinct_stats$sig_all <- c(abs(uncertainty_res$precinct_stats['precinct_mean']/uncertainty_res$precinct_stats['precinct_sd'])>=1.96)
            sim_hetero_precinct[[FigureName]] <- uncertainty_res$precinct_stats
            sim_hetero_stats2[[FigureName]] <- uncertainty_res$param_stats
            }
        } #end
      }
      }

  if (length(graph_results)>1){
    list_graphs  <-  lapply(graph_results, function(x) x[["drawfig"]])
    stats_table <- do.call(rbind.data.frame, lapply(names(graph_results), function(region) {
      x <- graph_results[[region]]
      x_clean <- x[setdiff(names(x), c("precinct_fraud", "drawfig"))]
      row <- as.data.frame(t(unlist(x_clean)))
      row$region <- region
      return(row)
    }))
    stats_table <- stats_table[, c("region", setdiff(names(stats_table), "region"))]
    }else{
      list_graphs <- graph_results$`Whole dataset`
    }

  if(!is.null(base_precinct)){
    mat_result <- base_precinct[['Whole dataset']]
    names(mat_result)[names(mat_result) == 'fraud.votes'] <- 'base.fraud.votes'
  }

  if(!is.null(base_precinct)&!is.null(sim_all_precinct)){
    d1 <- base_precinct[['Whole dataset']]; d2 <- sim_all_precinct[['Whole dataset']]
    names(d1)[names(d1) == 'fraud.votes'] <- 'base.fraud.votes'
    names(d2)[names(d2) == 'precinct_mean'] <- 'sim.precinct_mean'
    names(d2)[names(d2) == 'precinct_sd'] <- 'sim.precinct_sd'
    names(d2)[names(d2) == 'sig_all'] <- 'sim.sig_all'
    m1 <- merge(d1, d2, by="id")
    mat_result <- m1
    }

  if (length(base_sim_precinct)!=0){
    d3 <- do.call(rbind, lapply(names(base_sim_precinct), function(region) {
      df <- base_sim_precinct[[region]]
      data.frame(
        id = df$id,
        fraud.votes = df$fraud.votes
      )
    }))
    names(d3)[names(d3) == 'fraud.votes'] <- 'base.sim.fraud.votes'
    }

  if (length(sim_hetero_precinct) != 0) {
    d4 <- do.call(rbind, lapply(names(sim_hetero_precinct), function(region) {
      df <- sim_hetero_precinct[[region]]
      data.frame(
        #region = region,
        id = df$id,  # globally unique ID
        precinct_mean_hetero = df$precinct_mean,
        precinct_sd_hetero = df$precinct_sd,
        precinct_sd_hetero_sig = df$sig_all
      )
    }))
  }

  if (!is.null(d3) & !is.null(d4)){
    m2 <- merge(d3, d4, by="id")
    m3 <- merge(m1, m2, by="id")
    mat_result <- m3}

  if (!is.null(sim_hetero_stats1)){
    sim_hetero_stats_base <- do.call(rbind, sim_hetero_stats1)
    }

  if (!is.null(sim_hetero_stats2)){
    sim_hetero_stats_sims <- do.call(rbind, sim_hetero_stats2)
  }

  list_results <- list(list_graphs = list_graphs,
                       base_stats = base_stats,
                       sim_all_stats = sim_all_stats,
                       sim_hetero_stats_base = sim_hetero_stats_base,
                       sim_hetero_stats_sims = sim_hetero_stats_sims,
                       fraud_precinct_data = mat_result,
                       data = data,
                       Level = Level, creationdate = Sys.time())
  return(list_results)}
