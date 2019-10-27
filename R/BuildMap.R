#' @title Build Map function
#' @description This function constructs the map based on the results from BasicElectionForensics() object.
#' @usage BuildMap(eforensicsdata, geodata, Geoindex, Colorsig, xlab)
#' @param eforensicsdata  BasicElectionForensics() object
#' @param geodata shapefile object
#' @param Geoindex name of index variable from shapefile object used to merge geodata with BasicElectionForensics() object
#' @param Colorsig only statistically significant estimates are mapped (FALSE by default)
#' @param xlab X axis label
#' @export
#' @import foreign
#' @import rgdal
#' @import sp
#' @import RColorBrewer
#' @return A list containing the following parameters:
#' \itemize{
#'   \item figures - list of constructed figures
#'   \item Colorsig -  external parameter
#'   \item shpdata - updated shp file
#'   \item creationdate - date/time of analysis
#' }
#' @examples
#' library(EFToolkit)
#' library(rgdal)
#'
#' dat<-read.csv(system.file("Albania2013.csv", package="EFToolkit"))
#' eldata<-BasicElectionForensics(dat,
#'                               Candidates=c("C035", "C050"),
#'                               Level="Prefectures", TotalReg="Registered",
#'                               TotalVotes="Ballots",
#'                               Methods=c("P05s", "C05s"), R=100)
#' geodata<-readOGR(system.file("Albania2013_prefectures.shp",
#'                           package="EFToolkit"), verbose = FALSE)
#'
#' figures<-BuildMap(eforensicsdata=eldata, geodata=geodata, Geoindex="Level")

############################################################
##               Election Forensics Toolkit               ##
##  25sep2015, 24oct2019                                  ##
##  Kirill Kalinin and Walter R. Mebane, Jr               ##
############################################################

BuildMap<-function(eforensicsdata=eforensicsdata, geodata=geodata, Geoindex=Geoindex, Colorsig=FALSE, xlab=""){


  if(is.list(eforensicsdata)&any(names(eforensicsdata)=="table")){
    cdat<-eforensicsdata[[1]][seq(1,dim(eforensicsdata[[1]])[1],2),]
    mdat<-eforensicsdata[[4]][seq(1,dim(eforensicsdata[[4]])[1],2),]
    mdat<-apply(mdat, 2, function(x) as.numeric(gsub("color:red", 1, x)))

    mdat<-as.data.frame(mdat[,-c(1:2, dim(mdat)[2])])
    mdat$Level<-cdat$Level
  }
  #if(is.list(eforensicsdata)&any(names(eforensicsdata)=="stats_table")){
  #  eforensicsdata$stats_table$Level<-rownames(eforensicsdata$stats_table)
  #  cdat<-data.frame("Variable1", "Variable2", subset(eforensicsdata$stats_table, select=c("Level", "official_turnout", "real_turnout",
  #                                                  "fraud_measure1", "fraud_measure2", "votes_incumbent")))
  #  }

  CandidatesName<-unique(cdat[,2])
  CandidatesNum<-length(CandidatesName)
  MethodsNum<-dim(cdat)[2]-3
  MethodsName<-names(cdat)[3:(2+MethodsNum)]

  s_cdat<-split(cdat, cdat[,2]);
  s_mdat<-split(mdat, cdat[,2]);


  if(names(s_cdat)[1]==""){s_cdat[[1]]<-NULL}

  merger<-s_cdat[['Turnout']]
  merger_m<-s_mdat[['Turnout']];merger_m$Level<-merger$Level

  cnames <- unique(cdat[[2]])[-1]

  for(Cand in cnames){
    merger<-merge(merger,s_cdat[[Cand]], by = "Level")
    merger_m<-merge(merger_m,s_mdat[[Cand]], by = "Level")
  }

  whichCandidates<-which(grepl("Candidate",names(merger)));
  obs_num<-which(grepl("Obs",names(merger)))
  obs_no_nam<-merger[,-c(whichCandidates,obs_num)]
  obs_no_namN<-suppressWarnings(
    apply(obs_no_nam[2:dim(obs_no_nam)[2]], 2, function(x){as.numeric(x)}))
  merger_fin<-data.frame(merger[,1], obs_no_namN)

  n_names<-c("Level",
             paste(rep(substring(CandidatesName,1,4),each=MethodsNum),
                   rep(MethodsName,CandidatesNum), sep="_"))
  n_namesf<-c("Level",
              paste(rep(CandidatesName,each=MethodsNum),
                    rep(MethodsName,CandidatesNum), sep="_"))
  n_oldf<-names(geodata@data)
  colnames(merger_fin)<-n_names
  geoele<-merge(geodata@data, merger_fin, by.x = Geoindex, by.y = "Level", all.x=TRUE)
  names(geoele)<-make.unique(names(geoele))
  n_newf<-names(geoele)[!names(geoele)%in%n_oldf]

  #remove nas
  ndat <- subset(geoele, select=c(n_newf))
  remCol <- apply(ndat,2, function(x) sum(is.na(x))==length(x))
  remCol <- which(remCol==TRUE)
  if(length(remCol)>0) ndat<-ndat[,-remCol];


  if(Colorsig){
    if(length(remCol)>0) merger_m<-merger_m[,-c(remCol+1)];
    colnames(merger_m)[2:dim(merger_m)[2]]<-paste("S",names(ndat), sep="")
    geoele<-merge(geoele, merger_m, by.x = Geoindex, by.y = "Level", all.x=TRUE)
  }

  geodata@data<-geoele

  #Draw graph
  my.palette <- brewer.pal(n = 7, name = "OrRd")

  if(Colorsig){
    lmaps<-sapply(names(ndat), function(x){
      var_name<-x
      spplot(geodata, paste("S", var_name, sep=""), main = paste("Figure for ",  var_name, sep=""),
             sub = xlab,
             col.regions = my.palette,
             cuts = 1, col = "transparent")}, simplify = FALSE, USE.NAMES = TRUE)
  }else{
    lmaps<-sapply(names(ndat), function(x){
      var_name<-x
      spplot(geodata, var_name, main = paste("Figure for ",  var_name, sep=""),
             sub = xlab,
             col.regions = my.palette,
             cuts = 6, col = "transparent")}, simplify = FALSE, USE.NAMES = TRUE)}

  return_result=list(figures=lmaps, Colorsig=Colorsig, shpdata=geodata, creationdate=Sys.time())

  return(return_result)}
