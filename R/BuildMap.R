#' @title Build Map function
#' @description This function constructs the map based on the results from BasicElectionForensics() object.
#' @usage BuildMap(eforensicsdata, geodata, Geoindex, Colorsig, xlab)
#' @param eforensicsdata  BasicElectionForensics() object
#' @param geodata sf object
#' @param Geoindex name of index variable from sf object used to merge geodata with BasicElectionForensics() object
#' @param Colorsig only statistically significant estimates are mapped (FALSE by default)
#' @param xlab X axis label
#' @export
#' @import foreign
#' @import sf
#' @import RColorBrewer
#' @importFrom methods as
#' @importFrom sp spplot
#' @return A list containing the following parameters:
#' \itemize{
#'   \item figures - list of constructed figures
#'   \item Colorsig -  external parameter
#'   \item shpdata - updated sf object
#'   \item creationdate - date/time of analysis
#' }
#' @examples
#' library(EFToolkit)
#' library(sf)
#'
#' dat<-read.csv(system.file("extdata/Albania2013.csv", package="EFToolkit"))
#' eldata<-BasicElectionForensics(dat,
#'                               Candidates=c("C035", "C050"),
#'                               Level="Prefectures", TotalReg="Registered",
#'                               TotalVotes="Ballots",
#'                               Methods=c("P05s", "C05s"), R=100)
#' geodata<-st_read(system.file("extdata/Albania2013_prefectures.shp",
#'                           package="EFToolkit"), quiet = TRUE)
#'
#' figures<-BuildMap(eforensicsdata=eldata, geodata=geodata, Geoindex="Level")

BuildMap<-function(eforensicsdata=eforensicsdata, geodata=geodata, Geoindex=Geoindex, Colorsig=FALSE, xlab=""){

  if(is.list(eforensicsdata)&any(names(eforensicsdata)=="table")){
    cdat<-eforensicsdata[[1]][seq(1,dim(eforensicsdata[[1]])[1],2),]
    mdat<-eforensicsdata[[4]][seq(1,dim(eforensicsdata[[4]])[1],2),]
    mdat<-apply(mdat, 2, function(x) as.numeric(gsub("color:red", 1, x)))

    mdat<-as.data.frame(mdat[,-c(1:2, dim(mdat)[2])])
    mdat$Level<-cdat$Level
  }

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
  n_oldf<-names(geodata)
  colnames(merger_fin)<-n_names

  # Convert to data.frame for merging, then back to sf
  geodata_df <- st_drop_geometry(geodata)
  geoele <- merge(geodata_df, merger_fin, by.x = Geoindex, by.y = "Level", all.x=TRUE)
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

  # Reattach geometry and convert back to sf
  geodata <- st_as_sf(cbind(geoele, st_geometry(geodata)))

  #Draw graph
  my.palette <- brewer.pal(n = 7, name = "OrRd")

  # Convert to Spatial for spplot if needed, or use tmap for sf objects
  # Here I'll show both options - uncomment the one you prefer

  # Option 1: Keep using spplot by converting to Spatial (temporary)
  geodata_sp <- as(geodata, "Spatial")

  if(Colorsig){
    lmaps<-sapply(names(ndat), function(x){
      var_name<-x
      spplot(geodata_sp, paste("S", var_name, sep=""), main = paste("Figure for ",  var_name, sep=""),
             sub = xlab,
             col.regions = my.palette,
             cuts = 1, col = "transparent")}, simplify = FALSE, USE.NAMES = TRUE)
  }else{
    lmaps<-sapply(names(ndat), function(x){
      var_name<-x
      spplot(geodata_sp, var_name, main = paste("Figure for ",  var_name, sep=""),
             sub = xlab,
             col.regions = my.palette,
             cuts = 6, col = "transparent")}, simplify = FALSE, USE.NAMES = TRUE)}

  # Option 2: Use tmap for native sf plotting (recommended)
  # library(tmap)
  # if(Colorsig){
  #   lmaps<-sapply(names(ndat), function(x){
  #     var_name<-x
  #     tm_shape(geodata) +
  #       tm_fill(paste("S", var_name, sep=""),
  #               title = paste("Figure for ", var_name, sep=""),
  #               palette = my.palette,
  #               style = "fixed", breaks = c(0, 1)) +
  #       tm_layout(main.title = xlab)
  #   }, simplify = FALSE, USE.NAMES = TRUE)
  # }else{
  #   lmaps<-sapply(names(ndat), function(x){
  #     var_name<-x
  #     tm_shape(geodata) +
  #       tm_fill(var_name,
  #               title = paste("Figure for ", var_name, sep=""),
  #               palette = my.palette,
  #               style = "fixed", breaks = 6) +
  #       tm_layout(main.title = xlab)
  #   }, simplify = FALSE, USE.NAMES = TRUE)
  # }

  return_result=list(figures=lmaps, Colorsig=Colorsig, shpdata=geodata, creationdate=Sys.time())

  return(return_result)
}
