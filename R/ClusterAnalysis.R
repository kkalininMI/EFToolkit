#' @title Cluster Analysis
#' @description This function implements geographic clustering tests: Getis-Ord Gi and local Moran's I.
#' @usage ClusterAnalysis(geodata, Vars, IndexCL, cores)
#' @param geodata sf object within a list (the algorithm also takes both polygon and point sf objects within a list) or BuildMap object().
#' @param Vars variable names used for geographic clustering tests
#' @param IndexCL index variable name used for merging polygon and point sf objects
#' @param cores number of cores for parallel computing (2 cores by default)
#' @export
#' @import methods
#' @import digest
#' @import foreign
#' @import spdep
#' @import RANN
#' @import foreach
#' @import sf
#' @import ggplot2
#' @import Matrix
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach %dopar% %do%
#' @importFrom ggplot2 ggplot geom_sf geom_point scale_fill_identity
#' @importFrom ggplot2 scale_color_identity labs theme_void theme
#' @importFrom ggplot2 element_text element_rect aes
#' @importFrom stats var
#' @return A list of maps with results from cluster analysis:
#' \itemize{
#'   \item table - dataframe with results
#'   \item html -  html table with results
#'   \item tex - tex table with results
#'   \item sigMatrix - significance matrix
#' }
#' @examples
#' library(EFToolkit)
#' library(sf)
#'
#' dat<-read.csv(system.file("extdata/Albania2013.csv", package="EFToolkit"))
#' #Obtain election forensics estimates
#' eldata<-BasicElectionForensics(dat,
#'                               Candidates=c("C035", "C050"),
#'                               Level="Prefectures", TotalReg="Registered",
#'                               TotalVotes="Ballots",
#'                               Methods=c("P05s", "C05s"), R=100)
#'
#' #Create the map with results
#' geodata<-st_read(system.file(
#'            "extdata/Albania2013_prefectures.shp",
#'            package="EFToolkit"), quiet = TRUE)
#' figures<-BuildMap(eforensicsdata=eldata, geodata=geodata,
#'                   Geoindex="Level", Colorsig=FALSE)
#'
#' #Using the mapped results, implement cluster analysis
#' figureC<-ClusterAnalysis(figures, Vars=c("C050_P05s", "C050_C05s"))
#'

ClusterAnalysis <- function(geodata, Vars, IndexCL = NULL, cores = 2) {

  usecores <- cores

  uploadmap <- function(geodata) {
    if (length(geodata) == 1) {
      CLfile1 <- geodata
      Return_File <- list(CLfile1 = CLfile1)
    } else if (length(geodata) == 2) {
      CLfile1 <- geodata[[1]]
      CLfile2 <- geodata[[2]]

      # Attempt to detect the correct file order
      detectType <- ifelse(
        inherits(CLfile1, "sf") &&
          all(sf::st_geometry_type(CLfile1) %in% c("POLYGON", "MULTIPOLYGON")),
        2,
        1
      )

      if (detectType == 2) {
        # If CLfile1 is spatial, it's the second file
        temp <- CLfile1
        CLfile1 <- CLfile2
        CLfile2 <- temp
      }
      Return_File <- list(CLFile1 = CLfile1, CLFile2 = CLfile2)
    } else {
      stop("Unexpected geodata format.")
    }
    return(Return_File)
  }

  createfigures <- function(graph.object) {
    output <- list()
    for (i in seq(2, length(graph.object), 2)) {
      gr1 <- graph.object[[i]]$gr1
      gr2 <- graph.object[[i]]$gr2
      gr1$poshpZ$pcolor <- gr1$pcolorZ
      gr2$poshpZ$pcolor <- gr2$pcolorZ

      plotnameM <- paste("Moran's I for", graph.object[[i - 1]])
      plotnameG <- paste("Getis-Ord for", graph.object[[i - 1]])

      # Moran's I plot using ggplot2
      if (all(sf::st_geometry_type(gr1$poshpZ) %in% c("POLYGON", "MULTIPOLYGON"))) {
        output[[plotnameM]] <- ggplot(gr1$poshpZ) +
          geom_sf(aes(fill = pcolor), color = "white", size = 0.1) +
          scale_fill_identity() +
          labs(title = plotnameM) +
          theme_void() +
          theme(
            plot.title = element_text(hjust = 0.5, color = "black"),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA)
          )
      } else {
        # For point data
        coords <- st_coordinates(gr1$poshpZ)
        plot_data <- cbind(as.data.frame(gr1$poshpZ), coords)

        output[[plotnameM]] <- ggplot(plot_data) +
          geom_point(aes(x = X, y = Y, color = pcolor), size = 2) +
          scale_color_identity() +
          labs(title = plotnameM) +
          theme_void() +
          theme(
            plot.title = element_text(hjust = 0.5, color = "black"),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA)
          )
      }

      # Getis-Ord plot using ggplot2
      if (all(sf::st_geometry_type(gr2$poshpZ) %in% c("POLYGON", "MULTIPOLYGON"))) {
        output[[plotnameG]] <- ggplot(gr2$poshpZ) +
          geom_sf(aes(fill = pcolor), color = "white", size = 0.1) +
          scale_fill_identity() +
          labs(title = plotnameG) +
          theme_void() +
          theme(
            plot.title = element_text(hjust = 0.5, color = "black"),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA)
          )
      } else {
        # For point data
        coords <- st_coordinates(gr2$poshpZ)
        plot_data <- cbind(as.data.frame(gr2$poshpZ), coords)

        output[[plotnameG]] <- ggplot(plot_data) +
          geom_point(aes(x = X, y = Y, color = pcolor), size = 2) +
          scale_color_identity() +
          labs(title = plotnameG) +
          theme_void() +
          theme(
            plot.title = element_text(hjust = 0.5, color = "black"),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA)
          )
      }
    }
    return(output)
  }

  # Use the appropriate geodata structure
  if (is.list(geodata) && any(names(geodata) == "shpdata")) {
    dataCL <- list(list(CLFile1 = geodata[["shpdata"]]))
  } else {
    dataCL <- uploadmap(geodata)
  }

  # Call clusterdetector function
  result.clusters <- clusterdetector(
    dataCL = dataCL,
    IndexCL = IndexCL,
    Vars = Vars,
    usecores = usecores
  )

  result <- createfigures(result.clusters)
  return(result)
}

# Your existing clusterdetector function remains the same
clusterdetector <- function(dataCL=dataCL, Vars=Vars, IndexCL=IndexCL, usecores = usecores){
  asc <- function(z) { as.character(z) }

  fdr1 <- function(p, level=.05) {
    lenp <- length(p);
    test <- (1:lenp)/lenp*level > sort(p);
    nrejects <- sum(test);
    if (nrejects > 0) {
      oop <- order(order(p));
      pickrejects <- (1:lenp)[test[oop]];
      idx <- 1:max((1:lenp)[test]);
      nrejects <- length(idx);
      rejects <- sort(p)[idx];
    }
    else {
      pickrejects <- idx <- rejects <- NULL;
    }
    list(ntests=lenp, nrejects=nrejects, rejects=rejects,
         pickrejects=pickrejects);
  }

  coltable <- matrix(NA,5,3);
  dimnames(coltable) <- list(c("LL","HH","LH","HL","ns"),c("99","95","90"));
  coltable["LL",] <- c("#0000FF","#B4B4FF","#E6E6E6"); # blue
  coltable["HH",] <- c("#FF0000","#FFB4B4","#E6E6E6"); # red
  coltable["LH",] <- c("#32CD32","#C8FFC8","#E6E6E6"); # green
  coltable["HL",] <- c("#FF9900","#FFFFB4","#E6E6E6"); # orange
  coltable["ns",] <- c("#E6E6E6","#E6E6E6","#E6E6E6");

  permuteStart <- function(j,wshp,X,w,k,longlat) {
    x <- X[,j];
    x <- ifelse(x<0,NA,x);
    if (any(is.na(x))) {
      locm <- (1:length(x))[is.na(x)];
      wshpj <- wshp[-locm,];
      comatj <- sf::st_coordinates(sf::st_centroid(wshpj));
      NNj <- knearneigh(comatj, k=k, longlat=longlat);
      wj <- nb2listw(knn2nb(NNj),style="B");
      xj <- x[-locm];
      mj <- localmoran(xj, wj, zero.policy=FALSE);
      gj <- localG(xj, wj, zero.policy=FALSE);
      m <- matrix(NA,length(x),5);
      dimnames(m)[[2]] <- dimnames(mj)[[2]];
      attr(m,"class") <- attr(mj,"class");
      attr(m,"call") <- attr(mj,"call");
      m[-locm,] <- mj;
      g <- rep(NA,length(x));
      attr(g,"class") <- attr(gj,"class");
      attr(g,"call") <- attr(gj,"call");
      attr(g,"gstarti") <- attr(gj,"gstarti");
      g[-locm] <- gj;
    }
    else {
      NNj <- NULL;
      m <- localmoran(x, w, zero.policy=FALSE);
      g <- localG(x, w, zero.policy=FALSE);
    }
    list(m=m, g=g, NN=NNj);
  }

  permuteIter <- function(n,X,NN,NNlist) {
    s <- sample(n);
    NNs <- t(sapply(1:n,function(k){ s[NN[[1]][k,]]; }));
    loc <- sapply(1:n, function(k){ match(k,NNs[k,]); });
    for (ii in (1:n)[!is.na(loc)]) {
      NNs[ii,loc[!is.na(loc)][ii]] <- sample(s[!(s %in% NNs[ii,])],1);
    }
    NNi <- NN;
    NNi[[1]] <- NNs;
    wi <- nb2listw(knn2nb(NNi),style="B");
    J <- dim(X)[2];
    Garr <- Iarr <- array(NA,dim=c(J,n));
    for (j in 1:J) {
      x <- X[,j];
      x <- ifelse(x<0,NA,x);
      if (any(is.na(x))) {
        locm <- (1:length(x))[is.na(x)];
        xj <- x[-locm];
        NNj <- NNlist[[j]];
        sj <- sample(nj <- dim(NNj[[1]])[1]);
        NNjs <- t(sapply(1:nj,function(k){ sj[NNj[[1]][k,]]; }));
        locj <- sapply(1:nj, function(k){ match(k,NNjs[k,]); });
        for (ii in (1:nj)[!is.na(locj)]) {
          NNjs[ii,locj[!is.na(locj)][ii]] <- sample(sj[!(sj %in% NNjs[ii,])],1);
        }
        NNj[[1]] <- NNjs;
        wj <- nb2listw(knn2nb(NNj),style="B");
        Iarr[j,-locm] <- localmoran(xj, wj, zero.policy=FALSE)[,"Ii"];
        Garr[j,-locm] <- localG(xj, wj, zero.policy=FALSE);
      }
      else {
        Iarr[j,] <- localmoran(x, wi, zero.policy=FALSE)[,"Ii"];
        Garr[j,] <- localG(x, wi, zero.policy=FALSE);
      }
    }
    list(Iarr=Iarr, Garr=Garr);
  }

  permutePost <- function(j,n,R,mglist,Iarr,Garr) {
    simg <- simp <- rep(NA,n);
    for (i in 1:n) {
      lo <- sum(mglist[[j]]$m[i,1] < Iarr[j,i,], na.rm=TRUE);
      hi <- sum(mglist[[j]]$m[i,1] > Iarr[j,i,], na.rm=TRUE);
      if (all(is.na(Iarr[j,i,]))) { lo <- hi <- 1; }
      simp[i] <- min(lo,hi)/R;
      lo <- sum(mglist[[j]]$g[i] < Garr[j,i,], na.rm=TRUE);
      hi <- sum(mglist[[j]]$g[i] > Garr[j,i,], na.rm=TRUE);
      if (all(is.na(Garr[j,i,]))) { lo <- hi <- 1; }
      simg[i] <- min(lo,hi)/R;
    }
    list(mglistM=mglist[[j]]$m, simM=simp, mglistG=mglist[[j]]$g, simG=simg);
  }

  permuteMG <- function(X, NN, wshp, R=1000, k=8, longlat=TRUE) {
    w <- nb2listw(knn2nb(NN),style="B");
    J <- dim(X)[2];
    NNlist <- list();
    registerDoParallel(cores=usecores);
    mglist <- foreach(j=1:J, .export="permuteStart", .packages=c("spdep", "Matrix")) %dopar% permuteStart(j,wshp,X,w,k,longlat);
    stopImplicitCluster();
    for (j in 1:J) {
      NNlist[[j]] <- mglist[[j]]$NN;
      mglist[[j]]$NN <- list(NULL);
    }

    n <- dim(X)[1];
    registerDoParallel(cores=usecores);
    GIlist <- foreach(i=1:R, .export="permuteIter", .packages=c("spdep", "Matrix")) %dopar% permuteIter(n,X,NN,NNlist);
    stopImplicitCluster();
    Garr <- Iarr <- array(NA,dim=c(J,dim(X)[1],R));
    for (i in 1:R) {
      Garr[,,i] <- GIlist[[i]]$Garr;
      Iarr[,,i] <- GIlist[[i]]$Iarr;
    }
    registerDoParallel(cores=usecores);
    mgRlist <- foreach(j=1:J, .export="permutePost", .packages=c("Matrix")) %dopar% permutePost(j,n,R,mglist,Iarr,Garr);
    stopImplicitCluster();
    mgRlist;
  }

  dolocalIter <- function(j,X,NN,pMG) {
    x <- X[,j];
    x <- ifelse(x<0,NA,x);
    xm <- x-mean(x,na.rm=TRUE);
    Wx <- sapply(1:length(xm), function(i){ sum(xm[NN$nn[i,]],na.rm=TRUE) });
    widx <- Wx>0;
    rawtype <- ifelse(xm>0,ifelse(widx,"HH","HL"),ifelse(widx,"LH","LL"));
    zp <- pMG[[j]][["simM"]];
    f99 <- fdr1(zp, level=.01)$pickrejects;
    f95 <- fdr1(zp, level=.05)$pickrejects;
    f90 <- fdr1(zp, level=.1)$pickrejects;
    if (length(f95)>0) f90 <- f90[!(f90 %in% f95)];
    if (length(f99)>0) f95 <- f95[!(f95 %in% f99)];
    msig <- rep(0,length(x));
    msig[f99] <- 3;
    msig[f95] <- 2;
    msig[f90] <- 1;
    type <- ifelse(msig>0,rawtype,"");

    g <- pMG[[j]][["mglistG"]];
    zp <- pMG[[j]][["simG"]];
    f99 <- fdr1(zp, level=.01)$pickrejects;
    f95 <- fdr1(zp, level=.05)$pickrejects;
    f90 <- fdr1(zp, level=.1)$pickrejects;
    if (length(f95)>0) f90 <- f90[!(f90 %in% f95)];
    if (length(f99)>0) f95 <- f95[!(f95 %in% f99)];
    gsign <- ifelse(g>0,1,ifelse(g<0,-1,0));
    sig <- rep(0,length(g));
    sig[f99] <- ifelse(gsign[f99]==1,3,-3);
    sig[f95] <- ifelse(gsign[f95]==1,2,-2);
    sig[f90] <- ifelse(gsign[f90]==1,1,-1);

    list(localm=pMG[[j]][["mglistM"]],msig=msig,type=type,
         localg=pMG[[j]][["mglistG"]],gsig=sig,pMG=pMG[[j]]);
  }

  dolocal <- function(X, comat, k, wshp, R=1000, ll=TRUE) {
    NN <- knearneigh(comat, k=k, longlat=ll);
    pMG <- permuteMG(X, NN, wshp, R=R, k=k, longlat=ll);
    J <- dim(X)[2];
    registerDoParallel(cores=usecores);
    outlist <- foreach(j=1:J, .export=c("dolocalIter", "fdr1")) %dopar% dolocalIter(j,X,NN,pMG);
    stopImplicitCluster();
    outlist;
  }

  plocalM <- function(poshp, localobj) {
    z <- localobj$localm[,5];
    msig <- localobj$msig;
    sig <- rep("90",length(z));
    sig[msig==1] <- "90";
    sig[msig==2] <- "95";
    sig[msig==3] <- "99";
    COtype <- localobj$type;
    COtype <- ifelse(COtype=="","ns",as.character(COtype));
    COtype <- ifelse(is.na(COtype),"ns",COtype);
    pcolor <- sapply(1:length(sig),function(i){ coltable[COtype[i],sig[i]]; });
    oz <- order(abs(z),decreasing=FALSE);
    list(poshpZ=poshp[oz,], pcolorZ=pcolor[oz]);
  }

  plocalG <- function(poshp, localobj) {
    z <- localobj$localg;
    sig <- localobj$gsig;
    pcolor <-
      ifelse(sig==(-3),"#0000FF",
             ifelse(sig==(-2),"#2791C3",
                    ifelse(sig==(-1),"#CAE1FF",
                           ifelse(sig==(3),"#FF0000",
                                  ifelse(sig==(2),"#D08B6C",
                                         ifelse(sig==(1),"#FFAEB9","#E6E6E6"))))));
    oz <- order(abs(z),decreasing=FALSE);
    list(poshpZ=poshp[oz,], pcolorZ=pcolor[oz]);
  }

  mkpdf <- function(vnames,pdfname,wshp,wshpplot, ll=TRUE) {
    varID <- rep(NA,length(vnames));
    for (j in 1:length(vnames)) {
      varID[j] <-
        digest::digest(as.character(wshp[[vnames[j]]]), algo="sha1", serialize=FALSE);
    }
    comat <- sf::st_coordinates(sf::st_centroid(wshp));
    dimnames(comat)[[2]] <- c("long","lati");
    moutlist <- dolocal(as.data.frame(wshp)[,vnames,drop=FALSE], comat, 8, wshp, ll=ll);
    GraphStorage<-""
    for (j in 1:length(vnames)) {
      gr1 <- plocalM(wshpplot, moutlist[[j]]);
      gr2 <- plocalG(wshpplot, moutlist[[j]]);
      CandName<-vnames[j]
      GraphStorage<-c(GraphStorage,list(CandName=CandName, list(gr1=gr1,gr2=gr2)))
    }
    return(GraphStorage);
  }

  pdfbase<-names(dataCL)[1]
  MapsResultsCl <- list();

  if(length(dataCL)==2){
    ShpPoin<-dataCL[[1]]
    ShpPoly<-dataCL[[2]]
    loc <- match(asc(ShpPoin[[which(names(ShpPoin)%in%IndexCL)]]),
                 asc(ShpPoly[[which(names(ShpPoly)%in%IndexCL)]]));
    ShpPoly <- ShpPoly[loc,]
    idx <- sapply(Vars, function(v){ var(ShpPoin[[v]],na.rm=TRUE) })>0;
    if (any(idx)) {
      MapsResultsCl <- mkpdf(Vars[idx], pdfbase,ShpPoin,ShpPoly, ll=NULL);
    }
  }

  if(length(dataCL)==1){
    ShpData<-dataCL[[1]][1][[1]]
    idx <- sapply(Vars, function(v){ var(ShpData[[v]],na.rm=TRUE) })>0;
    if (any(idx)) {
      MapsResultsCl <- mkpdf(Vars, pdfbase,ShpData,ShpData, ll=NULL);
    }
  }
  MapsResultsCl[[1]] <- NULL

  return(MapsResultsCl)
}
