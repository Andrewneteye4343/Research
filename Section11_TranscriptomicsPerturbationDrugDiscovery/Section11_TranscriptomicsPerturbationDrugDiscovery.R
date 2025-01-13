#################################################################################################################
# R version 4.1.2
# Calculate the therapeutic score of compounds in d6 and d24 models with DEGs lists 
# DEGs lists only contain gene names, log2 fold change and p-value
DEBUG <- FALSE

# TO DO
PERT_TYPES               <- c("trt_cp__6h", "trt_cp__24h", "trt_sh__96h") # Vector of database
NCELL                    <- c("10C", "8C")[1]
REC.FDR                  <- 0.001
REC.SCORE.CUTOFF         <- c(5.622220831, 5.962561254, 5.465406168) ##### NOTION: REC scores for FDR = 0.001
names(REC.SCORE.CUTOFF)  <- PERT_TYPES
PERT_ABBR                <- c("d6", "d24", "sh96")
names(PERT_ABBR)         <- PERT_TYPES
#
z                        <- 2 # index for "term" and "pre" variable
term                     <- c("A", "B", "C")[z] # projects as elements
pre                      <- c("fct", "ft", "fc")[z] # The pre-name of input .txt data
post                     <- c(1,2,3) # The post-name of DEGs lists.txt
TUMOR_TYPES              <- paste0(pre, "_", post)
TUMOR_TYPES_DESC         <- paste0("Cluster", post)
names(TUMOR_TYPES_DESC)  <- TUMOR_TYPES
# N.CANCER.SIGNATURE       <- 1000
TUMOR.THRESHOLD.LOG2.FC  <- 0
TUMOR.THRESHOLD.LOG2.CPM <- 0
TUMOR.THRESHOLD.P.VALUE  <- 0.05
#
Nperm                    <- 1000 # number of permutations
Nexc                     <- 250 # number of exceedances used for fitting GPD
P.perm_threshold         <- 10/1000 # not necessarily having to conform to Nexc

# directory
myDelim <- function(names, delim = '\\', cr = T, after = F) {
  tmp <- ''; count <- 0
  for (name in names) {
    count <- count + 1
    if (count == 1) {tmp <- name} else {tmp <- paste(tmp, name, sep = delim)}
    if (cr == T) dir.create(tmp, showWarnings = F)
  }
  if (after) tmp <- paste(tmp, "\\", sep = "")
  return(tmp)
}
# P-value correction
AdjPvMtrx <- function(mx, method="BH") {
  tmp <- p.adjust(p=mx, method=method) # coerce the matrix col-by-col into a vector
  # NA allowed
  out <- matrix(tmp, nrow=nrow(mx), ncol=ncol(mx)) # coerce the vector back into a matrix
  dimnames(out) <- list(rownames(mx), colnames(mx))
  return(out)
}

# default writing and updating function (for internal use)
writeAndUpdate <- function(sc=0, pv=1) {
  # write therapeutic
  h5write(sc, oPath, "0/Sc" , index=list(H5_therap.current.row, H5_therap.current.col))
  h5write(pv, oPath, "0/Pv" , index=list(H5_therap.current.row, H5_therap.current.col))
  # update and write indices
  H5_therap.current.row <<- H5_therap.current.row + 1
  if (H5_therap.current.row > length(H5_therap.PERT.NAME)) {
    H5_therap.current.col <<- H5_therap.current.col + 1
    H5_therap.current.row <<- 1
  }
  h5write(H5_therap.current.row, oPath, "CURRENT.ROW")
  h5write(H5_therap.current.col, oPath, "CURRENT.COL")
}

# packages used
library(rhdf5)
library(stringr)
library(fExtremes)

# global directories
lincsDir  <- "D:\\REC\\LINCS"
sharedDir <- "D:\\REC"
jobDir    <- myDelim(c(sharedDir, "LabTask", "YourDirectory"))
shRefDir  <- myDelim(c(sharedDir, "LINCS"))
sDir      <- myDelim(c(jobDir, term))
oDir      <- sDir

# file names
oFN.h5    <- c() # to be determined for each perturbation type
shConvFN  <- "read.txt" # input file with Gene name (Column1) and log2FC (Column2, default should be 1)

# set paths
oPath      <- c() # to be determined for each perturbation type
shConvPath <- paste(shRefDir, shConvFN, sep="\\")

# global variables
TUMOR.DE.DATA.FC     <- NULL
TUMOR.DE.DATA.PV     <- NULL
TUMOR.DE.LIST        <- list()
TUMOR.DE.GENE.SPACE  <- c()
TUMOR.DE.GENE.SPACE.OVERLAPPING.WITH.REC.GENE.SPACE <- c() # for clarity
H5_rec.GIDS          <- c() # to be determined for each perturbation type
H5_rec.ALLC          <- c() # to be determined for each perturbation type
H5_rec.8C            <- c() # to be determined for each perturbation type
H5_rec.PERT          <- c() # to be determined for each perturbation type
H5_rec.G             <- c() # to be determined for each perturbation type
H5_therap.PERT.NAME  <- c() # to be determined for each therapeutic type
H5_therap.TUMOR.TYPE <- c() # to be determined for each therapeutic type

# get the gene space of the REC HDF5 data
tmpDir     <- myDelim(c(lincsDir, "trt_cp__6h"))
tmpPath    <- paste(tmpDir, "INFO.h5", sep="\\") # INFO.h5 path
rec.gids   <- h5read(tmpPath, "BACKGROUND.GIDS")
rec.gspace <- sapply(strsplit(rec.gids, "\\|"), "[")[1,] # gene symbols only

# get the overlapping genes within the tumor differential expression data
gat <- list()
de.gspace <- c()
de.gspace.overlapping.with.rec.gspace <- c()
for (tumor in TUMOR_TYPES) {
  sFN <- paste0(tumor, ".txt")
  sPath <- paste(sDir, sFN, sep="\\")
  data <- read.table(sPath, sep=" ", header=T, quote="\"")
  nms <- as.character(rownames(data))
  gat[[tumor]] <- nms
}

nms <- Reduce(union, gat)
nms <- nms[order(nms)]
nms <- c(nms, rec.gspace)
nms <- nms[!duplicated(nms)]
de.gspace <- nms

tmp1 <- de.gspace[de.gspace %in% rec.gspace]
tmp2 <- rec.gspace[rec.gspace %in% de.gspace]

if (TRUE) {
  stat <- table(c(tmp1, tmp2))
  bugs <- stat[stat!=2]
  dupl <- names(bugs)
}

de.gspace.overlapping.with.rec.gspace <- tmp2
TUMOR.DE.GENE.SPACE <- de.gspace
TUMOR.DE.GENE.SPACE.OVERLAPPING.WITH.REC.GENE.SPACE <- de.gspace.overlapping.with.rec.gspace
head(TUMOR.DE.GENE.SPACE)
head(TUMOR.DE.GENE.SPACE.OVERLAPPING.WITH.REC.GENE.SPACE)

# get the tumor differential expression data
gat <- list()
# fcs <- NULL
# pvs <- NULL
for (tumor in TUMOR_TYPES) {
  sFN <- paste0(tumor, ".txt")
  sPath <- paste(sDir, sFN, sep="\\")
  data <- read.table(sPath, sep=" ", header=T, quote="\"")
  tmp <- as.matrix(data)
  print(tmp)
  # dimnames(tmp) <- list(as.character(data[,1]), colnames(data[,-1]))
  
  x <- TUMOR.DE.GENE.SPACE.OVERLAPPING.WITH.REC.GENE.SPACE
  nms <- rownames(data)
  nms <- nms[nms %in% x]
  tmp <- tmp[nms,]
  idx <- which(abs(tmp[,1]) > TUMOR.THRESHOLD.LOG2.FC & tmp[,2] < TUMOR.THRESHOLD.P.VALUE)
  
  fc  <- tmp[idx,1]
  gat[[tumor]] <- fc
}
TUMOR.DE.LIST <- gat
sapply(TUMOR.DE.LIST, length)

##### iterate on PERT_TYPES
for (PERT_TYPE in PERT_TYPES[1:2]) {
  
  # get the abbreviation for the perturbation type
  pert.abbr <- PERT_ABBR[PERT_TYPE]
  
  # get paths
  recDir  <- myDelim(c(lincsDir, PERT_TYPE))
  recPath <- paste(recDir, "INFO.h5", sep="\\") # INFO.h5 path
  oFN.h5  <- paste(1, "way", pert.abbr, "therapeutic.h5", sep="-")
  oPath   <- paste(oDir, oFN.h5, sep="\\")
  
  # get REC HDF5 metadata
  h5ls(recPath)
  H5_rec.GIDS <- h5read(recPath, "BACKGROUND.GIDS")
  H5_rec.ALLC <- h5read(recPath, "CELL_TYPES"     )
  H5_rec.8C   <- h5read(recPath, "EIGHTCELLS"     )
  H5_rec.PERT <- h5read(recPath, "PERTURBAGENS"   )
  H5_rec.G    <- sapply(strsplit(H5_rec.GIDS, "\\|"), "[")[1,] # gene symbols only
  
  ### initializing THERAPEUTIC HDF5
  if (!file.exists(oPath)) {
    
    ### get informative perturbagens
    recName   <- ifelse(NCELL=="10C", "REC", "REC_8C")
    val       <- h5read(recPath, paste("0", recName, sep="/"), index=list(1, NULL))
    idx       <- which(!val %in% c("LTH", "NA"))
    pertNames <- H5_rec.PERT[idx]
    
    ### hdf5 file
    h5createFile(oPath)
    h5createGroup(oPath, "0")
    h5createDataset(oPath, "0/Sc", dims=c(length(pertNames), length(TUMOR_TYPES)),
                    storage.mode="double", chunk=c(1,1))
    h5createDataset(oPath, "0/Pv", dims=c(length(pertNames), length(TUMOR_TYPES)),
                    storage.mode="double", chunk=c(1,1))
    h5createGroup(oPath, "BH")
    h5createDataset(oPath, "BH/Pv", dims=c(length(pertNames), length(TUMOR_TYPES)),
                    storage.mode="double", chunk=c(1,1))
    h5write(pertNames, oPath, "PERTNAME")
    h5write(TUMOR_TYPES, oPath, "TUMORTYPE")
    h5write(1, oPath, "CURRENT.ROW")
    h5write(1, oPath, "CURRENT.COL")
  }
  
  # get therapeutic HDF5 metadata
  H5_therap.PERT.NAME   <- h5read(oPath, "PERTNAME")
  H5_therap.TUMOR.TYPE  <- h5read(oPath, "TUMORTYPE")
  H5_therap.current.row <- h5read(oPath, "CURRENT.ROW")
  H5_therap.current.col <- h5read(oPath, "CURRENT.COL")
  
  ### writing therapeutic HDF5
  LEN.COL <- length(H5_therap.TUMOR.TYPE)
  LEN.ROW <- length(H5_therap.PERT.NAME)
  while (H5_therap.current.col %in% 1:LEN.COL) {
    
    # get the duo of the perturgaben and the tissue descriptor
    pert.name  <- H5_therap.PERT.NAME[ H5_therap.current.row]
    tumor.type <- H5_therap.TUMOR.TYPE[H5_therap.current.col]
    cat(paste(as.character(Sys.time()), ": ",
              "Therapeutic: ", PERT_TYPE, " | ",
              "Tumor: ", tumor.type, " (", H5_therap.current.col, "/", LEN.COL, ") | ",
              "Pert: ", pert.name, " (", H5_therap.current.row, "/", LEN.ROW, ")",
              "\n", sep=""))
    
    # get the query 'signed' gene set of interest (QoI)
    fcs <- TUMOR.DE.LIST[[tumor.type]]
    tmp <- fcs
    # tmp <- sign(fcs)
    QoI <- tmp
    max.sc.for.QoI <- sum(abs(QoI)) * log2(REC.SCORE.CUTOFF[[PERT_TYPE]])
    
    # get REC data
    recName <- ifelse(NCELL=="10C", "REC", "REC_8C")
    idx.r   <- which(H5_rec.G %in% TUMOR.DE.GENE.SPACE)
    idx.c   <- which(H5_rec.PERT == pert.name)
    rec     <- h5read(recPath, paste("0" , recName, sep="/"), index=list(idx.r, idx.c))
    rec.bh  <- h5read(recPath, paste("BH", recName, sep="/"), index=list(idx.r, idx.c))
    rec     <- as.numeric(rec)
    rec.bh  <- as.numeric(rec.bh)
    rec.bh  <- 10^-abs(rec.bh)
    names(rec) <- names(rec.bh) <- H5_rec.G[idx.r]
    
    ### deriving therapeutic score
    valid.rec <- rec.bh < REC.FDR # noise reduction (REC FDR < 0.001)
    if (sum(valid.rec)==0) {
      
      writeAndUpdate()
      next
      
    } else {
      
      # compute the actual therapeutic score
      scs           <- log2(abs(rec)) * sign(rec) * valid.rec
      one.up.pos.sc <- {nms<-names(QoI)[QoI>0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp>0])),2)}
      one.up.neg.sc <- {nms<-names(QoI)[QoI>0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp<0])),2)}
      one.dn.pos.sc <- {nms<-names(QoI)[QoI<0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp<0])),2)}
      one.dn.neg.sc <- {nms<-names(QoI)[QoI<0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp>0])),2)}
      one.offset.sc <- one.up.neg.sc + one.dn.pos.sc - one.up.pos.sc - one.dn.neg.sc
      one.offset.sc <- one.offset.sc / max.sc.for.QoI ##### normalization
      # if the actual therapeutic score is zero
      if (one.offset.sc == 0 | is.na(one.offset.sc) == TRUE) {
        writeAndUpdate()
        next
      }
      
      # get permutation values
      permValues <- c()
      cnt <- 0
      while (cnt < Nperm) {
        names(scs) <- sample(names(scs))
        up.pos <- {nms<-names(QoI)[QoI>0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp>0])),2)}
        up.neg <- {nms<-names(QoI)[QoI>0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp<0])),2)}
        dn.pos <- {nms<-names(QoI)[QoI<0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp<0])),2)}
        dn.neg <- {nms<-names(QoI)[QoI<0]; tmp<-scs[nms]*QoI[nms]; round(abs(sum(tmp[tmp>0])),2)}
        # perm.v <- up.neg + dn.pos - up.pos - dn.neg
        perm.v <- up.neg + dn.pos - up.pos - dn.neg
        perm.v <- perm.v / max.sc.for.QoI 
        permValues <- c(permValues, perm.v)
        cnt <- cnt + 1
      }
      
      # get permutation P-value
      p.perm <- 1
      permValues <- abs(permValues) # using absolute values
      orderedPermValues <- permValues[order(permValues, decreasing=T)]
      x0 <- abs(one.offset.sc)
      p.perm.ecdf <- sum(permValues>=x0)/length(permValues) # without adding a pseudocount 1/Nperm
      if (p.perm.ecdf < P.perm_threshold) { # case where p.perm.gpd is computed
        
        # get one parameter
        mu <- mean(orderedPermValues[Nexc:(Nexc+1)]) # the location parameter ('t' in Bioinformatics, 2009)
        
        # model fitting
        dataForFit <- orderedPermValues[1:Nexc]
        gpd.object <- gpdFit(dataForFit, u=mu, type="pwm", information="expected")
        # gpd.object <- gpdFit(dataForFit, u=mu, type="mle", information="expected")
        
        # get the other two parameters
        xi <- gpd.object@fit$par.ests[1] # the shape parameter ('k' in Hosking and Wallis (1987))
        beta <- gpd.object@fit$par.ests[2] # the scale parameter ('a' in Hosking and Wallis (1987))
        
        # please check whether cdf.at.x0 == 1-(1+xi*(x0-mu)/beta)^(-1/xi)
        # note that pgpd.at.x0 = 1 - cdf.at.x0
        cdf.at.x0 <- pgpd(x0, xi, mu, beta) # at 'x0 - mu'
        p.perm.gpd <- as.numeric(Nexc/Nperm*(1-cdf.at.x0))
        
        # get permutation P-value
        p.perm <- p.perm.gpd
        
      } else { # case where p.perm.ecdf is sufficient
        
        # get permutation P-value
        p.perm <- p.perm.ecdf
        
      }
      
      # write and update
      writeAndUpdate(one.offset.sc, p.perm)
      
    } # end else
  } # end while
  
  # BH correction
  pvs    <- h5read(oPath, "0/Pv", index=list(NULL, NULL))
  pvs.bh <- AdjPvMtrx(pvs)
  h5write(pvs.bh, oPath, "BH/Pv", index=list(NULL, NULL))
}

# end PERT_TYPES
#################################################################################################################
# Output the results geneated by the above codes

DEBUG <- FALSE

# TO DO
PERT_TYPES               <- c("trt_cp__6h", "trt_cp__24h","trt_sh__96h")
PERT_ABBR                <- c("d6", "d24","sh96")
names(PERT_ABBR)         <- PERT_TYPES
#
z                        <- 2 # index for "term" and "pre" variable
term                     <- c("A", "B", "C")[z] # projects as elements
pre                      <- c("fct", "ft", "fc")[z] # The pre-name of input .txt data 
post                     <- c(1,2,3) # The post-name of DEGs lists.txt
TUMOR_TYPES              <- paste0(pre, "_", post)
TUMOR_TYPES_DESC         <- paste0("Cluster", post)
names(TUMOR_TYPES_DESC)  <- TUMOR_TYPES
#
THERAPEUTIC.P.CUTOFF     <- 0.05

# directory
myDelim <- function(names, delim = '\\', cr = T, after = F) {
  tmp <- ''; count <- 0
  for (name in names) {
    count <- count + 1
    if (count == 1) {tmp <- name} else {tmp <- paste(tmp, name, sep = delim)}
    if (cr == T) dir.create(tmp, showWarnings = F)
  }
  if (after) tmp <- paste(tmp, "\\", sep = "")
  return(tmp)
}

# P-value correction
AdjPvMtrx <- function(mx, method="BH") {
  tmp <- p.adjust(p=mx, method=method) # coerce the matrix col-by-col into a vector
  # NA allowed
  out <- matrix(tmp, nrow=nrow(mx), ncol=ncol(mx)) # coerce the vector back into a matrix
  dimnames(out) <- list(rownames(mx), colnames(mx))
  return(out)
}

# internal function
priorMxRowByCntThenMag <- function(Mx) {
  v <- apply(Mx, 1, function(X) sum(X != 0, na.rm=T)*1000000 + sum(abs(X), na.rm=T))
  o <- order(v, decreasing=T)
  out <- Mx[o,]
}

# internal function
priorMxRowByMagThenCnt <- function(Mx) {
  v <- apply(Mx, 1, function(X) sum(X != 0, na.rm=T)*1e-6 + sum(abs(X), na.rm=T))
  o <- order(v, decreasing=T)
  out <- Mx[o,]
}

# packages used
library(rhdf5)
library(stringr)

# global directories
sharedDir <- "D:\\REC"
jobDir    <- myDelim(c(sharedDir, "LabTask", "YourDirectory"))
sDir      <- myDelim(c(jobDir, term))
oDir      <- sDir

# file names
oFN.sc    <- "1-way-d6&d24-Table-Score.txt"
oFN.bh    <- "1-way-d6&d24-Table-Pval.txt"

# set paths
oPath.sc  <- paste(oDir, oFN.sc, sep="\\")
oPath.bh  <- paste(oDir, oFN.bh, sep="\\")

# global variables
THERAPEUTIC.SCORE    <- list()
THERAPEUTIC.PVALUE   <- list()
THERAPEUTIC.ALL.PERT <- list()

##### iterate on PERT_TYPES
for (PERT_TYPE in PERT_TYPES) {
  
  # get the abbreviation for the perturbation type
  pert.abbr <- PERT_ABBR[PERT_TYPE]
  
  
  # get paths
  sFN <- paste(1, "way", pert.abbr, "therapeutic.h5", sep="-")
  sPath <- paste(sDir, sFN, sep="\\")
  
  
  # get therapeutic HDF5 metadata
  H5_therap.PERT.NAME   <- h5read(sPath, "PERTNAME")
  H5_therap.TUMOR.TYPE  <- h5read(sPath, "TUMORTYPE")
  H5_therap.current.row <- h5read(sPath, "CURRENT.ROW")
  H5_therap.current.col <- h5read(sPath, "CURRENT.COL")
  
  
  # get therapeutic data
  scs <- h5read(sPath, "0/Sc", index=list(NULL, NULL))
  pvs.bh <- h5read(sPath, "BH/Pv", index=list(NULL, NULL))
  dimnames(pvs.bh) <- dimnames(scs) <- list(H5_therap.PERT.NAME, H5_therap.TUMOR.TYPE)
  
  
  # update
  THERAPEUTIC.SCORE[[pert.abbr]] <- scs
  THERAPEUTIC.PVALUE[[pert.abbr]] <- pvs.bh
  THERAPEUTIC.ALL.PERT[[pert.abbr]] <- H5_therap.PERT.NAME
}
# end PERT_TYPES

# preparing output matrix
mynms <- Reduce(union, THERAPEUTIC.ALL.PERT)
mynms <- mynms[order(mynms)]
ncol <- ncol(THERAPEUTIC.SCORE[[1]])*2
nrow <- length(mynms)
cnms <- rep(colnames(THERAPEUTIC.SCORE[[1]]), 2)

out.sc <- matrix(NA, nrow=nrow, ncol=ncol)
out.bh <- matrix(NA, nrow=nrow, ncol=ncol)
dimnames(out.sc) <- dimnames(out.bh) <- list(mynms, cnms)

for (i in 1:2) {
  nm <- rownames(THERAPEUTIC.SCORE[[i]])
  sel <- nm[nm %in% mynms]
  
  out.sc[sel, (1:3)+(i-1)*3] <- THERAPEUTIC.SCORE[[i]][sel,]
  out.bh[sel, (1:3)+(i-1)*3] <- THERAPEUTIC.PVALUE[[i]][sel,]
}

# reorder
e <- THERAPEUTIC.P.CUTOFF
tmp <- out.sc * (out.bh < e) * (out.sc > 0)
tmp <- priorMxRowByCntThenMag(tmp)
# tmp <- priorMxRowByMagThenCnt(tmp)
nms <- rownames(tmp)
out.sc <- out.sc[nms,]
out.bh <- out.bh[nms,]

# output the therapeutic score and BH.p-value results
write.table(out.sc, file=oPath.sc, sep="\t", quote=F, col.names=T, row.names=T)
write.table(out.bh, file=oPath.bh, sep="\t", quote=F, col.names=T, row.names=T)

#################################################################################################################
# Heatmap visualization

# Load packages
library(readxl)
library(xlsx)
library(ggplot2)
library(scales)
library(dichromat)
library(RColorBrewer)
library(dplyr)
library(circlize)
library(ComplexHeatmap)

# Load therapeutic score results
# The results can be acquired by the above codes
# We selected the results by considering the therapeutic score and significance
# and save the therapeutic score results as selected_results.xlsx, and BH.p-value results as select_pvalue.xlsx
data = read.xlsx2(file = "path/to/selected_results.xlsx", sheetIndex = 1, header = TRUE)

# Set the row name
# Remove the compounds with 0 sum value across datasets
rownames(data) = data[[1]]
data = data[,-1]
numeric_data <- apply(data, 2, as.numeric)
rowsum = rowSums(numeric_data, na.rm = TRUE)
list_name = c()
for (i in 1:length(rowsum)){
  if (rowsum[i] == 0){
    list_name = append(list_name, i)
  }
}
new_data = data[-list_name,]
my_data_frame <- as.data.frame(lapply(data, as.numeric))
rownames(my_data_frame) = rownames(data)

# Set the min and max value
min_value <- min(unlist(my_data_frame), na.rm = TRUE)
max_value <- max(unlist(my_data_frame), na.rm = TRUE)

# Load p-value results
pvalue_all = read.xlsx2(file = "path/to/select_pvalue.xlsx", sheetIndex = 1, header = TRUE)
rownames(pvalue_all) = pvalue_all[[1]]
pvalue_all = pvalue_all[,-1]
new_pvalue = pvalue_all[rownames(data), , drop = FALSE]

# Set the colors
col_fun = colorRamp2(c(min_value, max_value), c("white","#d6604d"))
col_fun(seq(-3, 3))

# Output the heatmap
pdf("path/to/output.pdf", width = 10, height = 8)
a = Heatmap(my_data_frame,
            rect_gp = gpar(col = "white", lwd = 2),
            col = col_fun,
            row_names_gp = gpar(fontsize = 14, fontfamily = "sans"),
            column_names_gp = gpar(fontsize = 14, fontfamily = "sans"),
            row_dend_width = unit(80, "mm"),
            name = "Therapeutic score",
            cell_fun = function(j, i, x, y, width, height, fill) {
              if(new_pvalue[i, j] < 0.05 & new_pvalue[i, j] >= 0.01)
                grid.text(sprintf("*"), x, y = y, just = "center", gp = gpar(fontsize = 10))
              if(new_pvalue[i, j] < 0.01 & new_pvalue[i, j] >= 0.001)
                grid.text(sprintf("**"), x, y = y, just = "center", gp = gpar(fontsize = 10))
              if(new_pvalue[i, j] < 0.001)
                grid.text(sprintf("***"), x, y = y, just = "center", gp = gpar(fontsize = 10))
            },
            column_title = "Selected inhibitory compounds",
            column_title_gp = gpar(fontsize = 20, fontface = "bold"),
            column_names_rot = 90,
            cluster_columns = FALSE)
draw(a)
dev.off()

# This is the end of the Transcriptomics Perturbation Drug Discovery
#################################################################################################################