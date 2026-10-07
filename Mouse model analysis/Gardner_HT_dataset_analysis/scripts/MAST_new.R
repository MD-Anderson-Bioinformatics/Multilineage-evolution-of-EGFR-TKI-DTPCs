require("MAST")
require("SingleCellExperiment", quietly=TRUE)
require("rhdf5")
require("data.table")
options(mc.cores=44)

library(Matrix)

# Define function to handle sparse input
# Function to import and format counts as sparse matrix (written with chatgpt-5)
csr_list_to_sparse <- function(lst, nrow = NULL, ncol = NULL) {
  x      <- as.numeric(lst[[1]])  # data
  j      <- as.integer(lst[[2]])  # column indices (0- or 1-based?)
  indptr <- as.integer(lst[[3]])  # row pointer, length = nrow + 1

  # infer nrow from indptr if not given
  if (is.null(nrow)) {
    nrow <- length(indptr) - 1L
  }

  # check if j is 0-based
  zero_based <- min(j) == 0L

  # build row indices from indptr
  # for each row r, entries are from indptr[r]:(indptr[r+1]-1)
  row_counts <- diff(indptr)
  i <- rep(seq_len(nrow) - 1L, row_counts)  # 0-based for now

  # shift to 1-based for R
  i <- i + 1L
  if (zero_based) {
    j <- j + 1L
  }

  # infer ncol if not provided
  if (is.null(ncol)) {
    ncol <- max(j)
  }

  as.matrix(sparseMatrix(i = i, j = j, x = x, dims = c(nrow, ncol)))
}

## Get infile from command line args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  # testing args
  args <- c("--infile", "/workdir/CIN/merged_MFP_pipeline_out/tmp_results/adata_cellTypes.h5ad", 
            "--outdir", "/workdir/CIN/merged_MFP_pipeline_out/Cancer_cells/highCIN_noCGAS_vs_others/MAST",  
            "--groups", "condition", "--comp-groups", "highCIN_noCGAS", "others", 
            "--subset", "Cell Type", "Cancer cells",
            "--numcores", "10")
}
cat("args = ", args, "\n")
w <- which(args == "--infile" | args == "-i")
if (length(w) == 1 && w < length(args)) {
    infile <- args[w+1]
} else {
    infile <- "/workdir/CIN/4T1_LUC911_pipeline_out/tmp_results/adata_clustered.h5ad"
}

w <- which(args == "--data_field")
if (length(w) == 1 && w < length(args)) {
    data_field <- args[w+1]
}

compGroups <- c()

# specify an extra string to append to end of output file
w <- which(args == "--outname-str")
outname_str <- ""
if (length(w) ==1 && w < length(args))
  outname_str <- args[w+1]

w <- which(args == "--numcores")
if (length(w) == 1 && w < length(args)) {
    numcores <- as.numeric(args[w+1])
    cat("Using ", numcores, " cores\n")
    options(mc.cores=numcores)
}


# if outfile already exists, just read it in and make sure NA values are fixed. If not,
# fix them. But if forceMast is TRUE, then always run MAST and overwrite old results.
forceMast <- FALSE
w <- which(args == "--force")
if (length(w) == 1) {
  forceMast <- TRUE
}

### Set output directory.
## note that R-studio cannot be run as root so when usin Rstudio, needs to be a writable path
## But when running from jupyter session within docker we will be root
w <- which(args == "--outdir" | args == "-o")
if (length(w) == 1 && w < length(args)) {
  outdir <- args[w+1]
} else {
  outdir <- sprintf("/opt/MAST/", sample)
}
system(sprintf("mkdir -p %s", outdir))
cat("infile =", infile, "\n")
cat("outdir = ", outdir, "\n")


w <- which(args == "--comp-groups")
if (length(w) == 1 && w < length(args)-1) {
    compGroups <- args[(w+1):(w+2)]
    cat("Comparing group ", compGroups[1], " to group ", compGroups[2], "\n")
} else {
    cat("Comparing each group to all others\n")
}

w <- which(args == "--comp-groups-name1")
cat("w=", w, " len(args)=", length(args), "\n")
if (length(w) ==1 && w < length(args)) {
    compGroups1Name <- args[w+1]
    cat("compGroups1Name = ", compGroups1Name, "\n")
} else if (length(compGroups) > 0) {
    cat("here\n")
    compGroups1Name <- compGroups[1]
}

# read "groups" from adata/obs/<group>
w <- which(args == "--groups" | args == "-g")
if (length(w) ==1 && w < length(args)) {
   groupName = args[w+1]
} else {
   groupName = "cluster"
}

subsetClass <- c()
w <- which(args == "--subset" | args == "-s")
if (length(w) ==1 && w < length(args)-1) {
    subsetClass <- c(args[w+1], args[w+2])
}

w <- which(args == "--hvg_filter")
if (length(w) == 1 && w < length(args)) {
  hvg_filter <- as.logical(args[w+1])
}

w <- which(args == "--sample_covariate")
if (length(w) == 1 && w < length(args)) {
  sample_covariate <- as.logical(args[w+1])
}

w <- which(args == "--covar_field")
if (length(w) == 1 && w < length(args)) {
  covar_field <- args[w+1]
}

w <- which(args == "--cell_filter")
if (length(w) == 1 && w < length(args)) {
  cell_filter <- args[w+1]
} 

w <- which(args == "--cell_filter_keep")
if (length(w) == 1 && w < length(args)) {
  cell_filter_keep <- args[w+1]
}

print(cell_filter)
print(cell_filter_keep)


## read in data
#h5ls(infile)

hvg <- as.logical(h5read(infile, "/var")$highly_variable)

counts <- h5read(infile, paste0("/layers/",data_field))#[hvg,]
counts <- csr_list_to_sparse(counts)
counts <- t(counts)

print(dim(counts))

cellNames <- as.character(h5read(infile, "/obs/_index"))
geneNames <- as.character(h5read(infile, "/var/_index"))
cellSample <- h5read(infile, sprintf("/obs/%s", covar_field))

if (cell_filter != FALSE) {
    cellFilterVAL <- h5read(infile, sprintf("/obs/%s", cell_filter))    
    cellFilterNames <- cellFilterVAL$categories
    cellFilterVAL <- cellFilterNames[as.numeric(cellFilterVAL$codes)+1]   
    cellFilterID <- (cellFilterVAL == cell_filter_keep)
}

#print(cellFilterVAL)
print(sum(cellFilterID))


tmp <- h5ls(infile)
w <- which(tmp$name == groupName)

if (length(w) == 1) {
   clusters <- h5read(infile, sprintf("/obs/%s", groupName))
} else if (length(w) == 2) {
   w2 <- which(tmp[w,]$dclass=="STRING")
   clusterNames <- h5read(infile, sprintf("%s/%s", tmp[w,][w2,"group"], tmp[w,][w2,"name"]))
   w2 <- which(tmp[w,]$dclass == "INTEGER")
   clusters <- clusterNames[as.numeric(h5read(infile, sprintf("%s/%s", tmp[w,][w2, "group"], tmp[w,][w2, "name"])))+1]
} else if (groupName == "geneExpression") {
    whichGene = which(geneNames == compGroups[1])
    if (length(whichGene) != 1) {
        if (compGroups[1] == "allGenes") {
            whichGene = 1
        } else {
            geneNamesSplit <- strsplit(compGroups[1], ",")[[1]]
            whichGene <- numeric(length(geneNamesSplit))
            for (i in 1:length(geneNamesSplit)) {
                tmp2 <- which(geneNames == geneNamesSplit[i])
                if (length(tmp2) != 1) {
                    stop("No gene ", geneNamesSplit[i], " in input data. Exiting")
                } else {
                    whichGene[i] <- tmp2
                }
            }
       }
     }
    cat("num gene = ", length(whichGene), "\n")
    cat("cutoff = ", compGroups[2], "\n")
    clusters <- ifelse(counts[whichGene[1],] > as.numeric(compGroups[2]), compGroups[1], sprintf("no %s", compGroups[1]))
    cat(whichGene[1], "\n")
    print(table(clusters))
    if (length(whichGene) > 1) {
        for (i in 2:length(whichGene)) {
            clusters <- ifelse(counts[whichGene[i],] > as.numeric(compGroups[2]), clusters, sprintf("no %s", compGroups[1]))
            cat(whichGene[i], "\n")
            print(table(clusters))
        }
    }
} else {
    cat("length(w)=", length(w), "\n")
    cat("groupName=", groupName, "\n")
    h5ls(infile)
    stop()
}

#print(clusters)
#print(cellNames)

if (length(clusters) != length(cellNames)) {
    
   clusterNames <- clusters$categories
   clusters <- clusterNames[as.numeric(clusters$codes)+1]
}

if (length(cellSample) != length(cellNames)) {
    
   sampleNames <- cellSample$categories
   cellSample <- sampleNames[as.numeric(cellSample$codes)+1]   
}

names(clusters) <- cellNames
names(cellSample) <- cellNames

if (length(compGroups) == 2) {
    if (compGroups[2] == 'others' || groupName == 'geneExpression') {
        keep <- rep(TRUE, length(cellNames))
    } else {
        keep <- (clusters == compGroups[1] | clusters == compGroups[2])
        clusters <- clusters[keep]
        counts <- counts[,keep]
        cellNames <- cellNames[keep]
        cellSample <- cellSample[keep]
        if (cell_filter != FALSE) {
            cellFilterID <- cellFilterID[keep]
        }
    }
    clusterList <- compGroups[1]
} else {
    clusterList <- unique(as.character(clusters))
    keep <- rep(TRUE, length(cellNames))
}

if (length(subsetClass) == 2) {
    w <- which(tmp$name == subsetClass[1])
    if (length(w) == 1) {
        groups <- h5read(infile, sprintf("/obs/%s", subsetClass[1]))
    } else if (length(w) == 2) {
        w2 <- which(tmp[w,]$dclass=="STRING")
        groupNames <- h5read(infile, sprintf("%s/%s", tmp[w,][w2,"group"], tmp[w,][w2,"name"]))
        w2 <- which(tmp[w,]$dclass == "INTEGER")
        groups <- groupNames[as.numeric(h5read(infile, sprintf("%s/%s", tmp[w,][w2, "group"], tmp[w,][w2, "name"])))+1]
    }
    groups <- groups[keep]
    keep2 <- (groups == subsetClass[2])
    clusters <- clusters[keep2]
    counts <- counts[,keep2]
    cellNames <- cellNames[keep2]
    cellSample <- cellSample[keep2]
    if (cell_filter != FALSE) {
            cellFilterID <- cellFilterID[keep2]
    }
    groups <- groups[keep2]
    cat("Clusters:\n")
    print(table(clusters))
    cat("Groups\n")
    print(table(groups))
}

if (groupName == "geneExpression" && compGroups[1] == "allGenes") {
    clusterList <- list()
    for (i in 1:length(geneNames)) {
        g <- geneNames[i]
        clusterList[[g]] <- ifelse(counts[i,] > as.numeric(compGroups[2]), g, sprintf("no %s", g))
    }
 }

if (hvg_filter == TRUE) {
    counts <- counts[hvg,]
    geneNames <- geneNames[hvg]
}

print(dim(counts))

if (cell_filter != FALSE) {
    counts <- counts[,cellFilterID]
    cellNames <- cellNames[cellFilterID]
    clusters <- clusters[cellFilterID]
    cellSample <- cellSample[cellFilterID]
}

print(dim(counts))

fixNAs <- function(fcHurdle, cluster, clusters, geneNames, counts) {
  pcol <- grep("Chisq", names(fcHurdle), value=TRUE)
  if (length(pcol) != 1) {
    cat("Error finding p-value column in fixNAs")
    return(fcHurdle)
  }
  fixNA <- (is.na(fcHurdle$coef) & fcHurdle[,pcol] < 1)
  coefRange <- range(fcHurdle$coef, na.rm=TRUE)
  geneOrder <- sapply(fcHurdle$primerid, function(x) {which(geneNames ==x)})
  fracIn <- apply(counts[geneOrder,clusters==cluster], 1, sum)/sum(clusters==cluster)
  fracOut <- apply(counts[geneOrder,clusters != cluster], 1, sum)/sum(clusters!=cluster)
  fcHurdle[fixNA & fracIn <= fracOut,"coef"] <- coefRange[1]
  fcHurdle[fixNA & fracIn > fracOut, "coef"] <- coefRange[2]
  fcHurdle
}


for (i in 1:length(clusterList)) {
    if (class(clusterList) == "list") {
        clusters <- clusterList[[i]]
        cluster <- names(clusterList)[i]
     } else {
         cluster <- clusterList[i]

     }
    clusterName <- sprintf("cluster%s", cluster)
  clusterAssign <- ifelse(clusters == cluster, cluster, "background")
  # get outfile name
  if (length(compGroups) == 2) {
    if (length(subsetClass) == 2) {
      if (groupName == 'geneExpression') {
        if (class(clusterList) == "list") {
          outfile <- sprintf("%s/%s_geneExp%s_%s%s.csv", outdir, clusterName, compGroups[2], gsub("/","_",subsetClass[2], fixed=TRUE), outname_str)
        } else {
          outfile <- sprintf("%s/%s_geneExp%s_%s%s.csv", outdir, compGroups1Name, compGroups[2], gsub("/","_",subsetClass[2], fixed=TRUE), outname_str)
        }
      } else {
        outfile <- sprintf("%s/%s_vs_%s_%s%s.csv", outdir, compGroups1Name, compGroups[2], gsub("/","_",subsetClass[2],fixed=TRUE), outname_str)
      }
    } else {
      outfile <- sprintf("%s/%s_vs_%s%s.csv", outdir, compGroups1Name, compGroups[2], outname_str)
    }
  } else {
    outfile <- sprintf("%s/%s%s.csv", outdir, clusterName, outname_str)
  }
  outfile <- gsub(" ", "_", outfile, fixed=TRUE)
  cat("outfile = \"", outfile, "\"\n")
  
  ## if outfile already exists, read it in and see if it has dealt with NA's properly
  
  if ((!forceMast) && file.exists(outfile)) {
    cat("outfile already exists, skipping\n")
    fcHurdle <- read.table(outfile, header=TRUE, sep=",", stringsAsFactors=FALSE)
    if (length(geneNames) == nrow(fcHurdle)) {
      fixNA <- (is.na(fcHurdle$coef) & fcHurdle$Pr..Chisq. < 1)
      if (sum(fixNA) ==0) next  ## this is already done
      fcHurdle <- fixNAs(fcHurdle, cluster, clusters, geneNames, counts)
      if (!file.exists(sprintf("%s.orig", outfile)))
        system(sprintf("mv %s %s.orig", outfile, outfile))
      write.csv(fcHurdle, outfile, quote=FALSE)
      next
    }
  }
  
  cat("Running MAST on ", clusterName, "\n")
    
  cat(length(cellNames), "\n")
  cat(length(clusterAssign), "\n")
  cat(length(cellSample), "\n")
    
  print('clusterAssign:')
  print(head(clusterAssign))
  cat("\n")
  print('covar assignments:')
  print(head(cellSample))
    
  sca <- FromMatrix(counts, 
                    cData=data.frame(wellKey=cellNames, barcode=cellNames, cluster=clusterAssign, sample=cellSample), 
                    fData=data.frame(primerid=geneNames, geneName=geneNames))
  tmpc <- factor(colData(sca)$cluster)
  tmpc <- relevel(tmpc, "background")
  colData(sca)$cluster <- tmpc
  colData(sca)$cngeneson <- scale(colSums(assay(sca) > 0))
  
  #print(colData(sca)$sample)
    
  if (sample_covariate == TRUE){
      cat("\n", "Running with", covar_field, "as Covariate", "\n")
      zlmResults <- zlm(~cluster + sample + cngeneson, sca)
  } else {
      zlmResults <- zlm(~cluster + cngeneson, sca)
  }
    
  cat("zlm function complete", "\n")
    
  zt <- summary(zlmResults, doLRT=clusterName)$datatable
  
  fcHurdle <- merge(zt[contrast==clusterName & component=='H',.(primerid, `Pr(>Chisq)`)],
                    zt[contrast==clusterName & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  setorder(fcHurdle, fdr)
  fcHurdle <- fixNAs(as.data.frame(fcHurdle), cluster, clusters, geneNames, counts)
  
  write.csv(fcHurdle, outfile, quote=FALSE)
}
