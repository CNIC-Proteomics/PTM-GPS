# import modules
library("limma")
library("logging")
library("yaml")

# Get get arguments
# configPath <- "D:\\ReportAnalysis\\2_ReportLimma_wo_GUI\\test\\test4\\limma_params_HFpEF_01.yaml"
args = commandArgs(trailingOnly=TRUE)
configPath <- args[1]


# set constants

obj <- list(
  "colNames"=c(), # First Line Header
  "colTypes"=c(), # Second Line Header
  
  "repIntegrations" = c(), # Which integrations are in the report
  "repSamples" = c(), # Which samples are in the report
  "Zcols" = c(), # Index of columns with Z values

  "repIntegrationSet" = c(), # set of integrations
  "repSampleSet" = c() # set of samples
  )


# define functions
readReport <- function(obj, config) {
  
  updatedObj = obj
  
  # read report data
  # updatedObj$repData <- read.csv(updatedObj$filePath, header=FALSE, sep="\t", skip=2)
  
  # read report metadata
  tmp <- read.csv(config$report_infile, header=FALSE, sep="\t", nrows=2)
  updatedObj$colNames <- unlist(tmp[1,])
  updatedObj$colTypes <- unlist(tmp[2,])
  
  # identify column indexes containing Z values
  index <- 1
  for (i in updatedObj$colNames) {
    if (i %in% config[['integrations']]) {
      updatedObj$Zcols <- c(updatedObj$Zcols, index)
      updatedObj$repIntegrations <- c(updatedObj$repIntegrations, i)
      updatedObj$repSamples <- c(updatedObj$repSamples, updatedObj$colTypes[index])
      
    }
    index <- index + 1 
  }
  # for (i in updatedObj$colTypes) {
  #   
  #   if (i != "LEVEL" & i != "REL" & i != "STATS" & i != "EXTRA") {
  #     updatedObj$Zcols <- append(updatedObj$Zcols, index)
  #     updatedObj$repIntegrations <- append(updatedObj$repIntegrations, 
  #                                       updatedObj$colNames[index])
  #     updatedObj$repSamples <- append(updatedObj$repSamples, i)
  #   }
  #  
  #   index <- index + 1 
  # }
  
  updatedObj$repIntegrationSet <- unique(updatedObj$repIntegrations)
  updatedObj$repSampleSet <- unique(updatedObj$repSamples)
  
  return (updatedObj)
}

LIMMA <- function(sampleGroups, Target, x, integration, eset, type) {
  
  f <- factor(Target, levels=names(sampleGroups))
  design <- model.matrix(~0+f)
  colnames(design) <- names(sampleGroups)
  
  if (!all(is.na(eset))) {
    
    fit <- lmFit(eset, design)
    
    #x <- gsub(" vs ", "-", obj$hypTesting)
    contrast.matrix <- makeContrasts(contrasts=x, levels=names(sampleGroups))
    
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    tmp <- as.data.frame(fit2$p.value)
    
    loginfo(paste0(integration, " - Prior Variance: ", fit2$s2.prior), logger="ReportLimma")
    loginfo(paste0(integration, " - Prior Degrees of Freedom: ", fit2$df.prior), logger="ReportLimma")
  
    } else {
      
    tmp <- data.frame(matrix(NA, nrow = nrow(eset), ncol = length(x)))
    colnames(tmp) <- x
    
  }
  
  newColname <- c()
  for (i in colnames(tmp)) {
    newColname <- append(newColname, paste(integration, type, i, sep="_"))
  }
  colnames(tmp) <- newColname
  

  
  return (tmp)
}

classicTTEST <- function(eset, Target, x, integration) {

  pvalues_ttest <- data.frame(row.names=1:nrow(eset))
  for (contrast in x) {
    g <- strsplit(contrast, '-')[[1]]
    g1_bool <- g[1] == Target
    g2_bool <- g[2] == Target
    pvalues <- apply(eset, 1, function(y) {
      if (sum(!is.na(y[g1_bool])) < 2 | sum(!is.na(y[g2_bool])) < 2) return (NA)
      return (t.test(x=y[g1_bool], y=y[g2_bool], alternative="two.sided", var.equal=TRUE)$p.value)
    })

    colname <- append(colnames(pvalues_ttest), paste(integration, "ttest", contrast, sep="_"))
    pvalues_ttest <- cbind(pvalues_ttest, pvalues)
    colnames(pvalues_ttest) <- colname
    
  }
  return(pvalues_ttest)
}

MeanDiff <- function(eset, Target, x, integration) {
  
  mean_diff <- data.frame(row.names=1:nrow(eset))
  for (contrast in x) {
    g <- strsplit(contrast, '-')[[1]]
    g1_bool <- g[1] == Target
    g2_bool <- g[2] == Target
    mean1 <- rowMeans(eset[, g1_bool], na.rm = TRUE)
    mean2 <- rowMeans(eset[, g2_bool], na.rm = TRUE)
    mean_diff_tmp <- mean1-mean2

    colname <- append(colnames(mean_diff), paste(integration, "dX", contrast, sep="_"))
    mean_diff <- cbind(mean_diff, mean_diff_tmp)
    colnames(mean_diff) <- colname
    
  }
  return(mean_diff)
}


calculatePvalues <- function(obj, config, sampleGroups) {
  
  reportData <- read.csv(config$report_infile, header=FALSE, sep="\t", skip=2)
  
  # data frame containing pvalues
  pvalues_df <- data.frame(row.names = 1:nrow(reportData))
  
  # second line header
  subHeader <- c()
  
  index <- 0
  for (integration in config$integrations) {
    index <- index +1 
    
    loginfo(paste0("Doing calculations: ", integration), logger="ReportLimma")
    
    # Get vector with low level
    lowLevel1 <- config[["lowLevel1"]][index]
    lowLevel2 <- config[["lowLevel2"]][index]
    # lowLevel <- strsplit(gsub("Z_", "", integration), "2")[[1]][1]
    lowLevelCol <- as.vector(
      unlist(reportData[obj$colTypes == lowLevel2 & obj$colNames == lowLevel1][1])
    )
    pvalues_df_tmp <- data.frame(row.names = 1:nrow(reportData))
    pvalues_df_tmp$low <- lowLevelCol
    
    # dataframe containing working Z
    eset <- data.frame(row.names = 1:nrow(reportData))
    
    Target <- c()
    
    for (group in names(sampleGroups)) {
      
      for (sample in sampleGroups[[group]]) {
        
        integrationBool <- integration == obj$repIntegrations
        sampleBool <- sample == obj$repSamples
        
        zColBool <- integrationBool & sampleBool
        zColIndex <- obj$Zcols[zColBool]
        
        eset <- cbind(eset, reportData[zColIndex])
        
        Target <- append(Target, group)
      }
      
    }
    

    # LIMMA
    x <- c()
    for (i in config$hypothesis_testing) {
      x <- c(x, paste0(i[1], '-', i[2]))
    }
    #x <- gsub(" vs ", "-", obj$hypTesting)
    
    #Mean Difference
    loginfo(paste0(integration, " - Calculating Mean difference"), logger="ReportLimma")

    mean_diff <- MeanDiff(eset, Target, x, integration)
    mean_sign <- ifelse(mean_diff > 0, 1, -1)
    subHeader <- c(subHeader, rep('dX', length(x)))
    
    
    if ('limma' %in% config$test_type) { # removing duplicates
      loginfo(paste0(integration, " - Applying LIMMA"), logger="ReportLimma")
      
      lowLevelSet_bool <- !duplicated(lowLevelCol)
      tmp <- LIMMA(sampleGroups, Target, x, integration, eset[lowLevelSet_bool,], "limma")
      tmp_log <- -log10(tmp)*mean_sign[lowLevelSet_bool,]
      names(tmp_log) <- gsub("_limma_", "_logLimma_", names(tmp))
      tmp <- cbind(tmp, tmp_log)
      
      lowLevelSet <- lowLevelCol[lowLevelSet_bool]
      tmp$low <- lowLevelSet
      
      # merge is changing the order... and we MUST preserve it
      pvalues_df_tmp$index <- 1:nrow(pvalues_df_tmp)
      pvalues_df_tmp <- merge(pvalues_df_tmp, tmp, by="low", sort=FALSE) # still changes order
      pvalues_df_tmp <- pvalues_df_tmp[order(pvalues_df_tmp$index),]
      pvalues_df_tmp <- pvalues_df_tmp[, !(names(pvalues_df_tmp) %in% c("index"))]
      
      
      subHeader <- c(subHeader, c(rep('pvalue', length(x)), rep('LPS', length(x))))
    }
    
    
    if ( 'limma_with_duplicates' %in% config$test_type) {
      loginfo(paste0(integration, " - Applying LIMMA with duplicates"), logger="ReportLimma")
      tmp <- LIMMA(sampleGroups, Target, x, integration, eset, "limmaDup")
      tmp_log <- -log10(tmp)*mean_sign
      names(tmp_log) <- gsub("_limmaDup_", "_logLimmaDup_", names(tmp))
      tmp <- cbind(tmp, tmp_log)
      pvalues_df_tmp <- cbind(pvalues_df_tmp, tmp)
      
      subHeader <- c(subHeader, c(rep('pvalue', length(x)), rep('LPS', length(x))))
      
    }
    
    #TTEST
    if ( 't-test' %in% config$test_type) {
      loginfo(paste0(integration, " - Calculating t-test"), logger="ReportLimma")
      pvalues_ttest <- classicTTEST(eset, Target, x, integration)
      tmp_log <- -log10(pvalues_ttest)*mean_sign
      names(tmp_log) <- gsub("_ttest_", "_logTtest_", names(pvalues_ttest))
      pvalues_ttest <- cbind(pvalues_ttest, tmp_log)
      pvalues_df_tmp <- cbind(pvalues_df_tmp, pvalues_ttest)
      
      subHeader <- c(subHeader, c(rep('pvalue', length(x)), rep('LPS', length(x))))
    }
    

    # Add data to output dataframe
    pvalues_df <- cbind(pvalues_df, mean_diff, pvalues_df_tmp[-1])
    
  }
  
  pvalues_df[sapply(pvalues_df, is.nan)] <- NA
  
  return (
    list(
      reportData=reportData, 
      pvalues_df=pvalues_df,
      subHeader=subHeader
    )
  )
}
  
writeOutputReport <- function (config, reportData, pvalues_df, subHeader) {
    
  loginfo("Generating output table...", logger="ReportLimma")

  header <- read.csv(config$report_infile, header=FALSE, sep="\t", nrows=2)
  
  mainHeader <- colnames(pvalues_df)
  
  for (i in 1:length(mainHeader)) { #colnames(pvalues_df)) {
    header <- cbind(header, c(mainHeader[i], subHeader[i]))
  }
  
  reportData <- cbind(reportData, pvalues_df)
  reportData <- data.frame(mapply('c', header, reportData))
  
  #outDir <- dirname(config$report_infile) # dirname changes \\ --> /, yielding error in \\tierra...
  #outDir <- gsub("/", "\\\\", outDir) # fix it //tierra...
  #outFile <- paste("LIMMA", basename(obj$filePath), sep="_")
  #outPath <- paste(outDir, outFile, sep="\\")
  
  write.table(reportData, file = config$outfile, quote = F, sep = "\t", row.names = F,
              col.names = F, na="")
  
  loginfo(paste0("Output table was written: ", config$outfile), logger="ReportLimma")
}


#
# MAIN
#

# Read YAML file
config <- read_yaml(configPath)

config[['integrations']] <- c()
config[['lowLevel1']] <- c()
config[['lowLevel2']] <- c()

for (i in config[['ColumnNames']]) {
  config[['integrations']] <- c(config[['integrations']], i[[2]])
  config[['lowLevel1']] <- c(config[['lowLevel1']], i[[1]][1])
  config[['lowLevel2']] <- c(config[['lowLevel2']], i[[1]][2])
}

# Read with samples and parse it
sampleTable <- read.csv(config$sample_table, sep='\t', colClasses="character")
sampleGroups <- as.list(sampleTable)

for (i in names(sampleGroups)) {
  g <- c()
  for (j in sampleGroups[[i]]) {
    if (j!="") {
      g <- c(g, j)
    }
  }
  sampleGroups[[i]] <- g
}

# Set Logging file
logFile <- paste0(strsplit(config$report_infile, split='[^.]+$', perl=T), 'log')
basicConfig()
addHandler(writeToFile, logger='ReportLimma', file=logFile)

# Get Header information from report
loginfo("Reading report header...", logger="ReportLimma")
obj = readReport(obj, config)


# Calculate Limma pvalues
outputList <- calculatePvalues(obj, config, sampleGroups)

# Write output
writeOutputReport(config, outputList[['reportData']], outputList[['pvalues_df']], outputList[['subHeader']])
