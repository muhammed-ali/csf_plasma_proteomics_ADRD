rm(list = ls())
library(data.table)
library(readxl)
library(dplyr)
library(EnhancedVolcano)
options(stringsAsFactors = FALSE)
options(width = 160)
set.seed(1)


## AD vs. CO ##

load("CSF_Complete_Expr_MAP_PPMI_SADRC_N4685.RData")
dim(CSF_complete_expr_pcs) # 4685 7663
CSF_complete_expr_pcs[1:4,1:12]
table(CSF_complete_expr_pcs$status_for_analysis)
#       AD    CO_All  CO_ATneg       DLB   Exclude       FTD       MCI    MCI-PD        PD Prodromal   Unknown
#     1157       818       726        37       799        46        21         3       739       332         1
AD_expr <- CSF_complete_expr_pcs[CSF_complete_expr_pcs$status_for_analysis %in% c("AD", "CO_ATneg"),]
dim(AD_expr) # 1883 7663
table(AD_expr$status_for_analysis)
#      AD CO_ATneg
#    1157      726
AD_expr$level <- ifelse(AD_expr$status_for_analysis == "AD", 1, 0)
table(AD_expr$level)
#   0    1
# 726 1157
round((nrow(AD_expr[AD_expr$Sex == "Female",])/nrow(AD_expr))*100, 1) # 51.8
round(mean(AD_expr$Age, na.rm=T), 1) # 71.3
round(sd(AD_expr$Age, na.rm=T), 1) # 8.8

# check discovery statistics
discovery_stats <- AD_expr[,c("UniquePhenoID", "Sex", "status_for_analysis", "Cohort")]

### get analytes
AD_expr <- AD_expr %>% select(UniquePhenoID, DrawDate, Cohort, Age, Sex, AT_class, status_for_analysis, level, PlateId, PC1, PC2, PC3, everything())
dim(AD_expr) # 1883 7664
AD_expr[1:4,c(1:13,7660:7664)]
# replace the NaN and Inf values in our data frame
AD_expr[is.na(AD_expr) | AD_expr == "Inf"] <- NA
# Remove features that do not have 2 values for "level" after NA removal
# Function to check if each feature column has at least two levels for each covariate
check_levels <- function(df, feature_col, covariates) {
  for (covariate in covariates) {
    non_na_data <- df[!is.na(df[[feature_col]]), c(covariate, feature_col)]
    if (length(unique(non_na_data[[covariate]])) < 2) {
      return(FALSE)
    }
  }
  return(TRUE)
}
covariates <- c("level", "Sex", "Cohort")
# Identify valid features that have at least two levels for each covariate
valid_features <- sapply(names(AD_expr)[13:ncol(AD_expr)], function(feature) {
  check_levels(AD_expr, feature, covariates)
})
table(valid_features)[2] # 7099
# Filter the dataframe to include only valid features and the covariate columns
AD_expr <- AD_expr[, c(names(AD_expr)[1:12], names(AD_expr)[13:ncol(AD_expr)][valid_features])]
dim(AD_expr) # 1883 7111

analyte <- as.data.frame(AD_expr[,13:ncol(AD_expr)])
dim(analyte) #  1883 7099
analyte_name <- as.data.frame(colnames(analyte))
head(analyte_name)

AD_vs_CO_daa <- data.frame()
# run linear regression model for AD vs CO
for (i in 1:ncol(analyte)) {
  analyte_ID <- analyte_name[i,1]
  protein <- analyte[,i]
  model <- lm(protein ~ as.factor(level) + as.numeric(Age) + as.factor(Sex) + as.factor(PlateId) + as.numeric(PC1) + as.numeric(PC2), data = AD_expr)
  output <- as.data.frame(summary(model)[4])
  result <- cbind(as.character(analyte_ID), output[2,1], output[2,2], output[2,4])
  print(i)
  AD_vs_CO_daa <- rbind(AD_vs_CO_daa, result) 
}
colnames(AD_vs_CO_daa) <- c("Analyte","Estimate", "Standard_error","Pvalue")
dim(AD_vs_CO_daa) # 7099    4
any(is.na(AD_vs_CO_daa)) # FALSE
head(AD_vs_CO_daa, 2)
#            Analyte          Estimate     Standard_error               Pvalue
# 1 X10000.28_P43320 0.032988262698921 0.0521020761507962    0.526721420165982
# 2  X10001.7_P04049 0.324627170610582 0.0291977922094886 8.13080212750864e-28
save(AD_expr, AD_vs_CO_daa, file="AD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")

# multiple test correction
cols <- c(2,3,4)    
AD_vs_CO_daa[,cols] <- apply(AD_vs_CO_daa[,cols], 2, function(x) as.numeric(as.character(x)))
AD_vs_CO_daa$FDR <- p.adjust(AD_vs_CO_daa$Pvalue, method = "fdr", n = nrow(AD_vs_CO_daa))
AD_vs_CO_daa$Bonferroni <- p.adjust(AD_vs_CO_daa$Pvalue, method = "bonferroni", n = nrow(AD_vs_CO_daa))
nrow(AD_vs_CO_daa[AD_vs_CO_daa$Pvalue < 0.05,]) # 4638
nrow(AD_vs_CO_daa[AD_vs_CO_daa$FDR < 0.05,]) # 4425
nrow(AD_vs_CO_daa[AD_vs_CO_daa$Bonferroni < 0.05,]) # 2278

# Analyte annotation
load("CSF_Soma7K_LongsGarfield_Expr_Pheno_Annot.RData")
dim(longsGar_annot) # 7291   16
longsGar_annot <- longsGar_annot[,c("Analytes", "Symbol", "UniProt", "TargetFullName", "External_ID")]
longsGar_annot$Symbol <- sapply(strsplit(as.character(longsGar_annot$Symbol), "\\."), '[', 1)
table(AD_vs_CO_daa$Analyte %in% longsGar_annot$External_ID)
# FALSE  TRUE
#     1  7098
AD_vs_CO_MAP <- inner_join(AD_vs_CO_daa, longsGar_annot, by=c("Analyte"="External_ID"))
dim(AD_vs_CO_MAP) # 7098   10
AD_vs_CO_daa <- anti_join(AD_vs_CO_daa, longsGar_annot, by=c("Analyte"="External_ID"))
dim(AD_vs_CO_daa) # 1 6
load("CSF_Soma5K_StanfordADRC_Expr_Pheno_Annot.RData")
dim(sadrc_annot) # 4735   11
sadrc_annot <- sadrc_annot[,c("Analytes", "Symbol", "UniProt", "TargetFullName", "External_ID")]
sadrc_annot$Symbol <- sapply(strsplit(as.character(sadrc_annot$Symbol), "\\."), '[', 1)
table(AD_vs_CO_daa$Analyte %in% sadrc_annot$External_ID)
# TRUE
#    1
AD_vs_CO_SADRC <- inner_join(AD_vs_CO_daa, sadrc_annot, by=c("Analyte"="External_ID"))
dim(AD_vs_CO_SADRC) # 1 10
AD_vs_CO_daa <- rbind(AD_vs_CO_MAP, AD_vs_CO_SADRC)
AD_vs_CO_daa <- AD_vs_CO_daa[order(AD_vs_CO_daa$Pvalue),]
dim(AD_vs_CO_daa) # 7099   10
table(AD_vs_CO_daa$Analyte %in% names(AD_expr))
# TRUE
# 7099
AD_vs_CO_daa <- AD_vs_CO_daa[,c("Analyte", "UniProt", "Symbol", "TargetFullName", "Estimate", "Standard_error", "Pvalue", "FDR", "Bonferroni")]
colnames(AD_vs_CO_daa) <- c("Analyte", "UniProt", "Symbol", "Target", "Estimate", "SE", "P", "FDR", "Bonferroni")
write.table(AD_vs_CO_daa, file="CSF_LongsGar_PPMI_SADRC_Zscore_AD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", row.names=F, quote=F)

# volcano plot
AD <- AD_vs_CO_daa
AD <- AD[!duplicated(AD$Symbol),]
top_10_labels <- unique(c(AD[AD$Estimate > 0,]$Symbol[1:7], AD[AD$Estimate < 0,]$Symbol[1:5]))
png("Volcano_CSF_LongsGar_PPMI_SADRC_Zscore_AD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes_v2.png", units="mm", width=200, height=170, res=1000)
EnhancedVolcano(AD,
    lab = AD$Symbol,
    selectLab = top_10_labels,
    title= 'CSF AD (n=1,157) vs. CO (n=726)', # NULL
    subtitle = 'Covariates: Age, Sex, Plate, PC1, PC2', # NULL
    x = 'Estimate',
    y = 'FDR',
    pCutoff = 0.05,
    FCcutoff = 0, # 0.05
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2,
    labSize = 5,
    colAlpha = 1/5,
    labFace = 'bold',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = "", # to remove "total variables = xxxx" at the bottom of plot
    ylab = "FDR", xlab = "Estimate",
    legendLabels = c("NS", "NS", "FDR", "FDR"),
    ylim = c(0, max(-log10(AD[,"FDR"]), na.rm = TRUE)),
    xlim = c(min(AD$Estimate, na.rm = TRUE), max(AD$Estimate, na.rm = TRUE)))
dev.off()


## DLB vs. CO ##

DLB_expr <- CSF_complete_expr_pcs[CSF_complete_expr_pcs$status_for_analysis %in% c("DLB", "CO_ATneg"),]
dim(DLB_expr) # 763 7663
table(DLB_expr$status_for_analysis)
# CO_ATneg      DLB
#      726       37
DLB_expr$level <- ifelse(DLB_expr$status_for_analysis == "DLB", 1, 0)
table(DLB_expr$level)
#   0    1
# 726   37
round((nrow(DLB_expr[DLB_expr$Sex == "Female",])/nrow(DLB_expr))*100, 1) # 52.2
round(mean(DLB_expr$Age, na.rm=T), 1) # 67.7
round(sd(DLB_expr$Age, na.rm=T), 1) # 9.2
# check discovery statistics
discovery_stats <- DLB_expr[,c("UniquePhenoID", "Sex", "status_for_analysis", "Cohort")]
### get analytes
DLB_expr <- DLB_expr %>% select(UniquePhenoID, DrawDate, Cohort, Age, Sex, AT_class, status_for_analysis, level, PlateId, PC1, PC2, PC3, everything())
dim(DLB_expr) # 763 7664
DLB_expr[1:4,c(1:13,7660:7664)]
# replace the NaN and Inf values in our data frame
DLB_expr[is.na(DLB_expr) | DLB_expr == "Inf"] <- NA
# Remove features that do not have 2 values for "level" after NA removal
# Function to check if each feature column has at least two levels for each covariate
check_levels <- function(df, feature_col, covariates) {
  for (covariate in covariates) {
    non_na_data <- df[!is.na(df[[feature_col]]), c(covariate, feature_col)]
    if (length(unique(non_na_data[[covariate]])) < 2) {
      return(FALSE)
    }
  }
  return(TRUE)
}
covariates <- c("level", "Sex", "Cohort")
# Identify valid features that have at least two levels for each covariate
valid_features <- sapply(names(DLB_expr)[13:ncol(DLB_expr)], function(feature) {
  check_levels(DLB_expr, feature, covariates)
})
table(valid_features)[2] # 7099
# Filter the dataframe to include only valid features and the covariate columns
DLB_expr <- DLB_expr[, c(names(DLB_expr)[1:12], names(DLB_expr)[13:ncol(DLB_expr)][valid_features])]
dim(DLB_expr) # 763 7111
which(names(DLB_expr) %in% "XIL23A.10365_P29460, Q9NPF7")
analyte <- as.data.frame(DLB_expr [,13:ncol(DLB_expr)])
dim(analyte) #  763 7099
analyte_name <- as.data.frame(colnames(analyte))
head(analyte_name)

DLB_vs_CO_daa <- data.frame()
# run model (for DLB vs CO)
for (i in 1:ncol(analyte)) {
  analyte_ID <- analyte_name[i,1]
  protein <- analyte[,i]
  model <- lm(protein ~ as.factor(level) + as.numeric(Age) + as.factor(Sex) + as.factor(PlateId) + as.numeric(PC1) + as.numeric(PC2), data = DLB_expr)
  output <- as.data.frame(summary(model)[4])
  result <- cbind(as.character(analyte_ID), output[2,1], output[2,2], output[2,4])
  print(i)
  DLB_vs_CO_daa <- rbind(DLB_vs_CO_daa, result) 
}
colnames(DLB_vs_CO_daa) <- c("Analyte","Estimate", "Standard_error","Pvalue")
dim(DLB_vs_CO_daa) # 7099    4
any(is.na(DLB_vs_CO_daa)) # FALSE
head(DLB_vs_CO_daa, 2)
#            Analyte            Estimate    Standard_error            Pvalue
# 1 X10000.28_P43320  -0.343834043928747 0.333117682072407 0.302363857725853
# 2  X10001.7_P04049 -0.0625789289094047 0.192674413017493 0.745435282131122
save(DLB_expr, DLB_vs_CO_daa, file="DLB_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")

# multiple test correction
cols <- c(2,3,4)    
DLB_vs_CO_daa[,cols] <- apply(DLB_vs_CO_daa[,cols], 2, function(x) as.numeric(as.character(x)))
DLB_vs_CO_daa$FDR <- p.adjust(DLB_vs_CO_daa$Pvalue, method = "fdr", n = nrow(DLB_vs_CO_daa))
DLB_vs_CO_daa$Bonferroni <- p.adjust(DLB_vs_CO_daa$Pvalue, method = "bonferroni", n = nrow(DLB_vs_CO_daa))
nrow(DLB_vs_CO_daa[DLB_vs_CO_daa$Pvalue < 0.05,]) # 3396
nrow(DLB_vs_CO_daa[DLB_vs_CO_daa$FDR < 0.05,]) # 2862
nrow(DLB_vs_CO_daa[DLB_vs_CO_daa$Bonferroni < 0.05,]) # 790

# Analyte annotation
load("CSF_Soma7K_LongsGarfield_Expr_Pheno_Annot.RData")
dim(longsGar_annot) # 7291   16
longsGar_annot <- longsGar_annot[,c("Analytes", "Symbol", "UniProt", "TargetFullName", "External_ID")]
longsGar_annot$Symbol <- sapply(strsplit(as.character(longsGar_annot$Symbol), "\\."), '[', 1)
table(DLB_vs_CO_daa$Analyte %in% longsGar_annot$External_ID)
# FALSE  TRUE
#     1  7098
DLB_vs_CO_MAP <- inner_join(DLB_vs_CO_daa, longsGar_annot, by=c("Analyte"="External_ID"))
dim(DLB_vs_CO_MAP) # 7098   10
DLB_vs_CO_daa <- anti_join(DLB_vs_CO_daa, longsGar_annot, by=c("Analyte"="External_ID"))
dim(DLB_vs_CO_daa) # 1 6
load("CSF_Soma5K_StanfordADRC_Expr_Pheno_Annot.RData")
dim(sadrc_annot) # 4735   11
sadrc_annot <- sadrc_annot[,c("Analytes", "Symbol", "UniProt", "TargetFullName", "External_ID")]
sadrc_annot$Symbol <- sapply(strsplit(as.character(sadrc_annot$Symbol), "\\."), '[', 1)
table(DLB_vs_CO_daa$Analyte %in% sadrc_annot$External_ID)
# TRUE
#    1
DLB_vs_CO_SADRC <- inner_join(DLB_vs_CO_daa, sadrc_annot, by=c("Analyte"="External_ID"))
dim(DLB_vs_CO_SADRC) # 1 10
DLB_vs_CO_daa <- rbind(DLB_vs_CO_MAP, DLB_vs_CO_SADRC)
DLB_vs_CO_daa <- DLB_vs_CO_daa[order(DLB_vs_CO_daa$Pvalue),]
dim(DLB_vs_CO_daa) # 7099   10
table(DLB_vs_CO_daa$Analyte %in% names(DLB_expr))
# TRUE
# 7099
DLB_vs_CO_daa <- DLB_vs_CO_daa[,c("Analyte", "UniProt", "Symbol", "TargetFullName", "Estimate", "Standard_error", "Pvalue", "FDR", "Bonferroni")]
colnames(DLB_vs_CO_daa) <- c("Analyte", "UniProt", "Symbol", "Target", "Estimate", "SE", "P", "FDR", "Bonferroni")
write.table(DLB_vs_CO_daa, file="CSF_LongsGar_PPMI_SADRC_Zscore_DLB_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", row.names=F, quote=F)

# volcano plot
DLB_vs_CO_daa <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_DLB_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, quote="")
DLB <- DLB_vs_CO_daa
DLB <- DLB[!duplicated(DLB$Symbol),]
top_10_labels <- unique(c(DLB[DLB$Estimate > 0,]$Symbol[1:7], DLB[DLB$Estimate < 0,]$Symbol[1:5]))
png("Volcano_CSF_LongsGar_PPMI_SADRC_Zscore_DLB_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes_v2.png", units="mm", width=200, height=170, res=1000)
EnhancedVolcano(DLB,
    lab = DLB$Symbol,
    selectLab = top_10_labels,
    title= 'CSF DLB (n=37) vs. CO (n=726)', # NULL
    subtitle = 'Covariates: Age, Sex, Plate, PC1, PC2', # NULL
    x = 'Estimate',
    y = 'FDR',
    pCutoff = 0.05,
    FCcutoff = 0, # 0.05
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2,
    labSize = 5,
    colAlpha = 1/5,
    labFace = 'bold',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = "", # to remove "total variables = xxxx" at the bottom of plot
    ylab = "FDR", xlab = "Estimate",
    legendLabels = c("NS", "NS", "FDR", "FDR"),
    ylim = c(0, max(-log10(DLB[,"FDR"]), na.rm = TRUE)),
    xlim = c(min(DLB$Estimate, na.rm = TRUE), max(DLB$Estimate, na.rm = TRUE)))
dev.off()


## FTD vs. CO ##

FTD_expr <- CSF_complete_expr_pcs[CSF_complete_expr_pcs$status_for_analysis %in% c("FTD", "CO_ATneg"),]
dim(FTD_expr) # 772 7663
table(FTD_expr$status_for_analysis)
# CO_ATneg      FTD
#      726       46
FTD_expr$level <- ifelse(FTD_expr$status_for_analysis == "FTD", 1, 0)
table(FTD_expr$level)
#   0    1
# 726   46
round((nrow(FTD_expr[FTD_expr$Sex == "Female",])/nrow(FTD_expr))*100, 1) # 52.5
round(mean(FTD_expr$Age, na.rm=T), 1) # 67.5
round(sd(FTD_expr$Age, na.rm=T), 1) # 9.2
# check discovery statistics
discovery_stats <- FTD_expr[,c("UniquePhenoID", "Sex", "status_for_analysis", "Cohort")]
### get analytes
FTD_expr <- FTD_expr %>% select(UniquePhenoID, DrawDate, Cohort, Age, Sex, AT_class, status_for_analysis, level, PlateId, PC1, PC2, PC3, everything())
dim(FTD_expr) # 772 7664
FTD_expr[1:4,c(1:13,7660:7664)]
# replace the NaN and Inf values in our data frame
FTD_expr[is.na(FTD_expr) | FTD_expr == "Inf"] <- NA
# Remove features that do not have 2 values for "level" after NA removal
# Function to check if each feature column has at least two levels for each covariate
check_levels <- function(df, feature_col, covariates) {
  for (covariate in covariates) {
    non_na_data <- df[!is.na(df[[feature_col]]), c(covariate, feature_col)]
    if (length(unique(non_na_data[[covariate]])) < 2) {
      return(FALSE)
    }
  }
  return(TRUE)
}
covariates <- c("level", "Sex", "Cohort")
# Identify valid features that have at least two levels for each covariate
valid_features <- sapply(names(FTD_expr)[13:ncol(FTD_expr)], function(feature) {
  check_levels(FTD_expr, feature, covariates)
})
table(valid_features)[2] # 7008
# Filter the dataframe to include only valid features and the covariate columns
FTD_expr <- FTD_expr[, c(names(FTD_expr)[1:12], names(FTD_expr)[13:ncol(FTD_expr)][valid_features])]
dim(FTD_expr) # 772 7020
which(names(FTD_expr) %in% "X9921.14_Q9UKA2")
analyte <- as.data.frame(FTD_expr[,13:ncol(FTD_expr)])
dim(analyte) #  772 7008
analyte_name <- as.data.frame(colnames(analyte))
head(analyte_name)

FTD_vs_CO_daa <- data.frame()
# run model (for FTD vs CO)
for (i in 1:ncol(analyte)) {
  analyte_ID <- analyte_name[i,1]
  protein <- analyte[,i]
  model <- lm(protein ~ as.factor(level) + as.numeric(Age) + as.factor(Sex) + as.factor(PlateId) + as.numeric(PC1) + as.numeric(PC2), data = FTD_expr)
  output <- as.data.frame(summary(model)[4])
  result <- cbind(as.character(analyte_ID), output[2,1], output[2,2], output[2,4])
  print(i)
  FTD_vs_CO_daa <- rbind(FTD_vs_CO_daa, result) 
}
colnames(FTD_vs_CO_daa) <- c("Analyte","Estimate", "Standard_error","Pvalue")
dim(FTD_vs_CO_daa) # 7008    4
any(is.na(FTD_vs_CO_daa)) # FALSE
head(FTD_vs_CO_daa, 2)
#            Analyte            Estimate    Standard_error             Pvalue
# 1 X10000.28_P43320  -0.668716651432551 0.312095302046033 0.0324930750241944
# 2  X10001.7_P04049 -0.0839726733793854  0.15683733118891  0.592529954737935
save(FTD_expr, FTD_vs_CO_daa, file="FTD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")

# multiple test correction
cols <- c(2,3,4)    
FTD_vs_CO_daa[,cols] <- apply(FTD_vs_CO_daa[,cols], 2, function(x) as.numeric(as.character(x)))
FTD_vs_CO_daa$FDR <- p.adjust(FTD_vs_CO_daa$Pvalue, method = "fdr", n = nrow(FTD_vs_CO_daa))
FTD_vs_CO_daa$Bonferroni <- p.adjust(FTD_vs_CO_daa$Pvalue, method = "bonferroni", n = nrow(FTD_vs_CO_daa))
nrow(FTD_vs_CO_daa[FTD_vs_CO_daa$Pvalue < 0.05,]) # 4194
nrow(FTD_vs_CO_daa[FTD_vs_CO_daa$FDR < 0.05,]) # 3873
nrow(FTD_vs_CO_daa[FTD_vs_CO_daa$Bonferroni < 0.05,]) # 1568

# Analyte annotation
load("CSF_Soma7K_LongsGarfield_Expr_Pheno_Annot.RData")
dim(longsGar_annot) # 7291   16
longsGar_annot <- longsGar_annot[,c("Analytes", "Symbol", "UniProt", "TargetFullName", "External_ID")]
longsGar_annot$Symbol <- sapply(strsplit(as.character(longsGar_annot$Symbol), "\\."), '[', 1)
table(FTD_vs_CO_daa$Analyte %in% longsGar_annot$External_ID)
# TRUE
# 7008
FTD_vs_CO_daa <- inner_join(FTD_vs_CO_daa, longsGar_annot, by=c("Analyte"="External_ID"))
dim(FTD_vs_CO_daa) # 7008   10
FTD_vs_CO_daa <- FTD_vs_CO_daa[order(FTD_vs_CO_daa$Pvalue),]
table(FTD_vs_CO_daa$Analyte %in% names(FTD_expr))
# TRUE
# 7008
FTD_vs_CO_daa <- FTD_vs_CO_daa[,c("Analyte", "UniProt", "Symbol", "TargetFullName", "Estimate", "Standard_error", "Pvalue", "FDR", "Bonferroni")]
colnames(FTD_vs_CO_daa) <- c("Analyte", "UniProt", "Symbol", "Target", "Estimate", "SE", "P", "FDR", "Bonferroni")
write.table(FTD_vs_CO_daa, file="CSF_LongsGar_PPMI_SADRC_Zscore_FTD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", row.names=F, quote=F)

# volcano plot
FTD_vs_CO_daa <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_FTD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, quote="")
FTD <- FTD_vs_CO_daa
FTD <- FTD[!duplicated(FTD$Symbol),]
top_10_labels <- unique(c(FTD[FTD$Estimate > 0,]$Symbol[1:7], FTD[FTD$Estimate < 0,]$Symbol[1:5]))
png("Volcano_CSF_LongsGar_PPMI_SADRC_Zscore_FTD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes_v2.png", units="mm", width=200, height=170, res=1000)
EnhancedVolcano(FTD,
    lab = FTD$Symbol,
    selectLab = top_10_labels,
    title= 'CSF FTD (n=46) vs. CO (n=726)', # NULL
    subtitle = 'Covariates: Age, Sex, Plate, PC1, PC2', # NULL
    x = 'Estimate',
    y = 'FDR',
    pCutoff = 0.05,
    FCcutoff = 0, # 0.05
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2,
    labSize = 4,
    colAlpha = 1/5,
    labFace = 'bold',
    colConnectors = 'black',
    caption = "", # to remove "total variables = xxxx" at the bottom of plot
    ylab = "FDR", xlab = "Estimate",
    legendLabels = c("NS", "NS", "FDR", "FDR"),
    ylim = c(0, max(-log10(FTD[,"FDR"]), na.rm = TRUE)),
    xlim = c(min(FTD$Estimate, na.rm = TRUE), max(FTD$Estimate, na.rm = TRUE)))
dev.off()


## PD vs. CO ##

PD_expr <- CSF_complete_expr_pcs[CSF_complete_expr_pcs$status_for_analysis %in% c("PD", "CO_ATneg"),]
dim(PD_expr) # 1465 7663
table(PD_expr$status_for_analysis)
# CO_ATneg       PD
#      726      739
PD_expr$level <- ifelse(PD_expr$status_for_analysis == "PD", 1, 0)
table(PD_expr$level)
#   0     1
# 726   739
round((nrow(PD_expr[PD_expr$Sex == "Female",])/nrow(PD_expr))*100, 1) # 47.2
round(mean(PD_expr$Age, na.rm=T), 1) # 65.1
round(sd(PD_expr$Age, na.rm=T), 1) # 9.6
# check discovery statistics
discovery_stats <- PD_expr[,c("UniquePhenoID", "Sex", "status_for_analysis", "Cohort")]
### get analytes
PD_expr <- PD_expr %>% select(UniquePhenoID, DrawDate, Cohort, Age, Sex, AT_class, status_for_analysis, level, PlateId, PC1, PC2, PC3, everything())
dim(PD_expr) # 1465 7664
PD_expr[1:4,c(1:13,7660:7664)]
# replace the NaN and Inf values in our data frame
PD_expr[is.na(PD_expr) | PD_expr == "Inf"] <- NA
# Remove features that do not have 2 values for "level" after NA removal
# Function to check if each feature column has at least two levels for each covariate
check_levels <- function(df, feature_col, covariates) {
  for (covariate in covariates) {
    non_na_data <- df[!is.na(df[[feature_col]]), c(covariate, feature_col)]
    if (length(unique(non_na_data[[covariate]])) < 2) {
      return(FALSE)
    }
  }
  return(TRUE)
}
covariates <- c("level", "Sex", "Cohort")
# Identify valid features that have at least two levels for each covariate
valid_features <- sapply(names(PD_expr)[13:ncol(PD_expr)], function(feature) {
  check_levels(PD_expr, feature, covariates)
})
table(valid_features)[2] # 7047
# Filter the dataframe to include only valid features and the covariate columns
PD_expr <- PD_expr[, c(names(PD_expr)[1:12], names(PD_expr)[13:ncol(PD_expr)][valid_features])]
dim(PD_expr) # 1465 7059
which(names(PD_expr) %in% "XIL23A.10365_P29460, Q9NPF7")
analyte <- as.data.frame(PD_expr [,13:ncol(PD_expr)])
dim(analyte) #  1465 7047
analyte_name <- as.data.frame(colnames(analyte))
head(analyte_name)

PD_vs_CO_daa <- data.frame()
# run model (for PD vs CO)
for (i in 1:ncol(analyte)) {
  analyte_ID <- analyte_name[i,1]
  protein <- analyte[,i]
  model <- lm(protein ~ as.factor(level) + as.numeric(Age) + as.factor(Sex) + as.factor(PlateId) + as.numeric(PC1) + as.numeric(PC2), data = PD_expr)
  output <- as.data.frame(summary(model)[4])
  result <- cbind(as.character(analyte_ID), output[2,1], output[2,2], output[2,4])
  print(i)
  PD_vs_CO_daa <- rbind(PD_vs_CO_daa, result) 
}
colnames(PD_vs_CO_daa) <- c("Analyte","Estimate", "Standard_error","Pvalue")
dim(PD_vs_CO_daa) # 7047    4
any(is.na(PD_vs_CO_daa)) # FALSE
head(PD_vs_CO_daa, 2)
#            Analyte             Estimate     Standard_error             Pvalue
# 1 X10000.28_P43320 -0.00948461956707186 0.0995835777370492  0.924136156053325
# 2  X10001.7_P04049   -0.128465417687495 0.0572635114757984 0.0250256961960041
save(PD_expr, PD_vs_CO_daa, file="PD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")

# multiple test correction
cols <- c(2,3,4)    
PD_vs_CO_daa[,cols] <- apply(PD_vs_CO_daa[,cols], 2, function(x) as.numeric(as.character(x)))
PD_vs_CO_daa$FDR <- p.adjust(PD_vs_CO_daa$Pvalue, method = "fdr", n = nrow(PD_vs_CO_daa))
PD_vs_CO_daa$Bonferroni <- p.adjust(PD_vs_CO_daa$Pvalue, method = "bonferroni", n = nrow(PD_vs_CO_daa))
nrow(PD_vs_CO_daa[PD_vs_CO_daa$Pvalue < 0.05,]) # 1291
nrow(PD_vs_CO_daa[PD_vs_CO_daa$FDR < 0.05,]) # 370
nrow(PD_vs_CO_daa[PD_vs_CO_daa$Bonferroni < 0.05,]) # 70

# Analyte annotation
load("CSF_Soma7K_LongsGarfield_Expr_Pheno_Annot.RData")
dim(longsGar_annot) # 7291   16
longsGar_annot <- longsGar_annot[,c("Analytes", "Symbol", "UniProt", "TargetFullName", "External_ID")]
longsGar_annot$Symbol <- sapply(strsplit(as.character(longsGar_annot$Symbol), "\\."), '[', 1)
table(PD_vs_CO_daa$Analyte %in% longsGar_annot$External_ID)
# FALSE  TRUE
#     1  7046
PD_vs_CO_MAP <- inner_join(PD_vs_CO_daa, longsGar_annot, by=c("Analyte"="External_ID"))
dim(PD_vs_CO_MAP) # 7046   10
PD_vs_CO_daa <- anti_join(PD_vs_CO_daa, longsGar_annot, by=c("Analyte"="External_ID"))
dim(PD_vs_CO_daa) # 1 6
load("CSF_Soma5K_StanfordADRC_Expr_Pheno_Annot.RData")
dim(sadrc_annot) # 4735   11
sadrc_annot <- sadrc_annot[,c("Analytes", "Symbol", "UniProt", "TargetFullName", "External_ID")]
sadrc_annot$Symbol <- sapply(strsplit(as.character(sadrc_annot$Symbol), "\\."), '[', 1)
table(PD_vs_CO_daa$Analyte %in% sadrc_annot$External_ID)
# TRUE
#    1
PD_vs_CO_SADRC <- inner_join(PD_vs_CO_daa, sadrc_annot, by=c("Analyte"="External_ID"))
dim(PD_vs_CO_SADRC) # 1 10
PD_vs_CO_daa <- rbind(PD_vs_CO_MAP, PD_vs_CO_SADRC)
PD_vs_CO_daa <- PD_vs_CO_daa[order(PD_vs_CO_daa$Pvalue),]
table(PD_vs_CO_daa$Analyte %in% names(PD_expr))
# TRUE
# 7047
PD_vs_CO_daa <- PD_vs_CO_daa[,c("Analyte", "UniProt", "Symbol", "TargetFullName", "Estimate", "Standard_error", "Pvalue", "FDR", "Bonferroni")]
colnames(PD_vs_CO_daa) <- c("Analyte", "UniProt", "Symbol", "Target", "Estimate", "SE", "P", "FDR", "Bonferroni")
write.table(PD_vs_CO_daa, file="CSF_LongsGar_PPMI_SADRC_Zscore_PD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", row.names=F, quote=F)

# volcano plot
PD_vs_CO_daa <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_PD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, quote="")
PD <- PD_vs_CO_daa
PD <- PD[!duplicated(PD$Symbol),]
top_10_labels <- unique(c(PD[PD$Estimate > 0,]$Symbol[1:7], PD[PD$Estimate < 0,]$Symbol[1:5]))
png("Volcano_CSF_LongsGar_PPMI_SADRC_Zscore_PD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes_v2.png", units="mm", width=200, height=170, res=1000)
EnhancedVolcano(PD,
    lab = PD$Symbol,
    selectLab = top_10_labels,
    title= 'CSF PD (n=739) vs. CO (n=726)', # NULL
    subtitle = 'Covariates: Age, Sex, Plate, PC1, PC2', # NULL
    x = 'Estimate',
    y = 'FDR',
    pCutoff = 0.05,
    FCcutoff = 0, # 0.05
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    pointSize = 2,
    labSize = 4,
    colAlpha = 1/5,
    labFace = 'bold',
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = "", # to remove "total variables = xxxx" at the bottom of plot
    ylab = "FDR", xlab = "Estimate",
    legendLabels = c("NS", "NS", "FDR", "FDR"),
    ylim = c(0, max(-log10(PD[,"FDR"]), na.rm = TRUE)),
    xlim = c(min(PD$Estimate, na.rm = TRUE), max(PD$Estimate, na.rm = TRUE)))
dev.off()


## Combined CSF Differential Abundance Analyses (DAA) Results

rm(list = ls())
options(stringsAsFactors = FALSE)
library(dplyr)
options(width = 160)
set.seed(1)

AD <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_AD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, quote="")
dim(AD) # 7099   10

DLB <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_DLB_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, quote="")
dim(DLB) # 7099   10

FTD <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_FTD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, quote="")
dim(FTD) # 7008   10

PD <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_PD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, quote="")
dim(PD) # 7047    9

save(AD, DLB, FTD, PD, file="CSF_AD_DLB_FTD_PD_DAA_Sumstats_Tables.RData")

DLB <- DLB[,c("Analyte", "Estimate", "SE", "P", "FDR", "Bonferroni")]
FTD <- FTD[,c("Analyte", "Estimate", "SE", "P", "FDR", "Bonferroni")]
PD <- PD[,c("Analyte", "Estimate", "SE", "P", "FDR", "Bonferroni")]
AD_DLB <- full_join(AD, DLB, by="Analyte")
dim(AD_DLB) # 7099   14
AD_DLB_FTD <- full_join(AD_DLB, FTD, by="Analyte")
dim(AD_DLB_FTD) # 7099   19
AD_DLB_FTD_PD <- full_join(AD_DLB_FTD, PD, by="Analyte")
dim(AD_DLB_FTD_PD) # 7099   24
AD_DLB_FTD_PD[is.na(AD_DLB_FTD_PD)] <- ""
write.table(AD_DLB_FTD_PD, file="CSF_AD_DLB_FTD_PD_DAA_Sumstats_SupplTable.txt", sep="\t", row.names=F, quote=F)



## Plasma DAA ##

rm(list = ls())
# --- Libraries ---
library(readr)
library(dplyr)
library(broom)
library(ggplot2)
library(openxlsx)

# --- Parameters ---
target_disease <- "AD"  # Change this to "AD", "PD", "FTD", or "DLB"
input_file <- paste0("processed_data_with_zscores_CO_", target_disease, ".csv")
output_file <- paste0("significance_results_CO_", target_disease, ".xlsx")

# --- Load data ---
data <- read_csv(input_file)

# --- Standardize Sex variable ---
data$Sex <- recode(data$Sex, "F" = "Female", "M" = "Male")

# --- Create binary outcome: 0 = Control (CO), 1 = Disease ---
data$Status <- ifelse(data$status_at_draw_mapping == "CO", 0, 1)

# --- Extract protein names (assuming starts from column 4 to 6597) ---
protein_cols <- colnames(data)[4:6597]

# --- Function to run linear model on one protein ---
run_lm_on_protein <- function(protein_name, data) {
  formula <- as.formula(paste0(protein_name, " ~ Status + Age_at_draw + Sex + Project_x + PC1 + PC2"))
  model <- lm(formula, data = data)
  tidy_model <- broom::tidy(model)
  status_row <- tidy_model %>% filter(grepl("Status", term))
  if (nrow(status_row) > 0) {
    return(data.frame(
      Protein = protein_name,
      EffectSize = status_row$estimate,
      SE = status_row$std.error,
      PValue = status_row$p.value
    ))
  } else {
    return(NULL)
  }
}

# --- Run regression on all proteins ---
results_list <- lapply(protein_cols, run_lm_on_protein, data = data)
significance_results <- bind_rows(results_list)

# --- Multiple testing correction ---
significance_results$FDR_BH_PValue <- p.adjust(significance_results$PValue, method = "BH")
significance_results$Bonferroni_PValue <- p.adjust(significance_results$PValue, method = "bonferroni")

# --- Output ---
write.xlsx(significance_results, output_file, rowNames = FALSE)
print(paste0("Results saved to: ", output_file))


