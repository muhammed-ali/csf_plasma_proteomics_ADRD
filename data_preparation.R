rm(list = ls())
library(data.table)
library(readxl)
library(dplyr)
library(factoextra)
library(ggfortify)
options(stringsAsFactors = FALSE)
options(width = 160)
set.seed(1)

## Knight-ADRC Cohort ##

# Phenotype data
longsGarPheno <- data.table(read_excel("Phenotype_CSF_Longs_Garfield_03272024.xlsx", sheet = "Sheet1"))
dim(longsGarPheno) # 3446   24
length(unique(longsGarPheno$UniquePhenoID)) # 3346
# Remove 100 longitudinal samples
longsGarPheno <- longsGarPheno[order(longsGarPheno$UniquePhenoID, longsGarPheno$DrawDate),]
head(longsGarPheno[duplicated(longsGarPheno$UniquePhenoID),], 6)
longsGarPheno <- longsGarPheno[!duplicated(longsGarPheno$UniquePhenoID),]
dim(longsGarPheno) # 3346   24
length(unique(longsGarPheno$UniquePhenoID)) # 3346

# Filter for individuals to be kept
as.data.frame(table(longsGarPheno$Status_at_draw))
longsGarPheno$status_for_analysis <- ifelse(longsGarPheno$Status_at_draw == "AD", "AD",  
  ifelse(longsGarPheno$Status_at_draw== "CO", "CO_All",
  ifelse(longsGarPheno$Status_at_draw== "DLB", "DLB",
  ifelse(longsGarPheno$Status_at_draw== "FTD", "FTD",
  ifelse(longsGarPheno$Status_at_draw== "C9ORF72_Affected_carrier", "FTD",
  ifelse(longsGarPheno$Status_at_draw== "PGRN_Unaffected_carrier_Unknown", "FTD",
  ifelse(longsGarPheno$Status_at_draw== "MAPT_Affected_carrier_Unknown", "FTD",
  ifelse(longsGarPheno$Status_at_draw== "MAPT_Affected_carrier", "FTD",
  ifelse(longsGarPheno$Status_at_draw== "MAPT_Presymptomatic_carrier", "FTD",
  ifelse(longsGarPheno$Status_at_draw== "C9ORF72_Presymptomatic_carrier", "FTD",
  ifelse(longsGarPheno$Status_at_draw== "PD", "PD", "Exclude")))))))))))
table(longsGarPheno$status_for_analysis)
#     AD  CO_All     DLB Exclude     FTD      PD
#   1138    1175      35     799      46     147
# Filter for CO that are not A-T-
table(longsGarPheno[longsGarPheno$status_for_analysis == "CO_All",]$AT_class)
#   A- A-T- A-T+   A+ A+T- A+T+
#   12  620   59    1  219   97

longsGarPheno <- longsGarPheno %>% mutate(status_for_analysis = case_when(
  status_for_analysis == "CO_All" & AT_class == "A-T-" ~ "CO_ATneg", TRUE~status_for_analysis))
table(longsGarPheno$status_for_analysis)
#      AD   CO_All CO_ATneg      DLB  Exclude      FTD       PD
#    1138      555      620       35      799       46      147
dim(longsGarPheno) # 3346   25
longsGarPheno_subset <- longsGarPheno[,c("UniquePhenoID", "DrawDate", "Project", "Age_at_draw", "Sex", "AT_class", "status_for_analysis")]
colnames(longsGarPheno_subset) <- c("UniquePhenoID", "DrawDate", "Cohort", "Age", "Sex", "AT_class", "status_for_analysis")
dim(longsGarPheno_subset) # 3346    7

# Expression data
longsGar_expr <- read.csv("Zscore_Somascan7k_CSF_Longs_Garfield_protein_matrix.csv")
dim(longsGar_expr) # 3446 7082
table(longsGar_expr$Project)
# Garfield    Longs
#      465     2981

# Analyte annotation data 
longsGar_annot <- read.table("CSF_SOMAscan7k_analyte_information.tsv", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(longsGar_annot) # 7291   14
length(unique(longsGar_annot$Analytes)) # 7291
length(unique(longsGar_annot$UniProt)) # 6404
length(unique(longsGar_annot$EntrezGeneSymbol)) # 6388
any(is.na(longsGar_annot$UniProt)) # FALSE
longsGar_annot <- longsGar_annot[order(longsGar_annot$Analytes),]
longsGar_annot$External_ID <- paste0(longsGar_annot$Analytes, "_", longsGar_annot$UniProt, sep="")
length(unique(longsGar_annot$External_ID)) # 7291
longsGar_annot$Symbol <- longsGar_annot$EntrezGeneSymbol
longsGar_annot$Symbol <- make.unique(as.character(longsGar_annot$EntrezGeneSymbol), sep = ".")
length(unique(longsGar_annot$Symbol)) # 7291
# Rename colnames to match AnlyteID_Uniprot format for mapping across Soma 7K and 5K
table(names(longsGar_expr)[8:ncol(longsGar_expr)] %in% longsGar_annot$Analytes)
# TRUE
# 7075
colnames(longsGar_expr)[8:ncol(longsGar_expr)] <- longsGar_annot$External_ID[match(colnames(longsGar_expr)[8:ncol(longsGar_expr)], longsGar_annot$Analytes)]
longsGar_expr[ ,c("Project", "Tissue", "ExtIdentifier", "Somalogic_Barcode")] <- list(NULL)

# Merge pheno with expression data
dim(longsGar_expr) # 3446 7078
dim(longsGarPheno_subset) # 3346    7
table(longsGar_expr$UniquePhenoID %in% longsGarPheno_subset$UniquePhenoID)
# TRUE
# 3446
longsGarPheno_subset$DrawDate <- as.character(longsGarPheno_subset$DrawDate)
longsGar_expr_pheno <- inner_join(longsGarPheno_subset, longsGar_expr, by=c("UniquePhenoID", "DrawDate"))
dim(longsGar_expr_pheno) # 3346 7083
table(longsGar_expr_pheno$status_for_analysis)
#      AD   CO_All CO_ATneg      DLB  Exclude      FTD       PD
#    1138      555      620       35      799       46      147
table(names(longsGar_expr_pheno)[9:ncol(longsGar_expr_pheno)] %in% longsGar_annot$External_ID)
# TRUE
# 7075
table(longsGar_expr_pheno$status_for_analysis, longsGar_expr_pheno$AT_class)
#            A- A-T- A-T+  A+ A+T- A+T+
#  AD         7  118   37   6  206  678
#  CO_All    12    0   59   1  219   97
#  CO_ATneg   0  620    0   0    0    0
#  DLB        0    1    0   0    1    2
#  Exclude    1  382   35   1  114  125
#  FTD        1    7    0   0    1    0
#  PD         0    0    0   0    0    0

save(longsGar_expr_pheno, longsGar_annot, longsGarPheno, file="CSF_Soma7K_LongsGarfield_Expr_Pheno_Annot.RData")


## PPMI ##

ppmiPheno <- read.csv("PPMI_clinical_information_202403.csv")
dim(ppmiPheno) # 1158    8
length(unique(ppmiPheno$PATNO)) # 1158
table(ppmiPheno$Diagnosis)
#    Healthy Control Parkinson's Disease           Prodromal
#                186                 617                 355
table(ppmiPheno$Diagnosis, ppmiPheno$AT_class)
#                      A-T- A-T+ A+T- A+T+
#  Healthy Control       53   82   16   11
#  Parkinson's Disease  181  158   82   33
#  Prodromal             26   33   12    2
ppmiPheno$Diagnosis <- ifelse(ppmiPheno$Diagnosis == "Healthy Control", "CO_All", ppmiPheno$Diagnosis)
ppmiPheno$Diagnosis <- ifelse(ppmiPheno$Diagnosis == "Parkinson's Disease", "PD", ppmiPheno$Diagnosis)
ppmiPheno <- ppmiPheno %>% mutate(status_for_analysis = case_when(
  Diagnosis == "CO_All" & AT_class == "A-T-" ~ "CO_ATneg", TRUE~Diagnosis))
table(ppmiPheno$status_for_analysis)
#   CO_All  CO_ATneg        PD Prodromal
#      133        53       617       355
ppmiPheno$Diagnosis <- NULL
colnames(ppmiPheno)[c(1,4,5)] <- c("UniquePhenoID", "Age", "Sex")
ppmiPheno$Age <- round(ppmiPheno$Age, 0)
# A subset of pheno for merging with expression
ppmiPheno_subset <- ppmiPheno[,c("UniquePhenoID", "DrawDate", "Age", "Sex", "AT_class", "status_for_analysis")]
ppmiPheno_subset$Cohort <- "PPMI"
ppmiPheno_subset$PlateId <- "None"
dim(ppmiPheno_subset) # 1158    8


# Expression data
ppmi_expr <- read.csv("Zscores_PPMI_CSF_Soma5k_human_only_protein_matrix.csv")
dim(ppmi_expr) # 1075 4778
length(unique(names(ppmi_expr))) # 4778
ppmi_expr[1:4,1:10]


# Analyte Annotation
ppmi_annot <- read.csv("PPMI_soma_annotation_by_Chengran.csv")
dim(ppmi_annot) # 4942    7
length(unique(ppmi_annot$SeqId)) # 4782
length(unique(ppmi_annot$UniProt)) # 4130
ppmi_annot$Analytes <- gsub('-', '.', ppmi_annot$SeqId, fixed = TRUE)
ppmi_annot$Analytes <- paste0("X", ppmi_annot$Analytes, sep="")
table(names(ppmi_expr) %in% ppmi_annot$Analytes)
# FALSE  TRUE
#     1  4777
# Fix the duplicated Analytes (duplication is due to one analyte being associated with different TARGET_GENE_ID)
ppmi_annot <- ppmi_annot[order(ppmi_annot$Analytes),]
head(ppmi_annot[duplicated(ppmi_annot$Analytes),], 6)
ppmi_annot[ppmi_annot$Analytes %in% c("X10365.132", "X10367.62", "X10387.1"),]
ppmi_annot <- ppmi_annot[!duplicated(ppmi_annot$Analytes),]
ppmi_annot[ppmi_annot$Analytes %in% c("X10365.132", "X10367.62", "X10387.1"),]
length(unique(ppmi_annot$Analytes)) # 4782
length(unique(ppmi_annot$UniProt)) # 4125
ppmi_annot$External_ID <- paste0(ppmi_annot$Analytes, "_", ppmi_annot$UniProt, sep="")
length(unique(ppmi_annot$External_ID)) # 4782
length(unique(ppmi_annot$TARGET_GENE_SYMBOL)) # 4105
ppmi_annot$Symbol <- make.unique(as.character(ppmi_annot$TARGET_GENE_SYMBOL), sep = ".")
length(unique(ppmi_annot$Symbol)) # 4782

# Rename colnames to match AnlyteID_Uniprot format for mapping across Soma 7K and 5K
table(names(ppmi_expr)[2:ncol(ppmi_expr)] %in% ppmi_annot$Analytes)
# TRUE
# 4777
colnames(ppmi_expr)[2:ncol(ppmi_expr)] <- ppmi_annot$External_ID[match(colnames(ppmi_expr)[2:ncol(ppmi_expr)], ppmi_annot$Analytes)]
names(ppmi_expr)[1] <- "UniquePhenoID"

# Merge pheno with expression data
dim(ppmi_expr) # 1075 4778
dim(ppmiPheno_subset) # 1158    6
table(ppmi_expr$UniquePhenoID %in% ppmiPheno_subset$UniquePhenoID)
# TRUE
# 1075
ppmi_expr_pheno <- inner_join(ppmiPheno_subset, ppmi_expr, by="UniquePhenoID")
dim(ppmi_expr_pheno) # 1075 4785
table(ppmi_expr_pheno$status_for_analysis)
#   CO_All  CO_ATneg        PD Prodromal
#      116        52       575       332
table(names(ppmi_expr_pheno)[9:ncol(ppmi_expr_pheno)] %in% ppmi_annot$External_ID)
# TRUE
# 4777
table(ppmi_expr_pheno$status_for_analysis, ppmi_expr_pheno$AT_class)
#            A-T- A-T+ A+T- A+T+
#  CO_All       0   72   15    9
#  CO_ATneg    52    0    0    0
#  PD         174  141   78   33
#  Prodromal   24   31    9    2

save(ppmi_expr_pheno, ppmi_annot, ppmiPheno, file="CSF_Soma5K_PPMI_Expr_Pheno_Annot.RData")


## Stanford-ADRC cohort ##

# Expression data
sadrc_expr <- read.csv("Zscore_Stanford_CSF_SOMAscan5k_matrix_cleaned.csv")
dim(sadrc_expr) # 274 4742
length(unique(sadrc_expr$Barcode)) # 274
length(unique(sadrc_expr$PIDN)) # 264
table(sadrc_expr$Visit)
#   1   2
# 251  23
sadrc_expr <- sadrc_expr[order(sadrc_expr$PIDN, sadrc_expr$Visit),]
sadrc_expr[duplicated(sadrc_expr$PIDN),][1:6:1:4]
sadrc_expr[sadrc_expr$PIDN %in% c("611", "613"),][1:4,1:10]
sadrc_expr <- sadrc_expr[!duplicated(sadrc_expr$PIDN),]
dim(sadrc_expr) # 264 4742
sadrc_expr[sadrc_expr$PIDN %in% c("611", "613"),][1:4,1:10]
length(unique(sadrc_expr$PIDN)) # 264
sadrc_expr$Barcode <- NULL
sadrc_expr$Visit <- NULL
names(sadrc_expr)[1:4] <- c("UniquePhenoID", "Age", "Sex", "status_for_analysis")
sadrc_expr$Sex <- ifelse(sadrc_expr$Sex == "F", "Female", "Male")
table(sadrc_expr$status_for_analysis)
#     AD      HC     LBD     MCI  MCI-PD      PD Unknown
#     19     201       2      21       3      17       1

# Analyte annotation data
sadrc_annot <- read.csv("Stanford_CSF_SOMAscan5k_annotation_by_Muhammad.csv")
dim(sadrc_annot) # 4735   10
head(sadrc_annot, 3)
colnames(sadrc_annot)[9:10] <- c("Analytes","External_ID")
length(unique(sadrc_annot$Analytes)) # 4735
length(unique(sadrc_annot$External_ID)) # 4735
length(unique(sadrc_annot$EntrezGeneSymbol)) # 4488
sadrc_annot <- sadrc_annot[order(sadrc_annot$Analytes),]
sadrc_annot$Symbol <- make.unique(as.character(sadrc_annot$EntrezGeneSymbol), sep = ".")
length(unique(sadrc_annot$Symbol)) # 4735
table(names(sadrc_expr)[8:ncol(sadrc_expr)] %in% sadrc_annot$Analytes)
# FALSE  TRUE
#     2  4733
# Check the samples not present in the expression-annotation data
`%!in%` <- Negate(`%in%`)
which(names(sadrc_expr)[6:ncol(sadrc_expr)] %!in% sadrc_annot$Analytes) # 2485 3078
# While preparing the "Stanford_CSF_SOMAscan5k_annotation_by_Muhammad.csv" file, I changed the follwoing two analytes due to inconsistent names 
# analyte_info$AnalyteID <- ifelse(analyte_info$SOMAseqID == "FGA.FGB.FGG.4907.56.1", "X4907.56", analyte_info$AnalyteID)
# analyte_info$AnalyteID <- ifelse(analyte_info$SOMAseqID == "FGA.FGB.FGG.2796.62.2", "X2796.62", analyte_info$AnalyteID)
# So, I need to change back the analyteID in the expressin matrix as well
sadrc_expr[1:4,c(1:6,2490, 3083)]
names(sadrc_expr)[c(2490, 3083)] <- c("X4907.56","X2796.62")
table(names(sadrc_expr)[6:ncol(sadrc_expr)] %in% sadrc_annot$Analytes)
# TRUE
# 4735
# Also change the Symbol and External_ID in the annotation file
sadrc_annot[sadrc_annot$Analytes %in% c("X4907.56","X2796.62"),]
sadrc_annot[sadrc_annot$Analytes == "X2796.62",]$External_ID <- "X2796.62_P02671|P02675|P02679"
sadrc_annot[sadrc_annot$Analytes == "X4907.56",]$External_ID <- "X4907.56_P02671|P02675|P02679"
sadrc_annot[sadrc_annot$Analytes == "X2796.62",]$Symbol <- "FGA"
sadrc_annot[sadrc_annot$Analytes == "X4907.56",]$Symbol <- "FGA.1"
sadrc_annot[sadrc_annot$Analytes %in% c("X4907.56","X2796.62"),]


# AT_Status
sadrc_AT <- read.csv("Stanford_ADRC_csf_biomarkers.csv")
dim(sadrc_AT) # 145  10
length(unique(sadrc_AT$PIDN)) # 145
sadrc_AT$AT <- paste0(sadrc_AT$'AB.', sadrc_AT$'Tau.')
dim(sadrc_AT) # 145  11
table(sadrc_AT$AT)
# 00 01 10 11
# 86  6 25 28
sadrc_AT$AT_class <- ifelse(sadrc_AT$AT == "00", "A-T-",  
  ifelse(sadrc_AT$AT == "01", "A-T+",
  ifelse(sadrc_AT$AT == "10", "A+T-",
  ifelse(sadrc_AT$AT == "11", "A+T+", "Exclude"))))
table(sadrc_AT$AT_class)
# A-T- A-T+ A+T- A+T+
#   86    6   25   28
sadrc_AT <- sadrc_AT[,c("PIDN", "AT_class")]
names(sadrc_AT) <- c("UniquePhenoID", "AT_class")

# Merge AT status with expression and simplify status_for_analysis
dim(sadrc_AT) # 145   2
dim(sadrc_expr) # 264 4740
sadrc_expr_pheno <- left_join(sadrc_expr, sadrc_AT, by="UniquePhenoID")
dim(sadrc_expr_pheno) # 264 4741
sadrc_expr_pheno$DrawDate <- "None"
sadrc_expr_pheno$Cohort <- "SADRC"
sadrc_expr_pheno <- sadrc_expr_pheno %>% select(UniquePhenoID, Age, Sex, status_for_analysis, AT_class, PlateId, DrawDate, Cohort, everything())

sadrc_expr_pheno$status_for_analysis <- ifelse(sadrc_expr_pheno$status_for_analysis == "HC", "CO_All", sadrc_expr_pheno$status_for_analysis)
sadrc_expr_pheno <- sadrc_expr_pheno %>% mutate(Final_Status = case_when(
  status_for_analysis == "CO_All" & AT_class == "A-T-" ~ "CO_ATneg",
  status_for_analysis == "LBD" ~ "DLB", TRUE~status_for_analysis))
table(sadrc_expr_pheno$Final_Status)
#      AD   CO_All CO_ATneg      DLB      MCI   MCI-PD       PD  Unknown
#      19      147       54        2       21        3       17        1
sadrc_expr_pheno$status_for_analysis <- sadrc_expr_pheno$Final_Status
sadrc_expr_pheno$Final_Status <- NULL
table(sadrc_expr_pheno$status_for_analysis, sadrc_expr_pheno$AT_class)
#           A-T- A-T+ A+T- A+T+
#  AD          2    2    0   14
#  CO_All      0    4   17    8
#  CO_ATneg   54    0    0    0
#  DLB         0    0    2    0
#  MCI         6    0    4    5
#  MCI-PD      3    0    0    0
#  PD         15    0    2    0
#  Unknown     0    0    0    0
table(names(sadrc_expr_pheno)[9:ncol(sadrc_expr_pheno)] %in% sadrc_annot$Analytes)
# TRUE
# 4735
colnames(sadrc_expr_pheno)[9:ncol(sadrc_expr_pheno)] <- sadrc_annot$External_ID[match(colnames(sadrc_expr_pheno)[9:ncol(sadrc_expr_pheno)], sadrc_annot$Analytes)]
save(sadrc_expr_pheno, sadrc_annot, sadrc_AT, file="CSF_Soma5K_StanfordADRC_Expr_Pheno_Annot.RData")


## Principal Compnent Analysis - Create Overlapping and Combined Matrix ##

overlapping_prot <- Reduce(intersect, list(names(longsGar_expr_pheno), names(ppmi_expr_pheno), names(sadrc_expr_pheno)))
length(overlapping_prot) # 3639
longsGar_expr_pheno <- as.data.frame(longsGar_expr_pheno)
longsGar_overlap <- longsGar_expr_pheno[, overlapping_prot]
ppmi_overlap <- ppmi_expr_pheno[, overlapping_prot]
sadrc_overlap <- sadrc_expr_pheno[, overlapping_prot]
PCA_expr_matrix <- rbind(longsGar_overlap, ppmi_overlap, sadrc_overlap)
dim(PCA_expr_matrix) # 4685 3639
save(PCA_expr_matrix, file="PCA_expr_matrix.RData")

# Check Analyte with NA values in the overlapping analyte expression matrix
abc <- PCA_expr_matrix[9:ncol(PCA_expr_matrix)]
na_count <-sapply(abc, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count$External_ID <- row.names(na_count)
na_count <- na_count[order(na_count$na_count, decreasing=T),]
nrow(na_count[na_count$na_count > 2500,]) # 31 analytes with more than 2500 NA in the expression values out of total 3639 analytes
most_na_columns <- c(names(PCA_expr_matrix)[1:8], na_count[na_count$na_count > 2000,]$External_ID)
length(most_na_columns) # 39
most_na_matrix <- PCA_expr_matrix[,most_na_columns]
dim(most_na_matrix) # 4685   39
write.csv(most_na_matrix, file="most_na_matrix.csv")
# Upon checking the excel file, I noted that these analytes are only present in some of the MAP samples, PD_MARS, PPMI, and SADRC. Not present in ADNI, PAU, or RUIZ.
png("NA_Distribution_MAP_PPMI_SADRC_OverlappingProteins3639.png", units="mm", width=200, height=170, res=1000)
hist(na_count$na_count)
dev.off()

# Data imputation (Random replacemenet with bootstraping before PCA) for the overallping analyte matrix
row.names(PCA_expr_matrix) <- PCA_expr_matrix$UniquePhenoID
expr_matrix <- PCA_expr_matrix[,9:ncol(PCA_expr_matrix)]
dim(expr_matrix) # 4685 3631
expr_matrix[1:4,1:10]
any(is.na(expr_matrix)) # TRUE
sum(is.na(expr_matrix)) # 415113
sum(!is.na(expr_matrix)) # 16596122
(sum(is.na(expr_matrix))/sum(!is.na(expr_matrix)))*100 # 2.5% missing data
expr_matrix_nonNA <- as.data.frame(as.matrix(sapply(expr_matrix, as.numeric)))
count.na = apply(is.na(expr_matrix_nonNA), 2, sum)
set.seed(1)
for (i in which(count.na!=0)) { # bootstrapping with replacement
  index = is.na(expr_matrix_nonNA[,i])
  expr_matrix_nonNA[index, i] = sample(expr_matrix_nonNA[!index,i], sum(index), replace = TRUE)
}
data_imputed <- as.data.frame(cbind(PCA_expr_matrix[,1:8], expr_matrix_nonNA))
dim(data_imputed) # 4685 3639
any(is.na(data_imputed[,9:ncol(data_imputed)])) # FALSE
pca_input <- as.data.frame(as.matrix(sapply(data_imputed[,9:ncol(data_imputed)], as.numeric)))

# PCA for the overallping analyte matrix
data_imputed_pca <- prcomp(pca_input, center=TRUE, scale=TRUE)
png("PCA_plot_Random_Replacemnt_MAP_PPMI_SADRC_OverlappingProteins3639.png", units="mm", width=200, height=170, res=1000)
fviz_eig(data_imputed_pca)
dev.off()

# how many PCs explain 90% of the variation?
summ <- summary(data_imputed_pca)
head(summ$importance[2,])
sum(head(summ$importance[2,], 1040)) # 0.90106 (1040 PCs explain 90% of variation)

PCs <- as.data.frame(data_imputed_pca$x)
PCs$UniquePhenoID <- PCA_expr_matrix$UniquePhenoID
PCs <- PCs %>% select(UniquePhenoID, everything())
dim(PCs) # 4685 3632
write.table(PCs, file="CSF_Soma7K_MAP_5K_PPMI_SADRC_PCA.txt", sep="\t", row.names=F)
save(PCs, data_imputed_pca, file="CSF_Soma7K_MAP_5K_PPMI_SADRC_PCA.RData")
complete_expr_df_pc <- inner_join(PCA_expr_matrix, PCs[,1:4], by="UniquePhenoID")
complete_expr_df_pc <- complete_expr_df_pc %>% select(UniquePhenoID, DrawDate, Cohort, Age, Sex, AT_class, status_for_analysis, PlateId, PC1, PC2, PC3, everything())
dim(complete_expr_df_pc) # 4685 3642
complete_expr_df_pc_subset <- complete_expr_df_pc[complete_expr_df_pc$status_for_analysis %in% c("AD", "CO_ATneg", "DLB", "FTD", "PD"),]
dim(complete_expr_df_pc_subset) # 2705 3642
# https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
png("PCA_plot_Random_Replacemnt_PCA_MAP_PPMI_SADRC_OverlappingProteins3639.png", units="mm", width=200, height=170, res=1000)
ggplot(complete_expr_df_pc_subset, aes(x=PC1, y=PC2, col=status_for_analysis)) + geom_point(alpha = 0.3)
# autoplot(data_imputed_pca, data = complete_expr_df_pc, colour = 'Final_Status', alpha = 0.3)
dev.off()
# if you want to do autoplot of PC1 and PC3
png("PCA_plot_Random_Replacemnt_PCA13_MAP_PPMI_SADRC_OverlappingProteins3639.png", units="mm", width=200, height=170, res=1000)
ggplot(complete_expr_df_pc_subset, aes(x=PC1, y=PC3, col=status_for_analysis)) + geom_point(alpha = 0.3)
# autoplot(data_imputed_pca, x=1, y=3, data = complete_expr_df_pc, colour = 'Final_Status', alpha = 0.3)
dev.off()

save(PCs, complete_expr_df_pc, file="PCA_plot_Random_Replacemnt_MAP_PPMI_SADRC_OverlappingProteins3639.RData")


# Combine all cohorts: Create a complete matrix of 3 studies by filling in non-overlapping columns with NAs
longsGar_expr_pheno[setdiff(names(sadrc_expr_pheno), names(longsGar_expr_pheno))] <- NA
longsGar_expr_pheno[setdiff(names(ppmi_expr_pheno), names(longsGar_expr_pheno))] <- NA
sadrc_expr_pheno[setdiff(names(longsGar_expr_pheno), names(sadrc_expr_pheno))] <- NA
sadrc_expr_pheno[setdiff(names(ppmi_expr_pheno), names(sadrc_expr_pheno))] <- NA
ppmi_expr_pheno[setdiff(names(longsGar_expr_pheno), names(ppmi_expr_pheno))] <- NA
ppmi_expr_pheno[setdiff(names(sadrc_expr_pheno), names(ppmi_expr_pheno))] <- NA
dim(longsGar_expr_pheno) # 3346 7660
dim(sadrc_expr_pheno) # 264 7660
dim(ppmi_expr_pheno) # 1075 7660
CSF_complete_expr <- rbind(longsGar_expr_pheno, sadrc_expr_pheno, ppmi_expr_pheno)
dim(CSF_complete_expr) # 4685 7660
CSF_complete_expr_pcs <- inner_join(CSF_complete_expr, PCs[,1:11], by="UniquePhenoID")
dim(CSF_complete_expr_pcs) # 4685 7670
CSF_complete_expr_pcs <- CSF_complete_expr_pcs %>% select(UniquePhenoID, DrawDate, Cohort, Age, Sex, AT_class, status_for_analysis, PlateId, PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10, everything())
CSF_complete_expr_pcs[1:4,c(1:25, 7668:7670)]
7670-18 # 7652
save(CSF_complete_expr_pcs, file="CSF_Complete_Expr_MAP_PPMI_SADRC_N4685_PC10.RData")
write.table(CSF_complete_expr_pcs, file="CSF_Complete_Expr_MAP_PPMI_SADRC_N4685_PC10.txt", sep="\t", row.names=F, quote=F)

# Check Analyte with NA values in the complete analyte expression matrix
abc <- CSF_complete_expr[9:ncol(CSF_complete_expr)]
na_count <-sapply(abc, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
na_count$External_ID <- row.names(na_count)
na_count <- na_count[order(na_count$na_count, decreasing=T),]
nrow(na_count[na_count$na_count > 4000,]) # 83 analytes with more than 4000 NA in the expression values out of total 7K analytes
nrow(na_count[na_count$na_count > 2500,]) # 31 analytes with more than 2500 NA in the expression values out of total 7K analytes
most_na_columns <- c(names(CSF_complete_expr)[1:8], na_count[na_count$na_count > 2000,]$External_ID)
length(most_na_columns) # 652
most_na_matrix <- CSF_complete_expr[,most_na_columns]
dim(most_na_matrix) # 4685  652
write.csv(most_na_matrix, file="most_na_matrix_7K_analytes.csv")
# Upon checking the excel file, I noted that these analytes are only present in some of the MAP samples, PD_MARS, PPMI, and SADRC. Not present in ADNI, PAU, or RUIZ.
png("NA_Distribution_MAP_PPMI_SADRC_7k_Analytes.png", units="mm", width=200, height=170, res=1000)
hist(na_count$na_count)
dev.off()