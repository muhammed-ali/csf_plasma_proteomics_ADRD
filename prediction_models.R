rm(list = ls())
library(dplyr) # dplyr_1.1.4
library(pROC) # pROC_1.18.5
library(ROCR) # ROCR_1.0-11
library(caret) # caret_6.0-94
library(glmnet) # glmnet_4.1-8
library(ggplot2) # ggplot2_3.5.1
options(stringsAsFactors = FALSE)
options(width = 160)
set.seed(1)


## AD ##

# Read in the expression matrices for all disease
load("AD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
dim(AD_expr) # 1883 7111
load("PD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
PD_expr <- PD_expr[PD_expr$status_for_analysis == "PD",]
dim(PD_expr) # 739 7059 
load("FTD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
FTD_expr <- FTD_expr[FTD_expr$status_for_analysis == "FTD",]
dim(FTD_expr) #  46 7020
load("DLB_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
DLB_expr <- DLB_expr[DLB_expr$status_for_analysis == "DLB",]
dim(DLB_expr) # 37 7111

# Read in the FDR proteins for each disease that are overlapping across 3 cohorts (Knight-ADRC, PPMI, SADRC) and are non-correlated
load("CSF_disease_specific_All_FDR_proteins_AD_PD_DLB_FTD_Overlapping_NonCorr_AgeSexPlatePCs.RData")
length(AD_FDR_nonCorr) # 1958
length(DLB_FDR_nonCorr) # 1267
length(FTD_FDR_nonCorr) # 1715
length(PD_FDR_nonCorr) # 85

novel_proteins <- na.omit(unique(c(AD_FDR_nonCorr, DLB_FDR_nonCorr, FTD_FDR_nonCorr, PD_FDR_nonCorr)))
length(novel_proteins) # 2638

AD_expr <- AD_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(AD_expr) # 1883 2643
row.names(AD_expr) <- AD_expr$UniquePhenoID

DLB_expr <- DLB_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(DLB_expr) # 37 2643
row.names(DLB_expr) <- DLB_expr$UniquePhenoID

FTD_expr <- FTD_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(FTD_expr) # 46 2643
row.names(FTD_expr) <- FTD_expr$UniquePhenoID

PD_expr <- PD_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(PD_expr) # 739 2643
row.names(PD_expr) <- PD_expr$UniquePhenoID

complete_expr <- rbind(AD_expr, DLB_expr, FTD_expr,  PD_expr)
dim(complete_expr) # 2705 2643

# Impute the missing values
data_imputation <- function(expression_data) {
 expr_matrix <- expression_data[,6:ncol(expression_data)]
 expr_matrix_nonNA <- as.data.frame(as.matrix(sapply(expr_matrix, as.numeric)))
 count.na = apply(is.na(expr_matrix), 2, sum)
 set.seed(1)
 for (i in which(count.na!=0)) { # bootstrapping with replacement
  index = is.na(expr_matrix_nonNA[,i])
  expr_matrix_nonNA[index, i] = sample(expr_matrix_nonNA[!index,i], sum(index), replace=T)
 }
 expression_data <- as.data.frame(cbind(expression_data[,1:5], expr_matrix_nonNA))
 print("any NA in the imputed data:")
 print(sum(is.na(expression_data)[,6:ncol(expression_data)]))
 return(expression_data)
}
complete_expr_imputed <- data_imputation(complete_expr)
dim(complete_expr_imputed) # 2705 2643
any(is.na(complete_expr_imputed)) # TRUE (One sample has Age missing, other has Sex missing)
complete_expr_imputed <- na.omit(complete_expr_imputed)
dim(complete_expr_imputed) # 2703 2643
save(complete_expr_imputed, file="CSF_AD_DLB_FTD_PD_Imputed_Expression_AgeSexPlate2PCs.RData")

# subset for AD novel proteins
AD_expr_imputed <- complete_expr_imputed[complete_expr_imputed$status_for_analysis %in% c("AD", "CO_ATneg"),]
table(AD_expr_imputed$level)
#    0    1
#  726 1156
table(AD_expr_imputed$Sex)
# Female   Male
#    974    908
selected_cols <- c("UniquePhenoID", "level", "Age", "Sex", AD_FDR_nonCorr)
length(selected_cols) # 1962
AD_expr_subset <- AD_expr_imputed[, selected_cols]
dim(AD_expr_subset) # 1882 1962
row.names(AD_expr_subset) <- AD_expr_subset$UniquePhenoID


## Lasso Model Generation (Iterative; n=50)

# comparison <- "ADvsCO"
iterative_lasso_modeling <- function(expr_data, comparison) {
#Define a "not in" opperator
`%!in%` <- Negate(`%in%`)
  
# Expression matrix
quantMatrix <- expr_data[,5:ncol(expr_data)]
dim(quantMatrix) # 1762  362
quantMatrix <- as.matrix(quantMatrix)
quantMatrix[1:4,1:4]

# Phenotype
pheno <- expr_data[,1:4]
colnames(pheno)[1] <- "sample_ID"
dim(pheno) # 1763    4
head(pheno)

#List the proteins
proteins <- colnames(quantMatrix)
length(proteins) # 362
prot <- as.data.frame(proteins)
names(prot) <- "Protein"
protWeight_df <- as.data.frame(prot)
colnames(protWeight_df) <- c("Proteins")
dim(protWeight_df) # 362   1
aucdf <- data.frame()

for (i in 1:50)
  {
    set.seed(i)
    #Set up testing and training data
    trainSamples <- c(sample(pheno$sample_ID[pheno$level=='1'], round(length(pheno$level[pheno$level=='1'])*0.7)),
                    sample(pheno$sample_ID[pheno$level=='0'], round(length(pheno$level[pheno$level=='0'])*0.7)))
    trainX <- quantMatrix[rownames(quantMatrix)%in%trainSamples,]
    trainX <- as.matrix(trainX)
    trainY <- as.numeric(pheno$level[pheno$sample_ID %in% trainSamples])
    #Build a lasso model
    set.seed(567)
    cvOut <- cv.glmnet(trainX, trainY, family="binomial",alpha=1)
    lamMin <- cvOut$lambda.min
    #Rebuild the model
    lassoBest <- glmnet(trainX, trainY, family="binomial",alpha=1, lambda = lamMin)
    coefProt <- coef(lassoBest) %>% as.matrix() %>% as_tibble(rownames = "Proteins")
    sigProt <- coefProt[coefProt$s0 != 0,]
    protWeight_df <- merge(protWeight_df, sigProt, by="Proteins", all.x = T)
    testX <- quantMatrix[rownames(quantMatrix)%!in%trainSamples,]
    testY <- as.numeric(pheno$level[pheno$sample_ID %!in% trainSamples])
    preds <- as.data.frame(predict(lassoBest, newx = testX, type = "response"))
    t <- ci.auc(testY, as.numeric(preds$s0), conf.level = 0.9)
    t <- data.frame(low = round(t[1],3), AUC = round(t[2],3), high = round(t[3],3), Proteins = (nrow(sigProt)-1))
    aucdf <- rbind(aucdf,t)
    print(i)
}
# Save predictions and performance obtained from lasso model training
dim(protWeight_df) # 362  51
dim(aucdf) # 50  4
# hist(aucdf$Proteins)
print("summary of proteins # selected in each of the 50 iterations:")
print(summary(aucdf$Proteins))
write.csv(protWeight_df, paste0("ProteinWeights_50runs_", comparison, "_AgeSexPlatePCs.csv", sep=""), row.names = F)
write.csv(aucdf, paste0("ModelPerformance_50runs_", comparison, "_AgeSexPlatePCs.csv", sep=""), row.names = F)

col_names <- c("Proteins", seq(from=1,to=50))
colnames(protWeight_df) <- col_names
protWeight_df$count <- rowSums(!is.na(select(protWeight_df, 2:51)))
dim(protWeight_df) # 362  52
protWeight_df <- protWeight_df[order(protWeight_df$count, decreasing=T),]
protWeight_df <<- protWeight_df[order(protWeight_df$count, decreasing=T),]
print("summary of most important protein selection:")
print(summary(protWeight_df$count))
write.csv(protWeight_df, paste0("ProteinWeights_Sorted_50runs_", comparison, "_AgeSexPlatePCs.csv", sep=""), row.names = F)

# Carlos asked me to try the model with top 10 analytes as they appeared in 45 of the 50 iterations.
iterativeModel_10analytes <<- protWeight_df[1:10,]$Proteins
print("Top 10 proteins selected for final Lasso model:")
print(iterativeModel_10analytes)
save(iterativeModel_10analytes, file=paste0("iterativeModel_10analytes_", comparison,"_AgeSexPlatePCs.RData", sep=""))
}

iterative_lasso_modeling(AD_expr_subset, "ADvsCO")

iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
iterativeModel_10analytes
# [1] "X11634.32_O43665"  "X12853.112_Q9NZR1" "X13118.5_Q9H4F8"   "X14157.21_P62258"  "X6521.35_P47972"   "X8070.88_P78352"   "X9900.36_P12036"
# [8] "X12501.10_O75347"  "X13109.82_Q7Z3B1"  "X8358.30_P30048"
# Analyte annotation
annotation <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_AD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(annotation) # 7099    9
annotation[annotation$Analyte %in% iterativeModel_10analytes,]$Symbol
#  [1] "YWHAE" "TMOD2" "SMOC1" "DLG4"  "PRDX3" "TBCA"  "NPTX2" "NEFH"  "NEGR1" "RGS10"


## Iterative Class Balanced Training {70%} and Testing {30%}

load("CSF_Complete_DiseaseSpecific_Imputed_Expression_DFs_For_CSFPredModel_AgeSexPlatePCs.RData")
protWeight_df <- read.csv("ProteinWeights_Sorted_50runs_ADvsCO_AgeSexPlatePCs.csv")
protWeight_df <- protWeight_df[order(protWeight_df$count, decreasing=T),]
iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
selected_cols <- c("level", "Age", "Sex", iterativeModel_10analytes)
length(selected_cols) # 13
AD_expr_subset <- AD_expr_imputed[, selected_cols]
dim(AD_expr_subset) # 1882 13
DLB_expr_subset <- DLB_expr_imputed[, selected_cols]
dim(DLB_expr_subset) # 763  13
FTD_expr_subset <- FTD_expr_imputed[, selected_cols]
dim(FTD_expr_subset) # 772  13
PD_expr_subset <- PD_expr_imputed[, selected_cols]
dim(PD_expr_subset) # 1464   13

# for saving performance info from the model
perf_df_ad <- c()
perf_df_ad_base <- c()
# for saving predictions from the model
iter_pred_ad_base <- list()
iter_pred_ad <- list()
iter_y_ad <- list()

# other dementias
perf_df_DLB <- c()
iter_pred_DLB <- list()
iter_y_DLB <- list()
perf_df_FTD <- c()
iter_pred_FTD <- list()
iter_y_FTD <- list()
perf_df_PD <- c()
iter_pred_PD <- list()
iter_y_PD <- list()

# round(table(AD_expr_subset$level)[1]*0.70, 0) # 508
# round(table(AD_expr_subset$level)[1]*0.30, 0) # 218
table(AD_expr_subset$level)
#   0    1
# 726 1156
min(table(AD_expr_subset$level)) # 726

for (i in 1:100) {
  set.seed(i)
  print(i)  
  # Separate data by levels # In each iteration 70% data is training (equal number of AD and CO) and 30% data is testing (equal number of AD and CO)
  level_0 <- AD_expr_subset[AD_expr_subset$level == 0, ]
  level_1 <- AD_expr_subset[AD_expr_subset$level == 1, ]
  # Sample rows for level = 1
  train_level_1 <- level_1[sample(1:nrow(level_1), size = round(0.7 * min(table(AD_expr_subset$level)))), ]
  test_level_1 <- level_1[!rownames(level_1) %in% rownames(train_level_1), ] # maker sure test has no train samples
  test_level_1 <- test_level_1[sample(1:nrow(test_level_1), size = round(0.3 * min(table(AD_expr_subset$level)))), ] # downsample test dataset
  # Sample rows for level = 0
  train_level_0 <- level_0[sample(1:nrow(level_0), size = nrow(train_level_1)), ]
  test_level_0 <- level_0[!rownames(level_0) %in% rownames(train_level_0), ] # maker sure test has no train samples
  test_level_0 <- test_level_0[sample(1:nrow(test_level_0), size = nrow(test_level_1)), ] # downsample test dataset
  # Combine training and testing sets
  Train_set <- rbind(train_level_1, train_level_0)
  Test_set <- rbind(test_level_1, test_level_0)
  
  # Baseline model
  bl_disc <- glm(as.formula(level ~ Age + Sex), data = Train_set, family = 'binomial')
  base_pred_val_ad <- predict(bl_disc, newdata = Test_set, type='response')
  ROC_base_val_ad <- roc(as.formula(Test_set$level ~ base_pred_val_ad), plot = FALSE, print.auc = FALSE)
  
  # 10-protein prediction model
  all_disc <<- glm(level ~ ., data = Train_set, family = 'binomial')
  pred_val_glm_ad <- predict(all_disc, newdata = Test_set, type='response')
  ROC_val_ad <- roc(as.formula(Test_set$level ~ pred_val_glm_ad), plot = FALSE, print.auc = FALSE)

  iter_pred_ad_base[[i]] <- as.numeric(base_pred_val_ad)
  iter_pred_ad[[i]] <- as.numeric(pred_val_glm_ad)
  iter_y_ad[[i]] <- Test_set$level

  # each performance result for iterations, if there is any ties, get the first one only
  perf_df_ad <- rbind(perf_df_ad, 
                       coords(ROC_val_ad, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_ad$auc, .before = 1) %>% 
                         summarise_all(first))
  perf_df_ad_base <- rbind(perf_df_ad_base, 
                            coords(ROC_base_val_ad, x = "best", best.method="youden", input = "threshold", 
                                   ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                              mutate(iteration = i, auc = ROC_base_val_ad$auc, .before = 1) %>% 
                              summarise_all(first))
  
  # Predictions in other dementia
  # DLB
  pred_val_glm_DLB <- predict(all_disc, newdata = DLB_expr_subset, type='response') # Replace "DLB_expr_subset" with "DLB_test" for class-balanced test data 
  ROC_val_DLB <- roc(as.formula(DLB_expr_subset$level ~ pred_val_glm_DLB), plot = FALSE, print.auc = FALSE)
  iter_pred_DLB[[i]] <- as.numeric(pred_val_glm_DLB)
  iter_y_DLB[[i]] <- DLB_expr_subset$level
  perf_df_DLB <- rbind(perf_df_DLB, 
                       coords(ROC_val_DLB, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_DLB$auc, .before = 1) %>%  summarise_all(first))
  # FTD
  pred_val_glm_FTD <- predict(all_disc, newdata = FTD_expr_subset, type='response')
  ROC_val_FTD <- roc(as.formula(FTD_expr_subset$level ~ pred_val_glm_FTD), plot = FALSE, print.auc = FALSE)
  iter_pred_FTD[[i]] <- as.numeric(pred_val_glm_FTD)
  iter_y_FTD[[i]] <- FTD_expr_subset$level
  perf_df_FTD <- rbind(perf_df_FTD, 
                       coords(ROC_val_FTD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_FTD$auc, .before = 1) %>%  summarise_all(first))
  # PD
  pred_val_glm_PD <- predict(all_disc, newdata = PD_expr_subset, type='response')
  ROC_val_PD <- roc(as.formula(PD_expr_subset$level ~ pred_val_glm_PD), plot = FALSE, print.auc = FALSE)
  iter_pred_PD[[i]] <- as.numeric(pred_val_glm_PD)
  iter_y_PD[[i]] <- PD_expr_subset$level
  perf_df_PD <- rbind(perf_df_PD, 
                       coords(ROC_val_PD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_PD$auc, .before = 1) %>% summarise_all(first))
}
roc.test(ROC_base_val_ad, ROC_val_ad, method="delong") # Z = -10.214, p-value < 2.2e-16, 95 percent confidence interval: -0.3065268 -0.2078239
roc.test(ROC_val_ad, ROC_val_DLB, method="delong") # p-value = 0.03652
roc.test(ROC_val_ad, ROC_val_FTD, method="delong") # p-value = 3.834e-07
roc.test(ROC_val_ad, ROC_val_PD, method="delong") # p-value < 2.2e-16

pred_ad <- prediction(iter_pred_ad, iter_y_ad)
pred_ad_base <- prediction(iter_pred_ad_base, iter_y_ad)
perf_ad <- performance(pred_ad, "tpr", "fpr")
perf_ad_base <- performance(pred_ad_base, "tpr", "fpr")
mean(perf_df_ad$auc) # 0.9478238
sd(perf_df_ad$auc) # 0.008406664
mean(perf_df_ad_base$auc) # 0.6964152
sd(perf_df_ad_base$auc) # 0.02268007
summary(perf_df_ad$auc)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.9173  0.9429  0.9483  0.9478  0.9530  0.9652
dim(perf_df_ad) # 100   8
dim(perf_df_ad_base) # 100   8
# get AD probablity for all individuals
out_prob_AD <- predict(all_disc, newdata = AD_expr_subset, type='response')
out_prob_AD <- as.data.frame(out_prob_AD)
out_prob_AD$UniquePhenoID <- row.names(out_prob_AD)
dim(out_prob_AD) # 2066    2
out_prob_AD_mean <- inner_join(out_prob_AD, AD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_AD, UniquePhenoID, status_for_analysis)
dim(out_prob_AD_mean) # 2066    3

png("AUC_Lasso_10_IterativeModel_ADvsCO_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
plot(perf_ad, avg = "vertical", 
     spread.estimate="stderror", lwd=6, lty=1, plotCI.col = alpha("blue", 0.5), col = alpha("blue", 0.5), main = c("ADvsCO"))
text(0.4, 0.25, paste("Test Data (30%): ",  round(mean(perf_df_ad$auc), 2), sep=""), col=("blue"), cex=0.85, pos=4)
plot(perf_ad_base, avg = "vertical", 
     spread.estimate="stderror", lty=2, add = T, plotCI.col = alpha("blue", 0.5), col = alpha("blue", 0.5))
text(0.4, 0.20, paste("Baseline: ", round(mean(perf_df_ad_base$auc), 2), sep=""), col=alpha("blue", 0.5), cex=0.75, pos=4)
dev.off()
save(perf_ad, perf_df_ad, perf_ad_base, perf_df_ad_base, 
    file="AUC_Lasso_10_IterativeModel_ADvsCO_AgeSexProjectPCs_Iterative_ClassBalanced.RData")

# prediction/performance in other dementias
# DLB
pred_DLB <- prediction(iter_pred_DLB, iter_y_DLB)
perf_DLB <- performance(pred_DLB, "tpr", "fpr")
mean(perf_df_DLB$auc) # 0.8512322
sd(perf_df_DLB$auc) # 0.009232644
dim(perf_df_DLB) # 100   8
out_prob_DLB_mean <- sapply(1:length(iter_pred_DLB[[1]]), function(i) {
  mean(sapply(iter_pred_DLB, function(x) x[i]))
})
out_prob_DLB_mean <- as.data.frame(cbind(row.names(DLB_expr_subset), out_prob_DLB_mean))
colnames(out_prob_DLB_mean) <- c("UniquePhenoID", "out_prob_DLB")
out_prob_DLB_mean <- inner_join(out_prob_DLB_mean, DLB_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_DLB, UniquePhenoID, status_for_analysis)

# FTD
pred_FTD <- prediction(iter_pred_FTD, iter_y_FTD)
perf_FTD <- performance(pred_FTD, "tpr", "fpr")
mean(perf_df_FTD$auc) # 0.6771574
sd(perf_df_FTD$auc) # 0.01392985
dim(perf_df_FTD) # 100   8
out_prob_FTD_mean <- sapply(1:length(iter_pred_FTD[[1]]), function(i) {
  mean(sapply(iter_pred_FTD, function(x) x[i]))
})
out_prob_FTD_mean <- as.data.frame(cbind(row.names(FTD_expr_subset), out_prob_FTD_mean))
colnames(out_prob_FTD_mean) <- c("UniquePhenoID", "out_prob_FTD")
out_prob_FTD_mean <- inner_join(out_prob_FTD_mean, FTD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_FTD, UniquePhenoID, status_for_analysis)

# PD
pred_PD <- prediction(iter_pred_PD, iter_y_PD)
perf_PD <- performance(pred_PD, "tpr", "fpr")
mean(perf_df_PD$auc) # 0.7585169
sd(perf_df_PD$auc) # 0.01154217
dim(perf_df_PD) # 100   8
out_prob_PD_mean <- sapply(1:length(iter_pred_PD[[1]]), function(i) {
  mean(sapply(iter_pred_PD, function(x) x[i]))
})
out_prob_PD_mean <- as.data.frame(cbind(row.names(PD_expr_subset), out_prob_PD_mean))
colnames(out_prob_PD_mean) <- c("UniquePhenoID", "out_prob_PD")
out_prob_PD_mean <- inner_join(out_prob_PD_mean, PD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_PD, UniquePhenoID, status_for_analysis)

png("AUC_Lasso_10_IterativeModel_ADvsCO_Testing_Other_NDs_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
# xlab='Average Specificity', ylab='Average Sensitivity', xlim=c(1,0), 
plot(perf_DLB, avg = "vertical",  
     spread.estimate="stderror", lwd=6, plotCI.col = alpha("orange", 0.5), col = alpha("orange", 0.5), main = c("AD model testing across other dementia"))
text(0.4, 0.05, paste("DLB: ",  round(mean(perf_df_DLB$auc), 2), sep=""), col=("orange"), cex=0.85, pos=4)
plot(perf_FTD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("green4", 0.5), col = alpha("green4", 0.5))
text(0.4, 0.10, paste("FTD: ",  round(mean(perf_df_FTD$auc), 2), sep=""), col=("green4"), cex=0.80, pos=4)
plot(perf_PD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("red", 0.5), col = alpha("red", 0.5))
text(0.4, 0.15, paste("PD: ",  round(mean(perf_df_PD$auc), 2), sep=""), col=("red"), cex=0.80, pos=4)
dev.off()
save(perf_DLB, perf_df_DLB, perf_FTD, perf_df_FTD, perf_PD, perf_df_PD, 
    file="AUC_Lasso_10_IterativeModel_ADvsCO_Testing_Other_NDs_AgeSexProjectPCs_Iterative_ClassBalanced.RData")

colnames(out_prob_DLB_mean)[1] <- "out_prob_AD"
colnames(out_prob_FTD_mean)[1] <- "out_prob_AD"
colnames(out_prob_PD_mean)[1] <- "out_prob_AD"
AD_Probablities <- rbind(out_prob_AD_mean, out_prob_DLB_mean, out_prob_FTD_mean, out_prob_PD_mean)
dim(AD_Probablities) # 4881    3
names(AD_Probablities)[3] <- "Final_Status"
length(unique(AD_Probablities$UniquePhenoID)) # 2703
AD_Probablities <- AD_Probablities[order(AD_Probablities$UniquePhenoID, AD_Probablities$out_prob_AD),]
AD_Probablities <- AD_Probablities[!duplicated(AD_Probablities$UniquePhenoID),]
AD_Probablities <- na.omit(AD_Probablities)
AD_Probablities$out_prob_AD <- as.numeric(as.character(AD_Probablities$out_prob_AD))
dim(AD_Probablities) # 2703    3
table(AD_Probablities$Final_Status)
#   AD  CO_ATneg     DLB      FTD       PD
# 1156       726      37       46      738
save(AD_Probablities, file="Plasma_AD_Probability_N3009_AgeSexPlatePCs_Iterative_ClassBalanced.RData")
png("Probability_Violinplot_for_Plasma_AD_Signature_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=100, height=100, res=1000)
ggplot(AD_Probablities, aes(x=Final_Status, y=out_prob_AD, color=Final_Status)) + geom_violin(trim = TRUE) + 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x ="Clinical Diagnosis", y = "Predicted AD Probability") + ggtitle("Probability Distribution for AD signature")
dev.off()

# Calculate mean and standard deviation for each column
calculate_mean_sd <- function(df) {
  df <- df[,-1]
  # Ensure only numeric columns are included
  numeric_cols <- df[sapply(df, is.numeric)]
  # Calculate mean and SD for each numeric column
  formatted_stats <- sapply(numeric_cols, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    paste0(round(mean_val, 2), "(", round(sd_val, 2), ")")
  })
  # Convert to a single-row dataframe
  result <- as.data.frame(t(formatted_stats))
  colnames(result) <- colnames(numeric_cols)
  return(result)
}
calculate_mean_sd(perf_df_ad)
#          auc  threshold   accuracy        npv        ppv sensitivity specificity
# 1 0.95(0.01) 0.57(0.09) 0.89(0.01) 0.87(0.02) 0.92(0.02)  0.87(0.03)  0.92(0.03)
calculate_mean_sd(perf_df_ad_base)
#         auc  threshold   accuracy        npv        ppv sensitivity specificity
# 1 0.7(0.02) 0.54(0.04) 0.66(0.02) 0.65(0.03) 0.68(0.04)  0.61(0.09)  0.71(0.09)
calculate_mean_sd(perf_df_DLB)
#       auc  threshold  accuracy     npv        ppv sensitivity specificity
# 1 0.85(0.01) 0.44(0.11) 0.86(0.04) 0.98(0) 0.24(0.06)  0.72(0.05)  0.87(0.05)
calculate_mean_sd(perf_df_FTD)
#          auc  threshold   accuracy     npv        ppv sensitivity specificity
# 1 0.68(0.01) 0.68(0.05) 0.92(0.01) 0.96(0) 0.36(0.07)   0.4(0.01)  0.95(0.01)
calculate_mean_sd(perf_df_PD)
#          auc  threshold   accuracy        npv        ppv sensitivity specificity
# 1 0.76(0.01) 0.16(0.03) 0.69(0.01) 0.7(0.02) 0.69(0.02)  0.71(0.04)  0.67(0.04)



## PD ##

## Lasso Model Generation (Iterative; n=50)

comparison <- "PDvsCO"
iterative_lasso_modeling(PD_expr_subset, "PDvsCO")

protWeight_df$count[1:15]
# [1] 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50
iterativeModel_10analytes <- protWeight_df$Proteins[11:20]
iterativeModel_10analytes
# [1] "X4413.3_P03973"  "X4801.13_P22079" "X6371.50_O43827" "X6496.60_P80370" "X6586.19_O75078" "X7200.4_Q9ULH4"  "X8358.30_P30048" "X8804.39_Q9P218"
# [9] "X8866.53_Q9NXS2" "X9025.5_Q9NZD4"
# Analyte annotation
annotation <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_PD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(annotation) # 7047    9
annotation[annotation$Analyte %in% iterativeModel_10analytes,]$Symbol
#  [1] "LPO"     "DLK1"    "ADAM11"  "PRDX3"   "ANGPTL7" "COL20A1" "LRFN2"   "QPCTL"   "SLPI"    "AHSP"


## Iterative Class Balanced Training {70%} and Testing {30%}

load("CSF_Complete_DiseaseSpecific_Imputed_Expression_DFs_For_CSFPredModel_AgeSexPlatePCs.RData")
protWeight_df <- read.csv("ProteinWeights_Sorted_50runs_PDvsCO_AgeSexPlatePCs.csv")
protWeight_df <- protWeight_df[order(protWeight_df$count, decreasing=T),]
protWeight_df$count[1:15]
# [1] 50 50 50 50 50 50 50 50 50 50 50 50 50 50 50
iterativeModel_10analytes <- protWeight_df$Proteins[11:20]
selected_cols <- c("level", "Age", "Sex", iterativeModel_10analytes)
length(selected_cols) # 13
AD_expr_subset <- AD_expr_imputed[, selected_cols]
dim(AD_expr_subset) # 1882 13
DLB_expr_subset <- DLB_expr_imputed[, selected_cols]
dim(DLB_expr_subset) # 763  13
FTD_expr_subset <- FTD_expr_imputed[, selected_cols]
dim(FTD_expr_subset) # 772  13
PD_expr_subset <- PD_expr_imputed[, selected_cols]
dim(PD_expr_subset) # 1464   13

# for saving performance info from the model
perf_df_pd <- c()
perf_df_pd_base <- c()
# for saving predictions from the model
iter_pred_pd_base <- list()
iter_pred_pd <- list()
iter_y_pd <- list()

# other dementias
perf_df_DLB <- c()
iter_pred_DLB <- list()
iter_y_DLB <- list()
perf_df_FTD <- c()
iter_pred_FTD <- list()
iter_y_FTD <- list()
perf_df_AD <- c()
iter_pred_AD <- list()
iter_y_AD <- list()

table(PD_expr_subset$level)
#   0    1
# 726  738
min(table(PD_expr_subset$level)) # 779

for (i in 1:100) {
  set.seed(i)
  print(i)  
  # Separate data by levels # In each iteration 70% data is training (equal number of AD and CO) and 30% data is testing (equal number of AD and CO)
  level_0 <- PD_expr_subset[PD_expr_subset$level == 0, ]
  level_1 <- PD_expr_subset[PD_expr_subset$level == 1, ]
  # Sample rows for level = 1
  train_level_1 <- level_1[sample(1:nrow(level_1), size = round(0.7 * min(table(PD_expr_subset$level)))), ]
  test_level_1 <- level_1[!rownames(level_1) %in% rownames(train_level_1), ] # maker sure test has no train samples
  test_level_1 <- test_level_1[sample(1:nrow(test_level_1), size = round(0.3 * min(table(PD_expr_subset$level)))), ] # downsample test dataset
  # Sample rows for level = 0
  train_level_0 <- level_0[sample(1:nrow(level_0), size = nrow(train_level_1)), ]
  test_level_0 <- level_0[!rownames(level_0) %in% rownames(train_level_0), ] # maker sure test has no train samples
  test_level_0 <- test_level_0[sample(1:nrow(test_level_0), size = nrow(test_level_1)), ] # downsample test dataset
  # Combine training and testing sets
  Train_set <- rbind(train_level_1, train_level_0)
  Test_set <- rbind(test_level_1, test_level_0)
  print(nrow(Train_set))
  print(nrow(Test_set))
  
  # Baseline model
  bl_disc <- glm(as.formula(level ~ Age + Sex), data = Train_set, family = 'binomial')
  base_pred_val_pd <- predict(bl_disc, newdata = Test_set, type='response')
  ROC_base_val_pd <- roc(as.formula(Test_set$level ~ base_pred_val_pd), plot = FALSE, print.auc = FALSE)
  
  # 10-protein prediction model
  all_disc <<- glm(level ~ ., data = Train_set, family = 'binomial')
  pred_val_glm_pd <- predict(all_disc, newdata = Test_set, type='response')
  ROC_val_pd <- roc(as.formula(Test_set$level ~ pred_val_glm_pd), plot = FALSE, print.auc = FALSE)

  iter_pred_pd_base[[i]] <- as.numeric(base_pred_val_pd)
  iter_pred_pd[[i]] <- as.numeric(pred_val_glm_pd)
  iter_y_pd[[i]] <- Test_set$level

  # each performance result for iterations, if there is any ties, get the first one only
  perf_df_pd <- rbind(perf_df_pd, 
                       coords(ROC_val_pd, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_pd$auc, .before = 1) %>% 
                         summarise_all(first))
  perf_df_pd_base <- rbind(perf_df_pd_base, 
                            coords(ROC_base_val_pd, x = "best", best.method="youden", input = "threshold", 
                                   ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                              mutate(iteration = i, auc = ROC_base_val_pd$auc, .before = 1) %>% 
                              summarise_all(first))
  
  # Predictions in other dementia
  # DLB
  pred_val_glm_DLB <- predict(all_disc, newdata = DLB_expr_subset, type='response')
  ROC_val_DLB <- roc(as.formula(DLB_expr_subset$level ~ pred_val_glm_DLB), plot = FALSE, print.auc = FALSE)
  iter_pred_DLB[[i]] <- as.numeric(pred_val_glm_DLB)
  iter_y_DLB[[i]] <- DLB_expr_subset$level
  perf_df_DLB <- rbind(perf_df_DLB, 
                       coords(ROC_val_DLB, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_DLB$auc, .before = 1) %>%  summarise_all(first))
  # FTD
  pred_val_glm_FTD <- predict(all_disc, newdata = FTD_expr_subset, type='response')
  ROC_val_FTD <- roc(as.formula(FTD_expr_subset$level ~ pred_val_glm_FTD), plot = FALSE, print.auc = FALSE)
  iter_pred_FTD[[i]] <- as.numeric(pred_val_glm_FTD)
  iter_y_FTD[[i]] <- FTD_expr_subset$level
  perf_df_FTD <- rbind(perf_df_FTD, 
                       coords(ROC_val_FTD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_FTD$auc, .before = 1) %>%  summarise_all(first))
  # AD
  pred_val_glm_AD <- predict(all_disc, newdata = AD_expr_subset, type='response')
  ROC_val_AD <- roc(as.formula(AD_expr_subset$level ~ pred_val_glm_AD), plot = FALSE, print.auc = FALSE)
  iter_pred_AD[[i]] <- as.numeric(pred_val_glm_AD)
  iter_y_AD[[i]] <- AD_expr_subset$level
  perf_df_AD <- rbind(perf_df_AD, 
                       coords(ROC_val_AD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_AD$auc, .before = 1) %>% summarise_all(first))
}
# Train_set = 1016 (508 PD, 508 CO)
# Test_set = 436 (218 PD, 218 CO)

pred_pd <- prediction(iter_pred_pd, iter_y_pd)
pred_pd_base <- prediction(iter_pred_pd_base, iter_y_pd)
perf_pd <- performance(pred_pd, "tpr", "fpr")
perf_pd_base <- performance(pred_pd_base, "tpr", "fpr")

mean(perf_df_pd$auc) # 0.806996
sd(perf_df_pd$auc) # 0.01753645
mean(perf_df_pd_base$auc) # 0.6630638
sd(perf_df_pd_base$auc) # 0.02174895
dim(perf_df_pd) # 100   8
dim(perf_df_pd_base) # 100   8
# get AD probablity for all individuals
out_prob_PD <- predict(all_disc, newdata = PD_expr_subset, type='response')
out_prob_PD <- as.data.frame(out_prob_PD)
out_prob_PD$UniquePhenoID <- row.names(out_prob_PD)
dim(out_prob_PD) # 1464    2
out_prob_PD_mean <- inner_join(out_prob_PD, PD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_PD, UniquePhenoID, status_for_analysis)
dim(out_prob_PD_mean) # 1464    3

png("AUC_Lasso_10_IterativeModel_PDvsCO_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
plot(perf_pd, avg = "vertical", 
     spread.estimate="stderror", lwd=6, lty=1, plotCI.col = alpha("red", 0.5), col = alpha("red", 0.5), main = c("PDvsCO"))
text(0.4, 0.25, paste("Test Data (30%): ",  round(mean(perf_df_pd$auc), 2), sep=""), col=("red"), cex=0.85, pos=4)
plot(perf_pd_base, avg = "vertical", 
     spread.estimate="stderror", lty=2, add = T, plotCI.col = alpha("red", 0.5), col = alpha("red", 0.5))
text(0.4, 0.20, paste("Baseline: ", round(mean(perf_df_pd_base$auc), 2), sep=""), col=alpha("red", 0.5), cex=0.75, pos=4)
dev.off()
save(perf_pd, perf_df_pd, perf_pd_base, perf_df_pd_base, 
    file="AUC_Lasso_10_IterativeModel_PDvsCO_AgeSexProjectPCs_Iterative_ClassBalanced.RData")

# prediction/performance in other dementias
# DLB
pred_DLB <- prediction(iter_pred_DLB, iter_y_DLB)
perf_DLB <- performance(pred_DLB, "tpr", "fpr")
mean(perf_df_DLB$auc) # 0.7331881
sd(perf_df_DLB$auc) # 0.009337895
dim(perf_df_DLB) # 100   8
out_prob_DLB_mean <- sapply(1:length(iter_pred_DLB[[1]]), function(i) {
  mean(sapply(iter_pred_DLB, function(x) x[i]))
})
out_prob_DLB_mean <- as.data.frame(cbind(row.names(DLB_expr_subset), out_prob_DLB_mean))
colnames(out_prob_DLB_mean) <- c("UniquePhenoID", "out_prob_DLB")
out_prob_DLB_mean <- inner_join(out_prob_DLB_mean, DLB_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_DLB, UniquePhenoID, status_for_analysis)

# FTD
pred_FTD <- prediction(iter_pred_FTD, iter_y_FTD)
perf_FTD <- performance(pred_FTD, "tpr", "fpr")
mean(perf_df_FTD$auc) # 0.7386798
sd(perf_df_FTD$auc) # 0.008761934
dim(perf_df_FTD) # 100   8
out_prob_FTD_mean <- sapply(1:length(iter_pred_FTD[[1]]), function(i) {
  mean(sapply(iter_pred_FTD, function(x) x[i]))
})
out_prob_FTD_mean <- as.data.frame(cbind(row.names(FTD_expr_subset), out_prob_FTD_mean))
colnames(out_prob_FTD_mean) <- c("UniquePhenoID", "out_prob_FTD")
out_prob_FTD_mean <- inner_join(out_prob_FTD_mean, FTD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_FTD, UniquePhenoID, status_for_analysis)

# AD
pred_AD <- prediction(iter_pred_AD, iter_y_AD)
perf_AD <- performance(pred_AD, "tpr", "fpr")
mean(perf_df_AD$auc) # 0.6896655
sd(perf_df_AD$auc) # 0.0101771
dim(perf_df_AD) # 100   8
out_prob_AD_mean <- sapply(1:length(iter_pred_AD[[1]]), function(i) {
  mean(sapply(iter_pred_AD, function(x) x[i]))
})
out_prob_AD_mean <- as.data.frame(cbind(row.names(AD_expr_subset), out_prob_AD_mean))
colnames(out_prob_AD_mean) <- c("UniquePhenoID", "out_prob_AD")
out_prob_AD_mean <- inner_join(out_prob_AD_mean, AD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_AD, UniquePhenoID, status_for_analysis)

png("AUC_Lasso_10_IterativeModel_PDvsCO_Testing_Other_NDs_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
# xlab='Average Specificity', ylab='Average Sensitivity', xlim=c(1,0), 
plot(perf_DLB, avg = "vertical",  
     spread.estimate="stderror", lwd=6, plotCI.col = alpha("orange", 0.5), col = alpha("orange", 0.5), main = c("PD model testing across other dementia"))
text(0.4, 0.05, paste("DLB: ",  round(mean(perf_df_DLB$auc), 2), sep=""), col=("orange"), cex=0.85, pos=4)
plot(perf_FTD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("green4", 0.5), col = alpha("green4", 0.5))
text(0.4, 0.10, paste("FTD: ",  round(mean(perf_df_FTD$auc), 2), sep=""), col=("green4"), cex=0.80, pos=4)
plot(perf_AD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("blue", 0.5), col = alpha("blue", 0.5))
text(0.4, 0.15, paste("AD: ",  round(mean(perf_df_AD$auc), 2), sep=""), col=("blue"), cex=0.80, pos=4)
dev.off()
save(perf_DLB, perf_df_DLB, perf_FTD, perf_df_FTD, perf_AD, perf_df_AD, 
    file="AUC_Lasso_10_IterativeModel_PDvsCO_Testing_Other_NDs_AgeSexProjectPCs_Iterative_ClassBalanced.RData")

colnames(out_prob_DLB_mean)[1] <- "out_prob_PD"
colnames(out_prob_FTD_mean)[1] <- "out_prob_PD"
colnames(out_prob_AD_mean)[1] <- "out_prob_PD"
PD_Probablities <- rbind(out_prob_PD_mean, out_prob_DLB_mean, out_prob_FTD_mean, out_prob_AD_mean)
dim(PD_Probablities) # 4881    3
names(PD_Probablities)[3] <- "Final_Status"
length(unique(PD_Probablities$UniquePhenoID)) # 2703
PD_Probablities <- PD_Probablities[order(PD_Probablities$UniquePhenoID, PD_Probablities$out_prob_PD),]
PD_Probablities <- PD_Probablities[!duplicated(PD_Probablities$UniquePhenoID),]
PD_Probablities <- na.omit(PD_Probablities)
PD_Probablities$out_prob_PD <- as.numeric(as.character(PD_Probablities$out_prob_PD))
dim(PD_Probablities) # 2703    3
table(PD_Probablities$Final_Status)
#      AD CO_ATneg      DLB      FTD       PD
#    1156      726       37       46      738
save(PD_Probablities, file="Plasma_PD_Probability_N2703_AgeSexPlatePCs_Iterative_ClassBalanced.RData")
png("Probability_Violinplot_for_Plasma_PD_Signature_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=100, height=100, res=1000)
ggplot(PD_Probablities, aes(x=Final_Status, y=out_prob_PD, color=Final_Status)) + geom_violin(trim = TRUE) + 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x ="Clinical Diagnosis", y = "Predicted PD Probability") + ggtitle("Probability Distribution for PD signature")
dev.off()

# Calculate mean and standard deviation for each column
calculate_mean_sd <- function(df) {
  df <- df[,-1]
  # Ensure only numeric columns are included
  numeric_cols <- df[sapply(df, is.numeric)]
  # Calculate mean and SD for each numeric column
  formatted_stats <- sapply(numeric_cols, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    paste0(round(mean_val, 2), "(", round(sd_val, 2), ")")
  })
  # Convert to a single-row dataframe
  result <- as.data.frame(t(formatted_stats))
  colnames(result) <- colnames(numeric_cols)
  return(result)
}
calculate_mean_sd(perf_df_pd)
calculate_mean_sd(perf_df_pd_base)
calculate_mean_sd(perf_df_DLB)
calculate_mean_sd(perf_df_FTD)
calculate_mean_sd(perf_df_AD)
#         auc  threshold   accuracy        npv        ppv sensitivity specificity
# 1 0.81(0.02) 0.45(0.06) 0.75(0.02) 0.78(0.04) 0.73(0.03)  0.79(0.06)  0.71(0.06)
#         auc  threshold   accuracy        npv        ppv sensitivity specificity
# 1 0.66(0.02) 0.49(0.04) 0.63(0.02) 0.64(0.03) 0.64(0.04)   0.64(0.1)  0.62(0.11)
#         auc  threshold   accuracy     npv        ppv sensitivity specificity
# 1 0.73(0.01) 0.44(0.04) 0.68(0.04) 0.98(0) 0.11(0.01)  0.75(0.04)  0.68(0.04)
#         auc  threshold   accuracy     npv        ppv sensitivity specificity
# 1 0.74(0.01) 0.51(0.06) 0.74(0.06) 0.97(0) 0.14(0.02)  0.65(0.06)  0.74(0.07)
#         auc  threshold   accuracy        npv        ppv sensitivity specificity
# 1 0.69(0.01) 0.41(0.03) 0.65(0.01) 0.54(0.01) 0.75(0.01)  0.65(0.03)  0.65(0.03)



## DLB ##

## Lasso Model Generation (Iterative; n=50)

comparison <- "DLBvsCO"
iterative_lasso_modeling(DLB_expr_subset, "DLBvsCO")
protWeight_df$count[1:15]
# [1] 26 25 25 24 24 23 23 21 20 20 19 19 18 18 15
iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
iterativeModel_10analytes
# [1] "X4272.46_P06744"   "X10082.251_P07196" "X9940.35_Q4G0W2"   "X11633.89_O95433"  "X9900.36_P12036"   "X3351.1_Q13554"    "X5354.11_P05783"
# [8] "X2973.15_P16671"   "X2578.67_P13500"   "X3291.30_P06734"

# Analyte annotation
annotation <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_DLB_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(annotation) # 7099   9
annotation[annotation$Analyte %in% iterativeModel_10analytes,]$Symbol
#  [1] "DUSP28" "GPI"    "CAMK2B" "NEFL"   "CCL2"   "AHSA1"  "FCER2"  "NEFH"   "KRT18"  "CD36"


## Iterative Class Balanced Training {70%} and Testing {30%}

load("CSF_Complete_DiseaseSpecific_Imputed_Expression_DFs_For_CSFPredModel_AgeSexPlatePCs.RData")
protWeight_df <- read.csv("ProteinWeights_Sorted_50runs_DLBvsCO_AgeSexProjectPCs.csv")
protWeight_df <- protWeight_df[order(protWeight_df$count, decreasing=T),]
protWeight_df$count[1:15]
# [1] 26 25 25 24 24 23 23 21 20 20 19 19 18 18 15
iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
selected_cols <- c("level", "Age", "Sex", iterativeModel_10analytes)
length(selected_cols) # 13
AD_expr_subset <- AD_expr_imputed[, selected_cols]
dim(AD_expr_subset) # 1882 13
DLB_expr_subset <- DLB_expr_imputed[, selected_cols]
dim(DLB_expr_subset) # 763  13
FTD_expr_subset <- FTD_expr_imputed[, selected_cols]
dim(FTD_expr_subset) # 772  13
PD_expr_subset <- PD_expr_imputed[, selected_cols]
dim(PD_expr_subset) # 1464   13

# for saving performance info from the model
perf_df_dlb <- c()
perf_df_dlb_base <- c()
# for saving predictions from the model
iter_pred_dlb_base <- list()
iter_pred_dlb <- list()
iter_y_dlb <- list()

# other dementias
perf_df_AD <- c()
iter_pred_AD <- list()
iter_y_AD <- list()
perf_df_FTD <- c()
iter_pred_FTD <- list()
iter_y_FTD <- list()
perf_df_PD <- c()
iter_pred_PD <- list()
iter_y_PD <- list()

table(DLB_expr_subset$level)
#   0    1
# 726   37
min(table(DLB_expr_subset$level)) # 37

for (i in 1:100) {
  set.seed(i)
  print(i)
  # Separate data by levels # In each iteration 70% data is training (equal number of AD and CO) and 30% data is testing (equal number of AD and CO)
  level_0 <- DLB_expr_subset[DLB_expr_subset$level == 0, ]
  level_1 <- DLB_expr_subset[DLB_expr_subset$level == 1, ]
  # Sample rows for level = 1
  train_level_1 <- level_1[sample(1:nrow(level_1), size = round(0.7 * min(table(DLB_expr_subset$level)))), ]
  test_level_1 <- level_1[!rownames(level_1) %in% rownames(train_level_1), ] # maker sure test has no train samples
  test_level_1 <- test_level_1[sample(1:nrow(test_level_1), size = round(0.3 * min(table(DLB_expr_subset$level)))), ] # downsample test dataset
  # Sample rows for level = 0
  train_level_0 <- level_0[sample(1:nrow(level_0), size = nrow(train_level_1)), ]
  test_level_0 <- level_0[!rownames(level_0) %in% rownames(train_level_0), ] # maker sure test has no train samples
  test_level_0 <- test_level_0[sample(1:nrow(test_level_0), size = nrow(test_level_1)), ] # downsample test dataset
  # Combine training and testing sets
  Train_set <- rbind(train_level_1, train_level_0)
  Test_set <- rbind(test_level_1, test_level_0)
  print(nrow(Train_set))
  print(nrow(Test_set))
  
  # Baseline model
  bl_disc <- glm(as.formula(level ~ Age + Sex), data = Train_set, family = 'binomial')
  base_pred_val_dlb <- predict(bl_disc, newdata = Test_set, type='response')
  ROC_base_val_dlb <- roc(as.formula(Test_set$level ~ base_pred_val_dlb), plot = FALSE, print.auc = FALSE)
  
  # 10-protein prediction model
  all_disc <<- glm(level ~ ., data = Train_set, family = 'binomial')
  pred_val_glm_dlb <- predict(all_disc, newdata = Test_set, type='response')
  ROC_val_dlb <- roc(as.formula(Test_set$level ~ pred_val_glm_dlb), plot = FALSE, print.auc = FALSE)

  iter_pred_dlb_base[[i]] <- as.numeric(base_pred_val_dlb)
  iter_pred_dlb[[i]] <- as.numeric(pred_val_glm_dlb)
  iter_y_dlb[[i]] <- Test_set$level

  # each performance result for iterations, if there is any ties, get the first one only
  perf_df_dlb <- rbind(perf_df_dlb, 
                       coords(ROC_val_dlb, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_dlb$auc, .before = 1) %>% 
                         summarise_all(first))
  perf_df_dlb_base <- rbind(perf_df_dlb_base, 
                            coords(ROC_base_val_dlb, x = "best", best.method="youden", input = "threshold", 
                                   ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                              mutate(iteration = i, auc = ROC_base_val_dlb$auc, .before = 1) %>% 
                              summarise_all(first))
  
  # Predictions in other dementia
  # AD
  pred_val_glm_AD <- predict(all_disc, newdata = AD_expr_subset, type='response')
  ROC_val_AD <- roc(as.formula(AD_expr_subset$level ~ pred_val_glm_AD), plot = FALSE, print.auc = FALSE)
  iter_pred_AD[[i]] <- as.numeric(pred_val_glm_AD)
  iter_y_AD[[i]] <- AD_expr_subset$level
  perf_df_AD <- rbind(perf_df_AD, 
                       coords(ROC_val_AD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_AD$auc, .before = 1) %>%  summarise_all(first))
  # FTD
  pred_val_glm_FTD <- predict(all_disc, newdata = FTD_expr_subset, type='response')
  ROC_val_FTD <- roc(as.formula(FTD_expr_subset$level ~ pred_val_glm_FTD), plot = FALSE, print.auc = FALSE)
  iter_pred_FTD[[i]] <- as.numeric(pred_val_glm_FTD)
  iter_y_FTD[[i]] <- FTD_expr_subset$level
  perf_df_FTD <- rbind(perf_df_FTD, 
                       coords(ROC_val_FTD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_FTD$auc, .before = 1) %>%  summarise_all(first))
  # PD
  pred_val_glm_PD <- predict(all_disc, newdata = PD_expr_subset, type='response')
  ROC_val_PD <- roc(as.formula(PD_expr_subset$level ~ pred_val_glm_PD), plot = FALSE, print.auc = FALSE)
  iter_pred_PD[[i]] <- as.numeric(pred_val_glm_PD)
  iter_y_PD[[i]] <- PD_expr_subset$level
  perf_df_PD <- rbind(perf_df_PD, 
                       coords(ROC_val_PD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_PD$auc, .before = 1) %>% summarise_all(first))
}
# Train_set = 52 (26 DLB, 26 CO)
# Test_set = 22 (11 DLB, 11 CO)

pred_dlb <- prediction(iter_pred_dlb, iter_y_dlb)
pred_dlb_base <- prediction(iter_pred_dlb_base, iter_y_dlb)
perf_dlb <- performance(pred_dlb, "tpr", "fpr")
perf_dlb_base <- performance(pred_dlb_base, "tpr", "fpr")

mean(perf_df_dlb$auc) # 0.9154132
sd(perf_df_dlb$auc) # 0.07834076
mean(perf_df_dlb_base$auc) # 0.6542562
sd(perf_df_dlb_base$auc) # 0.1073385
dim(perf_df_dlb) # 100   8
dim(perf_df_dlb_base) # 100   8
# get DLB probablity for all individuals
out_prob_DLB <- predict(all_disc, newdata = DLB_expr_subset, type='response')
out_prob_DLB <- as.data.frame(out_prob_DLB)
out_prob_DLB$UniquePhenoID <- row.names(out_prob_DLB)
dim(out_prob_DLB) # 763   2
out_prob_DLB_mean <- inner_join(out_prob_DLB, DLB_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_DLB, UniquePhenoID, status_for_analysis)
dim(out_prob_DLB_mean) # 763   3

png("AUC_Lasso_10_IterativeModel_DLBvsCO_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
plot(perf_dlb, avg = "vertical", 
     spread.estimate="stderror", lwd=6, lty=1, plotCI.col = alpha("orange", 0.5), col = alpha("orange", 0.5), main = c("DLBvsCO"))
text(0.4, 0.25, paste("Test Data (30%): ",  round(mean(perf_df_dlb$auc), 2), sep=""), col=("orange"), cex=0.85, pos=4)
plot(perf_dlb_base, avg = "vertical", 
     spread.estimate="stderror", lty=2, add = T, plotCI.col = alpha("orange", 0.5), col = alpha("orange", 0.5))
text(0.4, 0.20, paste("Baseline: ", round(mean(perf_df_dlb_base$auc), 2), sep=""), col=alpha("orange", 0.5), cex=0.75, pos=4)
dev.off()
save(perf_dlb, perf_df_dlb, perf_dlb_base, perf_df_dlb_base, 
    file="AUC_Lasso_10_IterativeModel_DLBvsCO_AgeSexProjectPCs_Iterative_ClassBalanced.RData")

# prediction/performance in other dementias
# PD
pred_PD <- prediction(iter_pred_PD, iter_y_PD)
perf_PD <- performance(pred_PD, "tpr", "fpr")
mean(perf_df_PD$auc) # 0.6416734
sd(perf_df_PD$auc) # 0.03296469
dim(perf_df_PD) # 100   8
out_prob_PD_mean <- sapply(1:length(iter_pred_PD[[1]]), function(i) {
  mean(sapply(iter_pred_PD, function(x) x[i]))
})
out_prob_PD_mean <- as.data.frame(cbind(row.names(PD_expr_subset), out_prob_PD_mean))
colnames(out_prob_PD_mean) <- c("UniquePhenoID", "out_prob_PD")
out_prob_PD_mean <- inner_join(out_prob_PD_mean, PD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_PD, UniquePhenoID, status_for_analysis)

# FTD
pred_FTD <- prediction(iter_pred_FTD, iter_y_FTD)
perf_FTD <- performance(pred_FTD, "tpr", "fpr")
mean(perf_df_FTD$auc) # 0.7382785
sd(perf_df_FTD$auc) # 0.03573307
dim(perf_df_FTD) # 100   8
out_prob_FTD_mean <- sapply(1:length(iter_pred_FTD[[1]]), function(i) {
  mean(sapply(iter_pred_FTD, function(x) x[i]))
})
out_prob_FTD_mean <- as.data.frame(cbind(row.names(FTD_expr_subset), out_prob_FTD_mean))
colnames(out_prob_FTD_mean) <- c("UniquePhenoID", "out_prob_FTD")
out_prob_FTD_mean <- inner_join(out_prob_FTD_mean, FTD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_FTD, UniquePhenoID, status_for_analysis)

# AD
pred_AD <- prediction(iter_pred_AD, iter_y_AD)
perf_AD <- performance(pred_AD, "tpr", "fpr")
mean(perf_df_AD$auc) # 0.7654994
sd(perf_df_AD$auc) # 0.03538492
dim(perf_df_AD) # 100   8
out_prob_AD_mean <- sapply(1:length(iter_pred_AD[[1]]), function(i) {
  mean(sapply(iter_pred_AD, function(x) x[i]))
})
out_prob_AD_mean <- as.data.frame(cbind(row.names(AD_expr_subset), out_prob_AD_mean))
colnames(out_prob_AD_mean) <- c("UniquePhenoID", "out_prob_AD")
out_prob_AD_mean <- inner_join(out_prob_AD_mean, AD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_AD, UniquePhenoID, status_for_analysis)

png("AUC_Lasso_10_IterativeModel_DLBvsCO_Testing_Other_NDs_AgeSexProjectPCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
# xlab='Average Specificity', ylab='Average Sensitivity', xlim=c(1,0), 
plot(perf_AD, avg = "vertical",  
     spread.estimate="stderror", lwd=6, plotCI.col = alpha("blue", 0.5), col = alpha("blue", 0.5), main = c("DLB model testing across other dementia"))
text(0.4, 0.15, paste("AD: ",  round(mean(perf_df_AD$auc), 2), sep=""), col=("blue"), cex=0.85, pos=4)
plot(perf_FTD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("green4", 0.5), col = alpha("green4", 0.5))
text(0.4, 0.10, paste("FTD: ",  round(mean(perf_df_FTD$auc), 2), sep=""), col=("green4"), cex=0.80, pos=4)
plot(perf_PD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("red", 0.5), col = alpha("red", 0.5))
text(0.4, 0.05, paste("PD: ",  round(mean(perf_df_PD$auc), 2), sep=""), col=("red"), cex=0.80, pos=4)
dev.off()
save(perf_AD, perf_df_AD, perf_FTD, perf_df_FTD, perf_PD, perf_df_PD, 
    file="AUC_Lasso_10_IterativeModel_DLBvsCO_Testing_Other_NDs_AgeSexProjectPCs_Iterative_ClassBalanced.RData")

colnames(out_prob_AD_mean)[1] <- "out_prob_DLB"
colnames(out_prob_FTD_mean)[1] <- "out_prob_DLB"
colnames(out_prob_PD_mean)[1] <- "out_prob_DLB"
DLB_Probablities <- rbind(out_prob_DLB_mean, out_prob_AD_mean, out_prob_FTD_mean, out_prob_PD_mean)
dim(DLB_Probablities) # 4881    3
names(DLB_Probablities)[3] <- "Final_Status"
length(unique(DLB_Probablities$UniquePhenoID)) # 2703
DLB_Probablities <- DLB_Probablities[order(DLB_Probablities$UniquePhenoID, DLB_Probablities$out_prob_DLB),]
DLB_Probablities <- DLB_Probablities[!duplicated(DLB_Probablities$UniquePhenoID),]
DLB_Probablities <- na.omit(DLB_Probablities)
DLB_Probablities$out_prob_DLB <- as.numeric(as.character(DLB_Probablities$out_prob_DLB))
dim(DLB_Probablities) # 2703    3
table(DLB_Probablities$Final_Status)
#      AD CO_ATneg      DLB      FTD       PD
#    1156      726       37       46      738
save(DLB_Probablities, file="CSF_DLB_Probability_N2703_AgeSexPlatePCs_Iterative_ClassBalanced.RData")
png("Probability_Violinplot_for_CSF_DLB_Signature_AgeSexPlatePCs_Iterative_ClassBalanced.png", units="mm", width=100, height=100, res=1000)
ggplot(DLB_Probablities, aes(x=Final_Status, y=out_prob_DLB, color=Final_Status)) + geom_violin(trim = TRUE) + 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x ="Clinical Diagnosis", y = "Predicted DLB Probability") + ggtitle("Probability Distribution for DLB signature")
dev.off()

# Calculate mean and standard deviation for each column
calculate_mean_sd <- function(df) {
  df <- df[,-1]
  # Ensure only numeric columns are included
  numeric_cols <- df[sapply(df, is.numeric)]
  # Calculate mean and SD for each numeric column
  formatted_stats <- sapply(numeric_cols, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    paste0(round(mean_val, 2), "(", round(sd_val, 2), ")")
  })
  # Convert to a single-row dataframe
  result <- as.data.frame(t(formatted_stats))
  colnames(result) <- colnames(numeric_cols)
  return(result)
}
calculate_mean_sd(perf_df_dlb)
calculate_mean_sd(perf_df_dlb_base)
calculate_mean_sd(perf_df_PD)
calculate_mean_sd(perf_df_FTD)
calculate_mean_sd(perf_df_AD)
#         auc  threshold  accuracy        npv        ppv sensitivity specificity
#1 0.92(0.08) 0.42(0.39) 0.9(0.07) 0.91(0.09) 0.91(0.09)   0.91(0.1)   0.9(0.12)
#         auc  threshold  accuracy        npv       ppv sensitivity specificity
#1 0.65(0.11) 0.47(0.12) 0.7(0.07) 0.77(0.14) 0.7(0.11)  0.76(0.19)  0.64(0.19)
#         auc  threshold   accuracy        npv        ppv sensitivity specificity
#1 0.64(0.03) 0.08(0.25) 0.63(0.03) 0.59(0.03) 0.72(0.04)  0.43(0.08)  0.83(0.06)
#         auc  threshold   accuracy     npv        ppv sensitivity specificity
#1 0.74(0.04) 0.35(0.44) 0.84(0.05) 0.97(0) 0.23(0.07)  0.59(0.07)  0.86(0.06)
#         auc  threshold   accuracy       npv        ppv sensitivity specificity
#1 0.77(0.04) 0.03(0.17) 0.72(0.03) 0.6(0.04) 0.86(0.03)  0.65(0.07)  0.82(0.05)



## FTD ##

## Lasso Model Generation (Iterative; n=50)

comparison <- "FTDvsCO"
iterative_lasso_modeling(FTD_expr_subset, "FTDvsCO")
protWeight_df$count[1:15]
# [1] 30 29 29 25 24 23 23 22 22 22 21 20 20 19 19
iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
iterativeModel_10analytes
# [1] "X2962.50_P12272"   "X10880.38_P0C2L3"  "X11330.15_P67870"  "X10082.251_P07196" "X13745.10_Q13563"  "X2864.2_Q02750"    "X7892.132_A2RU67"
# [8] "X14144.3_Q7L7L0"   "X6214.84_P32881"   "X8839.4_O75431"

# Analyte annotation
annotation <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_FTD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(annotation) # 7008   9
annotation[annotation$Analyte %in% iterativeModel_10analytes,]$Symbol
#  [1] "CSNK2B"  "H2AW"    "FAM163B" "MAP2K1"  "FAM234B" "NEFL"    "PKD2"    "PTHLH"   "IFNA8"   "MTX2"


## Iterative Class Balanced Training {70%} and Testing {30%}

load("CSF_Complete_DiseaseSpecific_Imputed_Expression_DFs_For_CSFPredModel_AgeSexPlatePCs.RData")
protWeight_df <- read.csv("ProteinWeights_Sorted_50runs_FTDvsCO_AgeSexProjectPCs.csv")
protWeight_df <- protWeight_df[order(protWeight_df$count, decreasing=T),]
protWeight_df$count[1:15]
# [1] 30 29 29 25 24 23 23 22 22 22 21 20 20 19 19
iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
selected_cols <- c("level", "Age", "Sex", iterativeModel_10analytes)
length(selected_cols) # 13
AD_expr_subset <- AD_expr_imputed[, selected_cols]
dim(AD_expr_subset) # 1882 13
DLB_expr_subset <- DLB_expr_imputed[, selected_cols]
dim(DLB_expr_subset) # 763  13
FTD_expr_subset <- FTD_expr_imputed[, selected_cols]
dim(FTD_expr_subset) # 772  13
PD_expr_subset <- PD_expr_imputed[, selected_cols]
dim(PD_expr_subset) # 1464   13

# for saving performance info from the model
perf_df_ftd <- c()
perf_df_ftd_base <- c()
# for saving predictions from the model
iter_pred_ftd_base <- list()
iter_pred_ftd <- list()
iter_y_ftd <- list()

# other dementias
perf_df_AD <- c()
iter_pred_AD <- list()
iter_y_AD <- list()
perf_df_DLB <- c()
iter_pred_DLB <- list()
iter_y_DLB <- list()
perf_df_PD <- c()
iter_pred_PD <- list()
iter_y_PD <- list()

table(FTD_expr_subset$level)
#   0    1
# 943   46
min(table(FTD_expr_subset$level)) # 46

for (i in 1:100) {
  set.seed(i)
  print(i)
  # Separate data by levels # In each iteration 70% data is training (equal number of AD and CO) and 30% data is testing (equal number of AD and CO)
  level_0 <- FTD_expr_subset[FTD_expr_subset$level == 0, ]
  level_1 <- FTD_expr_subset[FTD_expr_subset$level == 1, ]
  # Sample rows for level = 1
  train_level_1 <- level_1[sample(1:nrow(level_1), size = round(0.7 * min(table(FTD_expr_subset$level)))), ]
  test_level_1 <- level_1[!rownames(level_1) %in% rownames(train_level_1), ] # maker sure test has no train samples
  test_level_1 <- test_level_1[sample(1:nrow(test_level_1), size = round(0.3 * min(table(FTD_expr_subset$level)))), ] # downsample test dataset
  # Sample rows for level = 0
  train_level_0 <- level_0[sample(1:nrow(level_0), size = nrow(train_level_1)), ]
  test_level_0 <- level_0[!rownames(level_0) %in% rownames(train_level_0), ] # maker sure test has no train samples
  test_level_0 <- test_level_0[sample(1:nrow(test_level_0), size = nrow(test_level_1)), ] # downsample test dataset
  # Combine training and testing sets
  Train_set <- rbind(train_level_1, train_level_0)
  Test_set <- rbind(test_level_1, test_level_0)
  print(nrow(Train_set))
  print(nrow(Test_set))

  # Baseline model
  bl_disc <- glm(as.formula(level ~ Age + Sex), data = Train_set, family = 'binomial')
  base_pred_val_ftd <- predict(bl_disc, newdata = Test_set, type='response')
  ROC_base_val_ftd <- roc(as.formula(Test_set$level ~ base_pred_val_ftd), plot = FALSE, print.auc = FALSE)
  
  # 10-protein prediction model
  all_disc <<- glm(level ~ ., data = Train_set, family = 'binomial')
  pred_val_glm_ftd <- predict(all_disc, newdata = Test_set, type='response')
  ROC_val_ftd <- roc(as.formula(Test_set$level ~ pred_val_glm_ftd), plot = FALSE, print.auc = FALSE)

  iter_pred_ftd_base[[i]] <- as.numeric(base_pred_val_ftd)
  iter_pred_ftd[[i]] <- as.numeric(pred_val_glm_ftd)
  iter_y_ftd[[i]] <- Test_set$level

  # each performance result for iterations, if there is any ties, get the first one only
  perf_df_ftd <- rbind(perf_df_ftd, 
                       coords(ROC_val_ftd, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_ftd$auc, .before = 1) %>% 
                         summarise_all(first))
  perf_df_ftd_base <- rbind(perf_df_ftd_base, 
                            coords(ROC_base_val_ftd, x = "best", best.method="youden", input = "threshold", 
                                   ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                              mutate(iteration = i, auc = ROC_base_val_ftd$auc, .before = 1) %>% 
                              summarise_all(first))
  
  # Predictions in other dementia
  # AD
  pred_val_glm_AD <- predict(all_disc, newdata = AD_expr_subset, type='response')
  ROC_val_AD <- roc(as.formula(AD_expr_subset$level ~ pred_val_glm_AD), plot = FALSE, print.auc = FALSE)
  iter_pred_AD[[i]] <- as.numeric(pred_val_glm_AD)
  iter_y_AD[[i]] <- AD_expr_subset$level
  perf_df_AD <- rbind(perf_df_AD, 
                       coords(ROC_val_AD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_AD$auc, .before = 1) %>%  summarise_all(first))
  # DLB
  pred_val_glm_DLB <- predict(all_disc, newdata = DLB_expr_subset, type='response')
  ROC_val_DLB <- roc(as.formula(DLB_expr_subset$level ~ pred_val_glm_DLB), plot = FALSE, print.auc = FALSE)
  iter_pred_DLB[[i]] <- as.numeric(pred_val_glm_DLB)
  iter_y_DLB[[i]] <- DLB_expr_subset$level
  perf_df_DLB <- rbind(perf_df_DLB, 
                       coords(ROC_val_DLB, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_DLB$auc, .before = 1) %>%  summarise_all(first))
  # PD
  pred_val_glm_PD <- predict(all_disc, newdata = PD_expr_subset, type='response')
  ROC_val_PD <- roc(as.formula(PD_expr_subset$level ~ pred_val_glm_PD), plot = FALSE, print.auc = FALSE)
  iter_pred_PD[[i]] <- as.numeric(pred_val_glm_PD)
  iter_y_PD[[i]] <- PD_expr_subset$level
  perf_df_PD <- rbind(perf_df_PD, 
                       coords(ROC_val_PD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_PD$auc, .before = 1) %>% summarise_all(first))
}
# Train_set = 64 (32 FTD, 32 CO)
# Test_set = 28 (14 FTD, 14 CO)

pred_ftd <- prediction(iter_pred_ftd, iter_y_ftd)
pred_ftd_base <- prediction(iter_pred_ftd_base, iter_y_ftd)
perf_ftd <- performance(pred_ftd, "tpr", "fpr")
perf_ftd_base <- performance(pred_ftd_base, "tpr", "fpr")

mean(perf_df_ftd$auc) # 0.8951786
sd(perf_df_ftd$auc) # 0.06477022
mean(perf_df_ftd_base$auc) # 0.5932143
sd(perf_df_ftd_base$auc) # 0.07288022
dim(perf_df_ftd) # 100   8
dim(perf_df_ftd_base) # 100   8
# get FTD probablity for all individuals
out_prob_FTD <- predict(all_disc, newdata = FTD_expr_subset, type='response')
out_prob_FTD <- as.data.frame(out_prob_FTD)
out_prob_FTD$UniquePhenoID <- row.names(out_prob_FTD)
dim(out_prob_FTD) # 772    2
out_prob_FTD_mean <- inner_join(out_prob_FTD, FTD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_FTD, UniquePhenoID, status_for_analysis)
dim(out_prob_FTD_mean) # 772   3

png("AUC_Lasso_10_IterativeModel_FTDvsCO_AgeSexPlatePCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
plot(perf_ftd, avg = "vertical", 
     spread.estimate="stderror", lwd=6, lty=1, plotCI.col = alpha("green4", 0.5), col = alpha("green4", 0.5), main = c("FTDvsCO"))
text(0.4, 0.25, paste("Test Data (30%): ",  round(mean(perf_df_ftd$auc), 2), sep=""), col=("green4"), cex=0.85, pos=4)
plot(perf_ftd_base, avg = "vertical", 
     spread.estimate="stderror", lty=2, add = T, plotCI.col = alpha("green4", 0.5), col = alpha("green4", 0.5))
text(0.4, 0.20, paste("Baseline: ", round(mean(perf_df_ftd_base$auc), 2), sep=""), col=alpha("green4", 0.5), cex=0.75, pos=4)
dev.off()
save(perf_ftd, perf_df_ftd, perf_ftd_base, perf_df_ftd_base,
    file="AUC_Lasso_10_IterativeModel_FTDvsCO_AgeSexPlatePCs_Iterative_ClassBalanced.RData")

# prediction/performance in other dementias
# AD
pred_AD <- prediction(iter_pred_AD, iter_y_AD)
perf_AD <- performance(pred_AD, "tpr", "fpr")
mean(perf_df_AD$auc) # 0.5716528
sd(perf_df_AD$auc) # 0.03292279
dim(perf_df_AD) # 100   8
out_prob_AD_mean <- sapply(1:length(iter_pred_AD[[1]]), function(i) {
  mean(sapply(iter_pred_AD, function(x) x[i]))
})
out_prob_AD_mean <- as.data.frame(cbind(row.names(AD_expr_subset), out_prob_AD_mean))
colnames(out_prob_AD_mean) <- c("UniquePhenoID", "out_prob_AD")
out_prob_AD_mean <- inner_join(out_prob_AD_mean, AD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_AD, UniquePhenoID, status_for_analysis)
# DLB
pred_DLB <- prediction(iter_pred_DLB, iter_y_DLB)
perf_DLB <- performance(pred_DLB, "tpr", "fpr")
mean(perf_df_DLB$auc) # 0.7538672
sd(perf_df_DLB$auc) # 0.04991893
dim(perf_df_DLB) # 100   8
out_prob_DLB_mean <- sapply(1:length(iter_pred_DLB[[1]]), function(i) {
  mean(sapply(iter_pred_DLB, function(x) x[i]))
})
out_prob_DLB_mean <- as.data.frame(cbind(row.names(DLB_expr_subset), out_prob_DLB_mean))
colnames(out_prob_DLB_mean) <- c("UniquePhenoID", "out_prob_DLB")
out_prob_DLB_mean <- inner_join(out_prob_DLB_mean, DLB_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_DLB, UniquePhenoID, status_for_analysis)
# PD
pred_PD <- prediction(iter_pred_PD, iter_y_PD)
perf_PD <- performance(pred_PD, "tpr", "fpr")
mean(perf_df_PD$auc) # 0.6356982
sd(perf_df_PD$auc) # 0.04005721
dim(perf_df_PD) # 100   8
out_prob_PD_mean <- sapply(1:length(iter_pred_PD[[1]]), function(i) {
  mean(sapply(iter_pred_PD, function(x) x[i]))
})
out_prob_PD_mean <- as.data.frame(cbind(row.names(PD_expr_subset), out_prob_PD_mean))
colnames(out_prob_PD_mean) <- c("UniquePhenoID", "out_prob_PD")
out_prob_PD_mean <- inner_join(out_prob_PD_mean, PD_expr_imputed[,1:5], by="UniquePhenoID") %>% select(out_prob_PD, UniquePhenoID, status_for_analysis)

png("AUC_Lasso_10_IterativeModel_FTDvsCO_Testing_Other_NDs_AgeSexPlatePCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
# xlab='Average Specificity', ylab='Average Sensitivity', xlim=c(1,0), 
plot(perf_AD, avg = "vertical",  
     spread.estimate="stderror", lwd=6, plotCI.col = alpha("blue", 0.5), col = alpha("blue", 0.5), main = c("FTD model testing across other dementia"))
text(0.4, 0.15, paste("AD: ",  round(mean(perf_df_AD$auc), 2), sep=""), col=("blue"), cex=0.85, pos=4)
plot(perf_DLB, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("orange", 0.5), col = alpha("orange", 0.5))
text(0.4, 0.10, paste("DLB: ",  round(mean(perf_df_DLB$auc), 2), sep=""), col=("orange"), cex=0.80, pos=4)
plot(perf_PD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("red", 0.5), col = alpha("red", 0.5))
text(0.4, 0.05, paste("PD: ",  round(mean(perf_df_PD$auc), 2), sep=""), col=("red"), cex=0.80, pos=4)
dev.off()
save(perf_AD, perf_df_AD, perf_DLB, perf_df_DLB, perf_PD, perf_df_PD, 
    file="AUC_Lasso_10_IterativeModel_FTDvsCO_Testing_Other_NDs_AgeSexPlatePCs_Iterative_ClassBalanced.RData")

colnames(out_prob_AD_mean)[1] <- "out_prob_FTD"
colnames(out_prob_DLB_mean)[1] <- "out_prob_FTD"
colnames(out_prob_PD_mean)[1] <- "out_prob_FTD"
FTD_Probablities <- rbind(out_prob_FTD_mean, out_prob_AD_mean, out_prob_DLB_mean, out_prob_PD_mean)
dim(FTD_Probablities) # 4881    3
names(FTD_Probablities)[3] <- "Final_Status"
length(unique(FTD_Probablities$UniquePhenoID)) # 2703
FTD_Probablities <- FTD_Probablities[order(FTD_Probablities$UniquePhenoID, FTD_Probablities$out_prob_FTD),]
FTD_Probablities <- FTD_Probablities[!duplicated(FTD_Probablities$UniquePhenoID),]
FTD_Probablities <- na.omit(FTD_Probablities)
FTD_Probablities$out_prob_FTD <- as.numeric(as.character(FTD_Probablities$out_prob_FTD))
dim(FTD_Probablities) # 2703    3
table(FTD_Probablities$Final_Status)
#      AD CO_ATneg      DLB      FTD       PD
#    1156      726       37       46      738
save(FTD_Probablities, file="CSF_FTD_Probability_N2703_AgeSexPlatePCs_Iterative_ClassBalanced.RData")
png("Probability_Violinplot_for_CSF_FTD_Signature_AgeSexPlatePCs_Iterative_ClassBalanced.png", units="mm", width=100, height=100, res=1000)
ggplot(FTD_Probablities, aes(x=Final_Status, y=out_prob_FTD, color=Final_Status)) + geom_violin(trim = TRUE) + 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x ="Clinical Diagnosis", y = "Predicted FTD Probability") + ggtitle("Probability Distribution for FTD signature")
dev.off()

# Calculate mean and standard deviation for each column
calculate_mean_sd <- function(df) {
  df <- df[,-1]
  # Ensure only numeric columns are included
  numeric_cols <- df[sapply(df, is.numeric)]
  # Calculate mean and SD for each numeric column
  formatted_stats <- sapply(numeric_cols, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    paste0(round(mean_val, 2), "(", round(sd_val, 2), ")")
  })
  # Convert to a single-row dataframe
  result <- as.data.frame(t(formatted_stats))
  colnames(result) <- colnames(numeric_cols)
  return(result)
}
calculate_mean_sd(perf_df_ftd)
calculate_mean_sd(perf_df_ftd_base)
calculate_mean_sd(perf_df_PD)
calculate_mean_sd(perf_df_DLB)
calculate_mean_sd(perf_df_AD)
#        auc  threshold   accuracy       npv        ppv sensitivity specificity
#1 0.9(0.06) 0.42(0.33) 0.89(0.06) 0.9(0.08) 0.89(0.08)    0.9(0.1)  0.88(0.11)
#         auc  threshold   accuracy      npv       ppv sensitivity specificity
#1 0.59(0.07) 0.49(0.08) 0.65(0.05) 0.7(0.1) 0.66(0.1)   0.7(0.17)   0.59(0.2)
#         auc  threshold   accuracy       npv       ppv sensitivity specificity
#1 0.64(0.04) 0.04(0.17) 0.62(0.03) 0.6(0.04) 0.7(0.04)  0.46(0.13)  0.79(0.08)
#         auc  threshold   accuracy     npv       ppv sensitivity specificity
#1 0.75(0.05) 0.62(0.46) 0.85(0.05) 0.98(0) 0.2(0.05)  0.63(0.09)  0.86(0.05)
#         auc  threshold   accuracy        npv        ppv sensitivity specificity
#1 0.57(0.03) 0.23(0.41) 0.51(0.04) 0.43(0.02) 0.74(0.03)  0.31(0.09)  0.83(0.05)



## Common Prediction Model (for AD/PD/DLB/FTD) ##

# Read in the expression matrices for all disease
load("AD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
dim(AD_expr) # 1883 7111
load("PD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
PD_expr <- PD_expr[PD_expr$status_for_analysis == "PD",]
dim(PD_expr) # 739 7059
load("FTD_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
FTD_expr <- FTD_expr[FTD_expr$status_for_analysis == "FTD",]
dim(FTD_expr) #  46 7020
load("DLB_vs_CO_Zscore_DAA_AgeSexPlatePCs.RData")
DLB_expr <- DLB_expr[DLB_expr$status_for_analysis == "DLB",]
dim(DLB_expr) # 37 7111

# Read in the FDR proteins for each disease that are overlapping across 3 cohorts (Longs+Garfield, PPMI, SADRC) and are non-correlated
load("CSF_disease_specific_All_FDR_proteins_AD_PD_DLB_FTD_Overlapping_NonCorr_AgeSexPlatePCs.RData")
length(AD_FDR_nonCorr) # 1958
length(DLB_FDR_nonCorr) # 1267
length(FTD_FDR_nonCorr) # 1715
length(PD_FDR_nonCorr) # 85

# Novel overlapping FDR proteins across all disorders
novel_proteins <- Reduce(intersect, list(AD_FDR_nonCorr, DLB_FDR_nonCorr, FTD_FDR_nonCorr, PD_FDR_nonCorr))
length(novel_proteins) # 18

AD_expr <- AD_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(AD_expr) # 1883   23
row.names(AD_expr) <- AD_expr$UniquePhenoID

DLB_expr <- DLB_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(DLB_expr) # 37 23
row.names(DLB_expr) <- DLB_expr$UniquePhenoID

FTD_expr <- FTD_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(FTD_expr) # 46 23
row.names(FTD_expr) <- FTD_expr$UniquePhenoID

PD_expr <- PD_expr[,c("UniquePhenoID", "level", "Sex", "Age", "status_for_analysis", novel_proteins)]
dim(PD_expr) # 739 23
row.names(PD_expr) <- PD_expr$UniquePhenoID

complete_expr <- rbind(AD_expr, DLB_expr, FTD_expr,  PD_expr)
dim(complete_expr) # 2705   23

# Impute the missing values
data_imputation <- function(expression_data) {
 expr_matrix <- expression_data[,6:ncol(expression_data)]
 expr_matrix_nonNA <- as.data.frame(as.matrix(sapply(expr_matrix, as.numeric)))
 count.na = apply(is.na(expr_matrix), 2, sum)
 set.seed(1)
 for (i in which(count.na!=0)) { # bootstrapping with replacement
  index = is.na(expr_matrix_nonNA[,i])
  expr_matrix_nonNA[index, i] = sample(expr_matrix_nonNA[!index,i], sum(index), replace=T)
 }
 expression_data <- as.data.frame(cbind(expression_data[,1:5], expr_matrix_nonNA))
 print("any NA in the imputed data:")
 print(sum(is.na(expression_data)[,6:ncol(expression_data)]))
 return(expression_data)
}
complete_expr_imputed <- data_imputation(complete_expr)
dim(complete_expr_imputed) # 2705 23
any(is.na(complete_expr_imputed)) # TRUE (One sample has Age missing, other has Sex missing)
complete_expr_imputed <- na.omit(complete_expr_imputed)
dim(complete_expr_imputed) # 2703 23
save(complete_expr_imputed, file="CSF_AD_DLB_FTD_PD_Imputed_Expression_AgeSexPlate2PCs_CommonFDRProteins.RData")

selected_cols <- c("UniquePhenoID", "level", "Age", "Sex", novel_proteins)
length(selected_cols) # 22
complete_expr_subset <- complete_expr_imputed[, selected_cols]
dim(complete_expr_subset) # 2703   22
row.names(complete_expr_subset) <- complete_expr_subset$UniquePhenoID

## Lasso Model Generation (Iterative; n=50)

comparison <- "ALLvsCO"
iterative_lasso_modeling(complete_expr_subset, "ALLvsCO")
protWeight_df$count[1:15]
# [1] 50 50 50 50 50 50 50 50 50 50 50 49 45 29 28
iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
iterativeModel_10analytes
#  [1] "X12853.112_Q9NZR1" "X14283.12_P61106"  "X2597.8_P15692"    "X4140.3_P13232"    "X7178.59_Q9NTK1"   "X7200.4_Q9ULH4"    "X8002.27_Q9H3S3"
#  [8] "X8601.167_Q07954"  "X9025.5_Q9NZD4"    "X9202.309_P48047"

# Analyte annotation
AD_annot <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_AD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
DLB_annot <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_DLB_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
FTD_annot <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_FTD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
PD_annot <- read.table("CSF_LongsGar_PPMI_SADRC_Zscore_PD_vs_CO_DAA_AgeSexPlatePC12_AllAnalytes.txt", sep="\t", header=T, stringsAsFactors=F, quote="")
PD_annot[PD_annot$Analyte %in% iterativeModel_10analytes,]$Symbol
# [1] "VEGFA"   "LRP1"    "TMPRSS5" "IL7"     "LRFN2"   "TMOD2"   "RAB14"   "DEPP1"   "ATP5PO"  "AHSP"
# Check if selected proteins are in same direction for each disease
AD_annot[AD_annot$Analyte %in% iterativeModel_10analytes,][,c(3,5,7)]
DLB_annot[DLB_annot$Analyte %in% iterativeModel_10analytes,][,c(3,5,7)]
FTD_annot[FTD_annot$Analyte %in% iterativeModel_10analytes,][,c(3,5,7)]
PD_annot[PD_annot$Analyte %in% iterativeModel_10analytes,][,c(3,5,7)]
# consistent proteins: TMOD2 (+), LRFN2 (-), VEGFA (-),  RAB14 (-), ATP5PO (+), IL7 (-), TMPRSS5 (+)
# in-consistent proteins: LRP1, AHSP, DEPP1 (+ only in PD)


## Iterative Class Balanced Training {70%} and Testing {30%}

load("CSF_AD_DLB_FTD_PD_Imputed_Expression_AgeSexPlate2PCs_CommonFDRProteins.RData")
protWeight_df <- read.csv("ProteinWeights_Sorted_50runs_ALLvsCO_AgeSexPlatePCs.csv")
protWeight_df <- protWeight_df[order(protWeight_df$count, decreasing=T),]
protWeight_df$count[1:15]
# [1] 50 50 50 50 50 50 50 50 50 50 50 49 45 29 28
iterativeModel_10analytes <- protWeight_df$Proteins[1:10]
selected_cols <- c("level", "Age", "Sex", iterativeModel_10analytes)
length(selected_cols) # 13

complete_expr_subset <- complete_expr_imputed
complete_expr_imputed <- complete_expr_imputed[, selected_cols]
dim(complete_expr_imputed) # 2703   13

AD_expr_subset <- complete_expr_imputed[complete_expr_imputed$status_for_analysis %in% c("CO_ATneg", "AD"),]
AD_expr_subset <- AD_expr_subset[, selected_cols]
dim(AD_expr_subset) # 1882 13
DLB_expr_subset <- complete_expr_imputed[complete_expr_imputed$status_for_analysis %in% c("CO_ATneg", "DLB"),]
DLB_expr_subset <- DLB_expr_subset[, selected_cols]
dim(DLB_expr_subset) # 763  13
FTD_expr_imputed <- complete_expr_imputed[complete_expr_imputed$status_for_analysis %in% c("CO_ATneg", "FTD"),]
FTD_expr_subset <- FTD_expr_imputed[, selected_cols]
dim(FTD_expr_subset) # 772  13
PD_expr_imputed <- complete_expr_imputed[complete_expr_imputed$status_for_analysis %in% c("CO_ATneg", "PD"),]
PD_expr_subset <- PD_expr_imputed[, selected_cols]
dim(PD_expr_subset) # 1464   13

# for saving performance info from the model
perf_df_ftd <- c()
perf_df_ftd_base <- c()
# for saving predictions from the model
iter_pred_ftd_base <- list()
iter_pred_ftd <- list()
iter_y_ftd <- list()

# other dementias
perf_df_AD <- c()
iter_pred_AD <- list()
iter_y_AD <- list()
perf_df_DLB <- c()
iter_pred_DLB <- list()
iter_y_DLB <- list()
perf_df_PD <- c()
iter_pred_PD <- list()
iter_y_PD <- list()
perf_df_FTD <- c()
iter_pred_FTD <- list()
iter_y_FTD <- list()

table(complete_expr_imputed$level)
#   0    1
# 726 1977
min(table(complete_expr_imputed$level)) # 726

for (i in 1:100) {
  set.seed(i)
  print(i)
  # Separate data by levels # In each iteration 70% data is training (equal number of AD and CO) and 30% data is testing (equal number of AD and CO)
  level_0 <- complete_expr_imputed[complete_expr_imputed$level == 0, ]
  level_1 <- complete_expr_imputed[complete_expr_imputed$level == 1, ]
  # Sample rows for level = 1
  train_level_1 <- level_1[sample(1:nrow(level_1), size = round(0.7 * min(table(complete_expr_imputed$level)))), ]
  test_level_1 <- level_1[!rownames(level_1) %in% rownames(train_level_1), ] # maker sure test has no train samples
  test_level_1 <- test_level_1[sample(1:nrow(test_level_1), size = round(0.3 * min(table(complete_expr_imputed$level)))), ] # downsample test dataset
  # Sample rows for level = 0
  train_level_0 <- level_0[sample(1:nrow(level_0), size = nrow(train_level_1)), ]
  test_level_0 <- level_0[!rownames(level_0) %in% rownames(train_level_0), ] # maker sure test has no train samples
  test_level_0 <- test_level_0[sample(1:nrow(test_level_0), size = nrow(test_level_1)), ] # downsample test dataset
  # Combine training and testing sets
  Train_set <- rbind(train_level_1, train_level_0)
  Test_set <- rbind(test_level_1, test_level_0)
  print(nrow(Train_set))
  print(nrow(Test_set))
  # Baseline model
  bl_disc <- glm(as.formula(level ~ Age + Sex), data = Train_set, family = 'binomial')
  base_pred_val_ftd <- predict(bl_disc, newdata = Test_set, type='response')
  ROC_base_val_ftd <- roc(as.formula(Test_set$level ~ base_pred_val_ftd), plot = FALSE, print.auc = FALSE)
  
  # 10-protein prediction model
  all_disc <<- glm(level ~ ., data = Train_set, family = 'binomial')
  pred_val_glm_ftd <- predict(all_disc, newdata = Test_set, type='response')
  ROC_val_ftd <- roc(as.formula(Test_set$level ~ pred_val_glm_ftd), plot = FALSE, print.auc = FALSE)

  iter_pred_ftd_base[[i]] <- as.numeric(base_pred_val_ftd)
  iter_pred_ftd[[i]] <- as.numeric(pred_val_glm_ftd)
  iter_y_ftd[[i]] <- Test_set$level

  # each performance result for iterations, if there is any ties, get the first one only
  perf_df_ftd <- rbind(perf_df_ftd, 
                       coords(ROC_val_ftd, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_ftd$auc, .before = 1) %>% 
                         summarise_all(first))
  perf_df_ftd_base <- rbind(perf_df_ftd_base, 
                            coords(ROC_base_val_ftd, x = "best", best.method="youden", input = "threshold", 
                                   ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                              mutate(iteration = i, auc = ROC_base_val_ftd$auc, .before = 1) %>% 
                              summarise_all(first))
  
  # Predictions in other dementia
  # AD
  pred_val_glm_AD <- predict(all_disc, newdata = AD_expr_subset, type='response')
  ROC_val_AD <- roc(as.formula(AD_expr_subset$level ~ pred_val_glm_AD), plot = FALSE, print.auc = FALSE)
  iter_pred_AD[[i]] <- as.numeric(pred_val_glm_AD)
  iter_y_AD[[i]] <- AD_expr_subset$level
  perf_df_AD <- rbind(perf_df_AD, 
                       coords(ROC_val_AD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_AD$auc, .before = 1) %>%  summarise_all(first))
  # DLB
  pred_val_glm_DLB <- predict(all_disc, newdata = DLB_expr_subset, type='response')
  ROC_val_DLB <- roc(as.formula(DLB_expr_subset$level ~ pred_val_glm_DLB), plot = FALSE, print.auc = FALSE)
  iter_pred_DLB[[i]] <- as.numeric(pred_val_glm_DLB)
  iter_y_DLB[[i]] <- DLB_expr_subset$level
  perf_df_DLB <- rbind(perf_df_DLB, 
                       coords(ROC_val_DLB, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_DLB$auc, .before = 1) %>%  summarise_all(first))
  # PD
  pred_val_glm_PD <- predict(all_disc, newdata = PD_expr_subset, type='response')
  ROC_val_PD <- roc(as.formula(PD_expr_subset$level ~ pred_val_glm_PD), plot = FALSE, print.auc = FALSE)
  iter_pred_PD[[i]] <- as.numeric(pred_val_glm_PD)
  iter_y_PD[[i]] <- PD_expr_subset$level
  perf_df_PD <- rbind(perf_df_PD, 
                       coords(ROC_val_PD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_PD$auc, .before = 1) %>% summarise_all(first))
  # FTD
  pred_val_glm_FTD <- predict(all_disc, newdata = FTD_expr_subset, type='response')
  ROC_val_FTD <- roc(as.formula(FTD_expr_subset$level ~ pred_val_glm_FTD), plot = FALSE, print.auc = FALSE)
  iter_pred_FTD[[i]] <- as.numeric(pred_val_glm_FTD)
  iter_y_FTD[[i]] <- FTD_expr_subset$level
  perf_df_FTD <- rbind(perf_df_FTD, 
                       coords(ROC_val_FTD, x = "best", best.method="youden", input = "threshold", 
                              ret = c("threshold","auc", "acc","npv","ppv","sensitivity","specificity")) %>% 
                         mutate(iteration = i, auc = ROC_val_FTD$auc, .before = 1) %>%  summarise_all(first))
}
# Train_set = 1016 (508 CA, 508 CO)
# Test_set = 436 (218 CA, 218 CO)

pred_ftd <- prediction(iter_pred_ftd, iter_y_ftd)
pred_ftd_base <- prediction(iter_pred_ftd_base, iter_y_ftd)
perf_ftd <- performance(pred_ftd, "tpr", "fpr")
perf_ftd_base <- performance(pred_ftd_base, "tpr", "fpr")

mean(perf_df_ftd$auc) # 0.7805906
sd(perf_df_ftd$auc) # 0.02078624
mean(perf_df_ftd_base$auc) # 0.5656205
sd(perf_df_ftd_base$auc) # 0.02843641
dim(perf_df_ftd) # 100   8
dim(perf_df_ftd_base) # 100   8
# get ALL probablity for all individuals
out_prob_ALL <- predict(all_disc, newdata = complete_expr_imputed, type='response')
out_prob_ALL <- as.data.frame(out_prob_ALL)
out_prob_ALL$UniquePhenoID <- row.names(out_prob_ALL)
dim(out_prob_ALL) # 2703    2
out_prob_ALL_mean <- inner_join(out_prob_ALL, complete_expr_subset[,1:5], by="UniquePhenoID") %>% select(out_prob_ALL, UniquePhenoID, status_for_analysis)
dim(out_prob_ALL_mean) # 2703    3

png("AUC_Lasso_10_IterativeModel_ALLvsCO_AgeSexPlatePCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
plot(perf_ftd, avg = "vertical", 
     spread.estimate="stderror", lwd=6, lty=1, plotCI.col = alpha("brown", 0.5), col = alpha("brown", 0.5), main = c("ALLvsCO"))
text(0.4, 0.25, paste("Test Data (30%): ",  round(mean(perf_df_ftd$auc), 2), sep=""), col=("brown"), cex=0.85, pos=4)
plot(perf_ftd_base, avg = "vertical", 
     spread.estimate="stderror", lty=2, add = T, plotCI.col = alpha("brown", 0.5), col = alpha("brown", 0.5))
text(0.4, 0.20, paste("Baseline: ", round(mean(perf_df_ftd_base$auc), 2), sep=""), col=alpha("brown", 0.5), cex=0.75, pos=4)
dev.off()
save(perf_ftd, perf_df_ftd, perf_ftd_base, perf_df_ftd_base,
    file="AUC_Lasso_10_IterativeModel_ALLvsCO_AgeSexPlatePCs_Iterative_ClassBalanced.RData")

# prediction/performance in other dementias
# AD
pred_AD <- prediction(iter_pred_AD, iter_y_AD)
perf_AD <- performance(pred_AD, "tpr", "fpr")
mean(perf_df_AD$auc) # 0.8049084
sd(perf_df_AD$auc) # 0.005482738
dim(perf_df_AD) # 100   8
out_prob_AD_mean <- sapply(1:length(iter_pred_AD[[1]]), function(i) {
  mean(sapply(iter_pred_AD, function(x) x[i]))
})
out_prob_AD_mean <- as.data.frame(cbind(row.names(AD_expr_subset), out_prob_AD_mean))
colnames(out_prob_AD_mean) <- c("UniquePhenoID", "out_prob_AD")
out_prob_AD_mean <- inner_join(out_prob_AD_mean, complete_expr_subset[,1:5], by="UniquePhenoID") %>% select(out_prob_AD, UniquePhenoID, status_for_analysis)

# DLB
pred_DLB <- prediction(iter_pred_DLB, iter_y_DLB)
perf_DLB <- performance(pred_DLB, "tpr", "fpr")
mean(perf_df_DLB$auc) # 0.8110882
sd(perf_df_DLB$auc) # 0.01150495
dim(perf_df_DLB) # 100   8
out_prob_DLB_mean <- sapply(1:length(iter_pred_DLB[[1]]), function(i) {
  mean(sapply(iter_pred_DLB, function(x) x[i]))
})
out_prob_DLB_mean <- as.data.frame(cbind(row.names(DLB_expr_subset), out_prob_DLB_mean))
colnames(out_prob_DLB_mean) <- c("UniquePhenoID", "out_prob_DLB")
out_prob_DLB_mean <- inner_join(out_prob_DLB_mean, complete_expr_subset[,1:5], by="UniquePhenoID") %>% select(out_prob_DLB, UniquePhenoID, status_for_analysis)

# PD
pred_PD <- prediction(iter_pred_PD, iter_y_PD)
perf_PD <- performance(pred_PD, "tpr", "fpr")
mean(perf_df_PD$auc) # 0.747939
sd(perf_df_PD$auc) # 0.009526562
dim(perf_df_PD) # 100   8
out_prob_PD_mean <- sapply(1:length(iter_pred_PD[[1]]), function(i) {
  mean(sapply(iter_pred_PD, function(x) x[i]))
})
out_prob_PD_mean <- as.data.frame(cbind(row.names(PD_expr_subset), out_prob_PD_mean))
colnames(out_prob_PD_mean) <- c("UniquePhenoID", "out_prob_PD")
out_prob_PD_mean <- inner_join(out_prob_PD_mean, complete_expr_subset[,1:5], by="UniquePhenoID") %>% select(out_prob_PD, UniquePhenoID, status_for_analysis)

# FTD
pred_FTD <- prediction(iter_pred_FTD, iter_y_FTD)
perf_FTD <- performance(pred_FTD, "tpr", "fpr")
mean(perf_df_FTD$auc) # 0.7908507
sd(perf_df_FTD$auc) # 0.01293257
dim(perf_df_FTD) # 100   8
out_prob_FTD_mean <- sapply(1:length(iter_pred_FTD[[1]]), function(i) {
  mean(sapply(iter_pred_FTD, function(x) x[i]))
})
out_prob_FTD_mean <- as.data.frame(cbind(row.names(FTD_expr_subset), out_prob_FTD_mean))
colnames(out_prob_FTD_mean) <- c("UniquePhenoID", "out_prob_FTD")
out_prob_FTD_mean <- inner_join(out_prob_FTD_mean, complete_expr_subset[,1:5], by="UniquePhenoID") %>% select(out_prob_FTD, UniquePhenoID, status_for_analysis)


png("AUC_Lasso_10_IterativeModel_ALLvsCO_Testing_Other_NDs_AgeSexPlatePCs_Iterative_ClassBalanced.png", units="mm", width=120, height=120, res=1000)
# xlab='Average Specificity', ylab='Average Sensitivity', xlim=c(1,0), 
plot(perf_AD, avg = "vertical",  
     spread.estimate="stderror", lwd=6, plotCI.col = alpha("blue", 0.5), col = alpha("blue", 0.5), main = c("Generic model testing in all dementia"))
text(0.4, 0.20, paste("AD: ",  round(mean(perf_df_AD$auc), 2), sep=""), col=("blue"), cex=0.85, pos=4)
plot(perf_DLB, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("orange", 0.5), col = alpha("orange", 0.5))
text(0.4, 0.15, paste("DLB: ",  round(mean(perf_df_DLB$auc), 2), sep=""), col=("orange"), cex=0.80, pos=4)
plot(perf_FTD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("green4", 0.5), col = alpha("green4", 0.5))
text(0.4, 0.10, paste("FTD: ",  round(mean(perf_df_FTD$auc), 2), sep=""), col=("green4"), cex=0.80, pos=4)
plot(perf_PD, avg = "vertical", spread.estimate="stderror", lwd=6, add = T, plotCI.col = alpha("red", 0.5), col = alpha("red", 0.5))
text(0.4, 0.05, paste("PD: ",  round(mean(perf_df_PD$auc), 2), sep=""), col=("red"), cex=0.80, pos=4)
dev.off()
save(perf_AD, perf_df_AD, perf_DLB, perf_df_DLB, perf_PD, perf_df_PD, perf_FTD, perf_df_FTD, 
    file="AUC_Lasso_10_IterativeModel_ALLvsCO_Testing_Other_NDs_AgeSexPlatePCs_Iterative_ClassBalanced.RData")

colnames(out_prob_AD_mean)[1] <- "out_prob_ALL"
colnames(out_prob_DLB_mean)[1] <- "out_prob_ALL"
colnames(out_prob_PD_mean)[1] <- "out_prob_ALL"
colnames(out_prob_FTD_mean)[1] <- "out_prob_ALL"
ALL_Probablities <- rbind(out_prob_ALL_mean, out_prob_AD_mean, out_prob_DLB_mean, out_prob_PD_mean, out_prob_FTD_mean)
dim(ALL_Probablities) # 7584    3
names(ALL_Probablities)[3] <- "Final_Status"
length(unique(ALL_Probablities$UniquePhenoID)) # 2703
ALL_Probablities <- ALL_Probablities[order(ALL_Probablities$UniquePhenoID, ALL_Probablities$out_prob_ALL),]
ALL_Probablities <- ALL_Probablities[!duplicated(ALL_Probablities$UniquePhenoID),]
ALL_Probablities <- na.omit(ALL_Probablities)
ALL_Probablities$out_prob_ALL <- as.numeric(as.character(ALL_Probablities$out_prob_ALL))
dim(ALL_Probablities) # 2703    3
table(ALL_Probablities$Final_Status)
#      AD CO_ATneg      DLB      FTD       PD
#    1156      726       37       46      738
save(ALL_Probablities, file="CSF_ALL_Probability_N2703_AgeSexPlatePCs_Iterative_ClassBalanced.RData")
png("Probability_Violinplot_for_CSF_ALL_Signature_AgeSexPlatePCs_Iterative_ClassBalanced.png", units="mm", width=100, height=100, res=1000)
ggplot(ALL_Probablities, aes(x=Final_Status, y=out_prob_ALL, color=Final_Status)) + geom_violin(trim = TRUE) + 
  theme(legend.position = "none", panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x ="Clinical Diagnosis", y = "Predicted ALL Probability") + ggtitle("Probability Distribution for ALL signature")
dev.off()

# Calculate mean and standard deviation for each column
calculate_mean_sd <- function(df) {
  df <- df[,-1]
  # Ensure only numeric columns are included
  numeric_cols <- df[sapply(df, is.numeric)]
  # Calculate mean and SD for each numeric column
  formatted_stats <- sapply(numeric_cols, function(x) {
    mean_val <- mean(x, na.rm = TRUE)
    sd_val <- sd(x, na.rm = TRUE)
    paste0(round(mean_val, 2), "(", round(sd_val, 2), ")")
  })
  # Convert to a single-row dataframe
  result <- as.data.frame(t(formatted_stats))
  colnames(result) <- colnames(numeric_cols)
  return(result)
}
calculate_mean_sd(perf_df_ftd)
calculate_mean_sd(perf_df_ftd_base)
calculate_mean_sd(perf_df_PD)
calculate_mean_sd(perf_df_DLB)
calculate_mean_sd(perf_df_AD)
calculate_mean_sd(perf_df_FTD)
#         auc threshold   accuracy        npv        ppv sensitivity specificity
#1 0.78(0.02) 0.5(0.06) 0.72(0.02) 0.72(0.03) 0.72(0.04)  0.72(0.07)  0.72(0.07)
#         auc threshold   accuracy        npv        ppv sensitivity specificity
#1 0.57(0.03) 0.5(0.03) 0.57(0.02) 0.58(0.04) 0.58(0.05)   0.6(0.19)  0.53(0.19)
#         auc  threshold   accuracy        npv        ppv sensitivity specificity
#1 0.75(0.01) 0.45(0.03) 0.69(0.01) 0.71(0.02) 0.68(0.01)  0.73(0.03)  0.65(0.03)
#         auc  threshold   accuracy     npv        ppv sensitivity specificity
#1 0.81(0.01) 0.54(0.03) 0.76(0.04) 0.99(0) 0.14(0.02)  0.78(0.04)  0.76(0.04)
#        auc  threshold   accuracy        npv        ppv sensitivity specificity
#1 0.8(0.01) 0.56(0.02) 0.72(0.01) 0.61(0.01) 0.84(0.01)  0.68(0.03)  0.79(0.03)
#         auc threshold   accuracy        npv        ppv sensitivity specificity
#1 0.79(0.01) 0.5(0.08) 0.72(0.08) 0.98(0.01) 0.15(0.04)  0.75(0.09)  0.71(0.09)



## Heatmap of Lasso Proteins (from each disease) in CSF and Plasma

rm(list = ls())
library(dplyr) # dplyr_1.1.4
options(stringsAsFactors = FALSE)
options(width = 160)
set.seed(1)

load("CSF_AD_DLB_FTD_PD_DAA_Sumstats_Tables.RData")

AD_lasso <- read.csv("ProteinWeights_Sorted_50runs_ADvsCO_AgeSexPlatePCs.csv")
AD_lasso[1:10,c(1:4,50:52)]
AD_lasso_analytes <- AD_lasso$Proteins[1:10]
AD_lasso_proteins <- AD[AD$Analyte %in% AD_lasso_analytes,][,c("Analyte", "Symbol")]
AD_lasso_proteins$Disease <- "AD"

DLB_lasso <- read.csv("ProteinWeights_Sorted_50runs_DLBvsCO_AgeSexProjectPCs.csv")
DLB_lasso[1:10,c(1:4,50:52)]
DLB_lasso_analytes <- DLB_lasso$Proteins[1:10]
DLB_lasso_proteins <- DLB[DLB$Analyte %in% DLB_lasso_analytes,][,c("Analyte", "Symbol")]
DLB_lasso_proteins$Disease <- "DLB"

FTD_lasso <- read.csv("ProteinWeights_Sorted_50runs_FTDvsCO_AgeSexProjectPCs.csv")
FTD_lasso[1:10,c(1:4,50:52)]
FTD_lasso_analytes <- FTD_lasso$Proteins[1:10]
FTD_lasso_proteins <- FTD[FTD$Analyte %in% FTD_lasso_analytes,][,c("Analyte", "Symbol")]
FTD_lasso_proteins$Disease <- "FTD"

PD_lasso <- read.csv("ProteinWeights_Sorted_50runs_PDvsCO_AgeSexPlatePCs.csv")
PD_lasso[11:20,c(1:4,50:52)]
PD_lasso_analytes <- PD_lasso$Proteins[11:20]
PD_lasso_proteins <- PD[PD$Analyte %in% PD_lasso_analytes,][,c("Analyte", "Symbol")]
PD_lasso_proteins$Disease <- "PD"

ALL_lasso <- read.csv("ProteinWeights_Sorted_50runs_ALLvsCO_AgeSexPlatePCs.csv")
ALL_lasso[1:10,c(1:4,50:52)]
ALL_lasso_analytes <- ALL_lasso$Proteins[1:10]
ALL_lasso_proteins <- PD[PD$Analyte %in% ALL_lasso_analytes,][,c("Analyte", "Symbol")]
ALL_lasso_proteins$Disease <- "ALL"

uniq_lasso_analytes <- unique(c(AD_lasso_analytes, DLB_lasso_analytes, FTD_lasso_analytes, PD_lasso_analytes, ALL_lasso_analytes))
length(uniq_lasso_analytes) # 44

all_lasso_proteins <- rbind(AD_lasso_proteins, DLB_lasso_proteins, FTD_lasso_proteins, PD_lasso_proteins, ALL_lasso_proteins)
dim(all_lasso_proteins) # 50  3

# Check duplicated proteins
ind <- duplicated(all_lasso_proteins[,1:2])
all_lasso_proteins[ind,]
all_lasso_proteins <- all_lasso_proteins[!duplicated(all_lasso_proteins$Analyte),]
dim(all_lasso_proteins) # 44  3

all_lasso_proteins_AD <- inner_join(all_lasso_proteins, AD[,c(1,3,5,7)], by=c("Analyte", "Symbol"))
all_lasso_proteins_AD_DLB <- inner_join(all_lasso_proteins_AD, DLB[,c(1,3,5,7)], by=c("Analyte", "Symbol"))
all_lasso_proteins_AD_DLB_FTD <- inner_join(all_lasso_proteins_AD_DLB, FTD[,c(1,3,5,7)], by=c("Analyte", "Symbol"))
all_lasso_proteins_AD_DLB_FTD_PD <- inner_join(all_lasso_proteins_AD_DLB_FTD, PD[,c(1,3,5,7)], by=c("Analyte", "Symbol"))
names(all_lasso_proteins_AD_DLB_FTD_PD) <- c("Analyte", "Symbol", "Disease", "B_AD", "P_AD", "B_DLB", "P_DLB", "B_FTD", "P_FTD", "B_PD", "P_PD")
# sanity check
PD[PD$Analyte %in% PD_lasso_analytes,]
AD[AD$Analyte %in% PD_lasso_analytes,]
write.csv(all_lasso_proteins_AD_DLB_FTD_PD, file="Suppl_Table_CSF_Lasso_Proteins_Sumstats.csv", row.names=F, quote=F)

# Plasma
lasso_plasma <- read.csv("Suppl_Table_Plasma_Lasso_Proteins_Sumstats.csv")
lasso_plasma$Tissue <- "Plasma"
dim(lasso_plasma) # 41 13
load("Plasma_AD_DLB_FTD_PD_DAA_Sumstats_v3.RData")
AD_Plasma <- AD
DLB_Plasma <- DLB
FTD_Plasma <- FTD
PD_Plasma <- PD

common_analytes_plasma <- Reduce(intersect, list(AD_Plasma$Analyte, DLB_Plasma$Analyte, FTD_Plasma$Analyte, PD_Plasma$Analyte))
length(common_analytes_plasma) # 6594
table(lasso_plasma$Analyte %in% common_analytes_plasma)
# TRUE
#   41

# CSF
lasso_csf <- read.csv("Suppl_Table_CSF_Lasso_Proteins_Sumstats.csv")
lasso_csf$External_ID <- lasso_csf$Analyte
lasso_csf$Analyte <- sapply(strsplit(as.character(lasso_csf$External_ID), "_"), "[", 1)
lasso_csf <- lasso_csf[,c(1,2,12,3:11)]
lasso_csf$Tissue <- "CSF"
dim(lasso_csf) # 44 13

load("CSF_AD_DLB_FTD_PD_DAA_Sumstats_Tables.RData")
AD$External_ID <- AD$Analyte
AD$Analyte <- sapply(strsplit(as.character(AD$External_ID), "_"), "[", 1)
DLB$External_ID <- DLB$Analyte
DLB$Analyte <- sapply(strsplit(as.character(DLB$External_ID), "_"), "[", 1)
FTD$External_ID <- FTD$Analyte
FTD$Analyte <- sapply(strsplit(as.character(FTD$External_ID), "_"), "[", 1)
PD$External_ID <- PD$Analyte
PD$Analyte <- sapply(strsplit(as.character(PD$External_ID), "_"), "[", 1)

common_analytes_csf <- Reduce(intersect, list(AD$Analyte, DLB$Analyte, FTD$Analyte, PD$Analyte))
length(common_analytes_csf) # 6956
table(lasso_csf$Analyte %in% common_analytes_csf)
# TRUE
#   44
table(lasso_plasma$Analyte %in% common_analytes_csf)
# TRUE
#   41
table(lasso_csf$Analyte %in% common_analytes_plasma)
# FALSE  TRUE
#     4    40
# 4 analytes in CSF are not present in the plasma data
lasso_csf[!(lasso_csf$Analyte %in% common_analytes_plasma),]
#      Analyte Symbol      External_ID Disease       B_AD          P_AD      B_DLB        P_DLB       B_FTD        P_FTD        B_PD         P_PD
# 1  X14157.21  YWHAE X14157.21_P62258      AD  0.9801903 1.419266e-187  0.5757838 1.097736e-07  0.48088241 1.167462e-08  0.03705335 0.6806518331
# 7   X6521.35  NPTX2  X6521.35_P47972      AD -0.5192611  1.171276e-47 -1.2670447 4.208545e-23 -1.08599542 1.455752e-27 -0.21785169 0.0040043317
# 31  X6586.19 ADAM11  X6586.19_O75078      PD -0.1142535  3.944760e-05  0.0825013 4.690472e-01  0.04798633 5.672121e-01 -0.25654542 0.0001142414
# 37   X9025.5   AHSP   X9025.5_Q9NZD4      PD -0.2032967  2.189177e-07 -0.7299998 1.872562e-07 -0.63086415 1.628997e-09  0.26120221 0.0025314791

common_analytes_csf_plasma <- Reduce(intersect, list(common_analytes_csf, common_analytes_plasma))
length(common_analytes_csf_plasma) # 6418

lasso_csf_plasma <- rbind(lasso_csf, lasso_plasma)
dim(lasso_csf_plasma) # 85 13

# filter for only those analytes present in both tissues
lasso_csf_plasma <- lasso_csf_plasma[lasso_csf_plasma$Analyte %in% common_analytes_csf_plasma,]
dim(lasso_csf_plasma) # 81 13
length(unique(lasso_csf_plasma$Analyte)) # 80
length(unique(lasso_csf_plasma$Symbol)) # 79
ind <- duplicated(lasso_csf_plasma[,1:2])
lasso_csf_plasma[ind,]
# SMOC1 (X13118.5) = AD (in both CSF and Plasma)
uniq_lasso_analytes <- unique(lasso_csf_plasma$Analyte)
length(uniq_lasso_analytes) # 80


all_lasso_proteins <- lasso_csf_plasma[,c(1:4,13)]
dim(all_lasso_proteins) # 81  5
all_lasso_proteins <- all_lasso_proteins[!duplicated(all_lasso_proteins$Analyte),]
dim(all_lasso_proteins) # 80  5

sumstats_csf_AD <- inner_join(all_lasso_proteins, AD[,c(1,3,5,7)], by=c("Analyte", "Symbol")) # 80  7
sumstats_csf_AD_DLB <- inner_join(sumstats_csf_AD, DLB[,c(1,3,5,7)], by=c("Analyte", "Symbol")) # 80  9
sumstats_csf_AD_DLB_FTD <- inner_join(sumstats_csf_AD_DLB, FTD[,c(1,3,5,7)], by=c("Analyte", "Symbol")) # 80 11
sumstats_csf <- inner_join(sumstats_csf_AD_DLB_FTD, PD[,c(1,3,5,7)], by=c("Analyte", "Symbol")) # 80 13
names(sumstats_csf) <- c("Analyte", "Symbol", "External_ID", "Disease", "Tissue", "B_AD", "P_AD", "B_DLB", "P_DLB", "B_FTD", "P_FTD", "B_PD", "P_PD")


sumstats_plsm_AD <- inner_join(all_lasso_proteins, AD_Plasma[,c(1,2,6,8)], by=c("Analyte", "Symbol")) # 80  6
sumstats_plsm_AD_DLB <- inner_join(sumstats_plsm_AD, DLB_Plasma[,c(1,2,6,8)], by=c("Analyte", "Symbol")) # 80  8
sumstats_plsm_AD_DLB_FTD <- inner_join(sumstats_plsm_AD_DLB, FTD_Plasma[,c(1,2,6,8)], by=c("Analyte", "Symbol")) # 80 10
sumstats_plsm <- inner_join(sumstats_plsm_AD_DLB_FTD, PD_Plasma[,c(1,2,6,8)], by=c("Analyte", "Symbol")) # 80 12
names(sumstats_plsm) <- c("Analyte", "Symbol", "External_ID", "Disease", "Tissue", "B_AD", "P_AD", "B_DLB", "P_DLB", "B_FTD", "P_FTD", "B_PD", "P_PD")

final_sumstats_csf_plasma <- inner_join(sumstats_csf, sumstats_plsm, by=c("Analyte", "Symbol", "External_ID", "Disease"))
dim(final_sumstats_csf_plasma) # 80 22

sumstats_beta <- final_sumstats_csf_plasma[,c("Analyte", "Symbol", "Disease", "Tissue.x", "B_AD.x", "B_DLB.x", "B_FTD.x", "B_PD.x", "B_AD.y", "B_DLB.y", "B_FTD.y", "B_PD.y")]
names(sumstats_beta) <- c("Analyte", "Symbol", "Disease", "Tissue", "AD", "DLB", "FTD", "PD", "AD", "DLB", "FTD", "PD")
sumstats_beta$Disease <- factor(sumstats_beta$Disease, levels=c("AD", "DLB", "FTD", "PD", "ALL"))
sumstats_beta$Symbol <- make.unique(sumstats_beta$Symbol)
sumstats_beta <- sumstats_beta[order(sumstats_beta$Disease),]
# names(sumstats_beta) <- c("Analyte", "Symbol", "Disease", "Tissue", "AD_CSF", "DLB_CSF", "FTD_CSF", "PD_CSF", "AD_Plasma", "DLB_Plasma", "FTD_Plasma", "PD_Plasma")

library(ComplexHeatmap)
library(circlize)

heatmap_data <- sumstats_beta[, 5:ncol(sumstats_beta)]
rownames(heatmap_data) <- sumstats_beta$Symbol
heatmap_matrix <- as.matrix(heatmap_data)
# colnames(heatmap_matrix) <- c("AD_CSF", "DLB_CSF", "FTD_CSF", "PD_CSF", "AD_Plasma", "DLB_Plasma", "FTD_Plasma", "PD_Plasma")
colnames(heatmap_matrix) <- c("AD", "DLB", "FTD", "PD", "AD", "DLB", "FTD", "PD")

# Create row annotations for Disease and Tissue
# Color-blind friendly pallet (https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible)
row_annotation <- rowAnnotation(
  Disease = sumstats_beta$Disease,
  Tissue = sumstats_beta$Tissue,
  col = list(
    Disease = c("AD" = "#0072B2", "DLB" = "#E69F00", "FTD" = "#009E73", "PD" = "#DC3220", "ALL" = "#994F00"),
    Tissue = c("CSF" = "#56B4E9", "Plasma" = "#882255")
  )
)
png("CSF_Plasma_Lasso_Proteins_EffectSize_Heatmap.png", units="mm", width=180, height=300, res=1000)
Heatmap(
  heatmap_matrix,
  name = "Effect size",
  row_title = "LASSO Proteins",
  column_title = "Heatmap of LASSO protein effect size",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_split = sumstats_beta$Disease, # Split rows by Disease
  left_annotation = row_annotation,  # Add row annotation
  col = colorRamp2(c(-1, 0, 1), c("#6699CC", "white", "#D55E00")) # Color gradient for expression
)
dev.off()
save(heatmap_matrix, row_annotation, sumstats_beta, file="CSF_Plasma_Lasso_Proteins_EffectSize_Heatmap.RData", row.names=F, quote=F)