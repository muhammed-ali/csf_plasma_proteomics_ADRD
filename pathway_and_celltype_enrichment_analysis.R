suppressPackageStartupMessages({
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(tidyr)
	library(DOSE)
	library(ggnewscale)
	library(enrichplot)
	library(ReactomePA)
	library(multienrichjam) # for converting df to enrichResult object
	library(dplyr)
})
rm(list=ls())
set.seed(1)
options(width=130)


## REACTOME Pathways (FDR) for Unique DEPs ##

load("CSF_Nominal_FDR_AD_DLB_FTD_PD_DEPs_Sumstats.RData")

df_list <- list(AD = AD_FDR, PD = PD_FDR, FTD = FTD_FDR, DLB = DLB_FDR)
unique_DEPs_FDR <- lapply(1:4, function(i) setdiff(df_list[[i]]$Analyte, do.call(rbind, df_list[-i])$Analyte))
lengths(unique_DEPs_FDR) # 1376   60  551  132

AD_DAA <- AD_FDR[AD_FDR$Analyte %in% unique_DEPs_FDR[[1]],]
AD_DAA <- AD_DAA[order(AD_DAA$UniProt, AD_DAA$P),]
dim(AD_DAA) # 1376    9
AD_DAA_FDR <- AD_DAA[!duplicated(AD_DAA$UniProt),]
dim(AD_DAA_FDR) # 1317    9

DLB_DAA <- DLB_FDR[DLB_FDR$Analyte %in% unique_DEPs_FDR[[4]],]
DLB_DAA <- DLB_DAA[order(DLB_DAA$UniProt, DLB_DAA$P),]
dim(DLB_DAA) # 132   9
DLB_DAA_FDR <- DLB_DAA[!duplicated(DLB_DAA$UniProt),]
dim(DLB_DAA_FDR) # 131   9

FTD_DAA <- FTD_FDR[FTD_FDR$Analyte %in% unique_DEPs_FDR[[3]],]
FTD_DAA <- FTD_DAA[order(FTD_DAA$UniProt, FTD_DAA$P),]
dim(FTD_DAA) # 551   9
FTD_DAA_FDR <- FTD_DAA[!duplicated(FTD_DAA$UniProt),]
dim(FTD_DAA_FDR) # 543   9

PD_DAA <- PD_FDR[PD_FDR$Analyte %in% unique_DEPs_FDR[[2]],]
PD_DAA <- PD_DAA[order(PD_DAA$UniProt, PD_DAA$P),]
dim(PD_DAA) # 60  9
PD_DAA_FDR <- PD_DAA[!duplicated(PD_DAA$UniProt),]
dim(PD_DAA_FDR) # 59  9


## AD Reactome pathways

# All background proteins
background <- read.table("CSF_SOMAscan7k_analyte_information.tsv", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(background) # 7291   14
background <- as.data.frame(separate_rows(background, EntrezGeneID, UniProt, EntrezGeneSymbol, sep="\\|", convert = TRUE))
dim(background) # 7379   14
length(unique(background$EntrezGeneID)) # 6376
AD_FDR_entrez <- unique(background[background$UniProt %in% AD_DAA_FDR$UniProt,]$EntrezGeneID)
length(AD_FDR_entrez) # 1302
ADCO_Reactome <- enrichPathway(AD_FDR_entrez,
             organism = 'human',
             minGSSize = 3, 
             maxGSSize = 3000, 
             pvalueCutoff = 0.05,
             qvalueCutoff = 1,
             readable = TRUE,
             universe = unique(background$EntrezGeneID),
             pAdjustMethod = "fdr") 
dim(ADCO_Reactome) # 111   9
nrow(ADCO_Reactome@result[ADCO_Reactome@result$p.adjust < 0.05,]) # 111
nrow(ADCO_Reactome@result[ADCO_Reactome@result$pvalue < 0.05,]) # 367
nominal_pathways <- nrow(ADCO_Reactome@result[ADCO_Reactome@result$pvalue < 0.05,])
fdr_pathways <- nrow(ADCO_Reactome@result[ADCO_Reactome@result$p.adjust < 0.05,])
png("AD_vs_CO_enrichReactome_dotPlot.png", units="mm", width=150, height=180, res=1000)
dotplot(ADCO_Reactome, showCategory=10, font.size=16, title = paste0(" AD_vs_CO"))
dev.off()
write.csv(ADCO_Reactome, file="AD_vs_CO_enrichReactome_FDR.csv", row.names=F, quote=F)
save(ADCO_Reactome, AD_DAA_FDR, AD_FDR_entrez, file="AD_vs_CO_enrichReactome_FDR.RData")


## DLB Reactome pathways

DLB_FDR_entrez <- unique(background[background$UniProt %in% DLB_DAA_FDR$UniProt,]$EntrezGeneID)
length(DLB_FDR_entrez) # 127
DLBCO_Reactome <- enrichPathway(DLB_FDR_entrez,
             organism = 'human',
             minGSSize = 3, 
             maxGSSize = 3000, 
             pvalueCutoff = 0.05,
             qvalueCutoff = 1,
             readable = TRUE,
             universe = unique(background$EntrezGeneID),
             pAdjustMethod = "fdr") 
dim(DLBCO_Reactome) # 17  9
nrow(DLBCO_Reactome@result[DLBCO_Reactome@result$p.adjust < 0.05,]) # 17
nrow(DLBCO_Reactome@result[DLBCO_Reactome@result$pvalue < 0.05,]) # 108
nominal_pathways <- nrow(DLBCO_Reactome@result[DLBCO_Reactome@result$pvalue < 0.05,])
fdr_pathways <- nrow(DLBCO_Reactome@result[DLBCO_Reactome@result$p.adjust < 0.05,])
# Dotplot
write.csv(DLBCO_Reactome, file="DLB_vs_CO_enrichReactome_FDR.csv", row.names=F, quote=F)
save(DLBCO_Reactome, DLB_DAA_FDR, DLB_FDR_entrez, file="DLB_vs_CO_enrichReactome_FDR.RData")
DLBCO_Reactome_modified <- DLBCO_Reactome
DLBCO_Reactome_modified@result$Description[7] <- "Diseases of signal transduction .."
png("DLB_vs_CO_enrichReactome_dotPlot.png", units="mm", width=150, height=180, res=1000)
dotplot(DLBCO_Reactome_modified, showCategory=10, font.size=16, title = paste0(" DLB_vs_CO"))
dev.off()


## FTD Reactome pathways

FTD_FDR_entrez <- unique(background[background$UniProt %in% FTD_DAA_FDR$UniProt,]$EntrezGeneID)
length(FTD_FDR_entrez) # 536
FTDCO_Reactome <- enrichPathway(FTD_FDR_entrez,
             organism = 'human',
             minGSSize = 3, 
             maxGSSize = 3000, 
             pvalueCutoff = 0.05,
             qvalueCutoff = 1,
             readable = TRUE,
             universe = unique(background$EntrezGeneID),
             pAdjustMethod = "fdr") 
dim(FTDCO_Reactome) # 111   9
nrow(FTDCO_Reactome@result[FTDCO_Reactome@result$p.adjust < 0.05,]) # 111
nrow(FTDCO_Reactome@result[FTDCO_Reactome@result$pvalue < 0.05,]) # 367
nominal_pathways <- nrow(FTDCO_Reactome@result[FTDCO_Reactome@result$pvalue < 0.05,])
fdr_pathways <- nrow(FTDCO_Reactome@result[FTDCO_Reactome@result$p.adjust < 0.05,])
png("FTD_vs_CO_enrichReactome_dotPlot.png", units="mm", width=150, height=180, res=1000)
dotplot(FTDCO_Reactome, showCategory=10, font.size=16, title = paste0(" FTD_vs_CO"))
dev.off()
write.csv(FTDCO_Reactome, file="FTD_vs_CO_enrichReactome_FDR.csv", row.names=F, quote=F)
save(FTDCO_Reactome, FTD_DAA_FDR, FTD_FDR_entrez, file="FTD_vs_CO_enrichReactome_FDR.RData")


## PD Reactome pathways

PD_FDR_entrez <- unique(background[background$UniProt %in% PD_DAA_FDR$UniProt,]$EntrezGeneID)
length(PD_FDR_entrez) # 57
PDCO_Reactome <- enrichPathway(PD_FDR_entrez,
             organism = 'human',
             minGSSize = 3, 
             maxGSSize = 3000, 
             pvalueCutoff = 0.05,
             qvalueCutoff = 1,
             readable = TRUE,
             universe = unique(background$EntrezGeneID),
             pAdjustMethod = "none") # NOT FDR
dim(PDCO_Reactome) # 39  9
nrow(PDCO_Reactome@result[PDCO_Reactome@result$p.adjust < 0.05,]) # 39
nrow(PDCO_Reactome@result[PDCO_Reactome@result$pvalue < 0.05,]) # 39
nominal_pathways <- nrow(PDCO_Reactome@result[PDCO_Reactome@result$pvalue < 0.05,])
fdr_pathways <- nrow(PDCO_Reactome@result[PDCO_Reactome@result$p.adjust < 0.05,])
png("PD_vs_CO_enrichReactome_dotPlot.png", units="mm", width=160, height=180, res=1000)
dotplot(PDCO_Reactome, showCategory=10, font.size=16, title = paste0(" PD_vs_CO"))
dev.off()
write.csv(PDCO_Reactome, file="PD_vs_CO_enrichReactome_FDR.csv", row.names=F, quote=F)
save(PDCO_Reactome, PD_DAA_FDR, PD_FDR_entrez, file="PD_vs_CO_enrichReactome_FDR.RData")



## REACTOME Pathways (FDR) for Common DEPs ##

## AD_DLB_FTD

background <- read.table("/03-DryLab/02-Data/04-Proteomics/03-AnalysisReady/202309_Soma_CSF_Complete-Matrix_/SomaScan7k_projects/CSF_SOMAscan7k_analyte_information.tsv", sep="\t", header=T, stringsAsFactors=F, quote="")
dim(background) # 7291   14
background <- as.data.frame(separate_rows(background, EntrezGeneID, UniProt, EntrezGeneSymbol, sep="\\|", convert = TRUE))
dim(background) # 7379   14
length(unique(background$EntrezGeneID)) # 6376
load("CSF_Nominal_FDR_AD_DLB_FTD_PD_DEPs_Sumstats.RData")
AD_FDR <- as.data.frame(separate_rows(AD_FDR, UniProt, Symbol, sep="\\|", convert = TRUE))
dim(AD_FDR) # 4476    9
DLB_FDR <- as.data.frame(separate_rows(DLB_FDR, UniProt, Symbol, sep="\\|", convert = TRUE))
dim(DLB_FDR) # 2894    9
FTD_FDR <- as.data.frame(separate_rows(FTD_FDR, UniProt, Symbol, sep="\\|", convert = TRUE))
dim(FTD_FDR) # 3912    9
PD_FDR <- as.data.frame(separate_rows(PD_FDR, UniProt, Symbol, sep="\\|", convert = TRUE))
dim(PD_FDR) # 379   9
AD_DLB <- AD_FDR[AD_FDR$UniProt %in% DLB_FDR$UniProt,]
dim(AD_DLB) # 2427    9
AD_DLB_FTD <- AD_DLB[AD_DLB$UniProt %in% FTD_FDR$UniProt,]
dim(AD_DLB_FTD) # 2184    9
length(unique(AD_DLB_FTD$UniProt)) # 1903
AD_DLB_FTD_entrez <- unique(background[background$UniProt %in% AD_DLB_FTD$UniProt,]$EntrezGeneID)
length(AD_DLB_FTD_entrez) # 1900
AD_DLB_FTD_Reactome <- enrichPathway(AD_DLB_FTD_entrez,
             organism = 'human',
             minGSSize = 3, 
             maxGSSize = 3000, 
             pvalueCutoff = 0.05,
             qvalueCutoff = 1,
             readable = TRUE,
             universe = unique(background$EntrezGeneID),
             pAdjustMethod = "fdr") 
dim(AD_DLB_FTD_Reactome) # 2 9
nrow(AD_DLB_FTD_Reactome@result[AD_DLB_FTD_Reactome@result$p.adjust < 0.05,]) # 2
nrow(AD_DLB_FTD_Reactome@result[AD_DLB_FTD_Reactome@result$pvalue < 0.05,]) # 162
nominal_pathways <- nrow(AD_DLB_FTD_Reactome@result[AD_DLB_FTD_Reactome@result$pvalue < 0.05,])
fdr_pathways <- nrow(AD_DLB_FTD_Reactome@result[AD_DLB_FTD_Reactome@result$p.adjust < 0.05,])

# Dotplot
png("AD_DLB_FTD_enrichReactome_dotPlot.png", units="mm", width=250, height=80, res=1000)
dotplot(AD_DLB_FTD_Reactome, showCategory=15, font.size=16, title = paste0(" AD_DLB_FTD enrichReactome \n Proteins = ", length(AD_DLB_FTD_entrez), ", # Pathways = ", fdr_pathways))
dev.off()
write.csv(AD_DLB_FTD_Reactome, file="AD_DLB_FTD_enrichReactome_FDR.csv", row.names=F, quote=F)
write.csv(AD_DLB_FTD_Reactome@result[AD_DLB_FTD_Reactome@result$pvalue < 0.05,], file="AD_DLB_FTD_enrichReactome_Nominal.csv", row.names=F, quote=F)
save(AD_DLB_FTD_Reactome, AD_DLB_FTD, AD_DLB_FTD_entrez, file="AD_DLB_FTD_enrichReactome_FDR.RData")

# Check Unique Pathways
AD_DLB_FTD_pathways <- AD_DLB_FTD_Reactome@result[AD_DLB_FTD_Reactome@result$p.adjust < 0.05,]
dim(AD_DLB_FTD_pathways) # 2 9
pathway_list <- list(AD_DLB_FTD = AD_DLB_FTD_pathways, AD = AD_pathways, DLB = DLB_pathways, FTD = FTD_pathways, PD = PD_pathways)
unique_pathways <- lapply(1:5, function(i) setdiff(pathway_list[[i]]$ID, do.call(rbind, pathway_list[-i])$ID))
lengths(unique_pathways) # 2 91  7  4 37
# Both of the identified pathways are unique to this intersection (AD_DLB_FTD) and do not overlap with disease-specific pathways above.
AD_DLB_FTD_Reactome@result[1:2,c("Description", "geneID")]

## AD_DLB_FTD_PD
Common4 <- AD_DLB_FTD[AD_DLB_FTD$UniProt %in% PD_FDR$UniProt,]
dim(Common4) # 185   9
length(unique(Common4$UniProt)) # 146
Common4_entrez <- unique(background[background$UniProt %in% Common4$UniProt,]$EntrezGeneID)
length(Common4_entrez) # 146
Common4_Reactome <- enrichPathway(Common4_entrez,
             organism = 'human',
             minGSSize = 3, 
             maxGSSize = 3000, 
             pvalueCutoff = 0.05,
             qvalueCutoff = 1,
             readable = TRUE,
             universe = unique(background$EntrezGeneID),
             pAdjustMethod = "none") # fdr return 0 pathways
dim(Common4_Reactome) # 42  9
nrow(Common4_Reactome@result[Common4_Reactome@result$p.adjust < 0.05,]) # 42
nrow(Common4_Reactome@result[Common4_Reactome@result$pvalue < 0.05,]) # 42
nominal_pathways <- nrow(Common4_Reactome@result[Common4_Reactome@result$pvalue < 0.05,])
fdr_pathways <- nrow(Common4_Reactome@result[Common4_Reactome@result$p.adjust < 0.05,])
# Dotplot
png("Common4_enrichReactome_dotPlot.png", units="mm", width=250, height=220, res=1000)
dotplot(Common4_Reactome, showCategory=15, font.size=16, title = paste0(" Common4 enrichReactome \n Proteins = ", length(Common4_entrez), ", # Pathways = ", fdr_pathways))
dev.off()
write.csv(Common4_Reactome, file="Common4_enrichReactome_Nominal.csv", row.names=F, quote=F)
save(Common4_Reactome, Common4, Common4_entrez, file="Common4_enrichReactome_Nominal.RData")
# Check Unique Pathways
Common4_pathways <- Common4_Reactome@result[Common4_Reactome@result$p.adjust < 0.05,]
dim(Common4_pathways) # 42  9
pathway_list <- list(Common4 = Common4_pathways, AD_DLB_FTD = AD_DLB_FTD_pathways, AD = AD_pathways, DLB = DLB_pathways, FTD = FTD_pathways, PD = PD_pathways)
unique_pathways <- lapply(1:6, function(i) setdiff(pathway_list[[i]]$ID, do.call(rbind, pathway_list[-i])$ID))
lengths(unique_pathways) # 35  0 89  6  4 35
Common4_pathways_unique <- Common4_Reactome[Common4_Reactome$ID %in% unique_pathways[[1]], asis=T]
dim(Common4_pathways_unique) # 35  9
png("Common4_unique_enrichReactome_treePlot.png", units="mm", width=350, height=220, res=1000)
edox2 <- pairwise_termsim(Common4_pathways_unique)
treeplot(edox2, showCategory = 35) # nWords=4 # hclust_method = "centroid"
dev.off()
save(Common4_pathways, AD_DLB_FTD_pathways, AD_pathways, DLB_pathways, FTD_pathways, PD_pathways, file="CSF_enrichReactome_Pathways_AD_DLB_FTD_PD.RData")



## Get cell type enrichment for proteins in the pathways ##

suppressPackageStartupMessages({
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(tidyr)
	library(DOSE)
	library(ggnewscale)
	library(enrichplot)
	library(ReactomePA)
	library(multienrichjam) # for converting df to enrichResult object
	library(dplyr)
	library(data.table)
	library(readxl)
	library(gtools)
})
rm(list=ls())
set.seed(1)
options(width=130)

load("CSF_Nominal_FDR_AD_DLB_FTD_PD_DEPs_Sumstats.RData")

load("Venn_CSF_AD_PD_DLB_FTD_Common_Pathways_from_unique_DEPs.RData")

ct_composition <- data.table(read_excel("brain_cell_type_specificity.xlsx", sheet = "cell_types_with_percentiles"))
dim(ct_composition) # 23225    52
ct_composition <- ct_composition[,c(1,52)]
ct_composition <- ct_composition[3:nrow(ct_composition),] # drop first 2 rows
colnames(ct_composition) <- c("Gene", "Cell_Type")
ct_composition[is.na(ct_composition)] <- "Unknown" # Replace missing gene annotation with "Unknown"
dim(ct_composition) # 23223     2
length(unique(ct_composition$Gene)) # 21661
ct_composition <- ct_composition[order(ct_composition$Cell_Type),]
ct_composition <- ct_composition[!duplicated(ct_composition$Gene),]
dim(ct_composition) # 21661     2
head(ct_composition, 2)
#    Gene         Cell_Type
# 1:  ACE Human Endothelial
# 2: HBA1 Human Endothelial

ADCO_significant <- left_join(AD_FDR, ct_composition, by=c("Symbol"="Gene"))
dim(ADCO_significant) # 4425   10
ADCO_significant$Cell_Type[is.na(ADCO_significant$Cell_Type)] <- "Unknown"
table(ADCO_significant$Cell_Type)
#         Human Endothelial    Human mature astrocytes Human Microglia/Macrophage              Human Neurons
#                       166                        231                        358                        441
#    Human Oligodendrocytes                    Unknown
#                        95                       3134

background <- length(unique(ADCO_significant$Symbol)) # 4047

CT_df <- as.data.frame(table(ADCO_significant$Cell_Type))
colnames(CT_df) <- c("cell_type", "CT_labels")
CT_df
#                   cell_type CT_labels
#1          Human Endothelial       166
#2    Human mature astrocytes       231
#3 Human Microglia/Macrophage       358
#4              Human Neurons       441
#5     Human Oligodendrocytes        95
#6                    Unknown      3134

pathway_stats <- function(enrichment_object, pathway_name, ADCO_significant, background){
 pathway_id <- enrichment_object[enrichment_object$Description == pathway_name,]$Description
 print(pathway_id)
 pathway_genes <- enrichment_object[enrichment_object$Description == pathway_name,]$geneID
 pathway_genes <- strsplit(pathway_genes, "/")
 pathway_genes <- pathway_genes[[1]]
 print("genes in the pathway:")
 print(length(unique(pathway_genes))) # 111

 print("CT of genes in pathway:") # Either use "ADCO_significant$Gene" or "ADCO_significant$Entrez" based on input
 pathway_CT_groups <- as.data.frame(table(ADCO_significant[ADCO_significant$Symbol %in% pathway_genes,]$Cell_Type)) # Entrez
 colnames(pathway_CT_groups) <- c("cell_type", "freq")
 pathway_CT_groups <- pathway_CT_groups[order(pathway_CT_groups$cell_type),]
 pathway_CT_groups <- inner_join(pathway_CT_groups, CT_df, by="cell_type")
 pathway_CT_groups$percent_global <- (pathway_CT_groups$CT_labels/background)*100
 pathway_CT_groups$percent_local <- (pathway_CT_groups$freq/length(unique(pathway_genes)))*100
 pathway_CT_groups$cellType_FC <- round(foldchange(pathway_CT_groups$percent_local, pathway_CT_groups$percent_global), 2)
 pathway_CT_groups <- pathway_CT_groups[pathway_CT_groups$cell_type != "Unknown",]
 pathway_CT_groups$enrichment_P <- phyper(pathway_CT_groups$freq, pathway_CT_groups$CT_labels, background-pathway_CT_groups$CT_labels, length(unique(pathway_genes)), lower.tail = FALSE)
 print(pathway_CT_groups)
 list_res <<- list(pathway_genes, pathway_CT_groups)
}

Metabolism <- pathway_stats(AD_pathways_common, "Metabolism", ADCO_significant, background)
# [1] "genes in the pathway:"
# [1] 211
# [1] "CT of genes in pathway:"
#                    cell_type freq CT_labels percent_global percent_local cellType_FC enrichment_P
# 1          Human Endothelial   10       166       4.101804      4.739336        1.16 2.463805e-01
# 2    Human mature astrocytes   33       231       5.707932     15.639810        2.74 1.523080e-08
# 3 Human Microglia/Macrophage   16       358       8.846059      7.582938       -1.17 6.973705e-01
# 4              Human Neurons   22       441      10.896961     10.426540       -1.05 5.339996e-01
# 5     Human Oligodendrocytes    6        95       2.347418      2.843602        1.21 2.249017e-01

