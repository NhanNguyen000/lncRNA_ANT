###  =============================================== Analysis ImpulseDE2 outcome:

# Load code & data --------------------------------------------------------
# R scripts function
get.path <- function(list_categories) {
  require("here")
  library(here)
  path <- here::here(list_categories)
  return(path)
}
source(get.path("R_functions_forAnthracyclines_ImpulseDE2.R"))
ls()
## Input RNA data & get Ensemble data
setwd(get.path("data"))
load("Con_DF2data_Cardiac_NN_20191112.RData")
#load("Con_UNTRdata_Cardiac_NN_2020Feb26.RData")
load("DOXdata_Cardiac_NN_20191112.RData")
load("EPIdata_Cardiac_NN_20191112.RData") 
load("IDAdata_Cardiac_NN_20191112.RData")
#load("Biopsies_ANT_Cardiac_NN_2020Feb19.RData")
load("Biopsies_ANT_Cardiac_NN_2020May25.RData")


Ensemble_database    <- Loaded_data(folder = get.path("data"), file_name = "Ensemble_mart_export_NN_20190815.txt")

# Figure 1 --------------------------------------------------------
# Preparation: get normalized data (read count) + expression in gene categorical in plot (Figure 1.A) ---------------------------
## Get normalized data (without Con_UNTR)
RNA_data_count  <- get.RNA_data_combine(list(Con_DF2, DOX, EPI, IDA), "expected_count")
RNA_data_count  <- RNA_data_count[, -grep("Con_DF2_000_", colnames(RNA_data_count))]
RNA_total_read_count <- get.total_reads_plot(RNA_data_count) 
RNA_total_read_count[which(RNA_total_read_count<5e+6)] # read count of DOX_Tox_240_2, IDA_The_240_3 < 5 milions (5e+6) read count

RNA_data_count  <- RNA_data_count[, -grep("DOX_Tox_240_2", colnames(RNA_data_count))]
RNA_data_count  <- RNA_data_count[, -grep("IDA_The_240_3", colnames(RNA_data_count))]
norm_data       <- get.norm_data(RNA_data_count) # by DEseq2
norm_data_v1    <- get.remove_gene_low_express(norm_data, 0)

# figure 1A: clustering the invitro RNseq profiles
library(tidyverse)
library(dendextend)
library("scales")
dend <- t(norm_data_v1) %>% scale %>% dist %>% hclust("ward.D") %>% as.dendrogram %>% 
  set("branches_lwd", 2) %>% set("labels_cex", 0.8)
groupCodes <- substr(colnames(norm_data_v1), 1, 7)
colorCodes <- hue_pal()(8)[2:8]
names(colorCodes) <- unique(groupCodes)
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
#pdf("Figure1A.pdf", height = 20, width = 8)
par(mar=c(5.1, 10, 4.1, 15))
plot(dend, horiz = TRUE)
dev.off()


# Figure 1 B: PCA of the human biopsies data:
# Clean biopsies metadata -------------------------------------------------
Biospies_info <- read.table(paste0(get.path("data"), 
                                   "/ClinicalParametersHeCaToS_Patients2018FINAL_summary_NN_2020Feb26.txt"), 
                            header = TRUE, sep = "\t", fill= TRUE, stringsAsFactors = FALSE)
Biospies_info <- Biospies_info[which(Biospies_info$SID != ""),]
Biospies_info$SID[Biospies_info$SID == "10198 (!)"] <- "10198"

Biospies_info2 <- read.csv(paste0(get.path("data"), 
                                  "/Copy of ClinicalParametersHeCaToS_SendPatients_MOGESV2_4Florian4R.csv"),
                           na.strings = "99999", stringsAsFactors = FALSE)
Biospies_info2 <- Biospies_info2[which(is.na(Biospies_info2$SID) == FALSE),]

Overlap_items <- c("Matched._to_SID", "Year_of_birth", "Sex", "Sex_numerical", 
                   "Lenght", "Weight", "BMI", "Cancer_type", "Genetic_mutation",
                   "Comorbidities" , "Death")
Biospies_info2 <- Biospies_info2[-which(names(Biospies_info2) %in% Overlap_items)]
Biospies_info_All <- merge(Biospies_info, Biospies_info2, by = "SID", all = TRUE)
SID_Biopsies_RNAseq <- matrix(unlist(strsplit(colnames(Biopsies$expected_count), "_")),
                              ncol = 2, byrow = TRUE)[,2]
Biospies_info_All<-Biospies_info_All[which(Biospies_info_All$SID %in% SID_Biopsies_RNAseq),]

# Supplementary figures: Biopsies metadata PCA & MCA ---------------------------------------------
library(factoextra)
Selected_items <- c("BMI..kg.mm2", "LVEF..ultrasound....", 
                    "Collagen.volume.fraction.biopsy...",
                    "Inflammatory.cells.biopsy..CD45...CD68...n.per.mm2",
                    "LVEDD", "LVESD", "LV_mass_index",  "Viral._opiesPerMcgDNA_absolute")

#Biospies_info_All <- Biospies_info_All[which(Biospies_info_All$RNA.seq.Batch!=3),]
Biopsies_metadata.pca <- data.matrix(Biospies_info_All[Selected_items])

Convert_names <- c("BMI", "LVEF", 
                   "Collagen", "Inflammatory", "LVEDD",
                   "LVESD", "LV_mass_index", "Virus")
#Convert_names <- c("Year.of.birth", Convert_names)
colnames(Biopsies_metadata.pca) <- Convert_names

Biopsies_sample_names <- paste0(ifelse(Biospies_info_All$Type == "Control patients", "Heart_failure", "Heart_failure_Cancer"),
                                "_", Biospies_info_All$Treat.with, "_",
                                Biospies_info_All$SID)
rownames(Biopsies_metadata.pca)<- Biopsies_sample_names 

get.pca(Biopsies_metadata.pca, 
        paste0(Biospies_info_All$Type, "_",Biospies_info_All$Treat.with),
        name = "Metadata of the biopsies samples")

get.pca(Biopsies_metadata.pca, 
        paste("Batch", Biospies_info_All$RNA.seq.Batch),
        name = "Metadata of the biopsies samples")

# make MCA
# Use "Inflammation.biopsy" has same values and have no NA (while"Inflammation_biopsy_YN" have NA)
# Virus.in.biopsy..copies.per.microgram.DNA. have more value than "Virus_in_biopsy" - do not know what is true
Selected_items2 <- c("Sex", "Cancer.type", "Treat.with", "Inflammation.biopsy", 
                     "Genetic.mutation", "Comorbidities", "Death")
Biopsies_metadata.quali <- Biospies_info_All[Selected_items2]
Biopsies_metadata.quali$Inflammation.biopsy <- ifelse(Biopsies_metadata.quali$Inflammation.biopsy == "Yes", "Yes", "No")
for (i in 1:ncol(Biopsies_metadata.quali)) Biopsies_metadata.quali[,i]=as.factor(Biopsies_metadata.quali[,i])
Biopsies_metadata.mca <-cbind.data.frame(Biopsies_metadata.pca, Biopsies_metadata.quali)


library(FactoMineR)
res.mca <- MCA(Biopsies_metadata.mca,
               quanti.sup = 1:8,
               quali.sup = 9:10, graph=FALSE) #  the best combination

fviz_mca_ind(res.mca, geom="point",  pointsize = 4, 
             label="none", habillage= as.factor(paste0(Biospies_info_All$Type, "_",Biospies_info_All$Treat.with)),
             addEllipses=TRUE, ellipse.level=0.95,  
             title = paste0("MCA - ", "Biopsies"))


fviz_mca_ind(res.mca, geom="point",  pointsize = 4, 
             label="none", habillage= as.factor(paste("Batch",Biospies_info_All$RNA.seq.Batch)),
             addEllipses=TRUE, ellipse.level=0.95,  
             title = paste0("MCA - ", "Biopsies"))

# Biopsies - RNAseq normalization ----------------------------------------
which(colSums(Biopsies$expected_count)<5000000) # no sample had raw read cound lower than 5 milions

# make the metadata
metadata   <- (Biospies_info_All[c("SID", "Type", "Treat.with", "RNA.seq.Batch", "Matched.to.SID")])
metadata$Biopsies_type <- get.biopsies_type(metadata)
rownames(metadata)   <- paste0("HumanBiopsy_", metadata$SID)
metadata$Biopsies_type <- as.factor(metadata$Biopsies_type)
metadata <- metadata[which(metadata$RNA.seq.Batch != 3),] # remove sample in batch 3
metadata_biopsies <- metadata

Biopsies_readcount <- Biopsies$expected_count[,row.names(metadata)]
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(Biopsies_readcount), 
                              colData = metadata, design = ~ Biopsies_type)
dds2 <- estimateSizeFactors(dds)
norm_data_biopsies <- counts(dds2,normalized=TRUE)
norm_data_biopsies_v1    <- get.remove_gene_low_express(norm_data_biopsies, 0)

# DE lncRNA in biopsies ---------------------------------------------------
dds3 <- DESeq(dds, quiet=TRUE)

# compare "LateCardiotoxicity_with_ANT" vs. "Control"
res <- results(dds3, contrast = c("Biopsies_type", "LateCardiotoxicity_with_ANT", "Control"), 
               pAdjustMethod = "BH")
res.full <- as.data.frame(res)
biopsies_ANTvsCon <-subset(res.full,res.full$padj <= 0.01) # 37 DE genes
biopsies_ANTvsCon_lncRNA <- get.lncRNA(biopsies_ANTvsCon) # 5 DE lncRNAs
unique(Ensemble_database$Gene.name[which(Ensemble_database$Gene.stable.ID %in% biopsies_ANTvsCon_lncRNA)])
#write.table(rownames(biopsies_ANTvsCon), "DEgenes_biopsies.txt",
#            row.names=F,col.names=F,sep="\t", quote=FALSE)

Overlap_invitro_biopsies <- intersect(row.names(biopsies_ANTvsCon), unlist(ANTs_venn))
Outcome <- Ensemble_database[which(Ensemble_database$Gene.stable.ID %in% Overlap_invitro_biopsies),]
unique(Outcome$Gene.name)
unique(Outcome$Gene.type)

write.table(Overlap_invitro_biopsies, "DEgenes_Overlap_invitro_biopsies.txt",
            row.names=F,col.names=F,sep="\t", quote=FALSE)

# identify the DE lncRNA biopsies is differential express in which in vitro condition 
for (i in names(ANTs_venn)) {
  if (length(which(ANTs_venn[[i]] %in% "ENSG00000233478")) >0) print(i)
}

# Figure 1 B: PCA plot for biopsies data: dendrgraom --------------
norm_data_biopsies_v1    <- get.remove_gene_low_express(norm_data_biopsies, 0)
get.pca(t(norm_data_biopsies_v1 ), 
        paste0(metadata$Type, "_",metadata$Treat.with),
        name = "Biopsies samples have RNAseq")
library(tidyverse)
library(dendextend)
library("scales")
clust_data <- merge(metadata, t(norm_data_biopsies_v1), by = "row.names")
rownames(clust_data) <- paste0(clust_data$Biopsies_type, substring(clust_data$Row.names, 12))
clust_data <- clust_data[, -c(1:4)]

dend <- clust_data[, -1] %>% scale %>% dist %>% hclust("ward.D") %>% as.dendrogram %>% 
  set("branches_lwd", 2) %>% set("labels_cex", 0.8)
groupCodes <- clust_data$Biopsies_type
colorCodes <- hue_pal()(5)[2:4] #viridis(5)[1:3]
names(colorCodes) <- unique(groupCodes)
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
pdf("Figure1B.pdf", width=7, height = 5)
par(mar=c(5.1, 10, 4.1, 15))
plot(dend, horiz = TRUE)
dev.off()

# Figure 1 C: Count for gene expression in each gene type in in vitro data
database_sum          <- get.summary_database(norm_data, norm_data_v1)
#database_sum_log      <- get.log(database_sum, log_base = 2)

gene_types <- database_sum[c("protein_coding", "processed_pseudogene", "lncRNA"),]
other_genes <- colSums(database_sum) - colSums(gene_types)
gene_types <- rbind(gene_types, "other_genes" = other_genes)
get.pie_chart(gene_types$`number of gene expression`, row.names(gene_types))

# Load ImpulseDE2 outcome + DE gene categorical in plot -------------------------------------------------
setwd(get.path("outcome"))
load("drug_outcome_NN_2020April02.RData")
names(drug_outcome)
drug_ImpulseDE2Result <- list()
for(drug in names(drug_outcome)) {
  outcome_temp <- drug_outcome[[drug]]@dfImpulseDE2Results
  outcome_temp <- outcome_temp[which(outcome_temp$allZero == FALSE),]
  outcome_temp <- outcome_temp[outcome_temp$padj<0.01,]
  #outcome_temp <- outcome_temp[outcome_temp$mean>100,]
  
  #write.table(rownames(outcome_temp), paste0("DEgenes_", drug, ".txt"),
   #           row.names=F,col.names=F,sep="\t", quote=FALSE)
  drug_ImpulseDE2Result[[drug]] <- outcome_temp
}

#write.table(Ensemble_database$Gene.stable.ID, "Genes_background.txt",
#            row.names=F,col.names=F,sep="\t", quote=FALSE)
# check the p-adj value:
which((drug_ImpulseDE2Result$DOX_The$padj<0.01) == F)
which((drug_ImpulseDE2Result$DOX_Tox$padj<0.01) == F)
which((drug_ImpulseDE2Result$EPI_The$padj<0.01) == F)
which((drug_ImpulseDE2Result$EPI_Tox$padj<0.01) == F)
which((drug_ImpulseDE2Result$IDA_The$padj<0.01) == F)
which((drug_ImpulseDE2Result$IDA_Tox$padj<0.01) == F)

# Figure 1D: Pie chart - overlapped differentially expressed genes --------------------------
ANTs_venn <-get.vennDiam(drug_ImpulseDE2Result, names(drug_ImpulseDE2Result))
#write.table(ANTs_venn$`DOX_The:EPI_The:IDA_The:DOX_Tox:EPI_Tox:IDA_Tox`, "DEgenes_ANTs.txt",
#            row.names=F,col.names=F,sep="\t", quote=FALSE)

Overlap_DE_database <- get.selected_data(ANTs_venn$`DOX_The:EPI_The:IDA_The:DOX_Tox:EPI_Tox:IDA_Tox`,
                                         Ensemble_database, "Gene.stable.ID")
Overlap_DE_summary <- get.summary_gene_type(Overlap_DE_database)

Overlap_DE_gene_type <- Overlap_DE_summary[c("protein_coding", "processed_pseudogene", "lncRNA"),]
other_genes <- sum(Overlap_DE_summary) - sum(Overlap_DE_gene_type)
Overlap_DE_gene_type <- c(Overlap_DE_gene_type, "other_genes" = other_genes)
get.pie_chart(Overlap_DE_gene_type, names(Overlap_DE_gene_type))

# Figure 1E: Group gene types (protein_coding & lncRNA & other RNAs) -----------------------
DE_database   <- list()
DE_summary    <- list()
DE_gene_types <- c()

for (condition in names(drug_ImpulseDE2Result)) {
  DE_database[[condition]] <- get.selected_data(drug_ImpulseDE2Result[[condition]]$Gene, Ensemble_database, "Gene.stable.ID")
  DE_summary[[condition]] <- get.summary_gene_type(DE_database[[condition]])
  data <- as.data.frame(t(DE_summary[[condition]]))
  
  DE_gene_types_tem    <- c()
  DE_gene_types_tem[1] <- data$protein_coding
  DE_gene_types_tem[2] <- data$lncRNA
  DE_gene_types_tem[3] <- data$processed_pseudogene
  DE_gene_types_tem[4] <- sum(data) - (data$protein_coding + data$lncRNA + data$processed_pseudogene)
  DE_gene_types        <- cbind(DE_gene_types, DE_gene_types_tem)
}

rownames(DE_gene_types) <- c("protein_coding", "lncRNA", "processed_pseudogene", "other_RNAs")
colnames(DE_gene_types) <- names(drug_ImpulseDE2Result)

# Pathway anaylsis ------------------------------------
# Figure 2 A: 
setwd("D:/TGX/GitHub/lncRNA-EPI/outcome")
ORA_outcome <- read.csv("ORA_results_ANTs.tab", sep = "\t")
top_pathways <- ORA_outcome[c(1:10),]

top_pathways$pathway <- gsub('.{23}$', '', top_pathways$pathway)
for (i in 1:nrow(top_pathways)) {
  top_pathways$genes[i] <- length(unlist(strsplit(top_pathways$members_input_overlap[i], ";")))
}

library(tidyverse)
ggplot(data = top_pathways) + geom_point(mapping = aes(x=-log(p.value,10), y=pathway, size = genes)) +
  xlab("-log10(p.value)") + xlim(0, 40)

# Figure 2.B: gene types in the pathway anayslis outcome: --> different result ----------------
gene_list<- unique(gsub(" ", "", unlist(strsplit(top_pathways$members_input_overlap, ";")),fixed = TRUE))
gene_list_database <- get.selected_data(gene_list, Ensemble_database, "Gene.stable.ID")
gene_list_summary <- get.summary_gene_type(gene_list_database)

other_genes <- sum(gene_list_summary) - gene_list_summary["protein_coding",]
names(other_genes) <- "other_genes"
gene_types <- c(gene_list_summary["protein_coding",], "lncRNA" = 0, other_genes)
get.pie_chart(gene_types, names(gene_types))

# Table 1: The DE lncRNAs in different ANT treatment conditions --------------------------------
ANTs_venn <-get.vennDiam(drug_ImpulseDE2Result, names(drug_ImpulseDE2Result))

conditions <- c("DOX_The:EPI_The:IDA_The", "DOX_Tox:EPI_Tox:IDA_Tox",
                "DOX_The:DOX_Tox", "EPI_The:EPI_Tox", "IDA_The:IDA_Tox",
                "DOX_The:EPI_The:IDA_The:DOX_Tox:EPI_Tox:IDA_Tox")
DE_lncRNAs_conditions <- list()
for (i in conditions) {
  data_tem <- as.data.frame(ANTs_venn[[i]])
  row.names(data_tem) <- data_tem[,1]
  lncRNA_tem <- get.lncRNA(data_tem)
  DE_lncRNAs_conditions[[i]] <- unique(Ensemble_database$Gene.name[which(Ensemble_database$Gene.stable.ID %in% lncRNA_tem)])
}

# Extract information from lncRNA databases ---------------------------------------------------
# read the outcome of LncTarD
LncTarD_DElncRNAs <- Loaded_data("D:/TGX/GitHub/lncRNA-EPI/outcome", "LncTarD.csv")
unique(LncTarD_DElncRNAs$DiseaseName)
heart_diseases <- c("Acute myocardial infarction", "Aortic valve disease", 
                    "Cardiac hypertrophy", "Heart failure", "Myocardial infarction",
                    "Aortic valve stenosis", "Cardiac fibrosis")
lncRNA_heart_diseases_LncTarD <- get.selected_data(heart_diseases, LncTarD_DElncRNAs, "DiseaseName")

write.xlsx(lncRNA_heart_diseases_LncTarD, file = "outcome.xlsx", sheetName = "lncRNA_heart_diseases_LncTarD", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

# read the outcome of LncRNADisease_v2.0
setwd("D:/TGX/GitHub/lncRNA-EPI/outcome/")
LncRNADisease <- rbind(read.table(file = "LncRNADisease_v2.0_results_heart.xls", sep = "\t", header=TRUE), 
                       read.table(file = "LncRNADisease_v2.0_results_cardi.xls", sep = "\t", header=TRUE))

DElncRNAs <- get.selected_data(unlist(ANTs_venn), Ensemble_database, "Gene.stable.ID")
DElncRNAs <- get.selected_data("lncRNA", DElncRNAs, "Gene.type")
lncRNA_heart_diseases_LncRNADisease <-get.selected_data(DElncRNAs$Gene.name, LncRNADisease, "ncRNA_Symbol")

write.xlsx(lncRNA_heart_diseases_LncRNADisease, file = "outcome.xlsx", sheetName = "lncRNA_heart_diseases_LncRNADisease", 
           col.names = TRUE, row.names = FALSE, append = TRUE)


# Figure 4: lncRNA expression in the in vitro and biopsies samples
# Figure 4: log 2FC of selected genes
# 6 lncRNA genes were differential expressed in all ANT conditions
# RMRP (ENSG00000269900)
# SNHG7 (ENSG00000233016)
# H19 (ENSG00000130600)
# BDNF-AS (ENSG00000245573)
# SNHG29 (ENSG00000175061)
# PCAT19 (ENSG00000267107)
Selected_lncRNAs_1 <- c("ENSG00000269900", "ENSG00000233016", "ENSG00000130600", 
                        "ENSG00000245573", "ENSG00000175061", "ENSG00000267107")

# 7 lncRNAs genes were differential expressed in paticular ANT-treated conditions
# FTX (ENSG00000230590) in DOX-treated samples
# XIST (ENSG00000229807) in EPI-treated samples
# KCNQ1OT1 (ENSG00000269821) in IDA-treated samples
# SNHG22 (ENSG00000267322) in IDA-treated samples
# NEAT1 (ENSG00000245532) in EPI and IDA-treated samples
# MALAT1 (ENSG00000251562) in DOX and IDA-treated samples
# MEG3 (ENSG00000214548) in DOX and EPI-treated samples

Selected_lncRNAs_2 <- c("ENSG00000230590", "ENSG00000229807", "ENSG00000269821",
                        "ENSG00000267322", "ENSG00000245532", "ENSG00000251562",
                        "ENSG00000214548")
Selected_lncRNA <- c(Selected_lncRNAs_1, Selected_lncRNAs_2)

# only use matched pairs biopsies samples
matched_pair_IDs <- intersect(metadata_biopsies$SID, metadata_biopsies$Matched.to.SID)
metadata_biopsies_filled <- metadata_biopsies[which(metadata_biopsies$SID %in% matched_pair_IDs),]

get.Log2FC(Selected_lncRNA, norm_data_invitro,
           norm_data_biopsies, metadata_biopsies_filled)


