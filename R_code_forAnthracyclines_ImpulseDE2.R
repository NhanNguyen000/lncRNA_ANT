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
library(dendextend)
library(viridis)
dend <- t(norm_data_v1) %>% dist %>% hclust("mcquitty") %>% as.dendrogram
groupCodes <- substr(colnames(norm_data_v1), 1, 7)
colorCodes <- viridis_pal()(7)
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
plot(dend, horiz = TRUE)

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
get.pca <- function(data, condition, name) {
  data[is.na(data)] <- 0
  res.pca <- prcomp(data, scale = TRUE)
  fviz_pca_ind(res.pca, geom="point",  pointsize = 2, 
               habillage=condition, addEllipses=TRUE, ellipse.level=0.95,  
               title = paste0("PCA - ", name))
}

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
metadata   <- (Biospies_info_All[c("SID", "Type", "Treat.with", "RNA.seq.Batch")])
metadata$Biopsies_type <- get.biopsies_type(metadata)
metadata$SID        <- paste0("HumanBiopsy_", metadata$SID)
metadata   <- get.rownames(metadata)
metadata$Biopsies_type <- as.factor(metadata$Biopsies_type)
metadata <- metadata[which(metadata$RNA.seq.Batch != 3),]

# remove sampel in batch 3
Biopsies_readcount <- Biopsies$expected_count[,row.names(metadata)]
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(Biopsies$expected_count), 
                              colData = metadata, design = ~ Biopsies_type)
dds2 <- estimateSizeFactors(dds)
norm_data_biopsies <- counts(dds2,normalized=TRUE)

# Figure 1 B: PCA plot for biopsies data: low percetage --> mor to dendrgraom? --------------
norm_data_biopsies_v1    <- get.remove_gene_low_express(norm_data_biopsies, 0)
get.pca(t(norm_data_biopsies_v1 ), 
        paste0(metadata$Type, "_",metadata$Treat.with),
        name = "Biopsies samples have RNAseq")

dend <- t(norm_data_biopsies_v1) %>% dist %>% hclust("mcquitty") %>% as.dendrogram 
plot(dend, horiz = TRUE)

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

# FIgure 1D: Pie chart - overlapped differentially expressed genes --------------------------
ANTs_venn <-get.vennDiam(drug_ImpulseDE2Result, names(drug_ImpulseDE2Result))
get.summary_gene_type(ANTs_venn$`DOX_The:EPI_The:IDA_The:DOX_Tox:EPI_Tox:IDA_Tox`)

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

# Figure 2.B: gene types in the pathway anayslis outcome: --> different result
gene_list<- unique(gsub(" ", "", unlist(strsplit(top_pathways$members_input_overlap, ";")),fixed = TRUE))
gene_list_database <- get.selected_data(gene_list, Ensemble_database, "Gene.stable.ID")
gene_list_summary <- get.summary_gene_type(gene_list_database)

gene_types <- gene_list_summary["protein_coding"]
other_genes <- sum(gene_list_summary) - gene_list_summary["protein_coding",]
gene_types <- c(gene_list_summary["protein_coding",], "lncRNA" = 0,
                          "other_genes" = other_genes)
get.pie_chart(gene_types, names(gene_types))


# in ANTs: DOX_The, EPI_The, IDA_The, DOX_Tox, EPI_Tox, IDA_Tox 
ANTs_venn <-get.vennDiam(drug_ImpulseDE2Result, names(drug_ImpulseDE2Result))
ANTs_venn <-get.geneTypes_overlapLncRNA(Ensemble_database, ANTs_venn[[length(ANTs_venn)]])


Overlap_DE_protein_lncRNA_gene <- ANTs_venn$Sum_GeneType[c("protein_coding", "lncRNA"),]
Overlap_Other_DE_genes <- sum(ANTs_venn$Sum_GeneType) - sum(Overlap_DE_protein_lncRNA_gene)
Overlap_DE_genes <- c(Overlap_DE_protein_lncRNA_gene, Overlap_Other_DE_genes)
DE_gene_types <- cbind(DE_gene_types, Overlap_DE_genes)

Total <- colSums(DE_gene_types)
DE_gene_types <- rbind(DE_gene_types, Total)
library("xlsx")
write.xlsx(DE_gene_types, file = "outcome.xlsx", sheetName = "DE_gene_types", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
# bar chart
#par(mar=c(5,5,5,1))
#barplot(DE_gene_types, las = 1, cex.names = .9, col = c("orange","light green", "blue"),
#        main = "Number of differential expressed genes",  ylab = "Number of genes", ylim =c(0,3500))
#legend(5, 3700, cex = 1, 
#       c("Protein_coding", "LncRNA", "Other_RNAs"),
#       fill = c("orange","light green", "blue"))

## analysis outcome of pathway analysis
pathways <- c("Cardiac muscle contraction - Homo sapiens (human)",
              "Hypertrophic cardiomyopathy (HCM) - Homo sapiens (human)", 
              "Dilated cardiomyopathy (DCM) - Homo sapiens (human)", 
              "Adrenergic signaling in cardiomyocytes - Homo sapiens (human)",
              "Arrhythmogenic right ventricular cardiomyopathy (ARVC) - Homo sapiens (human)",
              "Doxorubicin Pathway (Cardiomyocyte Cell), Pharmacodynamics")

setwd("D:/TGX/GitHub/lncRNA-EPI/outcome")
ORA_outcome <- read.csv("ORA_results_ANTs.tab", sep = "\t")
ORA_heart <- ORA_outcome[which(ORA_outcome$pathway %in% pathways),]
ORA_heart_genes <- c()
ORA_heart_genes_input <- c()
for(i in 1:nrow(ORA_heart)) {
  genes_temp <- trimws(unlist(strsplit(ORA_heart$external_id[i], "[;]")))
  genes_input_temp <- trimws(unlist(strsplit(ORA_heart$members_input_overlap[i], "[;]")))
  
  ORA_heart$member_input_number[i] <- length(genes_input_temp)
  
  ORA_heart_genes <- c(ORA_heart_genes, genes_temp)
  ORA_heart_genes_input <- c(ORA_heart_genes_input, genes_input_temp)
}

Pathway_analysis <-  cbind(ORA_heart$pathway, ORA_heart$source, 
                 ORA_heart$size, ORA_heart$member_input_number, 
                 round(ORA_heart$member_input_number/ORA_heart$size*100, 2), 
                 ORA_heart$p.value, ORA_heart$q.value)
colnames(Pathway_analysis) <- c("Pathway", "Source", "Set size", "Candidate", "Percentage", "p-value", "q-value")

write.xlsx(Pathway_analysis, file = "outcome.xlsx", sheetName = "Pathway_analysis", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

DEgenes_inPathways <- get.selected_data(ORA_heart_genes, Ensemble_database, "Gene.stable.ID")
DElncRNA_inPathways <- DEgenes_inPathways[which(DEgenes_inPathways$Gene.type == "lncRNA"),]
length(unique(DElncRNA_inPathways$Gene.stable.ID))  # no lncRNA

# get total gene in the heart disease pathways:
Convert_gene_data <- read.csv("Accession number mapping - background list_cpdb.csv")

Convert_gene_data2 <- matrix(NA, ncol = 2)
colnames(Convert_gene_data2) <- names(Convert_gene_data)
for (i in 1:nrow(Convert_gene_data)){
  if (grepl(",", Convert_gene_data$Entrez.gene.identifier[i], fixed = TRUE) == FALSE) {
    Convert_gene_data2 <- rbind(Convert_gene_data2, Convert_gene_data[i,])
  } else {
    Convert_gene_data2 <- rbind(Convert_gene_data2, split.Entrez_IDs(Convert_gene_data[i,])) # split gene have multiple Entrez IDs
  }
}
Convert_gene_data2 <-  Convert_gene_data2[-1,]
Convert_gene_data2$Provided.identifier... <- substr(Convert_gene_data2$Provided.identifier..., 1, 15)
Convert_gene_data2$Entrez.gene.identifier <- as.numeric(Convert_gene_data2$Entrez.gene.identifier)

# Check gene in selected pathways
genes_Ensemble_Entrez <- merge(Ensemble_database, Convert_gene_data2, 
                               by.x= "Gene.stable.ID", by.y = "Provided.identifier...", all = TRUE)
cpdb_pathways <- list()
for (pathway in pathways) {
  data_temp <-read.csv(paste0(pathway, ".csv"), header = FALSE)
  colnames(data_temp) <- c("Entrez_ID", "Gene_name")
  cpdb_pathways[[pathway]] <- data_temp
}

table_pathway <- matrix(NA, ncol = length(drug_ImpulseDE2Result), nrow = length(cpdb_pathways))
colnames(table_pathway) <- names(drug_ImpulseDE2Result)
rownames(table_pathway) <- names(cpdb_pathways)
for (condition in names(drug_ImpulseDE2Result)) {
 condition_temp <- get.selected_data(drug_ImpulseDE2Result[[condition]]$Gene, genes_Ensemble_Entrez, "Gene.stable.ID")
 for (pathway in names(cpdb_pathways)) {
   genes_input_temp <- get.selected_data(cpdb_pathways[[pathway]]$Entrez_ID, condition_temp, "Entrez.gene.identifier")
   table_pathway[pathway, condition] <- length(unique(genes_input_temp$Gene.stable.ID))
 } 
}

write.xlsx(table_pathway, file = "outcome.xlsx", sheetName = "Genes_in_pathway_in_condition", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

for (pathway in names(cpdb_pathways)) {
  pathway_temp <- get.selected_data(cpdb_pathways[[pathway]]$Entrez_ID, genes_Ensemble_Entrez, "Entrez.gene.identifier")
  print(unique(pathway_temp$Gene.type))
} 
# --> no lncRNA

# extract list of DE lncRNAs
a<-rbind(drug_ImpulseDE2Result$DOX_The, drug_ImpulseDE2Result$EPI_The, drug_ImpulseDE2Result$IDA_The,
         drug_ImpulseDE2Result$DOX_Tox, drug_ImpulseDE2Result$EPI_Tox, drug_ImpulseDE2Result$IDA_Tox)
b<- get.selected_data(a$Gene, Ensemble_database, "Gene.stable.ID")
DElncRNAs <- b[which(b$Gene.type == "lncRNA"),]
writeLines(unique(DElncRNAs$Gene.stable.ID), "DE_lncRNA.txt", sep = ", ")

# read the outcome of LncTarD
LncTarD_DElncRNAs <- Loaded_data("D:/TGX/GitHub/lncRNA-EPI/outcome", "LncTarD_DElncRNAs.csv")

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
lncRNA_heart_diseases_LncRNADisease <-get.selected_data(DElncRNAs$Gene.name, LncRNADisease, "ncRNA_Symbol")

write.xlsx(lncRNA_heart_diseases_LncRNADisease, file = "outcome.xlsx", sheetName = "lncRNA_heart_diseases_LncRNADisease", 
           col.names = TRUE, row.names = FALSE, append = TRUE)


## extract list DE lncRNA The vs. Tox:
ANT_The<-rbind(drug_ImpulseDE2Result$DOX_The, drug_ImpulseDE2Result$EPI_The, drug_ImpulseDE2Result$IDA_The)
ANT_Tox<-rbind(drug_ImpulseDE2Result$DOX_Tox, drug_ImpulseDE2Result$EPI_Tox, drug_ImpulseDE2Result$IDA_Tox)

DOX<-rbind(drug_ImpulseDE2Result$DOX_The, drug_ImpulseDE2Result$DOX_Tox)
EPI<-rbind(drug_ImpulseDE2Result$EPI_The, drug_ImpulseDE2Result$EPI_Tox)
IDA<-rbind(drug_ImpulseDE2Result$IDA_The, drug_ImpulseDE2Result$IDA_Tox)

a<-get.lncRNA(ANT_The)
b<-get.lncRNA(ANT_Tox)


## Venn diagram & overlap lncRNA ==============
# in ANTs: DOX_The, EPI_The, IDA_The, DOX_Tox, EPI_Tox, IDA_Tox 
ANTs_venn <-get.vennDiam(drug_ImpulseDE2Result, names(drug_ImpulseDE2Result))
#write.table(ANTs_venn$`DOX_The:EPI_The:IDA_The:DOX_Tox:EPI_Tox:IDA_Tox`, 
#            "DEgenes_ANTs.txt", row.names=F,col.names=F,sep="\t", quote=FALSE)
names(ANTs_venn)[length(ANTs_venn)]
ANTs_venn <-get.geneTypes_overlapLncRNA(Ensemble_database, ANTs_venn[[length(ANTs_venn)]])

par(mar=c(15,5,2,1))
get.pie_chart(ANTs_venn$Sum_GeneType, "DE genes in ANTs over time-course")

length(unique(ANTs_venn$Overlap_lncRNA$Gene.name)) #16 overlap lncRNA
#write.table(unique(ANTs_venn$Overlap_lncRNA$Gene.stable.ID),
#            "Ovelap_lncRNA_ANTs_NN_2020Mar09.txt", col.names = FALSE,row.names = FALSE)

# in each ANTs: DOX, EPI, IDA
drug_venn<-list()
drug_venn_lncRNA<-list()

DE_inDrugs <-as.data.frame(matrix(NA, ncol = 3, nrow=2))
names(DE_inDrugs) <- c("DOX", "EPI", "IDA")
rownames(DE_inDrugs) <- c("DE_genes", "DE_lncRNAs")

for(drug in names(DE_inDrugs)){
  conditions <- names(drug_ImpulseDE2Result)[grep(drug, names(drug_ImpulseDE2Result))]
  drug_venn[[drug]] <-get.vennDiam(drug_ImpulseDE2Result, conditions)
  DE_inDrugs["DE_genes", drug] <- length(drug_venn[[drug]][[3]])
  
 # write.table(drug_venn[[drug]][[3]], paste0("DEgenes_ANTs_", drug, ".txt"),
 #             row.names=F,col.names=F,sep="\t", quote=FALSE)
  
  drug_venn_lncRNA[[drug]] <-get.geneTypes_overlapLncRNA(Ensemble_database, drug_venn[[drug]][[3]])
  
  par(mar=c(15,5,2,1))
  get.pie_chart(drug_venn_lncRNA[[drug]]$Sum_GeneType, paste0("DE genes in ", drug ," over time-course"))
  DE_inDrugs["DE_lncRNAs", drug] <- length(unique(drug_venn_lncRNA[[drug]]$Overlap_lncRNA$Gene.name)) #overlap lncRNA
}

# in each dose: The vs. Tox
dose_venn<-list()
dose_venn_lncRNA<-list()

DE_inDose <-as.data.frame(matrix(NA, ncol = 2, nrow=2))
names(DE_inDose) <- c("The", "Tox")
rownames(DE_inDose) <- c("DE_genes", "DE_lncRNAs")

for(dose in names(DE_inDose)){
  conditions <- names(drug_ImpulseDE2Result)[grep(dose, names(drug_ImpulseDE2Result))]
  dose_venn[[dose]] <-get.vennDiam(drug_ImpulseDE2Result, conditions)
  DE_inDose["DE_genes", dose] <- length(dose_venn[[dose]][[7]])
  
  #write.table(dose_venn[[dose]][[7]], paste0("DEgenes_ANTs_", dose, ".txt"),
   #          row.names=F,col.names=F,sep="\t", quote=FALSE)
  
  dose_venn_lncRNA[[dose]] <-get.geneTypes_overlapLncRNA(Ensemble_database, dose_venn[[dose]][[7]])
  
  par(mar=c(15,5,2,1))
  get.pie_chart(dose_venn_lncRNA[[dose]]$Sum_GeneType, paste0("DE genes in ", dose ," dose over time-course"))
  
  DE_inDose["DE_lncRNAs", dose] <- length(unique(dose_venn_lncRNA[[dose]]$Overlap_lncRNA$Gene.name)) #overlap lncRNA
}


# upset graph
library(UpSetR)
i = 2 # lncRNA
b2<- c(DE_gene_types["lncRNA",],
       "DOX_The&DOX_Tox" = DE_inDrugs$DOX[i], "EPI_The&EPI_Tox" = DE_inDrugs$EPI[i], "IDA_The&IDA_Tox" = DE_inDrugs$IDA[i],
       "DOX_The&EPI_The&IDA_The" = DE_inDose$The[i], "DOX_Tox&EPI_Tox&IDA_Tox" = DE_inDose$Tox[i],
       "DOX_The&EPI_The&IDA_The&DOX_Tox&EPI_Tox&IDA_Tox" = ANTs_venn$Sum_GeneType["lncRNA",])

upset(fromExpression(b2), nsets = 6)

unique(dose_venn_lncRNA$The$Overlap_lncRNA$Gene.stable.ID)
unique(dose_venn_lncRNA$Tox$Overlap_lncRNA$Gene.stable.ID)


unique(dose_venn_lncRNA$The$Overlap_lncRNA$Gene.stable.ID) %in% ANT_Tox$Gene
unique(dose_venn_lncRNA$Tox$Overlap_lncRNA$Gene.stable.ID)[which((unique(dose_venn_lncRNA$Tox$Overlap_lncRNA$Gene.stable.ID) %in% ANT_The$Gene) == FALSE)]

DOX_EPI <- rbind(DOX, EPI)
unique(drug_venn_lncRNA$IDA$Overlap_lncRNA$Gene.stable.ID)[which(unique(drug_venn_lncRNA$IDA$Overlap_lncRNA$Gene.stable.ID) %in% DOX_EPI$Gene == FALSE)]


# Make table 1 & 2:  ----------------------------------------------------------
#DE_Dose_table <- data.frame(NA, ncol = 7, nrow = 7)
#colnames(DE_Dose_table)<- c(NA, rep("Each condition", 2), rep("Dose specific", 2), rep("Total overlap", 2))
#DE_Dose_table[1,] <-c(NA, rep(c("Number of differential expression genes", "Number of differential expression lncRNAs"), 3))

DE_Drug_table <- as.data.frame(matrix(NA, ncol = 7, nrow = 6))
colnames(DE_Drug_table)<- c(NA, "Number of DE genes in each condition", "Number of DE lncRNAs in each condition", 
                            "Number of DE genes in drug specific", "Number of DE lncRNAs in drug specific",
                            "Number of DE genes with total overlap", "Number of DE lncRNAs with total overlap")
DE_Drug_table[,1] <- sort(colnames(DE_gene_types))
for (i in 1:nrow(DE_Drug_table)) {
  DE_Drug_table[i,2] <- colSums(DE_gene_types)[DE_Drug_table[i,1]]
  DE_Drug_table[i,3] <- DE_gene_types["lncRNA", DE_Drug_table[i,1]]
  DE_Drug_table[i,4] <- DE_inDrugs["DE_gene", substring(DE_Drug_table[i,1], 1, 3)]
  DE_Drug_table[i,5] <- DE_inDrugs["DE_lncRNAs", substring(DE_Drug_table[i,1], 1, 3)]
  DE_Drug_table[i,6] <- sum(ANTs_venn$Sum_GeneType)
  DE_Drug_table[i,7] <- ANTs_venn$Sum_GeneType["lncRNA",]
  
}

DE_Dose_table <- as.data.frame(matrix(NA, ncol = 7, nrow = 6))
colnames(DE_Dose_table)<- c(NA, "Number of DE genes in each condition", "Number of DE lncRNAs in each condition", 
                            "Number of DE genes in dose specific", "Number of DE lncRNAs in dose specific",
                            "Number of DE genes with total overlap", "Number of DE lncRNAs with total overlap")
DE_Dose_table[,1] <- colnames(DE_gene_types)
for (i in 1:nrow(DE_Dose_table)) {
  DE_Dose_table[i,2] <- colSums(DE_gene_types)[DE_Dose_table[i,1]]
  DE_Dose_table[i,3] <- DE_gene_types["lncRNA", DE_Dose_table[i,1]]
  DE_Dose_table[i,4] <- DE_inDose["DE_gene", substring(DE_Dose_table[i,1], 5)]
  DE_Dose_table[i,5] <- DE_inDose["DE_lncRNAs", substring(DE_Dose_table[i,1], 5)]
  DE_Dose_table[i,6] <- sum(ANTs_venn$Sum_GeneType)
  DE_Dose_table[i,7] <- ANTs_venn$Sum_GeneType["lncRNA",]
}



# 16 overlaped lncRNA expression ------------------------------------------
# lncRNA expression
Selected_data     <- Select_if_in_specific_list(norm_data, ANTs_venn$Overlap_lncRNA$Gene.stable.ID)
#colnames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "ANT_lncRNA_expression_NN_2020Apr09.pdf")

# 4 lncRNAs related to heart failure (RMRP, SNHG7, H19, BDNF-AS) ------------------------------------------
# RMRP (ENSG00000269900)
# SNHG7 (ENSG00000233016)
# H19 (ENSG00000130600)
# BDNF-AS (ENSG00000245573)

Selected_lncRNA <- c("ENSG00000269900", "ENSG00000233016", "ENSG00000130600", "ENSG00000245573")
Selected_data     <- Select_if_in_specific_list(norm_data, Selected_lncRNA)
print.Log2geneExpression(Selected_data)

# 2 lncRNAs not related to heart failure but have improtant function (SNHG29, PCAT19) ------------------------------------------
# SNHG29 (ENSG00000175061)
# PCAT19 (ENSG00000267107)

Selected_lncRNA <- c("ENSG00000175061", "ENSG00000267107")
Selected_data     <- Select_if_in_specific_list(norm_data, Selected_lncRNA)
print.Log2geneExpression(Selected_data)


cor(ENSG00000130600)
norm_data_v1["ENSG00000130600",]
cor_lncRNA <- cor(t(norm_data_v1), norm_data_v1["ENSG00000130600",])
RNAcor_H19_vitro <- rownames(norm_data_v1)[which(abs(cor_lncRNA) > 0.7)]

for(lncRNA in rownames(Selected_data)) {
  get.lncRNA_corrRNA_his(norm_data_v1, lncRNA)
}
print.lncRNA_corrRNA_condition(rownames(Selected_data), norm_data_v1, "hist_selected_lncRNA_NN_2020Apr14.pdf")


## outcome for each ANTs and each dose-------------------------------------
DE_lncRNAs <- list()
DE_lncRNAs[["ANT"]] <- as.vector(unique(ANTs_venn$Overlap_lncRNA$Gene.name))
for (conditions in names(drug_venn_lncRNA)) {
  DE_lncRNAs[[conditions]] <- as.vector(unique(drug_venn_lncRNA[[conditions]]$Overlap_lncRNA$Gene.name))
}

DE_lncRNAs_exclude16lncRNA <- list()
for (conditions in names(drug_venn_lncRNA)) {
  data_temp <- DE_lncRNAs[[conditions]]
  DE_lncRNAs_exclude16lncRNA[[conditions]] <- data_temp[!data_temp %in% DE_lncRNAs[["ANT"]]]
}

n <- 0
for (conditions in names(DE_lncRNAs)) {
  n <- max(n, length(DE_lncRNAs[[conditions]]))
}

for (conditions in names(DE_lncRNAs)) {
 length(DE_lncRNAs[[conditions]]) <- n
}

write.csv(matrix(unlist(DE_lncRNAs), ncol = 6), file = "DE_lncRNAs_lists.csv", row.names = F, col.names = TRUE)

n <- 0
for (conditions in names(DE_lncRNAs_exclude16lncRNA)) {
  n <- max(n, length(DE_lncRNAs_exclude16lncRNA[[conditions]]))
}

for (conditions in names(DE_lncRNAs_exclude16lncRNA)) {
  length(DE_lncRNAs_exclude16lncRNA[[conditions]]) <- n
}

write.csv(matrix(unlist(DE_lncRNAs_exclude16lncRNA), ncol = 5), file = "DE_lncRNAs_exclude16lncRNA_lists.csv", row.names = F, col.names = TRUE)


# lncRNA expression _ chuyen sang gene ID roi convert sang gene name  &&&&&&&&&%%%%%%
DE_lncRNAs <- list()
DE_lncRNAs[["ANT"]] <- as.vector(unique(ANTs_venn$Overlap_lncRNA$Gene.stable.ID))
for (conditions in names(drug_venn_lncRNA)) {
  DE_lncRNAs[[conditions]] <- as.vector(unique(drug_venn_lncRNA[[conditions]]$Overlap_lncRNA$Gene.stable.ID))
}

DE_lncRNAs_exclude16lncRNA <- list()
for (conditions in names(drug_venn_lncRNA)) {
  data_temp <- DE_lncRNAs[[conditions]]
  DE_lncRNAs_exclude16lncRNA[[conditions]] <- data_temp[!data_temp %in% DE_lncRNAs[["ANT"]]]
}

Selected_data     <- Select_if_in_specific_list(norm_data, DE_lncRNAs_exclude16lncRNA$The)
#colnames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "DE_lncRNA_expression_ANT_The_NN_2020Jul14.pdf")

Selected_data     <- Select_if_in_specific_list(norm_data, DE_lncRNAs_exclude16lncRNA$Tox)
#colnames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "DE_lncRNA_expression_ANT_Tox_NN_2020Jul14.pdf")


Selected_data     <- Select_if_in_specific_list(norm_data, DE_lncRNAs_exclude16lncRNA$DOX)
#colnames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "DE_lncRNA_expression_ANT_DOX_NN_2020Jul14.pdf")


Selected_data     <- Select_if_in_specific_list(norm_data, DE_lncRNAs_exclude16lncRNA$EPI)
#colnames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "DE_lncRNA_expression_ANT_EPI_NN_2020Jul14.pdf")


Selected_data     <- Select_if_in_specific_list(norm_data, DE_lncRNAs_exclude16lncRNA$IDA)
#colnames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "DE_lncRNA_expression_ANT_IDA_NN_2020Jul14.pdf")

# Cytoscape & gene network -------------------------------------

library(igraph)


for (condition in names(drug_ImpulseDE2Result)) {
  print(condition)
  data_temp<-norm_data_v1[drug_ImpulseDE2Result[[condition]]$Gene,]
  corr_temp<-cor(t(data_temp)) 
  cor_g <- graph_from_adjacency_matrix(corr_temp, mode='undirected', weighted = 'correlation')
  cor_edge_list <- as_data_frame(cor_g, 'edges')
  cor_edge_list_v2 <- cor_edge_list[which(cor_edge_list$correlation!=1),]
  cor_edge_list_v2 <- cor_edge_list_v2[which(cor_edge_list_v2$correlation>0.7),]
  
  cor_edge_list_v2[,3] <- 1
  write.csv(cor_edge_list_v2, paste0(condition, "_test.csv"), row.names = FALSE)
}


Ensemble_database_v2 <- Ensemble_database[,c("Gene.stable.ID", "Gene.type", "Gene.name", "Gene.description")]
K <- merge(cor_edge_list_v2, Ensemble_database_v2, by.x = "from", by.y = "Gene.stable.ID")
write.csv(K, "Ensemble_test.csv", row.names = FALSE)

g<- K[ which(K$Gene.type == "lncRNA"), ]
write.csv(g, "Ensemble_test2.csv", row.names = FALSE)

# option 1
All_DE_lncRNAs <- unique(unlist(DE_lncRNAs))
All_DE_lncRNAs_ensemble <- get.selected_data(All_DE_lncRNAs, Ensemble_database, "Gene.stable.ID")
write.csv(All_DE_lncRNAs_ensemble, "All_DE_lncRNAs.csv", row.names = FALSE)

All_DE_lncRNAs_ensemble <- read.csv("All_DE_lncRNAs.csv")
# option 2
cutoff = 0.7
CorrDElncRNA_genes <- c()
for(DElncRNA in unique(ANTs_venn$Overlap_lncRNA$Gene.stable.ID)) {
  corr_temp <- array(cor(t(norm_data_v1), norm_data_v1[DElncRNA,]))
  genes <- rownames(norm_data_v1)[which(corr_temp!=1 & abs(corr_temp)>cutoff)]
  corr_genes <- corr_temp[which(corr_temp!=1 & abs(corr_temp)>cutoff)]
  if (length(genes) == 0) {
    outcome_temp <- cbind(DElncRNA, NA, NA)
  } else {
    outcome_temp <- cbind(DElncRNA, genes, corr_genes)
  }
  
  CorrDElncRNA_genes <- rbind(CorrDElncRNA_genes, outcome_temp)
}

which(CorrDElncRNA_genes[,3] ==1)
CorrDElncRNA_genes[226,] # DElncRNA = genes: do not know why code not work for this gene

which(CorrDElncRNA_genes[,3] < 0)

write.csv(CorrDElncRNA_genes, "CorrDElncRNA_genes.csv", row.names = FALSE)

# option 3:
DE_genes <-  get.selected_data(drug_ImpulseDE2Result$DOX_The$Gene, norm_data_v1, "rownames")
selected_gene <- "ENSG00000130600"
data_temp <- get.selected_data(selected_gene, norm_data_v1, "rownames")
a<-t(cor(data_temp , t(DE_genes)))
b<- a[abs(a)>0.7,]
b2<- b[b!=1]

g <-  get.selected_data(names(b2), Ensemble_database, "Gene.stable.ID")

a<-read.csv("C:/Users/nhan.nguyen/Downloads/LncTarD.csv")
intersect(a$Target,  unique(g$Gene.name))

DE_genes_Ensembl <-  get.selected_data(drug_ImpulseDE2Result$DOX_The$Gene, Ensemble_database, "Gene.stable.ID")
intersect(a$Target,  unique(DE_genes_Ensembl$Gene.name))

# Link to biopsies data ---------------------------------------------------
## RUn code for biopsies
# the number of biopsies sample name to run RSEM
for (batch in c("Batch1", "Batch2", "Batch3")) {
  setwd("/ngs-data/data/hecatos/Cardiac/HumanBiopsies/TotalRNA/Batch3/concatenated/")
  print(getwd())
  sample_batch <- list.files()
  sample_batch_name <- matrix(unlist(strsplit(sample_batch, "[_]")), ncol = 3, byrow = TRUE)
  sample_batch_name <- unique(sample_batch_name[,2])
  print(sample_batch_name)
}
---------------
setwd(get.path("data"))
#load("Biopsies_ANT_Cardiac_NN_2020Feb19.RData")
load("Biopsies_ANT_Cardiac_NN_2020May25.RData")
#which(colSums(Biopsies$expected_count)<5000000) # no sample had raw read cound lower than 5 milions
#biopsies_norm_data       <- get.norm_data(Biopsies$expected_count)# by DEseq2


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
                     #"LVEF", "Inflammation_biopsy_YN", "Inflammatory_cells._biopsy")
Biospies_info2 <- Biospies_info2[-which(names(Biospies_info2) %in% Overlap_items)]
Biospies_info_All <- merge(Biospies_info, Biospies_info2, by = "SID", all = TRUE)

SID_Biopsies_RNAseq <- matrix(unlist(strsplit(colnames(Biopsies$expected_count), "_")),ncol = 2, byrow = TRUE)[,2]

Biospies_info_All<-Biospies_info_All[which(Biospies_info_All$SID %in% SID_Biopsies_RNAseq),]

# Biopsies metadata PCA & MCA ---------------------------------------------
library(factoextra)
get.pca <- function(data, condition, name) {
  data[is.na(data)] <- 0
  res.pca <- prcomp(data, scale = TRUE)
  fviz_pca_ind(res.pca, geom="point",  pointsize = 2, 
               habillage=condition, addEllipses=TRUE, ellipse.level=0.95,  
               title = paste0("PCA - ", name))
}

# Use "LVEF..ultrasound...." has same values and have no NA (while "LVEF" have NA)
# Use "Inflammatory.cells.biopsy..CD45...CD68...n.per.mm2" has same values and less NA than "Inflammatory_cells._biopsy" 
# Use "Collagen.volume.fraction.biopsy..." has same values and less NA than "Collagen_volume_fraction_biopsy"

Selected_items <- c("BMI..kg.mm2", "LVEF..ultrasound....", 
                    "Collagen.volume.fraction.biopsy...",
                    "Inflammatory.cells.biopsy..CD45...CD68...n.per.mm2",
                    "LVEDD", "LVESD", "LV_mass_index",  "Viral._opiesPerMcgDNA_absolute")

#Selected_items <- c("Year.of.birth.SID", Selected_items)
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
        name = "Biopsies samples have RNAseq")

get.pca(Biopsies_metadata.pca, 
        paste("Batch", Biospies_info_All$RNA.seq.Batch),
        name = "Biopsies samples have RNAseq")


par(mar = c(15, 10, 10, 15))
heatmap(Biopsies_metadata.pca, cexRow = 0.75, cexCol = 0.8)

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

# correlation of heart failure items
summary(as.matrix(Biopsies_metadata.mca)) # summary information sample is in the ppt slides
round(cor(Biopsies_metadata.pca[,-(which(colnames(Biopsies_metadata.pca) %in% c("BMI", "Virus")))],
          use="complete.obs") ,2)

metadata   <- (Biospies_info_All[c("SID", "Type", "Treat.with", Selected_items)])
metadata$Biopsies_type <- get.biopsies_type(metadata)
metadata$SID        <- paste0("HumanBiopsy_", metadata$SID)

metadata   <- merge(metadata, 
                    as.data.frame(colnames(Biopsies$expected_count)), 
                    by.x = "SID", by.y = "colnames(Biopsies$expected_count)")
#identical(metadata$SID, colnames(Biopsies$expected_count))
metadata   <- get.rownames(metadata)
metadata$Biopsies_type <- as.factor(metadata$Biopsies_type)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(Biopsies$expected_count), 
                              colData = metadata, design = ~ Biopsies_type)
dds2 <- estimateSizeFactors(dds)
norm_data_biopsies <- counts(dds2,normalized=TRUE)
norm_data_biopsies_v1    <- get.remove_gene_low_express(norm_data_biopsies, 0)

# DE lncRNA in biopsies ---------------------------------------------------
dds3 <- DESeq(dds, quiet=TRUE)

# compare "LateCardiotoxicity_with_ANT" vs. "Control"
res <- results(dds3, contrast = c("Biopsies_type", "LateCardiotoxicity_with_ANT", "Control"), pAdjustMethod = "BH")
res.full <- as.data.frame(res)
biopsies_ANTvsCon <-subset(res.full,res.full$padj <= 0.01) # 81 DE genes
biopsies_ANTvsCon_lncRNA <- get.lncRNA(biopsies_ANTvsCon) # 23 DE lncRNAs
#write.csv(biopsies_ANTvsCon, "Biopsies_ANT_vs_Control_NN_2020May28.csv")

# compare "LateCardiotoxicity" vs. "Control"
res <- results(dds3, contrast = c("Biopsies_type", "LateCardiotoxicity", "Control"), pAdjustMethod = "BH")
res.full <- as.data.frame(res)
biopsies_noANTvsCon <-subset(res.full,res.full$padj <= 0.01) # 128 genes
biopsies_noANTvsCon_lncRNA <- get.lncRNA(biopsies_noANTvsCon) # 32 lncRNAs
#write.csv(biopsies_noANTvsCon, "Biopsies_withoutANT_vs_Control_NN_2020May28.csv")

intersect(biopsies_ANTvsCon_lncRNA, biopsies_noANTvsCon_lncRNA) # 11 lncRNAs overlap
length(intersect(biopsies_ANTvsCon_lncRNA, Selected_lncRNA)) # no overlap lncRNAs
length(intersect(biopsies_noANTvsCon_lncRNA, Selected_lncRNA)) # no overlap lncRNAs

# Biopsies lncRNA correlation ---------------------------------------------
Biopsies_lncRNA <- get.lncRNA(norm_data_biopsies_v1)
Biopsies_lncRNA_data <- norm_data_biopsies_v1[as.vector(Biopsies_lncRNA),]

heart_failure_items <-data.matrix(metadata[5:10])

Selected_lncRNA <- c("ENSG00000269968", "ENSG00000244239", "ENSG00000269900", 
                     "ENSG00000260941", "ENSG00000233016", "ENSG00000251095", 
                     "ENSG00000175061", "ENSG00000130600", "ENSG00000240801",
                     "ENSG00000261069", "ENSG00000233521", "ENSG00000268518",
                     "ENSG00000250159", "ENSG00000245573", "ENSG00000271401",
                     "ENSG00000267107") # 16 DElncRNAs
cor_16DE_lncRNAs <- get.cor_lncRNA_phenotype(heart_failure_items, 
                                             Biopsies_lncRNA_data[Selected_lncRNA,])
cor_all_lncRNAs <-get.cor_lncRNA_phenotype(heart_failure_items, Biopsies_lncRNA_data)

cor_16DE_lncRNAs_log <- get.cor_lncRNA_phenotype(heart_failure_items, 
                                                 log(Biopsies_lncRNA_data[Selected_lncRNA,] + 1))
cor_all_lncRNAs_log <-get.cor_lncRNA_phenotype(heart_failure_items, 
                                               log(Biopsies_lncRNA_data + 1))

cor_16DE_lncRNAs_TPM <- get.cor_lncRNA_phenotype(heart_failure_items, 
                                                 Biopsies$TPM[Selected_lncRNA,])

cor_all_lncRNAs_TPM <-get.cor_lncRNA_phenotype(heart_failure_items, 
                                               Biopsies$TPM[Biopsies_lncRNA,])

cor_16DE_lncRNAs_FPKM <- get.cor_lncRNA_phenotype(heart_failure_items, 
                                                 Biopsies$FPKM[Selected_lncRNA,])

cor_all_lncRNAs_FPKM <-get.cor_lncRNA_phenotype(heart_failure_items, 
                                               Biopsies$FPKM[Biopsies_lncRNA,])

## plot the data:
dat <- merge(t(Biopsies_lncRNA_data[Selected_lncRNA,]), 
             metadata,
             by="row.names")
plot(as.numeric(dat$ENSG00000251095), as.numeric(dat$Collagen.volume.fraction.biopsy...), 
     xlab="Read count of ENSG00000251095 (AC093866.1)", ylab="Collagen_volume_fraction_biopsy",
     pch=19, xlim=c(0,max(as.numeric(dat$ENSG00000251095))), ylim=c(0,60),
     col=ifelse(dat$Biopsies_type == "Control", "red", 
                ifelse(dat$Biopsies_type == "LateCardiotoxicity", "green", "blue")))

legend("right", legend=c("Control","LateCardiotoxicity", "LateCardiotoxicity_ANT"), 
       col=c("red","green","blue"),
       pch=c(19,19), inset=.02)


lncRNA <- "ENSG00000251095"
dat2 <- dat[c("Row.names", lncRNA, "Collagen.volume.fraction.biopsy...")]
dat2 <- get.rownames(dat2)
dat2 <- na.omit(data.matrix(dat2))
rownames(dat2)[which(dat2[,lncRNA] == max(dat2[,lncRNA]))] # "HumanBiopsy_10033"

# TPM
dat <- merge(t(Biopsies$TPM[Selected_lncRNA,]), 
             metadata,
             by="row.names")
plot(as.numeric(dat$ENSG00000233521), as.numeric(dat$Inflammatory.cells.biopsy..CD45...CD68...n.per.mm2), 
     xlab="Read count of ENSG00000233521 (LINC01638)", ylab="Inflammatory.cells.biopsy..CD45...CD68...n.per.mm2",
     pch=19, xlim=c(0,max(as.numeric(dat$ENSG00000233521))), ylim=c(0,70),
     col=ifelse(dat$Biopsies_type == "Control", "red", 
                ifelse(dat$Biopsies_type == "LateCardiotoxicity", "green", "blue")))

legend("right", legend=c("Control","LateCardiotoxicity", "LateCardiotoxicity_ANT"), 
       col=c("red","green","blue"),
       pch=c(19,19), inset=.02)

lncRNA <- "ENSG00000233521"
dat2 <- dat[c("Row.names", lncRNA, "Inflammatory.cells.biopsy..CD45...CD68...n.per.mm2")]
dat2 <- get.rownames(dat2)
dat2 <- na.omit(data.matrix(dat2))
rownames(dat2)[which(dat2[,lncRNA] == max(dat2[,lncRNA]))] #"HumanBiopsy_10033"

# FPKM
dat <- merge(t(Biopsies$FKPM[Selected_lncRNA,]), 
             metadata,
             by="row.names")
plot(as.numeric(dat$ENSG00000260941), as.numeric(dat$Collagen.volume.fraction.biopsy...), 
     xlab="Read count of ENSG00000260941 (LINC00622)", ylab="Collagen.volume.fraction.biopsy...",
     pch=19, xlim=c(0,max(as.numeric(dat$ENSG00000260941))), ylim=c(0,70),
     col=ifelse(dat$Biopsies_type == "Control", "red", 
                ifelse(dat$Biopsies_type == "LateCardiotoxicity", "green", "blue")))

legend("top", legend=c("Control","LateCardiotoxicity", "LateCardiotoxicity_ANT"), 
       col=c("red","green","blue"),
       pch=c(19,19), inset=.02)


lncRNA <- "ENSG00000260941"
dat2 <- dat[c("Row.names", lncRNA, "Collagen.volume.fraction.biopsy...")]
dat2 <- get.rownames(dat2)
dat2 <- na.omit(data.matrix(dat2))
rownames(dat2)[which(dat2[,lncRNA] == max(dat2[,lncRNA]))] # ""HumanBiopsy_10027"


# calculate the lncRNA in biospies corre with heart failure measurement:
which((cor_all_lncRNAs[,3]>0.8) == TRUE)
which((cor_all_lncRNAs_TPM[,3]>0.8) == TRUE)
which((cor_all_lncRNAs_FPKM[,3]>0.8) == TRUE)

## Code until here: 
norm_data_biopsies_v1["ENSG00000130600",]
cor_lncRNA2 <- cor(t(norm_data_biopsies_v1), norm_data_biopsies_v1["ENSG00000130600",])
RNAcor_H19_biopsies <- rownames(norm_data_biopsies_v1)[which(abs(cor_lncRNA2) > 0.7)]
get.lncRNA_corrRNA_his(norm_data_biopsies_v1 , "ENSG00000130600")
intersect(RNAcor_H19_biopsies, RNAcor_H19_vitro)
# the high correlation RNA with H!( in invitro == in biopsies)

for (i in RNAcor_H19_vitro) {
  a <- cor(norm_data_biopsies_v1["ENSG00000130600",], norm_data_biopsies_v1[i,])
  print(a)
}


for (lncRNA in Selected_lncRNA) {
  cor_lncRNA1 <- cor(t(norm_data_v1), norm_data_v1[lncRNA,])
  RNAcor_invitro <- rownames(norm_data_v1)[which(abs(cor_lncRNA1) > 0.7)]
  
  cor_lncRNA2 <- cor(t(norm_data_biopsies_v1), norm_data_biopsies_v1[lncRNA,])
  RNAcor_biopsies <- rownames(norm_data_biopsies_v1)[which(abs(cor_lncRNA2) > 0.7)]
  
  
  if (length(intersect(RNAcor_biopsies, RNAcor_invitro)) >1) {
    print(lncRNA)
    print(intersect(RNAcor_biopsies, RNAcor_invitro))
  }
}

# "ENSG00000269900" (RMRP) with 27-1 = 26 genes cor in invitro + Biopsies
"ENSG00000200463" -SNORD118
"ENSG00000200792" -SNORA80A
"ENSG00000201512" -SNORA71C
"ENSG00000206597" -SNORA57 
"ENSG00000206811" -SNORA10
"ENSG00000206838" -SNORA5A
"ENSG00000207088" -SNORA7B
"ENSG00000207405" -SNORA64
"ENSG00000207496" -SNORA7A
"ENSG00000209582" -SNORA48
[11] "ENSG00000222489" -SNORA79B
"ENSG00000235408" -SNORA71B
"ENSG00000238917" -SNORD10
"ENSG00000239039" -SNORD13
"ENSG00000239149" -SNORA59A
[16] "ENSG00000251791" -SCARNA6
"ENSG00000251898" -SCARNA11
"ENSG00000253190" -AC084082.1 (noval transcript)
"ENSG00000260035" -AC051619.7 (noval transcript)
"ENSG00000263934" -SNORD3A
[21] "ENSG00000264940" -SNORD3C
"ENSG00000266079" -SNORA59B
"ENSG00000269900" -RMRP
"ENSG00000271798" -SNORA51
"ENSG00000277184" -SNORA9
[26] "ENSG00000278274" -SNORA61
"ENSG00000281808" -SNORA17

#RNase MRP (also called RMRP) is an enzymatically active ribonucleoprotein with two distinct roles in eukaryotes. RNAse MRP stands for RNAse for mitochondrial RNA processing. In mitochondria it plays a direct role in the initiation of mitochondrial DNA replication. In the nucleus it is involved in precursor rRNA processing, where it cleaves the internal transcribed spacer 1 between 18S and 5.8S rRNAs.[1] Despite distinct functions, RNase MRP has been shown to be evolutionarily related to RNase P. Like eukaryotic RNase P, RNase MRP is not catalytically active without associated protein subunits.
#--> https://www.sciencedirect.com/science/article/pii/S0300908416303972?casa_token=-2l9E1tYC3UAAAAA:raCSN03jqehwh0R2KXnDZhyaVY9Q-GJ7wd33CpXQR9dsga7t-nLj0ISi3gnM7thx3IodJNvz
#---> https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-017-0149-x

# ENSG00000233016 (SNHG7) with 1 gene (ENSG00000199293 = SNORA21) cor in invitro + Biopsies
# ENSG00000233521 (LINC01638) with 1 gene (ENSG00000273008 = AC010864.1) cor in invitro + Biopsies





# LncRNA with Heart function ----------------------------------------------
Selected_lncRNA <- c("ENSG00000269900", "ENSG00000233016", "ENSG00000130600", "ENSG00000245573")
#Selected_lncRNA <- unique(ANTs_venn$Overlap_lncRNA$Gene.stable.ID)
Selected_lncRNA <- c("ENSG00000269968", "ENSG00000244239", "ENSG00000269900", 
                     "ENSG00000260941", "ENSG00000233016", "ENSG00000251095", 
                     "ENSG00000175061", "ENSG00000130600", "ENSG00000240801",
                     "ENSG00000261069", "ENSG00000233521", "ENSG00000268518",
                     "ENSG00000250159", "ENSG00000245573", "ENSG00000271401",
                     "ENSG00000267107")
Selected_data <- t(Select_if_in_specific_list(norm_data_biopsies, Selected_lncRNA))
Selected_data <- t(Select_if_in_specific_list(Biopsies$TPM, Selected_lncRNA))
Selected_data <- t(Select_if_in_specific_list(Biopsies$FPKM, Selected_lncRNA))


# "LVEF..ultrasound...."

# read count
selection_heartSample <- (names(Sample_info)[c(4, 6:10)])

for (heart_measure in selection_heartSample) {
  Sample_info_temp <-Sample_info[c("SID", heart_measure)]
  print(heart_measure)
  Sample_info_temp[,2] <- as.numeric(Sample_info_temp[,2])
  a<-merge(Selected_data, Sample_info_temp, by.x="row.names", by.y = "SID")
  a1<-get.rownames(a)
  for(i in Selected_lncRNA ) {
    k<-cor(a1[i], a1[heart_measure],  use="complete.obs") # not strong correlation [0-0.2] with log or nomarl read count
    if (abs(k)>0.6) {print(k)}
  }
  
}

#--> results
[1] "LVEF..ultrasound...."
[1] "LVEF"
[1] "LVEDD"
[1] "LVESD"
[1] "LV_mass_index"
[1] "Collagen_volume_fraction_biopsy"
Collagen_volume_fraction_biopsy
ENSG00000251095                       0.7441094 #AC093866.1

dat <- merge(Selected_data, Sample_info[c("SID", "Treat.with", "Collagen_volume_fraction_biopsy")], by.x="row.names", by.y = "SID")
plot(as.numeric(dat$ENSG00000251095), as.numeric(dat$Collagen_volume_fraction_biopsy), 
     xlab="Read count of ENSG00000251095 (AC093866.1)", ylab="Collagen_volume_fraction_biopsy",
     pch=19, xlim=c(0,max(as.numeric(dat$ENSG00000251095))), ylim=c(0,60),
     col=ifelse(dat$Treat.with == "ANT", "red", "blue"))
legend("right", legend=c("ANT","non-ANT"), col=c("red","blue"),
       pch=c(19,19), inset=.02)


## Log read count
selection_heartSample <- (names(Sample_info)[c(4, 6:10)])

for (heart_measure in selection_heartSample) {
  Sample_info_temp <-Sample_info[c("SID", heart_measure)]
  print(heart_measure)
  Sample_info_temp[,2] <- as.numeric(Sample_info_temp[,2])
  a<-merge(Selected_data, Sample_info_temp, by.x="row.names", by.y = "SID")
  a1<-get.rownames(a)
  for(i in Selected_lncRNA ) {
    k<-cor(log(a1[i] + 1), a1[heart_measure],  use="complete.obs") # not strong correlation [0-0.2] with log or nomarl read count
    if (abs(k)>0.6) {print(k)}
  }
  
}

# --> result or read count:
[1] "LVEF..ultrasound...."
[1] "LVEF"
[1] "LVEDD"
[1] "LVESD"
[1] "LV_mass_index"
[1] "Collagen_volume_fraction_biopsy"
Collagen_volume_fraction_biopsy
ENSG00000251095                       0.7088608 # AC093866.1
Collagen_volume_fraction_biopsy
ENSG00000130600                      -0.6869546 # H19


dat <- merge(Selected_data, Sample_info[c("SID", "Treat.with", "Collagen_volume_fraction_biopsy")], by.x="row.names", by.y = "SID")
plot(log(as.numeric(dat$ENSG00000251095)+1), as.numeric(dat$Collagen_volume_fraction_biopsy), 
     xlab="Log read count of ENSG00000251095 (AC093866.1)", ylab="Collagen_volume_fraction_biopsy",
     pch=19, xlim=c(0,7),ylim=c(0,60),
     col=ifelse(dat$Treat.with == "ANT", "red", "blue"))
legend("right", legend=c("ANT","non-ANT"), col=c("red","blue"),
       pch=c(19,19), inset=.02)

dat <- merge(Selected_data, Sample_info[c("SID", "Treat.with", "Collagen_volume_fraction_biopsy")], by.x="row.names", by.y = "SID")
plot(log(as.numeric(dat$ENSG00000130600 )+1), as.numeric(dat$Collagen_volume_fraction_biopsy), 
     xlab="Log read count of ENSG00000130600 (H19)", ylab="Collagen_volume_fraction_biopsy",
     pch=19, xlim=c(5,9),ylim=c(0,60),
     col=ifelse(dat$Treat.with == "ANT", "red", "blue"))
legend("topright", legend=c("ANT","non-ANT"), col=c("red","blue"),
       pch=c(19,19), inset=.02)


##--> result wit TPM:
selection_heartSample <- (names(Sample_info)[c(4, 6:10)])

for (heart_measure in selection_heartSample) {
  Sample_info_temp <-Sample_info[c("SID", heart_measure)]
  print(heart_measure)
  Sample_info_temp[,2] <- as.numeric(Sample_info_temp[,2])
  a<-merge(Selected_data, Sample_info_temp, by.x="row.names", by.y = "SID")
  a1<-get.rownames(a)
  for(i in Selected_lncRNA ) {
    k<-cor(a1[i], a1[heart_measure],  use="complete.obs") # not strong correlation [0-0.2] with log or nomarl read count
    if (abs(k)>0.6) {print(k)}
  }
  
}

[1] "LVEF..ultrasound...."
[1] "LVEF"
[1] "LVEDD"
[1] "LVESD"
[1] "LV_mass_index"
[1] "Collagen_volume_fraction_biopsy"
Collagen_volume_fraction_biopsy
ENSG00000251095                       0.7438361 #AC093866.1


dat <- merge(Selected_data, Sample_info[c("SID", "Treat.with", "Collagen_volume_fraction_biopsy")], by.x="row.names", by.y = "SID")
plot(as.numeric(dat$ENSG00000251095), as.numeric(dat$Collagen_volume_fraction_biopsy), 
     xlab="TPM value of ENSG00000251095 (0.7438361)", ylab="Collagen_volume_fraction_biopsy",
     pch=19, xlim=c(0,20),ylim=c(0,60),
     col=ifelse(dat$Treat.with == "ANT", "red", "blue"))
legend("right", legend=c("ANT","non-ANT"), col=c("red","blue"),
       pch=c(19,19), inset=.02)

## --> result with FPKM:
selection_heartSample <- (names(Sample_info)[c(4, 6:10)])

for (heart_measure in selection_heartSample) {
  Sample_info_temp <-Sample_info[c("SID", heart_measure)]
  print(heart_measure)
  Sample_info_temp[,2] <- as.numeric(Sample_info_temp[,2])
  a<-merge(Selected_data, Sample_info_temp, by.x="row.names", by.y = "SID")
  a1<-get.rownames(a)
  for(i in Selected_lncRNA ) {
    k<-cor(a1[i], a1[heart_measure],  use="complete.obs") # not strong correlation [0-0.2] with log or nomarl read count
    if (abs(k)>0.6) {print(k)}
  }
  
}

[1] "LVEF..ultrasound...."
[1] "LVEF"
[1] "LVEDD"
[1] "LVESD"
[1] "LV_mass_index"
[1] "Collagen_volume_fraction_biopsy"
Collagen_volume_fraction_biopsy
ENSG00000251095                       0.7440069

##
dat <- merge(Selected_data, Sample_info[c("SID", "Treat.with", "Collagen_volume_fraction_biopsy")], by.x="row.names", by.y = "SID")
plot(as.numeric(dat$ENSG00000251095), as.numeric(dat$Collagen_volume_fraction_biopsy), 
     xlab="FPKM value of ENSG00000251095 (0.7438361)", ylab="Collagen_volume_fraction_biopsy",
     pch=19, xlim=c(0,25),ylim=c(0,60),
     col=ifelse(dat$Treat.with == "ANT", "red", "blue"))
legend("right", legend=c("ANT","non-ANT"), col=c("red","blue"),
       pch=c(19,19), inset=.02)



# Check patient group with only ANT treatment??
#--------------------
Sample_info_temp <-Sample_info[c("SID", "LVEF..ultrasound....")]
Sample_info_temp$LVEF..ultrasound.... <- as.numeric(Sample_info_temp$LVEF..ultrasound....)
a<-merge(Selected_data, Sample_info_temp, by.x="row.names", by.y = "SID")
a1<-get.rownames(a)
for(i in Selected_lncRNA ) {
  k<-cor(log(a[[i]] + 1), a[,2],  use="complete.obs") # not strong correlation [0-0.2] with log or nomarl read count
  print(k)
  plot(log(a[[i]]+1), a$LVEF..ultrasound....)
}

# LV mass index
# NYHA class
# Fractional shorting /morphology

# other code --------------------------------------------------------------
#a<- cor(t(Selected_data))
#b<- sort(as.vector(a), decreasing = TRUE) # 2 lncRNA have r > 0.7

# histogram of correlation 
for(lncRNA in rownames(Selected_data)) {
  get.lncRNA_corrRNA_his(norm_data_v1, lncRNA)
}
print.lncRNA_corrRNA_condition(rownames(Selected_data), norm_data_v1, "hist_selected_lncRNA_NN_2020Apr14.pdf")




# Venn diagram for only lncRNA
lncRNA_database   <- get.selected_data("lncRNA", Ensemble_database, column_in_database = "Gene.type")
ImpulseDE2Result_lncRNA <- list()
for(condition in names(drug_ImpulseDE2Result)) {
  ImpulseDE2Result_lncRNA[[condition]] <- get.selected_data(unique(lncRNA_database$Gene.stable.ID), 
                                                            drug_ImpulseDE2Result[[condition]], 
                                                            column_in_database = "rownames")
}

ANTs_venn <-get.vennDiam(ImpulseDE2Result_lncRNA, names(ImpulseDE2Result_lncRNA))
names(ANTs_venn)[length(ANTs_venn)]

drug_venn<- list()
for(drug in c("DOX", "EPI", "IDA")){
  conditions <- names(ImpulseDE2Result_lncRNA)[grep(drug, names(ImpulseDE2Result_lncRNA))]
  drug_venn[[drug]] <-get.vennDiam(ImpulseDE2Result_lncRNA, conditions)
}

dose_venn <- list()
for(dose in c("The", "Tox")){
  conditions <- names(ImpulseDE2Result_lncRNA)[grep(dose, names(ImpulseDE2Result_lncRNA))]
  dose_venn[[dose]] <-get.vennDiam(ImpulseDE2Result_lncRNA, conditions)
}

lncRNA_forThe <- setdiff(dose_venn$The$`DOX_The:EPI_The:IDA_The`,
                         ANTs_venn$`DOX_The:EPI_The:IDA_The:DOX_Tox:EPI_Tox:IDA_Tox`)
Selected_data     <- Select_if_in_specific_list(norm_data, lncRNA_forThe)
#colnames(Selected_data); rownames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "ANT_The_lncRNA_expression_NN_2020Apr015.pdf")
unique(Ensemble_database$Gene.name[which(Ensemble_database$Gene.stable.ID %in% lncNRA_forThe)])

lncNRA_forTox <- setdiff(dose_venn$Tox$`DOX_Tox:EPI_Tox:IDA_Tox`, 
                         ANTs_venn$`DOX_The:EPI_The:IDA_The:DOX_Tox:EPI_Tox:IDA_Tox`)
Selected_data     <- Select_if_in_specific_list(norm_data, lncNRA_forTox)
#colnames(Selected_data); rownames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "ANT_Tox_lncRNA_expression_NN_2020Apr015.pdf")
unique(Ensemble_database$Gene.name[which(Ensemble_database$Gene.stable.ID %in% lncNRA_forTox)])

## other way which give other result:
get.geneName(ANTs_venn$`DOX_The:DOX_Tox`, Ensemble_database) # 6 lncRNAs

Selected_data     <- Select_if_in_specific_list(norm_data, ANTs_venn$`DOX_The:DOX_Tox`)
print.geneExpression(Selected_data , "ANT_Dox_lncRNA_expression_NN_2020Apr015.pdf")

get.geneName(ANTs_venn$`EPI_The:EPI_Tox`, Ensemble_database) # 13 lncRNAs
get.geneName(ANTs_venn$`IDA_The:IDA_Tox`, Ensemble_database) # 17 lncRNAs

get.geneName(ANTs_venn$`DOX_Tox:EPI_Tox:IDA_Tox`, Ensemble_database) # 2 lncRNAs
get.geneName(ANTs_venn$`DOX_The:EPI_The:IDA_The`, Ensemble_database) # 1 lncRNAs

Selected_data     <- Select_if_in_specific_list(norm_data, 
                                                c(ANTs_venn$`DOX_The:EPI_The:IDA_The`,
                                                  ANTs_venn$`DOX_Tox:EPI_Tox:IDA_Tox`))
print.geneExpression(Selected_data , "ANT_Tox_The_lncRNA_expression_NN_2020Apr015.pdf")

## RUn code for biopsies
setwd(get.path("data"))
load("Biopsies_ANT_Cardiac_NN_2020Feb19.RData")
#which(colSums(Biopsies$expected_count)<5000000) # no sample had raw read cound lower than 5 milions
#biopsies_norm_data       <- get.norm_data(Biopsies$expected_count)# by DEseq2

Biospies_info <- read.table(paste0(get.path("data"), 
                                   "/ClinicalParametersHeCaToS_Patients2018FINAL_summary_NN_2020Feb26.txt"), 
                            header = TRUE, sep = "\t", fill= TRUE)
Sample_info  <- as.data.frame(as.matrix(Biospies_info[, c("SID", "Type", "RNA.seq.Batch","Treat.with")]),  stringsAsFactors=FALSE)
Sample_info$Biopsies_type <- get.biopsies_type(Sample_info)
Sample_info$SID        <- paste0("HumanBiopsy_", Sample_info$SID)

metadata   <- merge(Sample_info, as.data.frame(colnames(Biopsies$expected_count)), 
                    by.x = "SID", by.y = "colnames(Biopsies$expected_count)")
metadata   <- get.rownames(metadata)
metadata$Biopsies_type <- as.factor(metadata$Biopsies_type)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(Biopsies$expected_count), 
                              colData = metadata, design = ~ Biopsies_type)
dds2 <- estimateSizeFactors(dds)
norm_data_biopsies <- counts(dds2,normalized=TRUE)
norm_data_v1    <- get.remove_gene_low_express(norm_data, 0)

dds3 <- DESeq(dds, quiet=TRUE) 

# compare "LateCardiotoxicity_with_ANT" vs. "Control"
res <- results(dds3, contrast = c("Biopsies_type", "LateCardiotoxicity_with_ANT", "Control"), pAdjustMethod = "BH")
res.full <- as.data.frame(res)
biopsies_ANTvsCon <-subset(res.full,res.full$padj <= 0.01) # 89 DE genes
biopsies_ANTvsCon_lncRNA <- get.lncRNA(biopsies_ANTvsCon) # 21 DE lncRNAs
write.csv(biopsies_ANTvsCon, "Biopsies_ANT_vs_Control_NN_2020Apr16.csv")

# compare "LateCardiotoxicity" vs. "Control"
res <- results(dds3, contrast = c("Biopsies_type", "LateCardiotoxicity", "Control"), pAdjustMethod = "BH")
res.full <- as.data.frame(res)
biopsies_noANTvsCon <-subset(res.full,res.full$padj <= 0.01) # 302 genes
biopsies_noANTvsCon <- get.lncRNA(biopsies_noANTvsCon)
write.csv(biopsies_noANTvsCon, "Biopsies_withoutANT_vs_Control_NN_2020Apr16.csv")

setwd(get.path("outcome"))
#biopsies_ANTvsCon <- read.csv("Biopsies_ANT_vs_Control_NN_2020Apr16.csv")
#biopsies_ANTvsCon <- get.rownames(biopsies_ANTvsCon)
#biopsies_noANTvsCon <- read.csv("Biopsies_withoutANT_vs_Control_NN_2020Apr16.csv")
#biopsies_noANTvsCon <- get.rownames(biopsies_noANTvsCon)

biopsies_CancervsCon <-intersect(rownames(biopsies_ANTvsCon), rownames(biopsies_noANTvsCon)) #41 DE genes overlap
#gene_information <- Ensemble_database[which(Ensemble_database$Gene.stable.ID %in% biopsies_CancervsCon),]
#lncRNA_data <- get.selected_data("lncRNA", gene_information, "Gene.type")
#output <- unique(lncRNA_data$Gene.stable.ID)

biopsies_ANTvsCon_gene <- merge(Ensemble_database, biopsies_ANTvsCon, by.x = "Gene.stable.ID", by.y = "row.names")
biopsies_noANTvsCon_gene <- merge(Ensemble_database, biopsies_noANTvsCon, by.x = "Gene.stable.ID", by.y = "row.names")

biopsies_only_in_ANTvsCon <- biopsies_ANTvsCon[-which(rownames(biopsies_ANTvsCon) %in% biopsies_CancervsCon),]
# 48 genes only DE in ANT samples (not cancer_noANT)

biopsies_only_in_ANTvsCon_lncRNA <- get.lncRNA(biopsies_only_in_ANTvsCon) 
# 13 lncRNA (out of 48 genes above only DE in ANT samples (not cancer_noANT))

setwd(get.path("outcome"))
Invitro_Overlap_lncRNA_ANTs <- read.csv("Ovelap_lncRNA_ANTs_NN_2020Mar09.txt", header = FALSE)
intersect(as.vector(biopsies_only_in_ANTvsCon_lncRNA), Invitro_Overlap_lncRNA_ANTs) # no overlap genes

intersect(unique(biopsies_ANTvsCon_gene$Gene.stable.ID), Invitro_Overlap_lncRNA_ANTs) # no overlap genes


## print gene expresison in vitro
Selected_data     <- Select_if_in_specific_list(norm_data, 
                                                as.vector(unique(biopsies_only_in_ANTvsCon_lncRNA$Gene.stable.ID)))
#colnames(Selected_data)
setwd(get.path("outcome"))
print.geneExpression(Selected_data , "ANT_DElncRNABiopsies_expression_NN_2020Apr16.pdf")

## print gene experssion in biopsies
### find the mean
H19_expression <- get.gene_expression(Biospies_data_log_withSampleInfo, selected_lncRNA)
Mean_tem <-get.average(H19_expression)

boxplot(ENSG00000130600~Biopsies_type, H19_expression,
        main = "boxplots for biopsies groups",
        names, las = 1, cex.axis = 0.7, outline = TRUE
)

boxplot(log_ENSG00000130600~Biopsies_type, H19_expression,
        main = "boxplots for biopsies groups (log value)",
        names, las = 1, cex.axis = 0.7, outline = TRUE
)

get.gene_expression <- function (data, selected_gene) {
  output <- as.data.frame(data[,c("Biopsies_type", selected_gene)])
  output[, paste0("log_", selected_gene)] <- log(output[,selected_gene], 2)
  return(output)
}

Biospies_data_log_withSampleInfo   <- merge(Sample_info, t(biopsies_norm_data), 
                                            by.x = "SID", by.y = "row.names")
Biospies_data_log_withSampleInfo   <- get.rownames(Biospies_data_log_withSampleInfo)

H19_expression <- get.gene_expression(Biospies_data_log_withSampleInfo, selected_lncRNA)


Selected_data     <- Select_if_in_specific_list(norm_data, "ENSG00000130600")

#code done until here
#######################################


picked_lncRNA <- "ENSG00000130600"
picked_lncRNA <- "ENSG00000245532"

cor_lncNRA <- cor(t(norm_data_v1), norm_data_v1[picked_lncRNA,])
correlated_RNA<- norm_data_v1[which(abs(cor_lncNRA)>0.8),]
outcome <- c()
for(i in rownames(correlated_RNA)) {
  outcome_tem<- cor.test(correlated_RNA[i, ], norm_data_v1[picked_lncRNA,])
  outcome <- rbind(outcome, c(outcome_tem$p.value, outcome_tem$estimate))
}
rownames(outcome) <- rownames(correlated_RNA)
colnames(outcome) <- c("p.value", "cor_estimate")
outcome <- as.data.frame(outcome)
write.csv(outcome, "outcome_H19_cutoff0.8_NN_2020Jan23.csv")
write.csv(outcome, "outcome_NEAT1_cutoff0.8_NN_2020Jan23.csv")


######
## for each samples:
DE_lncRNA <- list()
for(i in names(drug_ImpulseDE2Result)) {
  DE_genes_info <- Ensemble_database[which(Ensemble_database$Gene.stable.ID %in% drug_ImpulseDE2Result[[i]]$Gene),]
  DE_lncRNA_temp <- DE_genes_info[which(DE_genes_info$Gene.type == "lncRNA"),]
  DE_lncRNA[[i]] <- unique(DE_lncRNA_temp$Gene.stable.ID)
}

#write.table(DE_lncRNA$DOX_The, "test.txt", row.names = FALSE, col.names = FALSE)




## Fix code
###  =============================================== ImpulseDE2 run in ngs-cal:

# functions
get.RNA_data_combine <- function(list_dataset, data_type) {
  RNA_data_combine  <- merge(list_dataset[[1]][[data_type]], list_dataset[[2]][[data_type]], by = "row.names")
  RNA_data_combine  <- get.rownames(RNA_data_combine)
  if (length(list_dataset) > 2) {
    for (i in c(3:length(list_dataset))) {
      RNA_data_combine  <- merge(RNA_data_combine, list_dataset[[i]][[data_type]], by = "row.names")
      RNA_data_combine  <- get.rownames(RNA_data_combine)
    }
  }
  
  return(RNA_data_combine)
}


get.rownames <- function(data) {
  rownames(data) <- data[, 1]
  data  <- data[, -1]
  return(data)
}

get.metadata_for_ImpulseDe2 <- function(Input_data) {
  metadata <- as.data.frame(matrix(NA, ncol = 3, nrow = length(colnames(Input_data))))
  colnames(metadata) <- c("Sample", "Condition", "Time")
  rownames(metadata) <- colnames(Input_data)
  
  metadata[, "Sample"]   <- rownames(metadata)
  for (i in 1:nrow(metadata)) {
    if (grepl("Con", rownames(metadata)[i])) {
      metadata[, "Condition"][i] <- "control"
    } else {
      metadata[, "Condition"][i] <- "case"
    }
  }
  
  metadata[, "Time"]     <- as.integer(substr(rownames(metadata), 9, 11))
  return(metadata)
}

get.ImpulseDE2_outcome <- function(Input_data, metadata) {
  library(ImpulseDE2)
  raw_data <- matrix(as.integer(as.matrix(Input_data)), ncol = ncol(Input_data))
  colnames(raw_data) <- colnames(Input_data)
  rownames(raw_data) <- rownames(Input_data)
  objectImpulseDE2 <- runImpulseDE2(matCountData = raw_data, 
                                    dfAnnotation = metadata,
                                    boolCaseCtrl = TRUE,
                                    vecConfounders = NULL,
                                    scaNProc = 1)
  return(objectImpulseDE2)
}

# code
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Con_Flu_DMSO_Gal/TotalRNA/Genes/Con_DF2data_Cardiac_NN_20191112.RData")
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Doxorubicin/TotalRNA/Genes/DOXdata_Cardiac_NN_20191112.RData")
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Epirubicin/TotalRNA/Genes/EPIdata_Cardiac_NN_20191112.RData")
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Idarubicin/TotalRNA/Genes/IDAdata_Cardiac_NN_20191112.RData")


RNA_data_count  <- get.RNA_data_combine(list(Con_DF2, DOX, EPI, IDA), "expected_count")
RNA_data_count  <- RNA_data_count[, -grep("Con_DF2_000_", colnames(RNA_data_count))]
RNA_data_count  <- RNA_data_count[, -grep("DOX_Tox_240_2", colnames(RNA_data_count))]
RNA_data_count  <- RNA_data_count[, -grep("IDA_The_240_3", colnames(RNA_data_count))]

drug_The_list <- c("DOX_The", "EPI_The", "IDA_The")
drug_outcome <- list()
for(drug in drug_The_list) {
  Control_data <- RNA_data_count[,grepl("Con_DF2", colnames(RNA_data_count))]
  Drug_data    <- RNA_data_count[,grepl(drug, colnames(RNA_data_count))]
  Input_data <- cbind(Drug_data,Control_data)
  
  metadata <-get.metadata_for_ImpulseDe2(Input_data)
  drug_outcome[[drug]] <- get.ImpulseDE2_outcome(Input_data, metadata)
}

drug_Tox_list <- c("DOX_Tox", "EPI_Tox", "IDA_Tox")
for(drug in drug_Tox_list) {
  Control_data <- RNA_data_count[,grepl("Con_DF2", colnames(RNA_data_count))]
  Drug_data    <- RNA_data_count[,grepl(drug, colnames(RNA_data_count))]
  Input_data <- cbind(Drug_data,Control_data)
  Input_data <- Input_data[, -grep("Con_DF2_336_", colnames(Input_data))]
  
  metadata <-get.metadata_for_ImpulseDe2(Input_data)
  drug_outcome[[drug]] <- get.ImpulseDE2_outcome(Input_data, metadata)
}


save(drug_outcome, file = "drug_outcome_NN_2020April02.RData")








###  =============================================== ImpulseDE2 run in ngs-cal:

# functions
get.RNA_data_combine <- function(list_dataset, data_type) {
  RNA_data_combine  <- merge(list_dataset[[1]][[data_type]], list_dataset[[2]][[data_type]], by = "row.names")
  RNA_data_combine  <- get.rownames(RNA_data_combine)
  if (length(list_dataset) > 2) {
    for (i in c(3:length(list_dataset))) {
      RNA_data_combine  <- merge(RNA_data_combine, list_dataset[[i]][[data_type]], by = "row.names")
      RNA_data_combine  <- get.rownames(RNA_data_combine)
    }
  }
  
  return(RNA_data_combine)
}


get.rownames <- function(data) {
  rownames(data) <- data[, 1]
  data  <- data[, -1]
  return(data)
}

get.metadata_for_ImpulseDe2 <- function(Input_data) {
  metadata <- as.data.frame(matrix(NA, ncol = 3, nrow = length(colnames(Input_data))))
  colnames(metadata) <- c("Sample", "Condition", "Time")
  rownames(metadata) <- colnames(Input_data)
  
  metadata[, "Sample"]   <- rownames(metadata)
  metadata[, "Condition"]<- rep(c("case","control"), each = (length(colnames(Input_data))/2))
  metadata[, "Time"]     <- as.integer(substr(rownames(metadata), 9, 11))
  return(metadata)
}

get.ImpulseDE2_outcome <- function(Input_data, metadata) {
  library(ImpulseDE2)
  raw_data <- matrix(as.integer(as.matrix(Input_data)), ncol = ncol(Input_data))
  colnames(raw_data) <- colnames(Input_data)
  rownames(raw_data) <- rownames(Input_data)
  objectImpulseDE2 <- runImpulseDE2(matCountData = raw_data, 
                                    dfAnnotation = metadata,
                                    boolCaseCtrl = TRUE,
                                    vecConfounders = NULL,
                                    scaNProc = 1)
  return(objectImpulseDE2)
}

# code
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Con_Flu_DMSO_Gal/TotalRNA/Genes/Con_DF2data_Cardiac_NN_20191112.RData")
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Doxorubicin/TotalRNA/Genes/DOXdata_Cardiac_NN_20191112.RData")
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Epirubicin/TotalRNA/Genes/EPIdata_Cardiac_NN_20191112.RData")
load("/ngs-data/analysis/hecatos/NhanNguyen/RNAseq/Cardiac/Idarubicin/TotalRNA/Genes/IDAdata_Cardiac_NN_20191112.RData")


RNA_data_count  <- get.RNA_data_combine(list(Con_DF2, DOX, EPI, IDA), "expected_count")
RNA_data_count  <- RNA_data_count[, -grep("Con_DF2_000_", colnames(RNA_data_count))]
RNA_data_count  <- RNA_data_count[, -grep("DOX_Tox_240_2", colnames(RNA_data_count))]
RNA_data_count  <- RNA_data_count[, -grep("IDA_The_240_3", colnames(RNA_data_count))]

drug_The_list <- c("DOX_The", "EPI_The", "IDA_The")
drug_outcome <- list()
for(drug in drug_The_list) {
  Control_data <- RNA_data_count[,grepl("Con_DF2", colnames(RNA_data_count))]
  Drug_data    <- RNA_data_count[,grepl(drug, colnames(RNA_data_count))]
  Input_data <- cbind(Drug_data,Control_data)
  
  metadata <-get.metadata_for_ImpulseDe2(Input_data)
  drug_outcome[[drug]] <- get.ImpulseDE2_outcome(Input_data, metadata)
}

drug_Tox_list <- c("DOX_Tox", "EPI_Tox", "IDA_Tox")
for(drug in drug_Tox_list) {
  Control_data <- RNA_data_count[,grepl("Con_DF2", colnames(RNA_data_count))]
  Drug_data    <- RNA_data_count[,grepl(drug, colnames(RNA_data_count))]
  Input_data <- cbind(Drug_data,Control_data)
  Input_data <- Input_data[, -grep("Con_DF2_336_", colnames(Input_data))]
  
  metadata <-get.metadata_for_ImpulseDe2(Input_data)
  drug_outcome[[drug]] <- get.ImpulseDE2_outcome(Input_data, metadata)
}


save(drug_outcome, file = "drug_outcome_NN_2020Mar13.RData")
