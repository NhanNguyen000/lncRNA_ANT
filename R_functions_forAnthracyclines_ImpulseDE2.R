
Loaded_data <- function(folder, file_name) {
  setwd(folder)
  output<-read.csv(file_name) 
  return(output)
}

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


get.total_reads_plot <- function(RNA_data_count) {
  RNA_total_read_count <- c()
  for(sample in colnames(RNA_data_count)) {
    RNA_total_read_count[sample] <-sum(RNA_data_count[[sample]])
  }
  
  barplot(RNA_total_read_count, 
          main = "Total RNA read count", 
          xlab = "Time point", ylab = "Total read count", col=c("grey", "darkblue","red"),
          legend = rownames(RNA_total_read_count), beside = TRUE)
  abline(h=5e+6)
  
  return(RNA_total_read_count)
}

get.norm_data <- function(RNA_data_count) {
  metadata <- generate_metadata_for_sample(colnames(RNA_data_count))
  library("DESeq2")
  dds <- DESeqDataSetFromMatrix(countData = round(RNA_data_count), 
                                colData = metadata, design = ~ Condition)
  dds2 <- estimateSizeFactors(dds)
  norm_data <- counts(dds2,normalized=TRUE)
  return(norm_data)
}

generate_metadata_for_sample <- function(list_data) {
  output <- matrix(data=NA, ncol=4, nrow= length(list_data))
  colnames(output) <- c("Compound", "Dose", "Time", "Condition")
  rownames(output) <- list_data
  output[, "Compound"] <- Select_Part_Names(list = list_data, separator = "_", selected_position = 1)
  output[, "Dose"] <- Select_Part_Names(list = list_data, separator = "_", selected_position = 2)
  output[, "Time"] <- Select_Part_Names(list = list_data, separator = "_", selected_position = 3)
  output[, "Condition"] <- paste0(output[, "Compound"], "_", output[, "Dose"], "_", output[, "Time"])
  return(output)
}

Select_Part_Names <- function(list, separator, selected_position) {
  # eg. for separator = "|" --> "[|]"
  output<-c()
  for (i in 1: length(list)) output[i] <- strsplit(list[i], paste0("[", separator, "]"))[[1]][selected_position]
  return(output)
}

get.remove_gene_low_express <- function(data, cutoff) {
  output <- data[rowMeans(data> cutoff) >0, ]
  return(output)
}


get.summary_database <- function(data_before_filter, data_after_filter) {
  database_before_filter <- get.selected_data(rownames(data_before_filter), Ensemble_database, "Gene.stable.ID")
  database_after_filter <- get.selected_data(rownames(data_after_filter), Ensemble_database, "Gene.stable.ID")
  
  sum_database_before_filter <- get.summary_gene_type(database_before_filter)
  sum_database_after_filter  <- get.summary_gene_type(database_after_filter)
  
  summary <- merge(sum_database_before_filter, sum_database_after_filter, 
                   by = "row.names", all.x = TRUE, all.y = TRUE)
  summary <- get.rownames(summary)
  summary[is.na(summary)] <- 0
  summary[,3] <- summary[,1] - summary[,2]
  
  colnames(summary) <- c("total number", "number of gene expression", "number of gene no expression")
  return(summary)
}

get.selected_data <- function(selected_criteria, database, column_in_database) {
  if (column_in_database == "rownames") output <- database[rownames(database) %in% selected_criteria,]
  if (column_in_database != "rownames") output <- database[database[,column_in_database] %in% selected_criteria,]
  return(output)
}

get.summary_gene_type <- function(database) {
  gene_type_list <- as.vector(unique(database$Gene.type))
  number_of_genes <- c()
  for (i in 1:length(gene_type_list)) {
    data_tem        <- get.selected_data(gene_type_list[i], database, column_in_database = "Gene.type")
    number_of_genes[i] <- length(unique(data_tem$Gene.stable.ID))
  }
  summary <- as.matrix(number_of_genes)
  rownames(summary) <- gene_type_list
  return(summary)
}


get.log <- function(data, log_base) {
  data[data == 0]  <- 1
  output           <- log(data, log_base)
  return(output)
}

get.pie_chart <- function(data, names) {
  library(viridis)
  Pie_chart_info_percent  <- round(100*data/sum(data), 2)
  labls <- paste0(names, " (", Pie_chart_info_percent, "%",")")
  labls[data < sort(data, decreasing = TRUE)[4]] <- ""
  pie(Pie_chart_info_percent, labels = labls, 
      init.angle = 45, col = viridis_pal()(length(data)))
}


get.list_Element <- function(list, selected_Element) {
  output <- sapply(list, function(x) x[selected_Element])
  return(output)
}


get.vennDiam <- function(drug_ImpulseDE2Result, selected_list) {
  library("venn")
  List_temp  <-drug_ImpulseDE2Result[selected_list]
  Input_venn <- get.list_Element(List_temp, "Gene")
  names(Input_venn) <- get.list_Element(strsplit(names(Input_venn), "[.]"), 1)
  
  Venn_outcome <- attr(venn(Input_venn), "intersections")
  return(Venn_outcome)
}

get.geneTypes_overlapLncRNA <- function(Ensemble_database, selected_list) {
  Overlap_genes_info    <- Ensemble_database[which(Ensemble_database$Gene.stable.ID %in% selected_list),]
  database_sum          <- get.summary_gene_type(Overlap_genes_info)
  colnames(database_sum)<- "NumberOfGene"
  
  Overlap_lncRNA     <- Overlap_genes_info[which(Overlap_genes_info$Gene.type == "lncRNA"),]
  
  outcome_list <- list("Sum_GeneType" = database_sum, "Overlap_lncRNA" = Overlap_lncRNA)
  return(outcome_list)
}

Select_if_in_specific_list <- function (data, selected_list) {
  output <- data[rownames(data) %in% selected_list,]
  return(output)
}

print.geneExpression <- function(Selected_data, file_names) {
  library("plyr")						
  library("ggplot2")		
  library("gridExtra")
  library("magrittr")
  library("dplyr")
  
  Log_Selected_data <- get.log(Selected_data, 2)
  
  metadata <- data.frame(generate_metadata_for_sample(colnames(Selected_data)), stringsAsFactors=FALSE)
  metadata$Dose <- paste0(metadata$Compound, "_", metadata$Dose)
  
  pdf(file_names, onefile = TRUE, width = 15)
  
  for(i in rownames(Selected_data)) {
    expression_tem  <- as.data.frame(cbind(metadata, Selected_data[i,]))
    
    Mean_tem        <- ddply(expression_tem, ~Dose+Time, summarise, 
                             Mean = mean(`Selected_data[i, ]`, na.rm = TRUE), sd = sd(`Selected_data[i, ]`, na.rm = TRUE))
    
    Log_expression_tem <- as.data.frame(cbind(metadata, Log_Selected_data[i,]))
    Log_Mean_tem       <- ddply(Log_expression_tem, ~Dose+Time, summarise, 
                                Mean = mean(`Log_Selected_data[i, ]`, na.rm = TRUE), sd = sd(`Log_Selected_data[i, ]`, na.rm = TRUE))
    
    plist <- list()
    plist[[1]]<- ggplot(Mean_tem, aes(x = Time, y = Mean, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
      xlab("Time (hours)") + ylab("Expression value") + ggtitle(i) +  theme_bw()
    
    plist[[2]]<- ggplot(Log_Mean_tem, aes(x = Time, y = Mean, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
      xlab("Time (hours)") + ylab("Log2 expression value") + ggtitle(i) +  theme_bw()
    ml <- marrangeGrob(plist, nrow=2, ncol=2)
    print(ml)
  }
  
  dev.off()
}

make.order_sample_byIDs <- function(data, sample_type, id_type) {
  output <- data[which(data$Biopsies_type == sample_type),]
  output <- output[order(as.numeric(output[,id_type])),]
  return(output)
}

get.Log2FC <- function(Selected_lncRNA, norm_data_invitro, norm_data_biopsies, metadata_biopsies) {
  library("plyr")						
  library("ggplot2")		
  library("gridExtra")
  library("magrittr")
  library("dplyr")
  library("ggpubr")
  
  # in vitro
  Selected_data_invitro     <- Select_if_in_specific_list(norm_data_invitro, Selected_lncRNA)
  Log_invitro <- get.log(Selected_data_invitro, 2)
  metadata <- data.frame(generate_metadata_for_sample(colnames(Log_invitro)), stringsAsFactors=FALSE)
  metadata$Dose <- paste0(metadata$Compound, "_", metadata$Dose)
  
  # biopsies
  Selected_data_biopsies <- Select_if_in_specific_list(norm_data_biopsies, Selected_lncRNA)
  Log_biopsies <- get.log(Selected_data_biopsies, 2)
  
  # plot
  plist <- list()
  for(i in Selected_lncRNA) {
    # invitro
    Log_invitro_tem  <- as.data.frame(cbind(metadata, Log_invitro[i,]))
    Log_invitro_mean <- ddply(Log_invitro_tem, ~Dose+Time, summarise,
                              Mean = mean(`Log_invitro[i, ]`, na.rm = TRUE))
    
    Control_value <- rbind(Log_invitro_mean[c(1:7),], Log_invitro_mean[c(1:6),],
                           Log_invitro_mean[c(1:7),], Log_invitro_mean[c(1:6),],
                           Log_invitro_mean[c(1:7),], Log_invitro_mean[c(1:6),])
    Log2FC_invitro <- Log_invitro_mean %>% filter(Dose != "Con_DF2")
    Log2FC_invitro$Log2FC <- Log2FC_invitro$Mean - Control_value$Mean
    
    # biopsies
    Log_biopsies_tem  <- merge(metadata_biopsies, Log_biopsies[i,], by = "row.names")
    biopsy_controls <- make.order_sample_byIDs(Log_biopsies_tem, "Control", "SID")
    biopsy_ANTs <- make.order_sample_byIDs(Log_biopsies_tem, "LateCardiotoxicity_with_ANT", "Matched.to.SID")
    Log2FC <- as.data.frame(biopsy_ANTs[,8] - biopsy_controls[,8])
    names(Log2FC) <- "Log2FC"
    
    # arange plot
    gene_name <- unique(Ensemble_database$Gene.name[which(Ensemble_database$Gene.stable.ID==i)])
    
    plot_invitro <- ggplot(Log2FC_invitro, aes(x = Time, y = Log2FC, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + xlab("Time (hours)") + ylab("Log2FC of read counts") +
      geom_hline(yintercept=0, linetype="dashed", color = "black") + 
      ggtitle(gene_name) +  theme_bw()
    
    plot_biopsies <- ggplot(Log2FC, aes(x= "", y=Log2FC)) +  geom_boxplot() + ggtitle(gene_name) +
      theme_bw() + geom_hline(yintercept=0, linetype="dashed", color = "black")
    
    plist[[i]] <- ggarrange(plot_invitro, plot_biopsies, widths = c(8,1.5))
  }
  ml <- marrangeGrob(plist, nrow=4, ncol=4)
  pdf("Figure4_Log2FC.pdf", width=25,height=12)
  print(ml)
  dev.off()
}

get.geneName <- function(list, Ensemble) {
  output <- unique(Ensemble$Gene.name[which(Ensemble$Gene.stable.ID %in% list)])
  return(output)
}

get.biopsies_type <- function(Biospies_info) {
  Samples <- c()
  Samples$Patient[Biospies_info$Type == "Late onset cardiotoxicity"] <- "LateCardiotoxicity"
  Samples$Patient[Biospies_info$Type == "Acute onset cardiotoxicity"]<- "AcuteCardiotoxicity"
  Samples$Patient[Biospies_info$Type == "Control patients"]          <- "Control"
  
  Samples$Drug[Biospies_info$Treat.with == ""]           <- ""
  Samples$Drug[Biospies_info$Treat.with == "ANT"]        <- "_with_ANT"
  #Samples$Drug[which(grepl("rubicin", Biospies_info$Chemotherapeutic.agents))]        <- "_with_ANT"
  #Samples$Drug[which(grepl("Anthracycline", Biospies_info$Chemotherapeutic.agents))]  <- "_with_ANT"
  
  
  summary_sample           <- matrix(unlist(Samples), ncol =2, byrow = FALSE)
  output <- apply(summary_sample, 1, paste , collapse = "" )
  return(output)
}

get.lncRNA <- function(data) {
  gene_information <- merge(Ensemble_database, data, 
                            by.x = "Gene.stable.ID", by.y = "row.names")
  
  lncRNA_data <- get.selected_data("lncRNA", gene_information, "Gene.type")
  output <- unique(lncRNA_data$Gene.stable.ID)
  return(output)
}

get.pca <- function(data, condition, name) {
  library(factoextra)
  data[is.na(data)] <- 0
  res.pca <- prcomp(data, scale = TRUE)
  fviz_pca_ind(res.pca, geom="point",  pointsize = 2, 
               habillage=condition, addEllipses=TRUE, ellipse.level=0.95,  
               title = paste0("PCA - ", name))
}
