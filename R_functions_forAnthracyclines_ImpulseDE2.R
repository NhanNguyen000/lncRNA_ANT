
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

get.pie_chart <- function(data, title) {
  Pie_chart_info_percent  <- round(100*data/sum(data), 1)
  labls <- paste0(Pie_chart_info_percent, "% ", rownames(data))
  labls[data < sort(data, decreasing = TRUE)[4]] <- ""
  pie(Pie_chart_info_percent, labels = labls, main = title)
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

print.Log2geneExpression <- function(Selected_data) {
  library("plyr")						
  library("ggplot2")		
  library("gridExtra")
  library("magrittr")
  library("dplyr")
  
  Log_Selected_data <- get.log(Selected_data, 2)
  
  metadata <- data.frame(generate_metadata_for_sample(colnames(Selected_data)), stringsAsFactors=FALSE)
  metadata$Dose <- paste0(metadata$Compound, "_", metadata$Dose)
  
  plist <- list()
  for(i in rownames(Selected_data)) {
    expression_tem  <- as.data.frame(cbind(metadata, Selected_data[i,]))
    
    Mean_tem        <- ddply(expression_tem, ~Dose+Time, summarise, 
                             Mean = mean(`Selected_data[i, ]`, na.rm = TRUE), sd = sd(`Selected_data[i, ]`, na.rm = TRUE))
    
    Log_expression_tem <- as.data.frame(cbind(metadata, Log_Selected_data[i,]))
    Log_Mean_tem       <- ddply(Log_expression_tem, ~Dose+Time, summarise, 
                                Mean = mean(`Log_Selected_data[i, ]`, na.rm = TRUE), sd = sd(`Log_Selected_data[i, ]`, na.rm = TRUE))
    plist[[i]]<- ggplot(Log_Mean_tem, aes(x = Time, y = Mean, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
      xlab("Time (hours)") + ylab("Log2 expression value") + ggtitle(i) +  theme_bw()
  }
  ml <- marrangeGrob(plist, nrow=nrow(Selected_data)/2, ncol=2)
  print(ml)
}


get.lncRNA_corrRNA_his <- function(data, lncRNA, samples = NULL) {
  cor_lncRNA <- cor(t(data), data[lncRNA,])
  hist(cor_lncRNA, main= paste0("Histogram for ", lncRNA, " ", samples), 
       xlab="Correlation", ylab = "Number of RNA",
       border="blue", col="green", xlim=c(-1, 1), las=1)
}

print.lncRNA_corrRNA_condition <- function(Selected_list, norm_data, file_names) {
  library("gridExtra")
  library("magrittr")
  
  pdf(file_names, onefile = TRUE, width = 8)
  
  for(lncRNA in Selected_list) {
    get.lncRNA_corrRNA_his(norm_data, lncRNA, "totalANT")
    for(drug in c("DOX", "EPI", "IDA")) {
      norm_data_temp <- cbind(norm_data[,grep("Con", colnames(norm_data))],
                              norm_data[,grep(drug, colnames(norm_data))])
      get.lncRNA_corrRNA_his(norm_data_temp, lncRNA, paste0("in ", drug, "-Con"))
    }
  }
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

get.cor_lncRNA_phenotype <- function(heart_failure_items, Biopsies_lncRNA_data) {
  output <- c()
  for (i in colnames(heart_failure_items)) {
    for(lncRNA in rownames(Biopsies_lncRNA_data)) {
      cor_temp <- cor(as.numeric(Biopsies_lncRNA_data[lncRNA,]),
                      as.numeric(heart_failure_items[,i]), use=p)
      if (is.na(cor_temp) != TRUE & abs(cor_temp )>0.6) {
        output<- rbind(output, c(i, lncRNA, cor_temp))
        #print(i)
        #print(lncRNA)
        #print(cor_temp)
      }
    }
  }
  return(output)
}

get.pca <- function(data, condition, name) {
  data[is.na(data)] <- 0
  res.pca <- prcomp(data, scale = TRUE)
  fviz_pca_ind(res.pca, geom="point",  pointsize = 2, 
               habillage=condition, addEllipses=TRUE, ellipse.level=0.95,  
               title = paste0("PCA - ", name))
}
