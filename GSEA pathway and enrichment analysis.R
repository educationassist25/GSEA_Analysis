# Install and load necessary packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
install.packages("tidyverse")
BiocManager::install("airway")
BiocManager::install("ReactomePA",force = TRUE)
BiocManager::install("KEGGREST", force = TRUE)
BiocManager::install("org.Hs.eg.db", force = TRUE)
BiocManager::install("msigdbr", force = TRUE)
install.packages("ggalt")
.rs.restartR()
install.packages("gtable")
rm(list = ls())


if (!requireNamespace("openxlsx", quietly = TRUE)) 
  install.packages("openxlsx")

install.packages("magick")

BiocManager::install("ComplexHeatmap")
install.packages("circlize")

if (!requireNamespace("enrichplot", quietly = TRUE)) {
  install.packages("enrichplot")
}



# Load necessary libraries
library(DESeq2)
library(tidyverse)
library(airway)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ReactomePA)
library(KEGGREST)
library(readr)
library(ggalt)
library(openxlsx)
library(stats)
library(ggrepel)
library(RColorBrewer)
library(magick)
library(ComplexHeatmap)
library(circlize)
library(enrichplot)
library(readxl)


# differential analysis (results of P value, FDR, log2Fold changes_KD over WT)

result <- read.csv('Example_Data_GSEA.csv', row.names = 1)


View(result)

# Separate the data and group information
data_matrix <- as.matrix(result[, -ncol(result)])
group <- result$condition

# Log2 transform data first
data_matrix <- log2(data_matrix + 1)

view(data_matrix)

# Calculate log2 fold changes
log2FC <- apply(data_matrix, 2, function(x) {
  (mean(x[group == "KD"]) - mean(x[group == "WT"]))
})

linearFC <- 2^log2FC

# Perform T-tests
p_values <- apply(data_matrix, 2, function(x) {
  t.test(x[group == "WT"], x[group == "KD"])$p.value
})

# Calculate FDR
fdr <- p.adjust(p_values, method = "fdr")


# Combine results into a data frame
results <- data.frame(p_value = p_values,
                      FDR = fdr,
                      log2FC = log2FC,
                      linearFC = linearFC)


View(results)

# Save the results to an Excel file and CSV file
write.xlsx(results, file = "Report File_KD over WT.xlsx", rowNames = TRUE)
write.csv(results, file = "Report File_KD over WT.csv", row.names = TRUE)



# -------------GSEA analysis (FDR <0.25)----------------------------------------

data_3 <- read.csv("Report File_KD over WT.csv", row.names = 1)

view(data_3)

# -------- FILTER BY FDR --------
data_4 <- data_3[data_3$FDR < 0.25, ] # if you want filter with p value 0.05, please use instead FDR 0.25


View(data_4)
 str(data_4)

gene_list <- data_4$log2FC

view(gene_list)

str(gene_list)
names(gene_list) <- rownames(data_4)  # Assuming rownames contain gene names
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene names to Entrez IDs
# Check the mapping
print(head(gene_list))  # Print the first few entries to verify


#-------------GO Pathway Analysis----------------------------------------------

# Perform GO GSEA analysis
gsea_go <- gseGO(
  geneList = gene_list,
  ont = "ALL",
  keyType = "SYMBOL",  # Adjust this based on geneList names
  nPerm = 10000,
  minGSSize = 3,
  maxGSSize = 100,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  OrgDb = org.Hs.eg.db,
  pAdjustMethod = "none"  # Consider using 'BH' or 'fdr' for multiple testing correction
)

# Save results as CSV
write.csv(as.data.frame(gsea_go), 'GSEA_GO_results.csv')

# Create and save the dot plot_top 10 pathway
gsea_go_plot <- dotplot(gsea_go, showCategory = 10) + ggtitle("GSEA GO Dotplot")

show(gsea_go_plot)

save_gsea_plots <- function(plot, name) {
  formats <- c("jpg")
  for (fmt in formats) {
    ggsave(paste0(name, ".", fmt), plot = plot, device = fmt)
  }
}
save_gsea_plots(gsea_go_plot, "GSEA_GO_plot")


#-------GSEA enrichment plot---------

# Function to sanitize filenames by replacing invalid characters
sanitize_filename <- function(name) {
  gsub("[^[:alnum:]_]", "_", name)  # Replace any non-alphanumeric or underscore character with an underscore
}

# Function to generate and customize GSEA plots using gseaplot2
generate_specific_gsea_plot2 <- function(gsea_result, pathway_id, save_prefix) {
  sanitized_id <- sanitize_filename(pathway_id)  # Sanitize the pathway ID for filename
  
  # Check if the pathway_id exists in the result
  if (pathway_id %in% rownames(gsea_result@result)) {
    gsea_plot <- gseaplot2(
      gsea_result,
      geneSetID = pathway_id,
      title = gsea_result@result[pathway_id, "Description"],
      pvalue_table = TRUE  # Include p-value table for context
    )
    
    # Display and save the plot
    print(gsea_plot)
    save_gsea_plots(gsea_plot, paste0(save_prefix, "_", sanitized_id))
  } else {
    warning(paste("Pathway ID", pathway_id, "not found in the GSEA results."))
  }
}


# Extract the top 2 pathways based on p-value
top_2_pathways <- head(gsea_go@result$ID, 2)

# Generate and save plots for the top 10 pathways
for (pathway_id in top_2_pathways) {
  generate_specific_gsea_plot2(gsea_go, pathway_id, "GO_GSEA_plot")
}



#---------Hallmark Pathway Analysis---------------------------------------------


# Obtain hallmark gene sets from MSigDB
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")  # 'H' for hallmark pathways

# Convert to a format suitable for GSEA
hallmark_gsea <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],  # Use SYMBOL as key
  pvalueCutoff = 0.05,
  nPerm = 10000,
  minGSSize = 3,
  maxGSSize = 500,  # Adjust as necessary
  verbose = TRUE,
  pAdjustMethod = "none"  # Consider 'BH' for multiple testing correction
)

# Save GSEA results
write.csv(as.data.frame(hallmark_gsea), 'GSEA_Hallmark_results.csv')

# Create and save a dot plot for visualization
hallmark_gsea_plot <- dotplot(hallmark_gsea, showCategory = 10) + ggtitle("GSEA Hallmark Pathways Dotplot")

# Save GSEA plots
save_gsea_plots <- function(plot, name) {
  formats <- c("jpg")
  for (fmt in formats) {
    ggsave(paste0(name, ".", fmt), plot = plot, device = fmt)
  }
}
save_gsea_plots(hallmark_gsea_plot, "GSEA_Hallmark_plot")

# Generate GSEA plots for the top 2 pathways
generate_specific_gsea_plot <- function(hallmark_gsea, pathway_id, save_prefix) {
  gsea_plot <- gseaplot2(hallmark_gsea, geneSetID = pathway_id, title = hallmark_gsea@result[pathway_id, "Description"])
  save_gsea_plots(gsea_plot, save_prefix)
}

# Extract the top 2 pathways based on p-value
top_2_pathways <- head(hallmark_gsea@result$ID, 2)

# Generate and save plots for the top 10 pathways
for (pathway_id in top_2_pathways) {
  generate_specific_gsea_plot(hallmark_gsea, pathway_id, paste0("Hallmark_GSEA_plot_", pathway_id))
}



#---------- Reactome Pathway Analysis-------------------------------------------


# Obtain Reactome gene sets from MSigDB
reactome_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")

# Convert to a format suitable for GSEA
reactome_gsea <- GSEA(
  geneList = gene_list,
  TERM2GENE = reactome_sets[, c("gs_name", "gene_symbol")],  # Use SYMBOL as key
  pvalueCutoff = 0.05,
  nPerm = 10000,
  minGSSize = 3,
  maxGSSize = 500,  # Adjust as necessary
  verbose = TRUE,
  pAdjustMethod = "none"  # Consider 'BH' for multiple testing correction
)

# Save GSEA results
write.csv(as.data.frame(reactome_gsea), 'GSEA_Reactome_results_FDR0.25.csv')

# Create and save a dot plot for visualization
reactome_gsea_plot <- dotplot(reactome_gsea, showCategory = 10) + ggtitle("GSEA Reactome Pathways Dotplot")

# Save GSEA plots
save_gsea_plots <- function(plot, name) {
  formats <- c("jpg")
  for (fmt in formats) {
    ggsave(paste0(name, ".", fmt), plot = plot, device = fmt)
  }
}
save_gsea_plots(reactome_gsea_plot, "GSEA_Reactome_plot")


#-------GSEA enrichment plot---------

# Create and save a dot plot for visualization
reactome_gsea_plot <- dotplot(reactome_gsea, showCategory = 10) + ggtitle("GSEA Reactome Pathways Dotplot")

# Save GSEA plots
save_gsea_plots <- function(plot, name) {
  formats <- c("jpg")
  for (fmt in formats) {
    ggsave(paste0(name, ".", fmt), plot = plot, device = fmt)
  }
}
save_gsea_plots(reactome_gsea_plot, "GSEA_Reactome_plot")


# Generate GSEA plots for the top 2 pathways
generate_specific_gsea_plot <- function(reactome_gsea, pathway_id, save_prefix) {
  gsea_plot <- gseaplot2(reactome_gsea, geneSetID = pathway_id, title = reactome_gsea@result[pathway_id, "Description"])
  save_gsea_plots(gsea_plot, save_prefix)
}

# Extract the top 2 pathways based on p-value
top_2_pathways <- head(reactome_gsea@result$ID, 2)

# Generate and save plots for the top 10 pathways
for (pathway_id in top_2_pathways) {
  generate_specific_gsea_plot(reactome_gsea, pathway_id, paste0("Reactome_GSEA_plot_", pathway_id))
}

#-------------- KEGG Pathway Analysis-------------------------------------------

# Convert gene symbols to Entrez IDs
gene_list_entrez <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list <- gene_list[gene_list_entrez$SYMBOL]
names(gene_list) <- gene_list_entrez$ENTREZID

# Perform KEGG GSEA analysis
kegg_gsea <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",  # 'hsa' for Homo sapiens
  nPerm = 10000,
  minGSSize = 3,
  maxGSSize = 500,  # Adjust as necessary
  pvalueCutoff = 0.05,
  verbose = TRUE,
  pAdjustMethod = "none"  # Consider 'BH' for multiple testing correction
)

# Convert core enrichment Entrez IDs to gene symbols
convert_entrez_to_symbol <- function(entrez_ids) {
  symbols <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  return(symbols$SYMBOL)
}

# Add gene symbols to core enrichment column
kegg_gsea@result$core_enrichment <- sapply(
  strsplit(kegg_gsea@result$core_enrichment, "/"),
  function(genes) {
    paste(convert_entrez_to_symbol(genes), collapse = "/")
  }
)

# Save GSEA results with gene symbols
write.csv(as.data.frame(kegg_gsea), 'GSEA_KEGG_results_FDR0.25.csv')

# Create and save a dot plot for visualization
kegg_gsea_plot <- dotplot(kegg_gsea, showCategory = 10) + ggtitle("GSEA KEGG Pathways Dotplot")

# Save GSEA plots
save_gsea_plots <- function(plot, name) {
  formats <- c("jpg")
  for (fmt in formats) {
    ggsave(paste0(name, ".", fmt), plot = plot, device = fmt)
  }
}
save_gsea_plots(kegg_gsea_plot, "GSEA_KEGG_plot")

#-------GSEA enrichment plot---------

# Generate GSEA plots for the top 2 pathways
generate_specific_gsea_plot <- function(kegg_gsea, pathway_id, save_prefix) {
  gsea_plot <- gseaplot2(kegg_gsea, geneSetID = pathway_id, title = kegg_gsea@result[pathway_id, "Description"])
  save_gsea_plots(gsea_plot, save_prefix)
}

# Extract the top 2 pathways based on p-value
top_2_pathways <- head(kegg_gsea@result$ID, 2)

# Generate and save plots for the top 10 pathways
for (pathway_id in top_2_pathways) {
  generate_specific_gsea_plot(kegg_gsea, pathway_id, paste0("KEGG_GSEA_plot_", pathway_id))
}

