library(dplyr)
library(tidyr)

path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")

pathResults1 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_GFAP\\")

pathResults2 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_ALDH1L1\\")

GFAP_res <- read.csv(paste0(pathResults1, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"))  

ALDH1L1_res <- read.csv(paste0(pathResults2, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"))  


pathResults3 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_GFAP\\Enrichment\\")


pathResults4 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_ALDH1L1\\Enrichment\\")

GFAP_enriched_pathways <- read.delim(paste0(pathResults3, "GFAP_enrichment.all2.tsv"))
colnames(GFAP_enriched_pathways)[8] <- "protein_IDs"
colnames(GFAP_enriched_pathways)[9] <- "protein_labels"

#Subset to only include particular databases
GFAP_enriched_pathways <- subset(GFAP_enriched_pathways, X.category %in% c("GO Function", "GO Process", "KEGG", "Reactome", "STRING clusters"))

df1 <- GFAP_enriched_pathways

ALDH1L1_enriched_pathways <- read.delim(paste0(pathResults4, "ALDH1L1_enrichment.all2.tsv"))
colnames(ALDH1L1_enriched_pathways)[8] <- "protein_IDs"
colnames(ALDH1L1_enriched_pathways)[9] <- "protein_labels"

#Subset to only include particular databases
ALDH1L1_enriched_pathways <- subset(ALDH1L1_enriched_pathways, X.category %in% c("GO Function", "GO Process", "KEGG", "Reactome", "STRING clusters"))

df2 <- ALDH1L1_enriched_pathways

unique_GFAP <- df1$term.description[which(!df1$term.description %in% df2$term.description)]

unique_ALDH1L1 <- df2$term.description[which(!df2$term.description %in% df1$term.description)]

shared_GFAP_ALDH1L1 <- intersect(df1$term.description, df2$term.description)

Select_pathways <- shared_GFAP_ALDH1L1

# Initialize the list to store data frames for each pathway
Pathways <- vector('list', length(Select_pathways))

# Loop through each selected pathway
for (x in 1:length(Select_pathways)) {
  # Extract genes associated with the current pathway
  pathway_genes <- ALDH1L1_enriched_pathways$protein_labels[
    grep(Select_pathways[x], ALDH1L1_enriched_pathways$term.description, ignore.case=TRUE)
  ]

  # Concatenate all elements into one long string, then split into individual genes
  all_genes <- unlist(strsplit(paste(pathway_genes, collapse = ","), ","))

  # Remove whitespace and extract unique gene names
  unique_genes <- unique(trimws(all_genes))

  # Store the unique genes in the list
  Pathways[[x]] <- unique_genes
  
  # Print the pathway name and the dataframe of unique genes
  #print(Select_pathways[x])
  #print(data.frame(Gene = unique_genes))
}

# Find the maximum length among the gene lists
max_length <- max(sapply(Pathways, length))

# Make each list the same length and create a dataframe for each, then name the columns
for (x in 1:length(Pathways)) {
  current_length <- length(Pathways[[x]])
  # Fill the list to ensure it has the same length as the longest list
  if (current_length < max_length) {
    Pathways[[x]] <- c(Pathways[[x]], rep(NA, max_length - current_length))
  }
  # Convert the list to a dataframe and assign the pathway name as the column name
  Pathways[[x]] <- data.frame(setNames(Pathways[[x]], Select_pathways[x]))
}

# Combine all dataframes into a single dataframe
all_pathway_genes <- do.call(cbind, Pathways)

# If needed, here's a quick fix to replace NA with an empty string or any other placeholder
all_pathway_genes[is.na(all_pathway_genes)] <- ""

Select_pathways1 <- gsub(" ", "_", Select_pathways)
Select_pathways1 <- gsub(",", "", Select_pathways1)
Select_pathways1 <- gsub("/", "__", Select_pathways1)
colnames(all_pathway_genes) <- Select_pathways1


#Make heatmap 
ALDH1L1_res_expanded <- ALDH1L1_res %>%
  separate_rows(GeneNames, sep = ";") %>%
  rename(Genes = GeneNames) # Rename the column back to 'Genes' if desired

# Initialize the heatmap_matrix with appropriate dimensions and names
heatmap_matrix <- matrix(NA, nrow = ncol(all_pathway_genes), ncol = length(unique(ALDH1L1_res_expanded$Genes)), 
                         dimnames = list(colnames(all_pathway_genes), unique(ALDH1L1_res_expanded$Genes)))

# Step 1: Collect all unique genes listed in the Pathways dataframe
all_pathway_genes2 <- unique(unlist(all_pathway_genes))

# Step 2: Filter the heatmap_matrix to keep only columns for genes found in the pathways
heatmap_matrix <- heatmap_matrix[, colnames(heatmap_matrix) %in% all_pathway_genes2]


# Assuming colnames(Pathways) are used as row names in heatmap_matrix
pathway_names <- colnames(all_pathway_genes) # Correct if pathways are actually column names in Pathways

for (pathway_name in pathway_names) {
  pathway_genes <- unlist(all_pathway_genes[pathway_name])
  pathway_genes <- pathway_genes[!is.na(pathway_genes)] # Remove NAs
  
  for (gene in pathway_genes) {
    if (gene %in% ALDH1L1_res_expanded$Genes) {
      logFC_values <- ALDH1L1_res_expanded$logFC[ALDH1L1_res_expanded$Genes == gene]
      if (length(logFC_values) > 0) {
        # Ensure the gene is a column in heatmap_matrix, and pathway_name is a valid row name
        if(gene %in% colnames(heatmap_matrix) && pathway_name %in% rownames(heatmap_matrix)) {
          heatmap_matrix[pathway_name, gene] <- mean(logFC_values, na.rm = TRUE)
        }
      }
    }
  }
}

library(ComplexHeatmap)
library(circlize)
# Define color mapping for logFC: Blue to White to Red, with NA as white
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Simple Heatmap call without the problematic top_annotation
filename1 <- paste0(pathResults4, "Heatmap_all_shared_ALDH1L1_GFAP_pathways_selectdatabase2.pdf")

pdf(filename1, width = 25, height = 10)
Heatmap(heatmap_matrix, name = "logFC", col = col_fun, na_col = "white", 
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = TRUE, 
        heatmap_legend_param = list(title = "LogFC", title_position = "topcenter", 
                                    position = "bottomright", legend_direction = "vertical"))
dev.off()
