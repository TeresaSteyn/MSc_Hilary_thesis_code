#Make Venn diagram of sig DEGs in GFAP vs ALDH1L astrocyte data (subsetted before processing)

#Load libraries
library(ggplot2)
library(dplyr)
library(VennDiagram)

path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")

pathResults1 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_GFAP\\")

pathResults2 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_ALDH1L1\\")

GFAP_res <- read.csv(paste0(pathResults1, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"))  

ALDH1L1_res <- read.csv(paste0(pathResults2, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"))  

##Add .2 and .3 to the second occurrence or third occurence of a gene name
GFAP_res <- GFAP_res %>%
  group_by(Gene) %>%
  mutate(
    occurrence = row_number(), # Assign an occurrence number to each gene within its group
    Gene = case_when(
      occurrence == 2 ~ paste0(Gene, ".2"), # Append ".2" to the second occurrence
      occurrence == 3 ~ paste0(Gene, ".3"), # Append ".3" to the third occurrence
      TRUE ~ Gene # Keep the gene name as is for all other cases
    )
  ) %>%
  select(-occurrence) %>% # Remove the helper 'occurrence' column
  ungroup() # Remove the grouping

ALDH1L1_res <- ALDH1L1_res %>%
  group_by(Gene) %>%
  mutate(
    occurrence = row_number(), # Assign an occurrence number to each gene within its group
    Gene = case_when(
      occurrence == 2 ~ paste0(Gene, ".2"), # Append ".2" to the second occurrence
      occurrence == 3 ~ paste0(Gene, ".3"), # Append ".3" to the third occurrence
      TRUE ~ Gene # Keep the gene name as is for all other cases
    )
  ) %>%
  select(-occurrence) %>% # Remove the helper 'occurrence' column
  ungroup() # Remove the grouping



# Assuming GFAP_res and ALDH1L1_res have a column named 'Gene'
# Extract the lists of Genes from each dataframe
Genes_GFAP <- unique(GFAP_res$Gene)
Genes_ALDH1L1 <- unique(ALDH1L1_res$Gene)


# Create a list of these Gene sets
Gene_lists <- list(GFAP = Genes_GFAP, ALDH1L1 = Genes_ALDH1L1)

# Generate the Venn diagram
venn.plot <- venn.diagram(
  x = Gene_lists,  width = 3500,
  category.names = c("GFAP", "ALDH1L1"),
  output = TRUE,  # To plot directly to the R plotting window
  filename = paste0(pathResults2, "GFAP_ALDH1L1_sig_genes_venn_diagram.png"), # Use NULL to plot interactively, or specify a filename to save the output
  imagetype = "png", # Specify the type of image if you're saving the file (e.g., "png", "pdf")
  fill = c("cornflowerblue", "darkorange"), # Colors for the sets
  alpha = 0.5, # Transparency of colors
  label.col = "black", # Color of labels
  cex = 1.5, # Scaling factor for text labels
  fontfamily = "sans", # Font family for text
  cat.col = c("cornflowerblue", "darkorange"), # Color for category names
  cat.cex = 1.5, # Scaling factor for category names
  cat.fontfamily = "sans" # Font family for category names
)

