#Make Venn diagram of sig DEGs in GFAP vs ALDH1L astrocyte data

#Load libraries
library(ggplot2)
library(VennDiagram)

path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")

pathResults <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM2\\")

GFAP_res <- read.csv(paste0(pathResults, "sig_res_preadjusted_GFAP_strict.csv"))  

ALDH1L1_res <- read.csv(paste0(pathResults, "sig_res_preadjusted_ALDH1L1_strict.csv"))  

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
  output = NULL,  # To plot directly to the R plotting window
  filename = paste0(pathResults, "GFAP_ALDH1L1_sig_genes_venn_diagram.png"), # Use NULL to plot interactively, or specify a filename to save the output
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
