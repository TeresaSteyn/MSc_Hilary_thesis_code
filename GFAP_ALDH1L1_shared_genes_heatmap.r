library(dplyr)
library(tidyr)

path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")

pathResults1 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_GFAP\\")

pathResults2 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_ALDH1L1\\")


pathResults3 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM2\\")

GFAP_res <- read.csv(paste0(pathResults1, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"))  

ALDH1L1_res <- read.csv(paste0(pathResults2, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"))  

GFAP_res[which(GFAP_res$Gene %in% ALDH1L1_res$Gene), ]

shared_genes <- GFAP_res$Gene[which(GFAP_res$Gene %in% ALDH1L1_res$Gene)]
shared_genes <- unique(shared_genes)
write.csv(shared_genes, paste0(pathResults3, "shared_genes_GFAP_ALDH1L1.csv"), row.names = FALSE)


# Load the necessary package
library(ComplexHeatmap)
library(circlize)

# Assuming GFAP_res and ALDH1L1_res are already loaded and shared_genes is defined

# Filter and prepare logFC data
gfap_logFC <- GFAP_res[!duplicated(GFAP_res$Gene) & GFAP_res$Gene %in% shared_genes, c("Gene", "logFC")]
aldh1l1_logFC <- ALDH1L1_res[ALDH1L1_res$Gene %in% shared_genes, c("Gene", "logFC")]

# Combine into a matrix with conditions as rows and genes as columns
# Note: Direct assignment of row names from shared_genes ensures alignment and correct labeling
logFC_matrix <- rbind(GFAP = gfap_logFC$logFC[match(shared_genes, gfap_logFC$Gene)],
                      ALDH1L1 = aldh1l1_logFC$logFC[match(shared_genes, aldh1l1_logFC$Gene)])

colnames(logFC_matrix) <- shared_genes  # Ensuring column names are the gene names

# Adjusting the color mapping based on the range of your logFC values
# Determine the range of logFC values for more accurate color mapping
range_logFC <- range(c(logFC_matrix), na.rm = TRUE)

my_colors <- colorRamp2(c(range_logFC[1], 0, range_logFC[2]), c("blue", "white", "red"))

# Draw the heatmap
filename1 <- paste0(pathResults3, "Heatmap_shared_genes.pdf")

pdf(filename1, width = 15, height = 15)
Heatmap(logFC_matrix,
        name = "logFC",
        col = my_colors,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 10, fontface = "bold", angle = 45),
        cluster_rows = FALSE, # Do not cluster rows
        cluster_columns = FALSE) # Do not cluster columns
dev.off()

