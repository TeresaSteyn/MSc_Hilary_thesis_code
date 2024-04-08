#Make PCA plots of GFAP and ALDH1L astrocyte data

#Load libraries
library(ggplot2)

path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")


pathResults <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_GFAP\\")

#pathResults2 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_ALDH1L1\\Enrichment\\")

adjusted_data <- read.csv(paste0(pathResults, "Step5_adjusted_empricialBayesLM.csv"))

replicate_samples_to_remove <- c("MS053_A","MS053_B","MS053_C")

outlier_samples <- read.csv(paste0(pathResults, "outlier_samples.csv"))

sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"))

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]

sample_metadata <- sample_metadata[grep("GFAP", sample_metadata$Astrocyte), ]

head(sample_metadata)

GFAP_samples <- sample_metadata$Name[grep("GFAP", sample_metadata$Astrocyte)]

clinical_metadata <- read.csv(paste0(path, "clinical_metadata_all_samples.csv"), sep=',')

clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% replicate_samples_to_remove), ]

clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% outlier_samples$x), ]

clinical_metadata <- clinical_metadata[which(clinical_metadata$ID %in% GFAP_samples), ]

head(clinical_metadata)  


#Make the rownames of the protein intensity data the same as the sample names from the metadata
rownames(adjusted_data) <- sample_metadata$Name

adjusted_data <- adjusted_data[,!(names(adjusted_data) %in% "X")]



#Run PCA
pca_result <- prcomp(adjusted_data, scale. = TRUE)

#Reformat metadata
metadata <- cbind(sample_metadata$Status, clinical_metadata$Age, clinical_metadata$Sex, clinical_metadata$Postmortem.Interval)
metadata <- as.data.frame(metadata)
colnames(metadata) <- c("Status", "Age", "Sex", "PMI")

metadata$Sex[which(metadata$Sex == 1)] = "Male" 
metadata$Sex[which(metadata$Sex == 2)] = "Female"

metadata$Age_range <- NA
metadata$Age_range[which(metadata$Age < 60)] = "<60"
metadata$Age_range[which(metadata$Age >= 60)] = ">=60"

metadata$PMI_range <- NA
metadata$PMI_range[which(metadata$PMI < 20)] = "<20"
metadata$PMI_range[which(metadata$PMI >= 20)] = ">=20"

metadata$Astrocyte_subclass <- sample_metadata$Astrocyte

# Extract PCA scores
#Can specify here which PCs you want to plot
pca_scores <- pca_result$x[, 1:2]  # Extracts the scores for first two consecutive PCs
colnames(pca_scores) <- c("PC1", "PC2")

# Ensure that your sample names in metadata match the row names in pca_scores
metadata$Sample <- rownames(pca_scores)

pca_scores <- as.data.frame(pca_scores)
pca_scores$Sample <- rownames(pca_scores)

# Merge PCA scores with metadata
pca_data <- merge(metadata, pca_scores, by = "Sample")


# Create PCA plot coloured by disease status
ggplot(pca_data, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 4) + # Increase point size here
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 30), # Increase axis titles size         
        axis.text = element_text(size = 30), # Increase axis text size         
        legend.title = element_text(size = 30), # Increase legend title size         
        legend.text = element_text(size = 30), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(x = "PC 1",
       y = "PC 2",
       shape = "Disease Status", color = "Disease status")

ggsave(paste0(pathResults, "PCA_plot_astrocytes_disease_status.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_disease_status.png"), height = 10, width = 14, bg = "white")
