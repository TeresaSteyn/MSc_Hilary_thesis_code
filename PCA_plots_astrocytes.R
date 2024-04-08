#Make PCA plots of astrocyte data

#Load libraries
library(ggplot2)

path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")

pathResults <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM2\\")

adjusted_data <- read.csv(paste0(pathResults, "Step5_adjusted_empricialBayesLM.csv"))

replicate_samples_to_remove <- c("MS053_A","MS053_B","MS053_C")

outlier_samples <- read.csv(paste0(pathResults, "outlier_samples.csv"))

clinical_metadata <- read.csv(paste0(path, "clinical_metadata_all_samples.csv"), sep=',')

clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% replicate_samples_to_remove), ]

clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% outlier_samples$x), ]

head(clinical_metadata)  

sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"))

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]

head(sample_metadata)

#Make the rownames of the protein intensity data the same as the sample names from the metadata
rownames(adjusted_data) <- sample_metadata$Name

adjusted_data <- adjusted_data[,!(names(adjusted_data) %in% "X")]



#Run PCA
pca_result <- prcomp(adjusted_data, scale. = TRUE)

#Plot explained variance 

# Get summary of PCA
pca_summary <- summary(pca_result)

# Plotting the variance explained by each principal component
plot(pca_summary$importance[2,], xlab = "Principal Component", ylab = "Proportion of Variance Explained", type = 'b', pch = 19, col = "blue")
title(main = "Scree Plot")

#Run the bar plot of %variance explained for each PC 

# Calculate variance explained by each principal component
var_explained <- pca_summary$importance[2,] * 100 # Converting to percentage

# Create a more aesthetic bar plot
pdf(paste0(pathResults, "PCA_percent_variance_explained.pdf"))
barplot(var_explained, names.arg = 1:length(var_explained), 
        xlab = "Principal Component", ylab = "Percentage of Variance Explained", 
        main = "Variance Explained by Each Principal Component", 
        col = "steelblue", border = "darkblue", ylim = c(0, 20), las = 1)
dev.off()

write.csv(var_explained, paste0(pathResults, "PCA_variance_explained.csv"))

##PCA plot

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

pca_data$Astrocyte_subclass <- gsub("GFAP", "GFAP+", pca_data$Astrocyte_subclass)
pca_data$Astrocyte_subclass <- gsub("ALDH1L1", "ALDH1L1+", pca_data$Astrocyte_subclass)
pca_data$Astrocyte_subclass <- gsub("Both", "GFAP+/ALDH1L1+", pca_data$Astrocyte_subclass)

# Create PCA plot coloured by disease status
ggplot(pca_data, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 3) + # Increase point size here
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 25), # Increase axis titles size         
        axis.text = element_text(size = 25), # Increase axis text size         
        legend.title = element_text(size = 25), # Increase legend title size         
        legend.text = element_text(size = 25), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(x = "PC 1",
       y = "PC 2",
       shape = "Disease Status", color = "Disease status")

ggsave(paste0(pathResults, "PCA_plot_astrocytes_disease_status.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_disease_status.png"), height = 10, width = 14, bg = "white")

# Create PCA plot coloured by Astrocyte subclass

ggplot(pca_data, aes(x = PC1, y = PC2, color = Astrocyte_subclass)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(text = element_text(size = 12), # Base text size; affects title, labels, and legend         
        axis.title = element_text(size = 25), # Axis titles size         
        axis.text = element_text(size = 25), # Axis text size         
        legend.title = element_text(size = 25), # Legend title size         
        legend.text = element_text(size = 25), # Legend text size
        axis.line = element_line(colour = "black")) + # Axis lines
  guides(color = guide_legend(keyheight = unit(2.2, "lines"), keywidth = unit(1, "lines"))) +
  labs(x = "PC1",
       y = "PC2",
       color = "Astrocyte subclass") + 
        scale_color_manual(values = c("orange", "#0CB702", "purple"))

ggsave(paste0(pathResults, "PCA_plot_astrocytes_subclass.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_subclass.png"), height = 10, width =14, bg="white")



# Create PCA plot coloured by Age
ggplot(pca_data, aes(x = PC1, y = PC2, color = Age)) +
  geom_point(size = 3) +
   theme_minimal() +   
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 12), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 12), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA plot of samples coloured by Age",
       x = "PC 1",
       y = "PC 2",
       color = "Age")
ggsave(paste0(pathResults, "PCA_plot_astrocytes_age.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_age.png"), height = 10, width =10, bg="white")


# Create PCA plot coloured by Sex
ggplot(pca_data, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA Plot of Control Samples by Sex",
       x = "PC 1",
       y = "PC 2", color = "Sex")

ggsave(paste0(pathResults, "PCA_plot_astrocytes_sex.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_sex.png"), height = 10, width =10, bg="white")

# Create PCA plot coloured by PMI
ggplot(pca_data, aes(x = PC1, y = PC2, color = PMI)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA plot of samples coloured by PMI",
       x = "PC 1",
       y = "PC 2",
       color = "PMI")
ggsave(paste0(pathResults, "PCA_plot_astrocytes_PMI.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_PMI.png"), height = 10, width =10, bg="white")


# Create PCA plot coloured by Age_range
ggplot(pca_data, aes(x = PC1, y = PC2, color = Age_range)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA plot of samples coloured by Age_range",
       x = "PC 1",
       y = "PC 2",
       color = "Age range")
ggsave(paste0(pathResults, "PCA_plot_astrocytes_AgeRange.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_AgeRange.png"), height = 10, width =10, bg="white")


# Create PCA plot coloured by PMI_range
ggplot(pca_data, aes(x = PC1, y = PC2, color = PMI_range)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA plot of samples coloured by PMI_range",
       x = "PC 1",
       y = "PC 2",
       color = "PMI range")
ggsave(paste0(pathResults, "PCA_plot_astrocytes_PMIRange.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_PMIRange.png"), height = 10, width =10, bg="white")


ggplot(pca_data, aes(x = PC1, y = PC2, shape = Status, color = Astrocyte_subclass)) +
 geom_point(size = 3) +
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA plot of samples identified by disease status and astrocyte subclass",
       x = "PC 1",
       y = "PC 2",
       shape = "Disease status", color = "Astrocyte subclass")
ggsave(paste0(pathResults, "PCA_plot_astrocytes_StatusandAstrosubclass.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_StatusandAstrosubclass.png"), height = 10, width =10, bg="white")


# Subset for Control samples
control_data <- subset(pca_data, Status == "Control")

# PCA plot for Control group
ggplot(control_data, aes(x = PC1, y = PC2, shape = Status, color = Astrocyte_subclass)) +
  xlim(-50, 100)+
  ylim(-20, 60)+
  geom_point(size = 3, shape = 16) +
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA Plot of Control Samples by Astrocyte Subclass",
       x = "PC 1",
       y = "PC 2",
       shape = "Disease Status", color = "Astrocyte Subclass")
# Save the Control plot
ggsave(paste0(pathResults, "PCA_plot_Control_AstrocyteSubclass.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_Control_AstrocyteSubclass.png"), height = 10, width = 10, bg = "white")


# Subset for MS samples
ms_data <- subset(pca_data, Status == "MS")

# PCA plot for MS group
ggplot(ms_data, aes(x = PC1, y = PC2, shape = Status, color = Astrocyte_subclass)) +
  xlim(-50, 100)+
  ylim(-20, 60)+
  geom_point(size = 3, shape = 17) +  # Specify triangle shape for all points
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA Plot of MS Samples by Astrocyte Subclass",
       x = "PC 1",
       y = "PC 2",
       shape = "Disease Status", color = "Astrocyte Subclass")

# Save the MS plot
ggsave(paste0(pathResults, "PCA_plot_MS_AstrocyteSubclass.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_MS_AstrocyteSubclass.png"), height = 10, width = 10, bg = "white")


#GFAP+ samples only coloured by disease status: 
GFAP_data <- subset(pca_data, Astrocyte_subclass == "GFAP+")

ggplot(GFAP_data, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 3) + # Increase point size here
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 12), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 12), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA Plot of Control Samples by Astrocyte disease status",
       x = "PC 1",
       y = "PC 2",
       shape = "Disease Status", color = "Disease status")

ggsave(paste0(pathResults, "PCA_plot_GFAP_astrocytes_disease_status.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_GFAP_astrocytes_disease_status.png"), height = 10, width = 10, bg = "white")


#ALDH1L1+ samples only coloured by disease status: 
ALDH1L1_data <- subset(pca_data, Astrocyte_subclass == "ALDH1L1+")

ggplot(ALDH1L1_data, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 3) + # Increase point size here
  theme_minimal() +
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 12), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 12), # Increase legend text size
        axis.line = element_line(colour = "black")) + # Add axis lines
  labs(title = "PCA Plot of Control Samples by Astrocyte disease status",
       x = "PC 1",
       y = "PC 2",
       shape = "Disease Status", color = "Disease status")

ggsave(paste0(pathResults, "PCA_plot_ALDH1L1_astrocytes_disease_status.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_ALDH1L1_astrocytes_disease_status.png"), height = 10, width = 10, bg = "white")
