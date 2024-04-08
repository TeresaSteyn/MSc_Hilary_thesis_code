#Add outlier removal to the script and remove the MS053 replicate

library(proteus)
library(DEP)
library(WGCNA)
library(dplyr)
library(pmartR)
library(ggplot2)

#Specify paths to read from and save data into
path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")

pathResults <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM2\\")



##QUALITY CONTROL
ProtData <- read.csv(paste0(path, "MSQ2258_DIANN_report.pg_matrix.csv"), sep = ";")

data <- ProtData

Number_of_proteins1 <- nrow(data)

#Remove contaminants

#Obtain indices for proteins that are contaminants
contaminants <- grep("contam", data$Protein_Ids, ignore.case=TRUE)

#Explore proteins IDs and descriptions that have been labelled as contaminants 
data[contaminants, c("Protein_Ids", "Protein_Description")]

#Remove the contaminants from the data
data <- data[-contaminants, ]

#Additional contaminants to remove: Keratin, Trypsin, casein, actin
Keratin <- grep("Keratin", data$Protein_Description, ignore.case=TRUE)

data <- data[-Keratin, ]

#Keratin2 <- grep("Keratin", data$Genes, ignore.case=TRUE)

#data <- data[-Keratin2, ]


Trypsin <- grep("Trypsin", data$Protein_Description, ignore.case=TRUE)

data <- data[-Trypsin, ]

Casein <- grep("Casein", data$Protein_Description, ignore.case=TRUE)

data <- data[-Casein, ]

Actin <- grep("Actin", data$Protein_Description, ignore.case=TRUE)

data <- data[-Actin, ]

Albumin <- grep("Albumin", data$Protein_Description, ignore.case=TRUE)

data <- data[-Albumin, ]

Number_of_proteins2 <- nrow(data)
print(Number_of_proteins2)

#write.csv(data, paste0(pathResults,"MSQ2258_DIANN_report.pg_matrix_filtered.csv"))

#Remove replicate samples 
replicate_samples_to_remove <- c("MS053_A","MS053_B","MS053_C")

data = data[,!(names(data) %in% replicate_samples_to_remove)]


##MAKE SUMMARIZEDEXPERIMENT OBJECT

###Construct Experiment object
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep=',')

sample_metadata = sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]


head(sample_metadata)

Experiment <- as.data.frame(sample_metadata$Name)

colnames(Experiment)[1] <- "label" 

Experiment$condition <- sample_metadata$Status

replicate <- rep("1", nrow(Experiment))

Experiment <- cbind(Experiment, replicate)

####


#Reformat object using make_unique
ProtData_unique <- make_unique(data, "Genes", "Protein_Ids", delim = ";")

#Index the columns with the protein intensity data
other_cols <- c("Protein_Group", "Protein_Ids", "Protein_Names",  "Genes", "Protein_Description", "name", "ID")

intensity_index <- which(!colnames(ProtData_unique) %in% other_cols)

ProData_se <-  make_se_parse(ProtData_unique, intensity_index) 

ProData_se$label <- Experiment$label
ProData_se$condition <- Experiment$condition
ProData_se$replicate <- Experiment$replicate
ProData_se$ID <- Experiment$label

#FILTER PROTEINS BY ROW: https://rdrr.io/bioc/DEP/man/filter_proteins.html
#Set the threshold for the minimum fraction of valid values allowed for any protein at 0.7
#Can just use the logged data as it gives the same result
se <- ProData_se
se <- filter_proteins(se, type = c("fraction"),
                      thr = NULL, min = 0.7)
ProData_se <- se 


#Apply VSN NORMALIZATION 
#Ideally, normalisation is applied before imputation?

ProData_norm_se <- normalize_vsn(ProData_se) # VSN is performed on original scale
assay(ProData_norm_se) #Explore assay data


#IMPUTATION
ProtData_imp <- impute(ProData_norm_se, fun = "QRILC")
ProtDataImpFinal <- assay(ProtData_imp)




###IDENTIFY and REMOVE OUTLIERS 
#Process metadata
clinical_metadata <- read.csv(paste0(path, "clinical_metadata_all_samples.csv"), sep=',')

replicate_samples_to_remove <- c("MS053_A","MS053_B","MS053_C")

clinical_metadata = clinical_metadata[!(clinical_metadata$ID %in% replicate_samples_to_remove), ]

sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep=',')

sample_metadata = sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]

measure <- rep("Intensity", nrow(sample_metadata))

metadata <- cbind(measure, sample_metadata$Name)

metadata <- as.data.frame(metadata)

colnames(metadata)[2] <- "sample"

metadata$condition <- NA
metadata$condition[grep("^M", metadata$sample)] <- "1"
metadata$condition[grep("^C", metadata$sample)] <- "0"

metadata$Age <- clinical_metadata$Age

metadata$Sex <- clinical_metadata$Sex

metadata$PMI <- clinical_metadata$Postmortem.Interval

outlier_samples <- read.csv(paste0(pathResults, "outlier_samples.csv"))

metadata$outlier <- NA
metadata$outlier[which(metadata$sample %in% outlier_samples$x)] <- "Outlier"
metadata$outlier[which(!metadata$sample %in% outlier_samples$x)] <- "Non-outlier"

#Process protein data
Data <- as.data.frame(ProtDataImpFinal)

colnames(Data) <- gsub("__", "_", colnames(Data))

##Run PCA 
pca_result <- prcomp(t(Data), scale. = TRUE)


# Extract PCA scores
#Can specify here which PCs you want to plot
pca_scores <- pca_result$x[, 1:2]  # Extracts the scores for first two consecutive principal components
colnames(pca_scores) <- c("PC1", "PC2")

# Ensure that your sample names in metadata match the row names in pca_scores
metadata$Sample <- rownames(pca_scores)

pca_scores <- as.data.frame(pca_scores)
pca_scores$Sample <- rownames(pca_scores)

# Merge PCA scores with metadata
pca_data <- merge(metadata, pca_scores, by = "Sample")

# Create PCA plot with labels for outliers
ggplot(pca_data, aes(x = PC1, y = PC2, color = outlier)) +
  geom_point(size =3) + # Plot points
  theme_minimal() +
  theme(text = element_text(size = 12), # Text size adjustments         
        axis.title = element_text(size = 14),         
        axis.text = element_text(size = 12),         
        legend.title = element_text(size = 12),         
        legend.text = element_text(size = 10),
       # panel.grid.major = element_blank(), # Remove major grid lines
       # panel.grid.minor = element_blank(), # Remove minor grid lines
        axis.line = element_line(colour = "black")) # Add axis lines
  labs(title = "PCA plot of samples coloured by disease status",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Outlier")

ggsave(paste0(pathResults, "PCA_plot_astrocytes_outliers.pdf"), height = 7, width = 9)
ggsave(paste0(pathResults, "PCA_plot_astrocytes_outliers.png"), height = 10, width = 10, bg = "white")


###Calculate Mahalanobis distance
#Process protein data
library(SummarizedExperiment)
Data2 <- as.data.frame(2^assay(ProData_se))

Data2$Genes <- rownames(Data2)

colnames(Data2) <- gsub("__", "_", colnames(Data2))

##MAKE proData OBJECT
MyProteinData<-as.proData(Data2,metadata,e_meta=NULL,edata_cname=c("Genes"),fdata_cname = c("sample"))

MyProteinData<-group_designation(MyProteinData,main_effects="condition")

MyProteinData<-edata_transform(MyProteinData,data_scale = 'log2')

#Sample outlier detection 
myfilter <- rmd_filter(MyProteinData, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness","Kurtosis"))

df <- as.data.frame(myfilter)

# Add a new column for color based on pvalue condition
df$color <- ifelse(df$pvalue < 0.001, "Significant", "Not Significant")

df$Status <- NA 
df$Status[grep("C", df$sample)] <- "Control"
df$Status[grep("M", df$sample)] <- "MS"

# Categorizing samples based on your conditions
df <- df %>%
  mutate(category = case_when(
    pvalue < 0.001 & Status == "MS" ~ "MS outlier",
    pvalue < 0.001 & Status == "Control" ~ "Control outlier",
    TRUE ~ "Non-outlier"
  ))

# Plotting
ggplot(df, aes(x = seq_along(Log2.md), y = Log2.md, color = category)) +
  geom_point(size = 3) + # Add points
  scale_color_manual(values = c("MS outlier" = "red", "Control outlier" = "blue", "Non-outlier" = "grey")) + # Custom colors
  labs(x = "Sample index", 
       y = "Log2 Mahalanobis distance", 
       color = "Outlier samples") + # Custom labels
  theme_minimal() + # Use a minimal theme
  theme(text = element_text(size = 12), # Increase base text size; affects title, labels, and legend         
        axis.title = element_text(size = 14), # Increase axis titles size         
        axis.text = element_text(size = 12), # Increase axis text size         
        legend.title = element_text(size = 12), # Increase legend title size         
        legend.text = element_text(size = 10)) # Increase legend text size
ggsave(paste0(pathResults, "Mahalanobis_distance_plot.pdf"), width = 10, height = 6)
ggsave(paste0(pathResults, "Mahalanobis_distance_plot.png"), width = 10, height = 10, bg = "white")