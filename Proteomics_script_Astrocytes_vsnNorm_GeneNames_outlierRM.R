#Add outlier removal to the script and remove the MS053 replicate

library(proteus)
library(DEP)
library(WGCNA)
library(dplyr)
library(pmartR)
library(tidyr)


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

#Remove replicate samples and outlier samples
replicate_samples_to_remove <- c("MS053_A","MS053_B","MS053_C")

data = data[,!(names(data) %in% replicate_samples_to_remove)]

outlier_samples <- read.csv(paste0(pathResults, "outlier_samples.csv"))

data = data[,!(names(data) %in% outlier_samples$x)]

write.csv(data, paste0(pathResults,"MSQ2258_DIANN_report.pg_matrix_filtered.csv"))

data$Genes <- gsub(";", "___", data$Genes)

##MAKE SUMMARIZEDEXPERIMENT OBJECT

###Construct Experiment object
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep=',')

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]

head(sample_metadata)

Experiment <- as.data.frame(sample_metadata$Name)

colnames(Experiment)[1] <- "label" 

Experiment$condition <- sample_metadata$Status

replicate <- rep("1", nrow(Experiment))

Experiment <- cbind(Experiment, replicate)

write.csv(Experiment, paste0(pathResults, "Experiment2.csv"), row.names = FALSE)

####


#Reformat object using make_unique
ProtData_unique <- make_unique(data, "Genes", "Protein_Ids", delim = ";")

#Index the columns with the protein intensity data
other_cols <- c("Protein_Group", "Protein_Ids", "Protein_Names",  "Genes",  "Protein_Description", "name", "ID")

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

library(SummarizedExperiment)
Number_of_proteins3 <- nrow(assay(se))
ProData_se <- se 



#Save SummarizedExperiment assay data before normalisation: 
library(SummarizedExperiment)

unlogged_unnorm_assay <- 2^assay(ProData_se)
unlogged_unnorm_assay <- as.data.frame(unlogged_unnorm_assay)
unlogged_unnorm_assay  <- unlogged_unnorm_assay[,order(colnames(unlogged_unnorm_assay))]

write.csv(unlogged_unnorm_assay, paste0(pathResults, "Step1_filtered_unlogged_unnormalised.csv"))


logged_unnorm_assay <- assay(ProData_se)
logged_unnorm_assay <- as.data.frame(logged_unnorm_assay)
logged_unnorm_assay <- logged_unnorm_assay[,order(colnames(logged_unnorm_assay))]

write.csv(logged_unnorm_assay, paste0(pathResults, "Step2_filtered_logged_unnormalised.csv"))

#Apply VSN NORMALIZATION 
#Ideally, normalisation is applied before imputation?

ProData_norm_se <- normalize_vsn(ProData_se) # VSN is performed on original scale
assay(ProData_norm_se) #Explore assay data
meanSdPlot(ProData_norm_se)

#plot_normalization(ProData_se, ProData_norm_se)

logged_norm_assay <- assay(ProData_norm_se)
logged_norm_assay <- as.data.frame(logged_norm_assay)
logged_norm_assay  <- logged_norm_assay[,order(colnames(logged_norm_assay))]


write.csv(logged_norm_assay, paste0(pathResults, "Step3_filtered_logged_normalised.csv"))



#IMPUTATION

ProtData_imp <- impute(ProData_norm_se, fun = "QRILC")

ProtDataImpFinal <- assay(ProtData_imp)

logged_norm_imputed_assay <- ProtDataImpFinal
logged_norm_imputed_assay <- as.data.frame(logged_norm_imputed_assay)
logged_norm_imputed_assay <- logged_norm_imputed_assay[,order(colnames(logged_norm_imputed_assay))]


write.csv(logged_norm_imputed_assay, paste0(pathResults, "Step4_filtered_normalised_logged_imputed.csv"))


##USE empiricialBayesLM() to perform moderated adjustment for unwanted covariates

#Transpose the ProtDataImpFinal object
transposed_ProtDataImpFinal <- t(ProtDataImpFinal)

transposed_logged_norm_imputed_assay <- t(logged_norm_imputed_assay)

write.csv(transposed_logged_norm_imputed_assay, paste0(pathResults, "Step4_filtered_normalised_logged_transposed_imputed.csv"))


#Check out the clinical metadata and sample_metadata
clinical_metadata <- read.csv(paste0(path, "clinical_metadata_all_samples.csv"), sep=',')

clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% replicate_samples_to_remove), ]

clinical_metadata <- clinical_metadata[!(clinical_metadata$ID %in% outlier_samples$x), ]

head(clinical_metadata)  

sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"))

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]

head(sample_metadata)

#Process metadata to obtain covariates and retained_covariates objects 
covariates <- select(clinical_metadata, c("Age", "Sex", "Postmortem.Interval"))

colnames(covariates)[which(colnames(covariates) == "Postmortem.Interval")] <- "PMI"

covariates$PMI <- as.numeric(covariates$PMI)

covariates$PMI <- gsub(",", ".", covariates$PMI)

covariates$PMI <- as.numeric(covariates$PMI)

rownames(covariates) <- NULL

covariates$Age <- scale(covariates$Age)

covariates$PMI <- scale(covariates$PMI)

covariates$Sex <- as.factor(covariates$Sex)

retained <- select(sample_metadata, c("Name", "Status"))

rownames(retained) <- NULL

retained$Status <- as.factor(retained$Status)

#Regress out covariates whilst explictly modelling MS Status
data.eblm <- empiricalBayesLM(transposed_ProtDataImpFinal, removedCovariates = covariates, retainedCovariates=retained$Status)

#Look at a subset of the data.eblm$adjustedData counts
print(data.eblm$adjustedData[1:5, 1:5])

data.eblm_assay <- data.eblm$adjustedData
data.eblm_assay <- as.data.frame(data.eblm_assay)
data.eblm_assay <- data.eblm_assay[order(rownames(data.eblm_assay)),]

write.csv(data.eblm_assay, paste0(pathResults, "Step5_adjusted_empricialBayesLM.csv"))





##PROCESS DATA for readProteinGroups() function

#Process metadata
sample_metadata <- read.csv(paste0(path, "sample_metadata_all_samples.csv"), sep=',')

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% replicate_samples_to_remove), ]

sample_metadata <- sample_metadata[!(sample_metadata$Name %in% outlier_samples$x), ]


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

write.csv(metadata, paste0(pathResults, "metadata1.csv"), row.names = FALSE, quote = FALSE)


#Process Proteins file
#Take imputed data from ProtDataImpFinal??? 
#But not sure if I should read in adjusted value from data.eblm

Proteins_df <- as.data.frame(ProtDataImpFinal)

colnames(Proteins_df) <- gsub("^", "Intensity ", colnames(Proteins_df))

Proteins_df$Majority.protein.IDs <- rownames(Proteins_df)

Proteins_df$Reverse <- NA 

Proteins_df$Potential.contaminant <- NA

colnames(Proteins_df)[which(colnames(Proteins_df) == "Majority.protein.IDs")] <- "Majority protein IDs"
colnames(Proteins_df)[which(colnames(Proteins_df) == "Potential.contaminant")] <- "Potential contaminant"

rownames(Proteins_df) <- NULL

colnames(Proteins_df) <- gsub("__", "_", colnames(Proteins_df))

write.table(Proteins_df, paste0(pathResults, "Genes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

Genes <- read.delim(paste0(pathResults, "Genes.txt"))


##MAKE PROTEUSDATA OBJECT
prodat <- readProteinGroups(paste0(pathResults, "Genes.txt"), metadata)


##RUN SOME PLOTS
#plotSampleDistributions(prodat, title="Not normalized", fill="condition", method="violin")

#plotSampleDistributions(prodat_norm, title="vsn normalization", fill="condition", method="violin")#


##DIFFERENTIAL EXPRESSION ANALYSIS
##Note: since limma is going to log transform the data and it is already log transformed, we should inverse log it beforehand!!!

prodat$tab <- 2^(prodat$tab)

res <- limmaDE(prodat, ~condition, transform.fun = log2)  

colnames(res)[1] <- "Gene"

# Ordering the DataFrame by the 'Genes' column in alphabetical order
data <- data[order(data$Genes), ]

#Sanity check with the GeneIDs 
head(data$Genes[which(data$Genes %in% res$Gene)], 20)

head(res$Gene, 20)

res$ProteinIDs <- data$Protein_Ids[which(data$Genes %in% res$Gene)]

#Split the protein names to get first name in list for each Protein
maximum_number_of_splits = 2

Proteins_separated <- res %>%
  separate(col = ProteinIDs, into = paste("ProteinIDs", 1:maximum_number_of_splits, sep = "_"), sep = ";", remove = FALSE, extra = "merge", fill = "right")

res$Protein <- Proteins_separated$ProteinIDs_1

res$GeneNames <- res$Gene

res$GeneNames <- gsub("___", ";", res$GeneNames)

Genes_separated <- res %>%
  separate(col = Gene, into = paste("Gene", 1:maximum_number_of_splits, sep = "_"), sep = "___", remove = FALSE, extra = "merge", fill = "right")

res$Gene <- Genes_separated$Gene_1


# Assuming res is your DataFrame
ncols <- ncol(res)  # Number of columns in the DataFrame

# New order of columns
new_order <- c(1, (ncols-1), ncols, (ncols-2), 2:(ncols-3))

# Reorder the DataFrame columns
res <- res[, new_order]


write.csv(res, paste0(pathResults, "res_nonregressed_MS_vs_Ctrl.csv"), row.names = FALSE)


adj.P.Val_thresh <- 0.05

logFC_thresh <- 1

sig_res <- res %>% 
  filter(adj.P.Val <= adj.P.Val_thresh) 

write.csv(sig_res, paste0(pathResults, "sig_res_preadjusted_MS_vs_Ctrl.csv"), row.names = FALSE)


sig_res_strict <- res %>% 
  filter(adj.P.Val <= adj.P.Val_thresh & (logFC > logFC_thresh | logFC < -logFC_thresh))

write.csv(sig_res_strict, paste0(pathResults, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"), row.names = FALSE)



#Use simpler design matrix with adjusted data from empiricalBayes
#Can make the prodat assay equal to the adjust data.eblm data:
prodat_adj <- prodat
prodat_adj$tab <- t(data.eblm$adjustedData)

prodat_adj$tab <- 2^prodat_adj$tab

res <- limmaDE(prodat_adj, ~condition, transform.fun = log2)  


####ADJUST PROTEIN AND GENE NAME COLUMNS IN RES FILE 
colnames(res)[1] <- "Gene"

# Ordering the DataFrame by the 'Genes' column in alphabetical order
data <- data[order(data$Genes), ]

#Sanity check with the GeneIDs 
head(data$Genes[which(data$Genes %in% res$Gene)], 20)

head(res$Gene, 20)

res$ProteinIDs <- data$Protein_Ids[which(data$Genes %in% res$Gene)]

#Split the protein names to get first name in list for each Protein
maximum_number_of_splits = 2

Proteins_separated <- res %>%
  separate(col = ProteinIDs, into = paste("ProteinIDs", 1:maximum_number_of_splits, sep = "_"), sep = ";", remove = FALSE, extra = "merge", fill = "right")

res$Protein <- Proteins_separated$ProteinIDs_1

res$GeneNames <- res$Gene

res$GeneNames <- gsub("___", ";", res$GeneNames)

Genes_separated <- res %>%
  separate(col = Gene, into = paste("Gene", 1:maximum_number_of_splits, sep = "_"), sep = "___", remove = FALSE, extra = "merge", fill = "right")

res$Gene <- Genes_separated$Gene_1


# Assuming res is your DataFrame
ncols <- ncol(res)  # Number of columns in the DataFrame

# New order of columns
new_order <- c(1, (ncols-1), ncols, (ncols-2), 2:(ncols-3))

# Reorder the DataFrame columns
res <- res[, new_order]

write.csv(res, paste0(pathResults, "res_preadjusted_MS_vs_Ctrl.csv"), row.names = FALSE)

adj.P.Val_thresh <- 0.05

logFC_thresh <- 1

sig_res <- res %>% 
  filter(adj.P.Val <= adj.P.Val_thresh) 

write.csv(sig_res, paste0(pathResults, "sig_res_preadjusted_MS_vs_Ctrl.csv"), row.names = FALSE)


sig_res_strict <- res %>% 
  filter(adj.P.Val <= adj.P.Val_thresh & (logFC > logFC_thresh | logFC < -logFC_thresh))

write.csv(sig_res_strict, paste0(pathResults, "sig_res_preadjusted_MS_vs_Ctrl_strict.csv"), row.names = FALSE)
