#Add outlier removal to the script and remove the MS053 replicate

library(proteus)
library(DEP)
library(WGCNA)
library(dplyr)
library(pmartR)


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

write.csv(Experiment, paste0(pathResults, "Experiment.csv"), row.names = FALSE)

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

Number_of_proteins3 <- nrow(assay(se))
ProData_se <- se 

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

write.csv(metadata, paste0(pathResults, "metadata1.csv"), row.names = FALSE, quote = FALSE)

#Process protein data
Data <- as.data.frame(2^assay(ProData_se))

Data$Genes <- rownames(Data)

colnames(Data) <- gsub("__", "_", colnames(Data))

##MAKE proData OBJECT
MyProteinData<-as.proData(Data,metadata,e_meta=NULL,edata_cname=c("Genes"),fdata_cname = c("sample"))

MyProteinData<-group_designation(MyProteinData,main_effects="condition")

MyProteinData<-edata_transform(MyProteinData,data_scale = 'log2')

#Sample outlier detection 
myfilter <- rmd_filter(MyProteinData, metrics = c("Correlation", "Proportion_Missing", "MAD", "Skewness","Kurtosis"))

Summary_of_filtered_data <- summary(myfilter, pvalue_threshold = 0.001)

#Filter out the outlier samples

write.csv(Summary_of_filtered_data$filtered_samples, paste0(pathResults, "outlier_samples.csv"), row.names = FALSE)