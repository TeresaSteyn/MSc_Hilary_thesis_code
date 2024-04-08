##Extract metadata for Astrocyte proteomics data samples from All Cases file 

path <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\data\\")

pathResults <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\")

path2 <-  ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Metadata\\")

ProtData <- read.csv(paste0(path, "MSQ2258_DIANN_report.pg_matrix.csv"), sep = ";")

data <- ProtData

MS_metadata <- read.csv(paste0(path2,"MS_Clinical_Information.csv"), sep = ";")

Controls_metadata <- read.csv(paste0(path2,"Controls_Clinical_Information.csv"), sep = ";")

#Rearrange the column and sample names from the raw data and metadata to match each other

colnames(MS_metadata)[5] <- "Subtype"

All_metadata <- rbind(MS_metadata, Controls_metadata)

samples <- colnames(data)[6:83]

samples <- gsub("_.*", "", samples)

#samples <- gsub("M", "MS", samples)

#samples <- gsub("C(\\d{1,2})\\b", "C0\\1", samples)

#samples <- gsub("MS(\\d{1,2})\\b", "MS0\\1", samples)

#Find list of unique samples from the raw data
samples_list <- unique(samples)

length(samples_list)

#Check all the samples are present in the metadata and again fix any naming issues if not 
metadata_in_samples <- All_metadata$ID[which(All_metadata$ID %in% samples_list)]
length(metadata_in_samples)

samples_list[which(samples_list %in% All_metadata$ID)]

not_in_list <- samples_list[which(!samples_list %in% All_metadata$ID)]

find_in_list <- gsub("^C", "PDC", not_in_list)

All_metadata$ID[c(which(All_metadata$ID %in% find_in_list))] <- not_in_list

samples_in_metadata <- samples_list[which(samples_list %in% All_metadata$ID)]
print(samples_in_metadata)

length(samples_in_metadata)

#Extract the metadata for each of the samples in the raw data
clinical_metadata <- All_metadata[which(All_metadata$ID %in% samples_list), ] #Clinical metadata???

write.csv(clinical_metadata, paste0(path, "clinical_metadata.csv"), row.names = FALSE, quote = TRUE)



##CREATE FILE WITH CLINICAL METADATA FOR ALL SAMPLES 

#Add metadata that we see in clinical_metadata to data_imp
clinical_metadata <- read.csv(paste0(path,"clinical_metadata.csv"), sep = ",")

colnames(clinical_metadata)[grep("Age", colnames(clinical_metadata))] <- "Age"
colnames(clinical_metadata)[grep("Sex", colnames(clinical_metadata))] <- "Sex"


#Extract sample names from data
sample_names <- colnames(data)[6:ncol(data)]

sample_names2 <- gsub("_.*", "", sample_names)


#Map each sample in data_imp to its corresponding entry in clinical_metadata
# Assuming the first column of clinical_metadata contains the unique sample names
metadata_mapping <- clinical_metadata[match(sample_names2, clinical_metadata[,1]), ]

metadata_mapping$ID <- sample_names

metadata_mapping$Postmortem.Interval <- gsub(",", ".", metadata_mapping$Postmortem.Interval)

write.csv(metadata_mapping, paste0(path, "clinical_metadata_all_samples.csv"), row.names = FALSE, quote = TRUE)


##CREATE FILE WITH SAMPLE METADATA FOR ALL SAMPLES 

sample_metadata <- read.csv(paste0(path,"Description_IDs_astrocytes.csv"), sep = ";")

#Adjust the samples in sample_metadata$Name to match with new data_imp colnames
sample_metadata$Name <- gsub("_.*", "", sample_metadata$Name)

sample_metadata$Name <- sample_names
  
write.csv(sample_metadata, paste0(path, "sample_metadata_all_samples.csv"), row.names = FALSE, quote = TRUE)


  