library(dplyr)
library(tidyr)

pathResults <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM2\\")

DE_results <- read.csv(paste0(pathResults, "sig_res_preadjusted_ALDH1L1_strict.csv"))

DE_results$GeneNames

maximum_number_of_splits <- 6

Genes_separated <- DE_results%>%
  separate(col = GeneNames, into = paste("Genes", 1:maximum_number_of_splits, sep = "_"), sep = ";", remove = FALSE, extra = "merge", fill = "right")

all_DEGs <- unlist(Genes_separated[4:9], use.names = FALSE)

all_DEGs <- data.frame(all_DEGs)

all_DEGs <- na.omit(all_DEGs)

write.csv(all_DEGs, paste(pathResults, "all_sig_DEGs_ALDH1L1.csv"), row.names = FALSE, quote = FALSE)