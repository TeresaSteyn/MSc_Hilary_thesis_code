#Plot volcano plots coloured by upregulated, downregulated, and non-signficant and add labels for genes of a certain significance level

library(tidyverse)
library(cowplot)
library(dplyr)
library(tibble)
library(png)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(stringr)

pathResults <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_ALDH1L1\\")

#Import csv of all genes (i.e. res_tbl) 

res_tbl <- read.csv(paste0(pathResults, "res_preadjusted_MS_vs_Ctrl.csv"))  


#Can set this part individually for specific clusters depending on their sets of DEGs
#e.g. Example threshold setting: threshold = padj < 0.001 & (abs(log2FoldChange) >= 2.0 | log2FoldChange <= -2.5)



res_tbl <- res_tbl %>% 
       mutate(threshold = (adj.P.Val < 0.01 & (abs(logFC) > 3 | logFC < -1)) | (adj.P.Val < 0.001 & abs(logFC) > 1) | (adj.P.Val < 0.05 & logFC < -2))


#Plot all genes and sets threshold that they colour the genes by 
res_tbl <- na.omit(res_tbl)

direction = ifelse(res_tbl$adj.P.Val > 0.05 | abs(res_tbl$logFC) < 1, "Non-significant", 
                   ifelse(res_tbl$logFC > 0, "Upregulated", 
                          ifelse(res_tbl$logFC < 0, "Downregulated", "NA")))

res_tbl$direction_of_change <- direction


  ggplot(res_tbl) +
  geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = direction_of_change), size = 5) +
  ggtitle("Volcano plot MS vs control ALDH1L1 Astrocytes") +
  scale_color_manual(name = "direction_of_change", values = c("Upregulated" = "#D22B2B", "Downregulated" = "#0096FF", "Non-significant" = "grey")) +
  xlab("logFC") + 
  ylab("-log10 adjusted p-value") +
  ylim(0, 4) +
  scale_y_continuous(limits = c(0,05)) +
  theme(legend.position = "left",
        text = element_text(size = 25), # Increase base text size; affects title, labels, and legend
        axis.title = element_text(size = 25), # Increase axis titles size
        axis.text = element_text(size = 25), # Increase axis text size
        legend.title = element_text(size = 25), # Increase legend title size
        legend.text = element_text(size = 25), # Increase legend text size
        axis.line = element_line(colour = "black"), # Add axis lines
        plot.title = element_text(size = rel(1.0), hjust = 0.5)) +
  guides(col = guide_legend("Direction of change")) + 
  geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val), label = ifelse(threshold == TRUE, as.character(Protein), "")),
                  size = 7,  # Try smaller sizes if needed
                  point.padding = unit(0.3, "lines"),  # Adjust padding
                  min.segment.length = 0.1,  # Reduce segment length to avoid extensive repulsion calculations
                  max.overlaps = 100)
	      ggsave(paste0(pathResults, "MS_vs_Control_GFAP_volcano_plot.pdf"), width = 17, height = 17)
	      ggsave(paste0(pathResults, "MS_vs_Control_GFAP_volcano_plot.png"), height = 14, width =15)
	      
	      


