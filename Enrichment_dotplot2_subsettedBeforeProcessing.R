library(ggplot2)
library(dplyr)


pathResults1 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_GFAP\\Enrichment\\")


pathResults2 <- ("C:\\Users\\teres\\OneDrive - Nexus365\\MSc Oxford\\Hilary Project\\Proteomics\\Astrocytes\\results\\vsnNorm\\outlierRM_ALDH1L1\\Enrichment\\")


#GFAP
GFAP_enriched_pathways <- read.delim(paste0(pathResults1, "GFAP_enrichment.all2.tsv"))

GFAP_enriched_pathways$GeneRatio <- GFAP_enriched_pathways$observed.gene.count/GFAP_enriched_pathways$background.gene.count

# New order of columns
new_order <- c(1, 2, 3, 4, 5, 10, 6, 7, 8, 9) 

GFAP_enriched_pathways <- GFAP_enriched_pathways[ ,new_order]

write.csv(GFAP_enriched_pathways, paste0(pathResults1, "GFAP_all_enriched_pathways2.csv"), row.names = FALSE)

df1 <- GFAP_enriched_pathways


#ALDH1L1
ALDH1L1_enriched_pathways <- read.delim(paste0(pathResults2, "ALDH1L1_enrichment.all2.tsv"))

ALDH1L1_enriched_pathways$GeneRatio <- ALDH1L1_enriched_pathways$observed.gene.count/ALDH1L1_enriched_pathways$background.gene.count

# New order of columns
new_order <- c(1, 2, 3, 4, 5, 10, 6, 7, 8, 9) 

ALDH1L1_enriched_pathways <- ALDH1L1_enriched_pathways[ ,new_order]

write.csv(ALDH1L1_enriched_pathways, paste0(pathResults2, "ALDH1L1_all_enriched_pathways2.csv"), row.names = FALSE)

df2 <- ALDH1L1_enriched_pathways


#Unique pathways GFAP: 
unique_GFAP <- df1$term.description[which(!df1$term.description %in% df2$term.description)]

unique_ALDH1L1 <- df2$term.description[which(!df2$term.description %in% df1$term.description)]


GFAP_selected_pathways <- c("Fibrinogen complex", "Hemostasis, and Dissolution of Fibrin Clot", "Regulation of TLR by endogenous ligand",
                                "Platelet Aggregation (Plug Formation)", "Common Pathway of Fibrin Clot Formation", 
                                "GRB2:SOS provides linkage to MAPK signaling for Integrins", "p130Cas linkage to MAPK signaling for integrins",
                                "MyD88:MAL(TIRAP) cascade initiated on plasma membrane", "Toll-like Receptor Cascades", "Blood clotting cascade",
                                "Fibrin complement receptor 3 signaling pathway", "MyD88:MAL(TIRAP) cascade initiated on plasma membrane", "VEGFA-VEGFR2 signaling")


# ALDH1L1 selected pathways 
ALDH1L1_selected_pathways <-  c("S100/CaBP-9k-type, calcium binding, subdomain, and Annexin", "Mixed, incl. S100/CaBP-9k-type, calcium binding, subdomain, and Cystatin superfamily", 
                                       "Proteinase inhibitor I25, cystatin, conserved site, and Endolysosome lumen",
                                       "Sulfuric ester hydrolase activity, and Other glycan degradation",
                                      "Positive regulation of plasma membrane repair, and EF hand",
                                      "Neutrophil degranulation",
                                      "Apoptotic cleavage of cell adhesion  proteins",
                                      "Dissolution of Fibrin Clot",
                                      "SUMOylation of DNA damage response and repair proteins",
                                      "Apoptotic cleavage of cellular proteins",
                                      "Interferon Signaling",
                                      "Adherens junctions interactions",
                                      "AnxA2-p11 complex", "Disorders of transmembrane transporters")

#Plot GFAP_selected_pathways 
selected_pathways <- GFAP_selected_pathways


# Filter the dataframe for only the selected pathways
filtered_df <- df1[df1$term.description %in% selected_pathways, ]

# Order the pathways by GeneRatio before plotting. We use dplyr for convenience
filtered_df <- filtered_df %>%
  group_by(term.description) %>%
  summarize(GeneRatio = mean(GeneRatio), # Assuming you want to order by the mean GeneRatio if there are duplicates
            false.discovery.rate = dplyr::first(false.discovery.rate), # Just take the first if there are duplicates
            observed.gene.count = dplyr::first(observed.gene.count)) %>% # Same for observed.gene.count
  ungroup() %>%
  arrange(GeneRatio) %>%
  mutate(term.description = factor(term.description, levels = unique(term.description)))


# Create the dotplot
ggplot(filtered_df, aes(x = GeneRatio, y = term.description)) +
  geom_point(aes(size = observed.gene.count, color = false.discovery.rate)) +
  scale_color_gradient(low = "blue", high = "red") + # Adjust colors as needed
  theme_minimal() +
  labs(title = "Pathway Analysis",
       x = "Gene Ratio",
       y = "Pathway",
       color = "False Discovery Rate",
       size = "Observed Protein Count") +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(pathResults1, "Enrichment_GFAP_astros_selected2.pdf"), height = 9, width = 12)
ggsave(paste0(pathResults1, "Enrichment_GFAP_astros_selected2.png"), height = 10, width =13, bg="white")


##Plot the ALDH1L1 selected pathways

selected_pathways <- ALDH1L1_selected_pathways


# Filter the dataframe for only the selected pathways
filtered_df <- df2[df2$term.description %in% selected_pathways, ]

# Order the pathways by GeneRatio before plotting. We use dplyr for convenience
filtered_df <- filtered_df %>%
  group_by(term.description) %>%
  summarize(GeneRatio = mean(GeneRatio), # Assuming you want to order by the mean GeneRatio if there are duplicates
            false.discovery.rate = dplyr::first(false.discovery.rate), # Just take the first if there are duplicates
            observed.gene.count = dplyr::first(observed.gene.count)) %>% # Same for observed.gene.count
  ungroup() %>%
  arrange(GeneRatio) %>%
  mutate(term.description = factor(term.description, levels = unique(term.description)))


# Create the dotplot
ggplot(filtered_df, aes(x = GeneRatio, y = term.description)) +
  geom_point(aes(size = observed.gene.count, color = false.discovery.rate)) +
  scale_color_gradient(low = "blue", high = "red") + # Adjust colors as needed
  theme_minimal() +
  labs(title = "Pathway Analysis",
       x = "Gene Ratio",
       y = "Pathway",
       color = "False Discovery Rate",
       size = "Observed Protein Count") +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0(pathResults2, "Enrichment_ALDH1L1_astros_selected2.pdf"), height = 9, width = 12)
ggsave(paste0(pathResults2, "Enrichment_ALDH1L1_astros_selected2.png"), height = 9, width =13, bg="white")





# #Select Pathways when all isoforms are included in analysis: 
# GFAP_selected_pathways <- c("Fibrinogen complex", "Hemostasis, and Dissolution of Fibrin Clot", "Immunoglobulin, and Immunoglobulin complex",
#                             "Fibrinogen, and Thrombophilia", "Complement and coagulation cascades, and Positive regulation of opsonization", 
#                             "Huntington disease", "Parkinson disease", "Alzheimer disease", "MHC class II antigen presentation", 
#                             "Neurotransmitter receptors and postsynaptic signal transmission", "Macroautophagy", 
#                             "Adaptive Immune System", "Regulation of TLR by endogenous ligand", "Platelet Aggregation (Plug Formation)", 
#                             "Common Pathway of Fibrin Clot Formation", "GRB2:SOS provides linkage to MAPK signaling for Integrins", 
#                             "p130Cas linkage to MAPK signaling for integrins", "Axon guidance","Signaling by Rho GTPases",  
#                             "MyD88:MAL(TIRAP) cascade initiated on plasma membrane", "Toll-like Receptor Cascades", "MAP2K and MAPK activation", 
#                             "Binding and Uptake of Ligands by Scavenger Receptors", "Blood clotting cascade", "Fibrin complement receptor 3 signaling pathway", 
#                             "Abnormality of circulating fibrinogen", "Microtubule cytoskeleton", "Beta tubulin, autoregulation binding site", 
#                             "Alpha tubulin", "Potassium channel tetramerisation-type BTB domain")
                            
                            
# # ALDH1L1 selected pathways 
# ALDH1L1_selected_pathways <- c("S100/CaBP-9k-type, calcium binding, subdomain, and Annexin", 
#                                "Mixed, incl. S100/CaBP-9k-type, calcium binding, subdomain, and Cystatin superfamily", 
#                                "Proteinase inhibitor I25, cystatin, conserved site, and Endolysosome lumen",
#                                "Sulfuric ester hydrolase activity, and Other glycan degradation", 
#                                "Positive regulation of plasma membrane repair, and EF hand", 
#                                "Tight junction", "Neutrophil degranulation", 
#                                "Apoptotic cleavage of cell adhesion  proteins",  
#                                "Dissolution of Fibrin Clot", 
#                                "SUMOylation of DNA damage response and repair proteins", "Apoptotic cleavage of cellular proteins", 
#                                "Interferon Signaling", "Adherens junctions interactions", "AnxA2-p11 complex")
                             