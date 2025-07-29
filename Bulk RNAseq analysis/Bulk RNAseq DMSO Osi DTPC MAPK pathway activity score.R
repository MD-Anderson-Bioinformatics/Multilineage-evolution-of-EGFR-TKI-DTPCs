library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)


##### HCC827 DMSO vs Osi DTPC bulk RNAseq analysis #####

#Read in data
hcc827_mapk <- HCC827_DMSO_Osi_DTPC_bulk_RNAseq_MAPK_pathway_activity_score_gene_expression

#Remove first column names
hcc827_mapk <- hcc827_mapk[, -1]

#Z score
hcc827_mapk <- t(scale(t(hcc827_mapk)))

#Convert to dataframe
hcc827_mapk <- as.data.frame(hcc827_mapk)

#Add gene names back as first column
hcc827_mapk <- hcc827_mapk %>% 
  mutate(Gene = c("PHLDA1", "SPRY4", "DUSP4", "DUSP6", 
                  "CCND1", "EPHA2", "EPHA4", "ETV4", "ETV5")) %>%
  select(Gene, everything())

hcc827_mapk

#Calculate MAPK pathway activity score
hcc827_mapk_score <- hcc827_mapk %>% select(-Gene) %>%
  summarise(across(everything(), ~ sum(.x)/sqrt(9)))

#Note: using sqrt of 9 as SPRY2 not detected in HCC827 RNAseq

#Transpose data frame
hcc827_mapk_score_t <- as.data.frame(t(hcc827_mapk_score))

#Rename MAPK activity score column
colnames(hcc827_mapk_score_t) <- "MAPK_pathway_activity_score"

#Add back sample names 
hcc827_mapk_score_t$sample <- rownames(hcc827_mapk_score_t)

#Set blank rownames
rownames(hcc827_mapk_score_t) <- NULL

#Add condition column
hcc827_mapk_score_t <- hcc827_mapk_score_t %>% 
  mutate(Condition = c("DMSO", "DMSO", "DMSO", "Osi_DTPC", "Osi_DTPC", "Osi_DTPC"))


###MAPK pathway activity score boxplot plot
hcc827_mapk_score_bp <- ggplot(hcc827_mapk_score_t, aes(x = Condition, y = MAPK_pathway_activity_score,
                                                        fill = Condition)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test",
                     label.x = 1.1,
                     label.y = 2.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black")) +
  scale_fill_manual(values = c("DMSO" = "blue",
                               "Osi_DTPC" = "red")) +
  scale_y_continuous("MAPK Pathway Activity Score") +
  labs(title = "HCC827")



hcc827_mapk_score_bp #Plot: 330 x 400



##### HCC4006 DMSO vs Osi DTPC bulk RNAseq analysis #####

#Read in data
hcc4006_mapk <- HCC4006_DMSO_Osi_DTPC_bulk_RNAseq_MAPK_pathway_activity_score_gene_expression

#Remove first column names
hcc4006_mapk <- hcc4006_mapk[, -1]

#Z score
hcc4006_mapk <- t(scale(t(hcc4006_mapk)))

#Convert to dataframe
hcc4006_mapk <- as.data.frame(hcc4006_mapk)

#Add gene names back as first column
hcc4006_mapk <- hcc4006_mapk %>% 
  mutate(Gene = c("PHLDA1", "SPRY2", "SPRY4", "DUSP4", "DUSP6", 
                  "CCND1", "EPHA2", "EPHA4", "ETV4", "ETV5")) %>%
  select(Gene, everything())

hcc4006_mapk

#Calculate MAPK pathway activity score
hcc4006_mapk_score <- hcc4006_mapk %>% select(-Gene) %>%
  summarise(across(everything(), ~ sum(.x)/sqrt(10)))

#Transpose data frame
hcc4006_mapk_score_t <- as.data.frame(t(hcc4006_mapk_score))

#Rename MAPK activity score column
colnames(hcc4006_mapk_score_t) <- "MAPK_pathway_activity_score"

#Add back sample names 
hcc4006_mapk_score_t$sample <- rownames(hcc4006_mapk_score_t)

#Set blank rownames
rownames(hcc4006_mapk_score_t) <- NULL

#Add condition column
hcc4006_mapk_score_t <- hcc4006_mapk_score_t %>% 
  mutate(Condition = c("DMSO", "DMSO", "DMSO", "Osi_DTPC", "Osi_DTPC", "Osi_DTPC"))


###MAPK pathway activity score boxplot plot
hcc4006_mapk_score_bp <- ggplot(hcc4006_mapk_score_t, aes(x = Condition, y = MAPK_pathway_activity_score,
                                                          fill = Condition)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test",
                     label.x = 1.1,
                     label.y = 3.4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black")) +
  scale_fill_manual(values = c("DMSO" = "blue",
                               "Osi_DTPC" = "red")) +
  scale_y_continuous("MAPK Pathway Activity Score") +
  labs(title = "HCC4006")



hcc4006_mapk_score_bp #Plot: 330 x 400



##### H1975 DMSO vs Osi DTPC bulk RNAseq analysis #####

#Read in data
h1975_mapk <- H1975_DMSO_Osi_DTPC_bulk_RNAseq_MAPK_pathway_activity_score_gene_expression

#Remove first column names
h1975_mapk <- h1975_mapk[, -1]

#Z score
h1975_mapk <- t(scale(t(h1975_mapk)))

#Convert to dataframe
h1975_mapk <- as.data.frame(h1975_mapk)

#Add gene names back as first column
h1975_mapk <- h1975_mapk %>% 
  mutate(Gene = c("PHLDA1", "SPRY2", "SPRY4", "DUSP4", "DUSP6", 
                  "CCND1", "EPHA2", "EPHA4", "ETV4", "ETV5")) %>%
  select(Gene, everything())

h1975_mapk

#Calculate MAPK pathway activity score
h1975_mapk_score <- h1975_mapk %>% select(-Gene) %>%
  summarise(across(everything(), ~ sum(.x)/sqrt(10)))

#Transpose data frame
h1975_mapk_score_t <- as.data.frame(t(h1975_mapk_score))

#Rename MAPK activity score column
colnames(h1975_mapk_score_t) <- "MAPK_pathway_activity_score"

#Add back sample names 
h1975_mapk_score_t$sample <- rownames(h1975_mapk_score_t)

#Set blank rownames
rownames(h1975_mapk_score_t) <- NULL

#Add condition column
h1975_mapk_score_t <- h1975_mapk_score_t %>% 
  mutate(Condition = c("DMSO", "DMSO", "DMSO", "Osi_DTPC", "Osi_DTPC", "Osi_DTPC"))


###MAPK pathway activity score boxplot plot
h1975_mapk_score_bp <- ggplot(h1975_mapk_score_t, aes(x = Condition, y = MAPK_pathway_activity_score,
                                                      fill = Condition)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test",
                     label.x = 1.1,
                     label.y = 3.1) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black")) +
  scale_fill_manual(values = c("DMSO" = "blue",
                               "Osi_DTPC" = "red")) +
  scale_y_continuous("MAPK Pathway Activity Score") +
  labs(title = "H1975")



h1975_mapk_score_bp #Plot: 330 x 400











