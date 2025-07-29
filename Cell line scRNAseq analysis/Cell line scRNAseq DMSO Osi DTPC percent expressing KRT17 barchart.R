library(plyr)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(forcats)



##### HCC827 HCC4006 H1975 DMSO Osi DTPC KRT17 expression barchart #####

#Create dataframe
krt17_exp <- data.frame(Condition = c("HCC827_DMSO", "HCC827_Osi_DTPC",
                                      "HCC4006_DMSO", "HCC4006_Osi_DTPC",
                                      "H1975_DMSO", "H1975_Osi_DTPC"),
                        KRT17_pos = c(0.5, 8.2,
                                      1.8, 54.3,
                                      6.9, 66.2))

#Specify order for plotting
krt17_exp$Condition <- factor(krt17_exp$Condition, levels = c("HCC827_DMSO", "HCC827_Osi_DTPC",
                                                              "HCC4006_DMSO", "HCC4006_Osi_DTPC",
                                                              "H1975_DMSO", "H1975_Osi_DTPC"))

#Draw barchart
krt17_exp_bp <- ggplot(krt17_exp, aes(x = Condition, y = KRT17_pos, fill = Condition)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = KRT17_pos), vjust = -0.25) +
  scale_y_continuous("%KRT17+ cells", limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(values = c("HCC827_DMSO" = "lightblue",
                               "HCC827_Osi_DTPC" = "blue",
                               "HCC4006_DMSO" = "lightgreen",
                               "HCC4006_Osi_DTPC" = "green",
                               "H1975_DMSO" = "lightcoral",
                               "H1975_Osi_DTPC" = "red")) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none") +
  labs(title = "%KRT17+ cells")

krt17_exp_bp #275 x 500 (WXH)




