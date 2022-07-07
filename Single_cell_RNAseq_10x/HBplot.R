feats <- features
fc <- 0.3
#Modifications to SC table
v0 <- v0 %>% mutate(expression = ifelse(abs(avg_log2FC) < fc, 
                                            "No Significant",
                                        "Significant"
                                        ),
                    expression = ifelse(gene %in% feats,
                                        "Suppl table 1 gene list",
                                         expression),
                    expression = factor(expression,
                                        levels=c("No Significant", "Significant",
                                                 "Suppl table 1 gene list"))
                     
            )
      
set.seed(1)
p <- v0 %>% ggplot(aes(avg_log2FC, celltype)) +
geom_jitter(data = v0,
    size = 1.5, aes(color = expression),
    position = position_jitter(seed = 1))+
geom_label_repel(data = v0,
     aes(label = ifelse((gene %in% feats) | 
     (abs(avg_log2FC) > fc & !(gene %in% feats)),
                                          gene, NA),
                           color = expression),
                       size = 2.5,
                       force = 1,
                       vjust = -0.1,
                       max.overlaps = 7,
                      position = position_jitter(seed = 1),show.legend = F)+
      scale_x_continuous(trans = pseudolog10_trans,
                         limits = symmetric_limits) +
      scale_color_manual(values = c("grey", "darkgreen", "red")) +
      xlab("log2 Fold Change") +
  geom_vline(xintercept = c(-fc, fc), linetype = 2, size = 0.3, col = "grey20")+
      ggprism::theme_prism()+
  ggtitle("5 weeks vs 10 weeks")+
  theme(
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_blank(),
        legend.text = element_text(size=16),
        legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_blank(),
        plot.title = element_text(size=24, face="bold", hjust=0.5)
                  )+
  guides(color = guide_legend(override.aes = list(size=3)))
