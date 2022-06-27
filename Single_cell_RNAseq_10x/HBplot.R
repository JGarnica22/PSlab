feats <- read.table("data/list_TR1.txt")$V1

#Modifications to SC table
sc <- sc %>% mutate(expression = ifelse(abs(avg_log2FC) < 1, 
                                            "NS", "SIG")) %>% 
mutate(expression = ifelse(gene %in% feats,
                                            "107", expression))
      
set.seed(1)
p <- sc %>% ggplot(aes(avg_log2FC, cluster)) +
geom_jitter(data = sc,
    size = 1.5, aes(color = expression),
    position = position_jitter(seed = 1))+
geom_label_repel(data = sc,
     aes(label = ifelse((gene %in% feats) | 
     (abs(avg_log2FC) > 1.5 & !(gene %in% features)),
                                          gene, NA),
                           color = expression),
                       size = 2.5,
                       force = 1,
                       vjust = -0.1,
                       max.overlaps = 7,
                      position = position_jitter(seed = 1))+
      scale_x_continuous(trans = pseudolog10_trans,
                         limits = symmetric_limits) +
      scale_color_manual(labels = c("107 Gene List", "Non Significant",
                                    "Significant"),
                         values = c("red", "grey", "darkgreen")) +
      ylab("Cluster") +
      xlab("Avglog2FC DE genes") +
      theme_bw()+
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            axis.title.y = element_blank()
      )+
      labs(title = paste0("DE analysis: ", check[i], " vs ", check[j]))
