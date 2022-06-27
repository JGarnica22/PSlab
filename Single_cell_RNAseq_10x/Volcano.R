Vol <- FindMarkers(treg.combined2,
                    ident.1 = check[i], ident.2 = check[j],
                           test.use = "wilcox", logfc.threshold = 0.25,
                           min.pct = 0.05)

      Vol$gene <- row.names(Vol)
      Vol <- Vol[!grepl("mt", Vol$gene) & !grepl("Rp", Vol$gene),]
      Vol$sig <- -log10(Vol$p_val_adj)
      Vol$log2FC <- Vol$avg_log2FC
      Vol$S <- "0"
      Vol[which(Vol$sig<1.5), "S"] <- "NS"
      Vol[which(Vol$log2FC>0.5 & Vol$S!="NS"), "S"] <- "UP"
      Vol[which(Vol$log2FC<(-0.5) & Vol$S!="NS"), "S"] <- "DOWN"
      Vol[which(Vol$S=="0"), "S"] <- "NS"
      Vol$S <- as.factor(Vol$S)
      Vol$sig[Vol$sig == "Inf"] <- 320
          
    volcano <- ggplot(Vol, aes(x=log2FC, y=sig, color = S)) +
            geom_point (size = 1.75,
                        show.legend = F) +
            geom_vline(xintercept = c(-0.5, 0.5), linetype = 1, size = 0.3, col = "grey20") +
            geom_hline(yintercept = 2, linetype = 1, size = 0.3, col = "grey20") +
            scale_color_manual(values = c("skyblue2", "grey60", "green3")) +
            theme_light() +
            ggtitle(paste0("DE analysis: ", check[i], " vs ", check[j])) +
            xlab("avg_log2FC") + 
            ylab("-log10 pval_adj") +
            scale_y_continuous(trans = pseudolog10_trans) +
            scale_x_continuous(trans = pseudolog10_trans,
                           limits = symmetric_limits) +
            theme(axis.line = element_line(size = 0.3, colour = "grey20", linetype=1),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank())+
            geom_text_repel(data= subset(Vol, Vol$S != "NS"),
                            aes(label = ifelse((S != "NS"), gene, NA)),
                            color = "black",
                            size = 3,
                            force = 2,
                            vjust = 0.5, max.overlaps = 8)
