
# li is a rowbinded dataframe with all comparisons (used for HB plot)
# comparison col would be TFH_vs_TR1.2 etc
# comparison the type of pathway: cell cycle, mitosis, ribosomal...

fc <- 0.25 # log2FC threeshold
lol <- table(li2$avg_log2FC>=fc, li2$comparison, li2$pathway_type) %>%
  as.data.frame() %>% filter(Var1 == TRUE) %>% rename(Var1="UP")
lol2 <- table(li2$avg_log2FC <= -fc, li2$comparison, li2$pathway_type) %>%
  as.data.frame() %>% filter(Var1 == TRUE) %>% rename(Var1="DOWN")
lol3 <- rbind(lol2, lol)%>% arrange(Var3, desc(Var2))

#dataframe for positions
pt <- c()
l <- length(unique(lol3$comparison))
mi <- min(li2)*1.1
max <- max(li2)*1.1
for(i in unique(lol3$comparison)){
  pt <- append(pt,c(rep(i,l)))
}
dat_text <- data.frame(
  label= as.character(lol3$Freq),
  pathway_type= pt,
  x = rep(c(mi,ma),l),
  y = rep(seq_len(l), each=2),
  comparison = lol3$Var2
)

hb2 <- hb + geom_text(
                      data    = dat_text,
                      mapping = aes(x = x, y = y,
                                    label = label, fontface=2),
                      size= 8.5,
                      show.legend = F
                      )

