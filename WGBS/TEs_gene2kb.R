#plot number DM TE, TE with at least 1 gene in 5kb
#number of TEs with a DEG in 5kb

#Differentially methylated TEs
SQ_TE_count2<-n_distinct(SQ_TEs_ranges_2$TEgene)
FL_TE_count2<-n_distinct(FL_TEs_ranges_2$TEgene)
HL_TE_count2<-n_distinct(HL_TEs_ranges_2$TEgene)
LL_TE_count2<-n_distinct(LL_TEs_ranges_2$TEgene)

#TEs with gene in 5kb
SQ_TE_g2<-n_distinct(SQ_TE_2kbgenes_df$TEgene)
FL_TE_g2<-n_distinct(FL_TE_2kbgenes_df$TEgene)
HL_TE_g2<-n_distinct(HL_TE_2kbgenes_df$TEgene)
LL_TE_g2<-n_distinct(LL_TE_2kbgenes_df$TEgene)

#TEs with DEG in 5kb

SQ_d2<-merge(SQ_DEGs,SQ_TE_2kbgenes_df, by.x="X", by.y="ID")
SQ_degs2<-n_distinct(SQ_d2$X)
SQ_DEG_2<-n_distinct(SQ_d2$TEgene)
FL_d2<-merge(FL_DEGs,FL_TE_2kbgenes_df, by.x="X", by.y="ID")
FL_degs2<-n_distinct(FL_d2$X)
FL_DEG_2<-n_distinct(FL_d2$TEgene)
HL_d2<-merge(HL_DEGs,HL_TE_2kbgenes_df, by.x="X", by.y="ID")
HL_DEG_2<-n_distinct(HL_d2$TEgene)
HL_degs2<-n_distinct(HL_d2$X)
LL_d2<-merge(LL_DEGs,LL_TE_2kbgenes_df, by.x="X", by.y="ID")
LL_DEG_2<-n_distinct(LL_d2$TEgene)
LL_degs2<-n_distinct(LL_d2$X)

plot(SQ_d2$log2FoldChange, SQ_d2$methylLevel)
plot(FL_d2$log2FoldChange, FL_d2$methylLevel)
plot(HL_d2$log2FoldChange, HL_d2$methylLevel)
plot(LL_d2$log2FoldChange, LL_d2$methylLevel)

list_TE_2_plot<-list(
  "SQHvSQL, number of DM TEs"=c(SQ_TE_count2),
  "SQHvSQL, number of TEs with gene within 2kb"=c(SQ_TE_g2),
  "SQHvSQL, number of TEs with at least one DEG"=c(SQ_DEG_2),
  "FLHvFLL, number of DM TEs"=c(FL_TE_count2),
  "FLHvFLL, number of TEs with gene within 2kb"=c(FL_TE_g2),
  "FLHvFLL, number of TEs with at least one DEG"=c(FL_DEG_2),
  "SQHvFLH, number of DM TEs"=c(HL_TE_count2),
  "SQHvFLH, number of TEs with gene within 2kb"=c(HL_TE_g2),
  "SQHvFLH, number of TEs with at least one DEG"=c(HL_DEG_2),
  "SQLvFLL, number of DM TEs"=c(LL_TE_count2),
  "SQLvFLL, number of TEs with gene within 2kb"=c(LL_TE_g2),
  "SQLvFLL, number of TEs with at least one DEG"=c(LL_DEG_2)
)

sort_list_TE_2_plot<-list(
  "SQHvSQL, number of DM TEs"=c(SQ_TE_count2),
  "FLHvFLL, number of DM TEs"=c(FL_TE_count2),
  "SQHvFLH, number of DM TEs"=c(HL_TE_count2),
  "SQLvFLL, number of DM TEs"=c(LL_TE_count2),
  "SQHvSQL, number of TEs with gene within 2kb"=c(SQ_TE_g2),
  "FLHvFLL, number of TEs with gene within 2kb"=c(FL_TE_g2),
  "SQHvFLH, number of TEs with gene within 2kb"=c(HL_TE_g2),
  "SQLvFLL, number of TEs with gene within 2kb"=c(LL_TE_g2),
  "SQHvSQL, number of TEs with at least one DEG"=c(SQ_DEG_2),
  "FLHvFLL, number of TEs with at least one DEG"=c(FL_DEG_2),
  "SQHvFLH, number of TEs with at least one DEG"=c(HL_DEG_2),
  "SQLvFLL, number of TEs with at least one DEG"=c(LL_DEG_2)
)

No_DM_TEs<-list(
  "SQHvSQL, DM TEs"=c(SQ_TE_count2),
  "FLHvFLL, DM TEs"=c(FL_TE_count2),
  "SQHvFLH, DM TEs"=c(HL_TE_count2),
  "SQLvFLL, DM TEs"=c(LL_TE_count2)
)

No_TEs_2kb<-list(
  "SQHvSQL,  DM TEs with gene within 2kb"=c(SQ_TE_g2),
  "FLHvFLL,  DM TEs with gene within 2kb"=c(FL_TE_g2),
  "SQHvFLH, DM TEs with gene within 2kb"=c(HL_TE_g2),
  "SQLvFLL,  DM TEs with gene within 2kb"=c(LL_TE_g2)
)

No_TEs_DEG_2kb<-list(
  "SQHvSQL, DM TEs with at least one DEG"=c(SQ_DEG_2),
  "FLHvFLL, DM TEs with at least one DEG"=c(FL_DEG_2),
  "SQHvFLH, DM TEs with at least one DEG"=c(HL_DEG_2),
  "SQLvFLL, DM TEs with at least one DEG"=c(LL_DEG_2)
)

No_DEGs<-list(
  "SQHvSQL, DEGs within 2kb TE"=c(SQ_degs),
  "FLHvFLL, DEGs within 2kb TE"=c(FL_degs),
  "SQHvFLH, DEGs within 2kb TE"=c(HL_degs),
  "SQLvFLL, DEGs within 2kb TE"=c(LL_degs)
)



colPal<- c("#785EF0","#DC627F","#FE6100","#FFB000" )

pdf("TE_2kb_counts_grouped.pdf", height=7, width=10, pointsize=12)
par(las=2)
par(mar=c(14,5,4,1)+1)
barplot(unlist(sort_list_TE_2_plot), col=colPal[c(1,2,3,4)], ylim=c(0,3500))
dev.off()

colPal2<- c("#785EF0","#9D8CE8","#C0B7EA","#DC627F","#E24F97","#D278A4","#FE6100","#EC8342","#E29464","#FFB000","#ECBA4A","#ECC263" )

pdf("TE_2kb_counts_byRegime.pdf", height=7, width=10, pointsize=12)
par(las=2)
par(mar=c(14,5,4,1)+1)
barplot(unlist(list_TE_2_plot), col=colPal2[c(1,2,3,4,5,6,7,8,9,10,11,12)], ylim=c(0,3500))
dev.off()

pdf("TE_2kb_grouped_separated.pdf")
par(las=2)
par(mar=c(14,5,4,1)+1)
barplot(unlist(No_DM_TEs), col=colPal)
barplot(unlist(No_TEs_2kb), col=colPal)
barplot(unlist(No_TEs_DEG_2kb), col=colPal)
barplot(unlist(No_DEGs),col=colPal)
dev.off()

pdf("DEGs_w_TE_2kb_grouped.pdf")
par(las=2)
par(mar=c(14,5,4,1)+1)
barplot(unlist(No_DEGs),col=colPal)
dev.off()

SQ_DEGs_2kbTE<-as.data.frame(SQ_d)
FL_DEGs_2kbTE<-as.data.frame(FL_d)
HL_DEGs_2kbTE<-as.data.frame(HL_d)
LL_DEGs_2kbTE<-as.data.frame(LL_d)

write.csv(SQ_DEGs_2kbTE, file="SQ_DEGs_2kb_TE.csv")
write.csv(FL_DEGs_2kbTE, file="FL_DEGs_2kb_TE.csv")
write.csv(HL_DEGs_2kbTE, file="HL_DEGs_2kb_TE.csv")
write.csv(LL_DEGs_2kbTE, file="LL_DEGs_2kb_TE.csv")