###Robyn Emmerson RNASeq 180822###

library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(rtracklayer)
library(ggplot2)

##using p<=0.05 and log2FC of 0.5- tested different parameters and settled on this

#load data table from preprocessing
my_RNASeq<- read.table("FLRnaSeq.tab", row.names=1)

#ensuring correct data labels and removing columns not needed
my_RNASeq<- my_RNASeq[, -c(1:5)]
colnames(my_RNASeq)<-gsub(".bam", "", colnames(my_RNASeq))
cData<- data.frame(lightRegime=c("FLH","FLH","FLH", "FLL", "FLL", "FLL", "SQH", "SQH", "SQH", "SQL", "SQL", "SQL"))
rownames(cData) <- colnames(my_RNASeq)

#create DESeq2 object from data for further use
ddsMat<- DESeqDataSetFromMatrix(countData=my_RNASeq, colData=cData, design= ~lightRegime)
dds<-DESeq(ddsMat)
res <- results(dds)
#saving just in case
save(dds, file=lightRnaSeq.Rdata)

hits_over <- res[which ((res$padj<0.05) & (res$log2FoldChange >=0.5)),]
hits_over[order(hits_over$padj)]

vsd<- vst(dds, blind=FALSE)
###SQHvSQL

res_SQHvSQL <- results(dds, contrast= c("lightRegime", "SQH", "SQL"))
SQhits_over <- res_SQHvSQL[which((res_SQHvSQL$padj<0.05) & (res_SQHvSQL$log2FoldChange >=0.5)),]
write.table(SQhits_over, "SQHvSQLRnaSeq.tab")

SQ_DGE<- read.table("SQHvSQLRnaSeq.tab", header=TRUE, rownames=1)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/SQ_volcano_enhanced.pdf", width=10, height=6, pointsize=12)
#simple volcano plot
plot(SQ_DGE$log2FoldChange, -log10(SQ_DGE$padj), xlab="log2 Fold Change", ylab="-log10(q)", pch=16, col="black")
abline(h=-log10(0.05), v=c(-0.5,0.5), lty='dashed')
points(SQ_DGE_sig$log2FoldChange, -log10(SQ_DGE_sig$padj), pch=16, col=cbbPalette[6])
#enhanced volcano plot
EnhancedVolcano(res_SQHvSQL,
    lab = rownames(res_SQHvSQL),
    x = 'log2FoldChange',
    y = 'pvalue',
    title = 'SQH vs SQL',
    pCutoff = 0.05,
    FCcutoff = 0.5,
    #labSize = 4.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)

dev.off()
#enchanced volcano reduced gene labels
pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/SQ_volcano_enhanced.pdf", width=10, height=10, pointsize=12)
EnhancedVolcano(res_SQHvSQL,
    lab = rownames(res_SQHvSQL),
    x = 'log2FoldChange',
    y = 'padj',
     xlim=c(-12,12),
    title = 'SQH vs SQL',
    selectLab= NA,
    pCutoff = 0.05,
    FCcutoff = 0.5,
    #labSize = 4.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)



dev.off()

##FLHvFLL
res_FLHvFLL <- results(dds, contrast= c("lightRegime", "FLH", "FLL"))
FLhits_over <- res_FLHvFLL[which((res_FLHvFLL$padj<0.05) & (res_FLHvFLL$log2FoldChange >=0.5)),]
write.table(FLhits_over, "FLHvFLLRnaSeq.tab")

FL_DGE<- read.table("FLHvFLLRnaSeq.tab", header=TRUE, rownames=1)


##PCA

pcaData<- plotPCA(vsd, intgroup="lightRegime")
percentVar<- round(100*attr(pcaData, "percentVar"))
pcaData$lightRegime <- factor(pcaData$lightRegime , levels= c("SQH","SQL","FLH", "FLL"))

pcaPlot<- ggplot(pcaData, aes(PC1, PC2, col=lightRegime))+geom_point(size=3) + guides(fill="none") +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
coord_fixed() + scale_colour_manual(values = c( "dimgray","darkgoldenrod1", "cyan4","darkorange2")) +
theme(panel.background=element_blank(),axis.line=element_line(colour="black"))
+ theme(panel.grid.major = element_line(color = "gray48",linetype = 2))
+ theme(text=element_text(size=16))+
+ theme(legend.title=element_blank())



cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/FL_volcano.pdf", width=10, height=6, pointsize=12)
#simple volcano plot
plot(FL_DGE$log2FoldChange, -log10(FL_DGE$padj), xlab="log2 Fold Change", ylab="-log10(q)", pch=16, col="black")
abline(h=-log10(0.05), v=c(-0.5,0.5), lty='dashed')
points(FL_DGE_sig$log2FoldChange, -log10(FL_DGE_sig$padj), pch=16, col=cbbPalette[7])
#enhanced volcano plot

EnhancedVolcano(res_FLHvFLL,
    lab = rownames(res_FLHvFLL),
    x = 'log2FoldChange',
    y = 'padj',
    title = 'FLH vs FLL',
    pCutoff = 0.05,
    FCcutoff = 0.5,
     selectLab= NA,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)

dev.off()

pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/FL_volcano_enhanced.pdf", width=10, height=10, pointsize=12)


#enchanced volcano reduced gene labels
EnhancedVolcano(res_FLHvFLL,
    lab = rownames(res_FLHvFLL),
    x = 'log2FoldChange',
    y = 'padj',
     xlim=c(-12,12),
    title = 'FLH vs FLL',
     
    pCutoff = 0.05,
    FCcutoff = 0.5,
     selectLab= NA,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)

dev.off()


##SQHvFLH

res_SQHvFLH <- results(dds, contrast= c("lightRegime", "SQH", "FLH"))
HLhits_over <- res_SQHvFLH[which((res_SQHvFLH$padj<0.05) & (res_SQHvFLH$log2FoldChange >=0.5)),]
write.table(HLhits_over, "SQHvFLHRnaSeq.tab")

HL_DGE<- read.table("SQHvFLHRnaSeq.tab", header=TRUE, rownames=1)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/HL_volcano.pdf", width=15, height=6, pointsize=12)
#simple volcano plot
plot(HL_DGE$log2FoldChange, -log10(HL_DGE$padj), xlab="log2 Fold Change", ylab="-log10(q)", pch=16, col="black")
abline(h=-log10(0.05), v=c(-0.5,0.5), lty='dashed')
points(HL_DGE_sig$log2FoldChange, -log10(HL_DGE_sig$padj), pch=16, col=cbbPalette[3])

pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/HL_volcano_enhanced.pdf", width=16, height=6, pointsize=12)

#enhanced volcano plot
EnhancedVolcano(res_SQHvFLH,
    lab = rownames(res_SQHvFLH),
    x = 'log2FoldChange',
    y = 'padj',
    title = 'SQH vs FLH',
    pCutoff = 0.05,
    FCcutoff = 0.5,
    #labSize = 4.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)

dev.off()

pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/HL_volcano_enhanced.pdf", width=10, height=10, pointsize=12)


#enchanced volcano reduced gene labels
EnhancedVolcano(res_SQHvFLH,
    lab = rownames(res_SQHvFLH),
    x = 'log2FoldChange',
    y = 'pvalue',
     xlim=c(-12,12),
    title = 'SQH vs FLH',
 selectLab= NA,
    pCutoff = 0.05,
    FCcutoff = 0.5,
    #labSize = 4.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)

dev.off()


#SQLvFLL

res_SQLvFLL <- results(dds, contrast= c("lightRegime", "SQL", "FLL"))
LLhits_over <- res_SQLvFLL[which((res_SQLvFLL$padj<0.05) & (res_SQLvFLL$log2FoldChange >=0.5)),]
write.table(LLhits_over, "SQHvFLHRnaSeq.tab")

LL_DGE<- read.table("SQHvFLHRnaSeq.tab", header=TRUE, rownames=1)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/LL_volcano.pdf", width=15, height=6, pointsize=12)
#simple volcano plot
plot(LL_DGE$log2FoldChange, -log10(LL_DGE$padj), xlab="log2 Fold Change", ylab="-log10(q)", pch=16, col="black")
abline(h=-log10(0.05), v=c(-0.5,0.5), lty='dashed')
points(LL_DGE_sig$log2FoldChange, -log10(LL_DGE_sig$padj), pch=16, col=cbbPalette[2])
#enhanced volcano plot
pdf("/storage/projects/ZabetLab/BS_seq_TAIR10/Robyn_data/RNAseq/LL_volcano_enhanced.pdf", width=10, height=10, pointsize=12)

EnhancedVolcano(res_SQLvFLL,
    lab = rownames(res_SQLvFLL),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim=c(-12,12),
    title = 'SQL vs FLL',
    selectLab= NA,
    pCutoff = 0.05,
    FCcutoff = 0.5,
    #labSize = 4.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)

dev.off()



#enchanced volcano reduced gene labels
EnhancedVolcano(res_SQLvFLL,
    lab = rownames(res_SQLvFLL),
    x = 'log2FoldChange',
    y = 'pvalue',
    xlim=c(-12,12),
    title = 'SQL vs FLL',
     selectLab = c('AT2G41260','AT5G42900', 'AT4G27130', 'AT1G23205', 'AT3G15680', 'AT1G35290'),
    pCutoff = 0.05,
    FCcutoff = 0.5,
    #labSize = 4.0,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    legendPosition = 'top',
   legendLabSize = 10,
    legendIconSize = 3.0
)

dev.off()




#simple heatmap
sampleDists<-dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
heatmap(sampleDistMatrix)

#better heatmap

pheatmap(sampleDistMatrix)


res<-results(dds)
hits1 <- res[which((res$padj<0.05) & (abs(res$log2FoldChange) >=0.5)),]
hits2 <- res[which((res$padj<0.05) & (abs(res$log2FoldChange) <=0.5)),]
hits<- rbind(hits1, hits2)
hits_genes<-rownames(hits)
vsd_hits <- vsd[hits_genes, ]
pheatmap(as.matrix(assay(vsd_hits)))
hit_genes_sorted <- rownames(hits[order(hits$log2FoldChange),])
vsd_hits_sorted <- vsd[hit_genes_sorted, ]
pheatmap(as.matrix(assay(vsd_hits_sorted)), cluster_cols=FALSE, cluster_rows=FALSE)

SQ_DGE_sig <- SQ_DGE[which(SQ_DGE$padj<=0.05 & abs(SQ_DGE$log2FoldChange)>=0.5),]

FL_DGE_sig <- FL_DGE[which(FL_DGE$padj<=0.05 & abs(FL_DGE$log2FoldChange)>=0.5),]

HL_DGE_sig <- HL_DGE[which(HL_DGE$padj<=0.05 & abs(HL_DGE$log2FoldChange)>=0.5),]

LL_DGE_sig <- LL_DGE[which(LL_DGE$padj<=0.05 & abs(LL_DGE$log2FoldChange)>=0.5),]


annotation<-import.gff3("TAIR10_GFF3_genes_transposons.gff")
ranges<-ranges(annotation)


