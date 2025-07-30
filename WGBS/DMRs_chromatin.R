##data from Jamge et al., 2023 https://doi.org/10.7554/eLife.87714.3 

At_chromatin <- read.table("C:/Users/biol0346/Downloads/AT_26_ChromHMM_ChromHMM_July_AT_segments.bed", header=FALSE, sep="\t")

TE_gff<- readGFFAsGRanges("C:/Users/biol0346/Downloads/Araport11_GFF3_genes_transposons.Jun2016.gff/Araport11_GFF3_genes_transposons.201606.gff")
seqlevelsStyle(TE_gff)<-"Ensembl"
#convert to GRanges
colnames(At_chromatin)<- c("seqnames", "start", "end", "chromatin")
At_chromatin_gr<- makeGRangesFromDataFrame(At_chromatin, keep.extra.columns = TRUE)


##assign chromatin state to TE

TE_eu<- findOverlaps(TE_gff, At_chromatin_gr)
TE_gff$chromatin1<- ""
TE_gff$chromatin1[queryHits(TE_eu)]<-At_chromatin_gr$chromatin[subjectHits(TE_eu)]

map_to_paper = data.frame(
  from=c("E23","E25","E24","E22","E21","E26",
         "E1","E3","E2","E4","E5","E6",
         "E18","E19","E20",
         "E7","E8","E9","E10","E11","E12","E13","E14","E17","E15","E16"),
  to=factor(c("H1","H2","H3","H4","H5","H6",
              "F1","F2","F3","F4","F5","F6",
              "I1","I2","I3",
              "E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11"),
            levels = c("E11","E10","E9","E8","E7","E6","E5","E4","E3","E2","E1",
                       "I3","I2","I1",
                       "F6","F5","F4","F3","F2","F1",
                       "H6","H5","H4","H3","H2","H1")))

##get emission files chromatin code
mcols_df <- as.data.frame(mcols(TE_gff))

# Add rownames to keep track of ordering if needed
mcols_df$row_id <- seq_len(nrow(mcols_df))

# Merge data- map_to_paper (contains chromatin code) by chromatin1
merged_df <- merge(mcols_df, map_to_paper, by.x = "chromatin1", by.y = "from", all.x = TRUE)

# Restore original order 
merged_df <- merged_df[order(merged_df$row_id), ]
merged_df$row_id <- NULL  # drop

# Assign the merged metadata back to GRanges
mcols(TE_gff) <- merged_df

# Rename new column if needed
names(mcols(TE_gff))[names(mcols(TE_gff)) == "to"] <- "chromatin2"

#Add new column specifying chromatin
#heterochromatin (H1-6), facultative heterochromatin (F1-6), euchromatin (E1-11), and intergenic regions (I1-3)

classify_chromatin_state <- function(state_code) {
  if (grepl("^H[1-6]$", state_code)) {
    return("heterochromatin")
  } else if (grepl("^F[1-6]$", state_code)) {
    return("facultative heterochromatin")
  } else if (grepl("^E([1-9]|10|11)$", state_code)) {
    return("euchromatin")
  } else if (grepl("^I[1-3]$", state_code)) {
    return("intergenic regions")
  } else {
    return(NA_character_)  
  }
}

#apply to TE_gff$chromatin 2 where the state is defined
mcols(TE_gff)$chromatin_state <- vapply(
  mcols(TE_gff)$chromatin2,
  classify_chromatin_state,
  FUN.VALUE = character(1)
)

#returns NAs
table(mcols(TE_gff)$chromatin_state, useNA = "ifany")


##subset into TEs and genes
gene_chromatin <- TE_gff[mcols(TE_gff)$type == "gene"]

TE_chromatin<- TE_gff [mcols(TE_gff)$type =="transposable_element" | mcols(TE_gff)$type =="transposable_element_gene"]

##methylation changes in heterochromatic and euchromatic TEs and genes- PLOT EACH CONTEXT

#read in DMRs
library(readxl)
all_DMRs<-read_xlsx("D:/PhD/DNA METHYLATION PAPER- DATA SUBMISSION/All_DMRs.xlsx")
all_DMRs_gr<- makeGRangesFromDataFrame(all_DMRs, keep.extra.columns = TRUE)
all_DMRs_gr$proportionDifference<- all_DMRs_gr$proportion2 - all_DMRs_gr$proportion1
seqlevelsStyle(all_DMRs_gr)<-"NCBI"
seqlevelsStyle(TE_gff)<-"NCBI"


#get related chromatin state
allDMR_chromatin<- findOverlaps(all_DMRs_gr, TE_gff)
all_DMRs_gr$chromatin<- ""
all_DMRs_gr$chromatin[queryHits(allDMR_chromatin)]<-TE_gff$chromatin_state[subjectHits(allDMR_chromatin)]

#subset by regime
SQ_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$Regime == "SQHvSQL"]
FL_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$Regime == "FLHvFLL"]
HL_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$Regime == "SQHvFLH"]
LL_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$Regime == "SQLvFLL"]

#subset by chromatin status
eu_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$chromatin == "euchromatin"]
# make sure in the order i want
dmr_df$Regime <- factor(dmr_df$Regime, levels = c("SQHvSL", "FLHvFLL", "SQHvFLH", "SQLvFLL"))

het_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$chromatin == "heterochromatin"]
fachet_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$chromatin == "facultative heterochromatin"]
intergenic_DMRs<- all_DMRs_gr[mcols(all_DMRs_gr)$chromatin == "intergenic regions"]

##boxplot for each context and regime
col_regime<-  c(
  "SQHvSQL" = "cyan4",
  "FLHvFLL" = "darkorange2",
  "SQHvFLH" = "dimgray",
  "SQLvFLL" = "darkgoldenrod1"
)
  


library(ggplot2)
euDMR_df<- as.data.frame(eu_DMRs)
euDMR_df$Regime <- factor(euDMR_df$Regime, levels = c("SQHvSQL", "FLHvFLL", "SQHvFLH", "SQLvFLL"))

ggplot(euDMR_df, aes(x = context, y = proportionDifference, fill = Regime)) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(values = col_regime) +
  labs(title = "DMR Proportion Difference- Euchromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()

hetDMR_df<- as.data.frame(het_DMRs)
hetDMR_df$Regime <- factor(hetDMR_df$Regime, levels = c("SQHvSQL", "FLHvFLL", "SQHvFLH", "SQLvFLL"))

ggplot(hetDMR_df, aes(x = context, y = proportionDifference, fill = Regime)) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(values = col_regime) +
  labs(title = "DMR Proportion Difference- Constitutive Heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()

fachetDMR_df<- as.data.frame(fachet_DMRs)
fachetDMR_df$Regime <- factor(fachetDMR_df$Regime, levels = c("SQHvSQL", "FLHvFLL", "SQHvFLH", "SQLvFLL"))

ggplot(fachetDMR_df, aes(x = context, y = proportionDifference, fill = Regime)) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(values = col_regime) +
  labs(title = "DMR Proportion Difference- Facultative Heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()

igDMR_df<- as.data.frame(intergenic_DMRs)
igDMR_df$Regime <- factor(igDMR_df$Regime, levels = c("SQHvSQL", "FLHvFLL", "SQHvFLH", "SQLvFLL"))

ggplot(igDMR_df, aes(x = context, y = proportionDifference, fill = Regime)) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(values = col_regime) +
  labs(title = "DMR Proportion Difference- Intergenic Regions",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()


library(ggpubr)


# Create the boxplot with stats
ggplot(euDMR_df, aes(x = Regime, y = proportionDifference, fill = Regime)) +
  geom_boxplot() +
  facet_wrap(~ context, scales = "free_y") +
  scale_fill_manual(values = col_regime) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("SQHvSQL", "FLHvFLL"),
                       c("SQHvSQL", "SQHvFLH"),
                       c("SQHvSQL", "SQLvFLL"),
                       c("FLHvFLL", "SQHvFLH"),
                       c("FLHvFLL", "SQLvFLL"),
                       c("SQHvFLH", "SQLvFLL")
                     ),
                     label = "p.signif",
                     hide.ns = TRUE,
                     method.args = list(alternative = "two.sided")) +
  labs(title = "Proportion Difference by Regime and Context with Statistics",
       x = "Light Regime Comparison",
       y = "Proportion Difference") +
  theme_minimal()

ggplot(hetDMR_df, aes(x = Regime, y = proportionDifference, fill = Regime)) +
  geom_boxplot() +
  facet_wrap(~ context, scales = "free_y") +
  scale_fill_manual(values = col_regime) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("SQHvSQL", "FLHvFLL"),
                       c("SQHvSQL", "SQHvFLH"),
                       c("SQHvSQL", "SQLvFLL"),
                       c("FLHvFLL", "SQHvFLH"),
                       c("FLHvFLL", "SQLvFLL"),
                       c("SQHvFLH", "SQLvFLL")
                     ),
                     label = "p.signif",
                     hide.ns = TRUE,
                     method.args = list(alternative = "two.sided")) +
  labs(title = "Proportion Difference by Regime and Context with Statistics",
       x = "Light Regime Comparison",
       y = "Proportion Difference") +
  theme_minimal()

ggplot(fachetDMR_df, aes(x = Regime, y = proportionDifference, fill = Regime)) +
  geom_boxplot() +
  facet_wrap(~ context, scales = "free_y") +
  scale_fill_manual(values = col_regime) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("SQHvSQL", "FLHvFLL"),
                       c("SQHvSQL", "SQHvFLH"),
                       c("SQHvSQL", "SQLvFLL"),
                       c("FLHvFLL", "SQHvFLH"),
                       c("FLHvFLL", "SQLvFLL"),
                       c("SQHvFLH", "SQLvFLL")
                     ),
                     label = "p.signif",
                     hide.ns = TRUE,
                     method.args = list(alternative = "two.sided")) +
  labs(title = "Proportion Difference by Regime and Context with Statistics",
       x = "Light Regime Comparison",
       y = "Proportion Difference") +
  theme_minimal()
ggplot(igDMR_df, aes(x = Regime, y = proportionDifference, fill = Regime)) +
  geom_boxplot() +
  facet_wrap(~ context, scales = "free_y") +
  scale_fill_manual(values = col_regime) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("SQHvSQL", "FLHvFLL"),
                       c("SQHvSQL", "SQHvFLH"),
                       c("SQHvSQL", "SQLvFLL"),
                       c("FLHvFLL", "SQHvFLH"),
                       c("FLHvFLL", "SQLvFLL"),
                       c("SQHvFLH", "SQLvFLL")
                     ),
                     label = "p.signif",
                     hide.ns = TRUE,
                     method.args = list(alternative = "two.sided")) +
  labs(title = "Proportion Difference by Regime and Context with Statistics",
       x = "Light Regime Comparison",
       y = "Proportion Difference") +
  theme_minimal()



##separate out into genes and TEs?##


##correlation between gene expression and hetero/euchromatic DMR TEs
SQ_DEG_TE5<-read.csv("D:/PhD/DATA IN THESIS/DNA methylation/RAW MEHTYLATION/SQ_DEGs_5kb_TE.csv", header=TRUE)
FL_DEG_TE5<-read.csv("D:/PhD/DATA IN THESIS/DNA methylation/RAW MEHTYLATION/FL_DEGs_5kb_TE.csv", header=TRUE)
HL_DEG_TE5<-read.csv("D:/PhD/DATA IN THESIS/DNA methylation/RAW MEHTYLATION/HL_DEGs_5kb_TE.csv", header=TRUE)
LL_DEG_TE5<-read.csv("D:/PhD/DATA IN THESIS/DNA methylation/RAW MEHTYLATION/LL_DEGs_5kb_TE.csv", header=TRUE)
SQ_DEGS_TE5_gr<-makeGRangesFromDataFrame(SQ_DEG_TE5, keep.extra.columns = TRUE)
FL_DEGS_TE5_gr<-makeGRangesFromDataFrame(FL_DEG_TE5, keep.extra.columns = TRUE)
HL_DEGS_TE5_gr<-makeGRangesFromDataFrame(HL_DEG_TE5, keep.extra.columns = TRUE)
LL_DEGS_TE5_gr<-makeGRangesFromDataFrame(LL_DEG_TE5, keep.extra.columns = TRUE)
seqlevelsStyle(SQ_DEGS_TE5_gr)<-("NCBI")
seqlevelsStyle(FL_DEGS_TE5_gr)<-("NCBI")
seqlevelsStyle(HL_DEGS_TE5_gr)<-("NCBI")
seqlevelsStyle(LL_DEGS_TE5_gr)<-("NCBI")

# Convert to data.frame for easy matching
SQ_deg_df <- as.data.frame(SQ_DEGS_TE5_gr)
SQ_dmr_df <- as.data.frame(SQ_DMRs)
# Match by TEgene
SQ_idx <- match(SQ_deg_df$TEgene, SQ_dmr_df$TEgene)
# Assign values from DMRs to GRanges metadata
SQ_DEGS_TE5_gr$chromatin <- SQ_dmr_df$chromatin[SQ_idx]
SQ_DEGS_TE5_gr$methylationProportion <- SQ_dmr_df$proportionDifference[SQ_idx]
SQ_DEGS_TE5_gr
table(SQ_DEGS_TE5_gr$chromatin)
#separate by chromatin
SQ_fachet<-all_DMRs_gr[mcols(all_DMRs_gr)$chromatin == "facultative heterochromatin"]
SQ_fachet
#separate by chromatin
SQ_fachet<-SQ_DEGS_TE5_gr[mcols(SQ_DEGS_TE5_gr)$chromatin == "facultative heterochromatin"]
# Safe subsetting: exclude NAs first
SQ_fachet <- SQ_DEGS_TE5_gr[
  !is.na(mcols(SQ_DEGS_TE5_gr)$chromatin) &
    mcols(SQ_DEGS_TE5_gr)$chromatin == "facultative heterochromatin"
]
ggplot(as.data.frame(SQ_fachet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(SQ_fachet$log2FoldChange), label.y = max(SQ_fachet$methylationProportion)) +
  labs(title = "SQHvSQL- DEG with 5kb DM TE- facultative heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
SQ_conhet <- SQ_DEGS_TE5_gr[
  !is.na(mcols(SQ_DEGS_TE5_gr)$chromatin) &
    mcols(SQ_DEGS_TE5_gr)$chromatin == "heterochromatin"
]
ggplot(as.data.frame(SQ_conhet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(SQ_conhet$log2FoldChange), label.y = max(SQ_conhet$methylationProportion)) +
  labs(title = "SQHvSQL- DEG with 5kb DM TE- constitutive heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
SQ_igen <- SQ_DEGS_TE5_gr[
  !is.na(mcols(SQ_DEGS_TE5_gr)$chromatin) &
    mcols(SQ_DEGS_TE5_gr)$chromatin == "intergenic regions"
]
ggplot(as.data.frame(SQ_igen), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(SQ_igen$log2FoldChange), label.y = max(SQ_igen$methylationProportion)) +
  labs(title = "SQHvSQL- DEG with 5kb DM TE- constitutive heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
FL_deg_df <- as.data.frame(FL_DEGS_TE5_gr)
FL_dmr_df <- as.data.frame(FL_DMRs)
HL_deg_df <- as.data.frame(HL_DEGS_TE5_gr)
HL_dmr_df <- as.data.frame(HL_DMRs)
LL_deg_df <- as.data.frame(LL_DEGS_TE5_gr)
LL_dmr_df <- as.data.frame(LL_DMRs)
#match by TEgene
SQ_idx <- match(SQ_deg_df$TEgene, SQ_dmr_df$TEgene)
FL_idx <- match(FL_deg_df$TEgene, FL_dmr_df$TEgene)
HL_idx <- match(HL_deg_df$TEgene, HL_dmr_df$TEgene)
LL_idx <- match(LL_deg_df$TEgene, LL_dmr_df$TEgene)
FL_DEGS_TE5_gr$chromatin <- FL_dmr_df$chromatin[FL_idx]
FL_DEGS_TE5_gr$methylationProportion <- FL_dmr_df$proportionDifference[FL_idx]
HL_DEGS_TE5_gr$chromatin <- HL_dmr_df$chromatin[HL_idx]
HL_DEGS_TE5_gr$methylationProportion <- HL_dmr_df$proportionDifference[HL_idx]
LL_DEGS_TE5_gr$chromatin <- LL_dmr_df$chromatin[LL_idx]
LL_DEGS_TE5_gr$methylationProportion <- LL_dmr_df$proportionDifference[LL_idx]
FL_fachet <- FL_DEGS_TE5_gr[
  !is.na(mcols(FL_DEGS_TE5_gr)$chromatin) &
    mcols(FL_DEGS_TE5_gr)$chromatin == "facultative heterochromatin"
]
FL_conhet <- FL_DEGS_TE5_gr[
  !is.na(mcols(FL_DEGS_TE5_gr)$chromatin) &
    mcols(FL_DEGS_TE5_gr)$chromatin == "heterochromatin"
]
FL_igen <- FL_DEGS_TE5_gr[
  !is.na(mcols(FL_DEGS_TE5_gr)$chromatin) &
    mcols(FL_DEGS_TE5_gr)$chromatin == "intergenic regions"
]
HL_fachet <- HL_DEGS_TE5_gr[
  !is.na(mcols(HL_DEGS_TE5_gr)$chromatin) &
    mcols(HL_DEGS_TE5_gr)$chromatin == "facultative heterochromatin"
]
HL_conhet <- HL_DEGS_TE5_gr[
  !is.na(mcols(HL_DEGS_TE5_gr)$chromatin) &
    mcols(HL_DEGS_TE5_gr)$chromatin == "heterochromatin"
]
HL_igen <- HL_DEGS_TE5_gr[
  !is.na(mcols(HL_DEGS_TE5_gr)$chromatin) &
    mcols(HL_DEGS_TE5_gr)$chromatin == "intergenic regions"
]
LL_fachet <- LL_DEGS_TE5_gr[
  !is.na(mcols(LL_DEGS_TE5_gr)$chromatin) &
    mcols(LL_DEGS_TE5_gr)$chromatin == "facultative heterochromatin"
]
LL_conhet <- LL_DEGS_TE5_gr[
  !is.na(mcols(LL_DEGS_TE5_gr)$chromatin) &
    mcols(LL_DEGS_TE5_gr)$chromatin == "heterochromatin"
]
LL_igen <- LL_DEGS_TE5_gr[
  !is.na(mcols(LL_DEGS_TE5_gr)$chromatin) &
    mcols(LL_DEGS_TE5_gr)$chromatin == "intergenic regions"
]




#plot
pdf("D:/PhD/DNA METHYLATION PAPER- DATA SUBMISSION/DEG_DMTE_5kb_facheterochromatin.pdf", width=8, height=5)
ggplot(as.data.frame(SQ_fachet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(SQ_fachet$log2FoldChange), label.y = max(SQ_fachet$methylationProportion)) +
  labs(title = "SQHvSQL- DEG with 5kb DM TE- facultative heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
ggplot(as.data.frame(FL_fachet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(FL_fachet$log2FoldChange), label.y = max(FL_fachet$methylationProportion)) +
  labs(title = "FLHvFLL- DEG with 5kb DM TE- facultative heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
ggplot(as.data.frame(HL_fachet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(HL_fachet$log2FoldChange), label.y = max(HL_fachet$methylationProportion)) +
  labs(title = "SQHvFLH- DEG with 5kb DM TE- facultative heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
ggplot(as.data.frame(LL_fachet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(LL_fachet$log2FoldChange), label.y = max(LL_fachet$methylationProportion)) +
  labs(title = "SQLvFLL- DEG with 5kb DM TE- facultative heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
dev.off()

LL_eu <- LL_DEGS_TE5_gr[
  !is.na(mcols(LL_DEGS_TE5_gr)$chromatin) &
    mcols(LL_DEGS_TE5_gr)$chromatin == "euchromatin"
]
HL_eu <- HL_DEGS_TE5_gr[
  !is.na(mcols(HL_DEGS_TE5_gr)$chromatin) &
    mcols(HL_DEGS_TE5_gr)$chromatin == "euchromatin"
]
FL_eu <- FL_DEGS_TE5_gr[
  !is.na(mcols(FL_DEGS_TE5_gr)$chromatin) &
    mcols(FL_DEGS_TE5_gr)$chromatin == "euchromatin"
]

pdf("D:/PhD/DNA METHYLATION PAPER- DATA SUBMISSION/DEG_DMTE_5kb_conheterochromatin.pdf", width=8, height=5)
ggplot(as.data.frame(SQ_conhet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(SQ_conhet$log2FoldChange), label.y = max(SQ_conhet$methylationProportion)) +
  labs(title = "SQHvSQL- DEG with 5kb DM TE- constitutive heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
ggplot(as.data.frame(FL_conhet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(FL_conhet$log2FoldChange), label.y = max(FL_conhet$methylationProportion)) +
  labs(title = "FLHvFLL- DEG with 5kb DM TE- constitutive heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
ggplot(as.data.frame(HL_conhet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(HL_conhet$log2FoldChange), label.y = max(HL_conhet$methylationProportion)) +
  labs(title = "SQHvFLH- DEG with 5kb DM TE- constitutive heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
ggplot(as.data.frame(LL_conhet), aes(x=log2FoldChange, y=methylationProportion))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  stat_cor(method = "pearson", label.x = min(LL_conhet$log2FoldChange), label.y = max(LL_conhet$methylationProportion)) +
  labs(title = "SQLvFLL- DEG with 5kb DM TE- constitutive heterochromatin",
       x = "Methylation Context",
       y = "Proportion Difference") +
  theme_minimal()
dev.off()


##by class##
TE_anno<-read.table("D:/PhD/DATA IN THESIS/DNA methylation/TAIR10_Transposable_Elements.txt", header=TRUE)
TE_anno
TE_gr<-makeGRangesFromDataFrame(TE_anno, keep.extra.columns = TRUE)
TE_gr
seqlevelsStyle(TE_gr)<-("NCBI")


SQ_TE_annotated<-merge(SQ_dmr_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
SQ_DEG_TE_annotated<-merge(SQ_deg_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
#subset to keep: TEgene, proportionDifference, chromatin, Transposon_Super_Family
SQ_chr_TE<- SQ_TE_annotated %>% select(TEgene, proportionDifference, chromatin, Transposon_Super_Family)
SQ_DEG_TE<- SQ_DEG_TE_annotated %>% select(TEgene, seqnames.x, start.x, end.x, X, log2FoldChange, padj)
SQ_chr_DEG_TE<-merge(SQ_chr_TE, SQ_DEG_TE, by="TEgene")
#subset by TE family
library(stringr)
SQ_TE_grouped <- SQ_chr_DEG_TE %>%
  mutate(Group = case_when(
    str_detect(Transposon_Super_Family, "^DNA") ~ "DNA",
    str_detect(Transposon_Super_Family, "^MITE") ~ "MITE",
    str_detect(Transposon_Super_Family, "^LINE") ~ "LINE",  
    str_detect(Transposon_Super_Family, "^RC") ~ "Helitron",
    str_detect(Transposon_Super_Family, "^LTR") ~ "LTR",
    str_detect(Transposon_Super_Family, "^RathE1") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE2") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE3") ~ "RathE",
    TRUE ~ "Other"  # Catch any unclassified entries
  ))
head(SQ_TE_grouped)

pdf("D:/PhD/DNA METHYLATION PAPER- DATA SUBMISSION/DEG_DMTE_correlationByFamily.pdf", width=8, height=5)


#plot
ggplot(as.data.frame(SQ_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  stat_cor(method = "pearson", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x.npc = "left", 
           label.y.npc = c(0.9, 0.8, 0.7, 0.6, 0.5)) +  # Adjust these values if you have more/fewer groups
  labs(
    title = "SQHvSQL - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"  ) +
  theme_minimal()

ggplot(as.data.frame(SQ_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  labs(
    title = "SQHvSQL - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"  ) +
  theme_minimal()

#FLHvFLL
FL_TE_annotated<-merge(FL_dmr_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
FL_DEG_TE_annotated<-merge(FL_deg_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
#subset to keep: TEgene, proportionDifference, chromatin, Transposon_Super_Family
FL_chr_TE<- FL_TE_annotated %>% select(TEgene, proportionDifference, chromatin, Transposon_Super_Family)
FL_DEG_TE<- FL_DEG_TE_annotated %>% select(TEgene, seqnames.x, start.x, end.x, X, log2FoldChange, padj)
FL_chr_DEG_TE<-merge(FL_chr_TE, FL_DEG_TE, by="TEgene")
table(FL_chr_DEG_TE$Transposon_Super_Family)
#subset by TE family

FL_TE_grouped <- FL_chr_DEG_TE %>%
  mutate(Group = case_when(
    str_detect(Transposon_Super_Family, "^DNA") ~ "DNA",
    str_detect(Transposon_Super_Family, "^MITE") ~ "MITE",
    str_detect(Transposon_Super_Family, "^LINE") ~ "LINE",  
    str_detect(Transposon_Super_Family, "^RC") ~ "Helitron",
    str_detect(Transposon_Super_Family, "^LTR") ~ "LTR",
    str_detect(Transposon_Super_Family, "^RathE1") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE2") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE3") ~ "RathE",
    str_detect(Transposon_Super_Family, "^SINE") ~ "SINE",
    TRUE ~ "Unassigned"  # Catch any unclassified entries
  ))
head(FL_TE_grouped)
table(FL_TE_grouped$Group)
#plot
ggplot(as.data.frame(FL_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  stat_cor(method = "pearson", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x.npc = "left", 
           label.y.npc = c(0.9, 0.8, 0.7, 0.6, 0.5)) +  # Adjust these values if you have more/fewer groups
  labs(
    title = "FLHvFLL - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"  ) +
  theme_minimal()

ggplot(as.data.frame(FL_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  labs(
    title = "FLHvFLL - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"  ) +
  theme_minimal()

#SQHvFLH
HL_TE_annotated<-merge(HL_dmr_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
HL_DEG_TE_annotated<-merge(HL_deg_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
#subset to keep: TEgene, proportionDifference, chromatin, Transposon_Super_Family
HL_chr_TE<- HL_TE_annotated %>% select(TEgene, proportionDifference, chromatin, Transposon_Super_Family)
HL_DEG_TE<- HL_DEG_TE_annotated %>% select(TEgene, seqnames.x, start.x, end.x, X, log2FoldChange, padj)
HL_chr_DEG_TE<-merge(HL_chr_TE, HL_DEG_TE, by="TEgene")
table(HL_chr_DEG_TE$Transposon_Super_Family)

HL_TE_grouped <- HL_chr_DEG_TE %>%
  mutate(Group = case_when(
    str_detect(Transposon_Super_Family, "^DNA") ~ "DNA",
    str_detect(Transposon_Super_Family, "^MITE") ~ "MITE",
    str_detect(Transposon_Super_Family, "^LINE") ~ "LINE",  
    str_detect(Transposon_Super_Family, "^RC") ~ "Helitron",
    str_detect(Transposon_Super_Family, "^LTR") ~ "LTR",
    str_detect(Transposon_Super_Family, "^RathE1") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE2") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE3") ~ "RathE",
    str_detect(Transposon_Super_Family, "^SINE") ~ "SINE",
    TRUE ~ "Unassigned"  # Catch any unclassified entries
  ))
   
head(HL_TE_grouped)
table(HL_TE_grouped$Group)

#plot
ggplot(as.data.frame(HL_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  stat_cor(method = "pearson", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x.npc = "left", 
           label.y.npc = c(0.9, 0.8, 0.7, 0.6, 0.5)) +  # Adjust these values if you have more/fewer groups
  labs(
    title = "SQHvFLH - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"  ) +
  theme_minimal()

ggplot(as.data.frame(HL_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  labs(
    title = "SQHvFLH - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"  ) +
  theme_minimal()

#SQLvFLL
LL_TE_annotated<-merge(LL_dmr_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
LL_DEG_TE_annotated<-merge(LL_deg_df, TE_anno, by.x="TEgene", by.y="Transposon_Name")
#subset to keep: TEgene, proportionDifference, chromatin, Transposon_Super_Family
LL_chr_TE<- LL_TE_annotated %>% select(TEgene, proportionDifference, chromatin, Transposon_Super_Family)
LL_DEG_TE<- LL_DEG_TE_annotated %>% select(TEgene, seqnames.x, start.x, end.x, X, log2FoldChange, padj)
LL_chr_DEG_TE<-merge(LL_chr_TE, LL_DEG_TE, by="TEgene")
table(LL_chr_DEG_TE$Transposon_Super_Family)

LL_TE_grouped <- LL_chr_DEG_TE %>%
  mutate(Group = case_when(
    str_detect(Transposon_Super_Family, "^DNA") ~ "DNA",
    str_detect(Transposon_Super_Family, "^MITE") ~ "MITE",
    str_detect(Transposon_Super_Family, "^LINE") ~ "LINE",  
    str_detect(Transposon_Super_Family, "^RC") ~ "Helitron",
    str_detect(Transposon_Super_Family, "^LTR") ~ "LTR",
    str_detect(Transposon_Super_Family, "^RathE1") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE2") ~ "RathE",
    str_detect(Transposon_Super_Family, "^RathE3") ~ "RathE",
    str_detect(Transposon_Super_Family, "^SINE") ~ "SINE",
    TRUE ~ "Unassigned"  # Catch any unclassified entries
  ))

head(LL_TE_grouped)
table(LL_TE_grouped$Group)

#plot
ggplot(as.data.frame(LL_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  stat_cor(method = "pearson", 
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           label.y.npc = "bottom",
          label.x.npc = "left") +  # Adjust these values if you have more/fewer groups
  labs(
    title = "SQLvFLL - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"
  ) +
  theme_minimal()

ggplot(as.data.frame(LL_TE_grouped), aes(x = log2FoldChange, y = proportionDifference, col = Group)) +
  geom_point() +
  geom_smooth(method = "lm", fill = NA) +
  labs(
    title = "SQLvFLL - DEG with 5kb DM TE by class",
    x = "Log2Fold Change",
    y = "Proportion Difference Methylation"  ) +
  theme_minimal()
dev.off()
