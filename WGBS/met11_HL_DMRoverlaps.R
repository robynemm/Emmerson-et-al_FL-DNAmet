###overlaps in DMRs between met1-1 and SQHvFLH
library(GRanges)


setwd("E:/DNA METHYLATION PAPER- DATA SUBMISSION/")

met11_CG<- read.csv("DMRs_WT_to_met1_1.CG.csv")
met11_CHG<- read.csv("DMRs_WT_to_met1_1.CHG.csv")
met11_CHH<- read.csv("DMRs_WT_to_met1_1.CHH.csv")
SQHvFLH<- read.csv("SQHvFLH_DMRs.csv")

SQHvFLH_CG<- SQHvFLH %>%
  filter(context =="CG")
SQHvFLH_CHG<- SQHvFLH %>%
  filter(context =="CHG")
SQHvFLH_CHH<- SQHvFLH %>%
  filter(context =="CHH")

gr_met<-makeGRangesFromDataFrame(met11_CG)
gr_metCHG<-makeGRangesFromDataFrame(met11_CHG)
gr_metCHH<-makeGRangesFromDataFrame(met11_CHH)
gr_CG<-makeGRangesFromDataFrame(SQHvFLH_CG, keep.extra.columns = TRUE)
gr_CHG<-makeGRangesFromDataFrame(SQHvFLH_CHG, keep.extra.columns = TRUE)
gr_CHH<-makeGRangesFromDataFrame(SQHvFLH_CHH, keep.extra.columns = TRUE)

seqlevelsStyle(gr_met)<-"UCSC"
seqlevelsStyle(gr_metCHG)<-"UCSC"
seqlevelsStyle(gr_metCHH)<-"UCSC"

CG_overlaps<-findOverlaps(gr_met,gr_CG)
CHG_overlaps<-findOverlaps(gr_metCHG, gr_CHG)
CHH_overlaps<-findOverlaps(gr_metCHH, gr_CHH)

CG_hits<-subjectHits(CG_overlaps)
CG_list<-gr_CG[CG_hits]

CHG_hits<-subjectHits(CHG_overlaps)
CHG_list<-gr_CHG[CHG_hits]

CHH_hits<-subjectHits(CHH_overlaps)
CHH_list<-gr_CHH[CHH_hits]
