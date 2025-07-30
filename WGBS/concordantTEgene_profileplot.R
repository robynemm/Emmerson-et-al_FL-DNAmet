names(SQ_neg_cocordant_1)[1]<-paste("gene")
names(SQ_pos_cocordant_1)[1]<-paste("gene")
names(FL_neg_cocordant_1)[1]<-paste("gene")
names(FL_pos_cocordant_1)[1]<-paste("gene")
names(HL_neg_cocordant_1)[1]<-paste("gene")
names(HL_pos_cocordant_1)[1]<-paste("gene")
names(LL_neg_cocordant_1)[1]<-paste("gene")
names(LL_pos_cocordant_1)[1]<-paste("gene")

names(SQ_neg_cocordant_2)[1]<-paste("gene")
names(SQ_pos_cocordant_2)[1]<-paste("gene")
names(FL_neg_cocordant_2)[1]<-paste("gene")
names(FL_pos_cocordant_2)[1]<-paste("gene")
names(HL_neg_cocordant_2)[1]<-paste("gene")
names(HL_pos_cocordant_2)[1]<-paste("gene")
names(LL_neg_cocordant_2)[1]<-paste("gene")
names(LL_pos_cocordant_2)[1]<-paste("gene")

names(SQ_neg_cocordant_5)[1]<-paste("gene")
names(SQ_pos_cocordant_5)[1]<-paste("gene")
names(FL_neg_cocordant_5)[1]<-paste("gene")
names(FL_pos_cocordant_5)[1]<-paste("gene")
names(HL_neg_cocordant_5)[1]<-paste("gene")
names(HL_pos_cocordant_5)[1]<-paste("gene")
names(LL_neg_cocordant_5)[1]<-paste("gene")
names(LL_pos_cocordant_5)[1]<-paste("gene")




for (i in 1:nrow(LL_neg_cocordant_2)) {
  # Extract the current gene information
  gene <- LL_neg_cocordant_1[i, ]
  # Ensure gene$gene is a character and handle missing values
  gene_name <- as.character(gene$gene)
  if (is.na(gene_name)) {
    print(paste("Skipping row", i, "due to missing gene name."))
    next
  }
  # Define the region of interest (10kb upstream and downstream)
  region_start <- max(1, gene$start - 10000)  # Ensure start is not below 1
  region_end <- gene$end + 10000
  # Create a GRanges object for the region
  region <- GRanges(
    seqnames = gene$seqnames,
    ranges = IRanges(start = region_start, end = region_end)
  )
  # Informative message about the current gene being processed
  print(paste("Processing gene:", gene_name))
  # Create a unique file name for each gene
  pdf_filename <- paste0("Methylation_ProfileLL_", gene_name, ".pdf")
  # Open a PDF device for the current gene
  pdf(pdf_filename, width = 9, height = 3.5, pointsize=12)
  # Try to plot the local methylation profile
  tryCatch({
    plotLocalMethylationProfile(
      SQL_all, FLL_all,
      region = region,
      DMRs = list("CG"=LL_DMR_CG), 
      context = "CG",  # Adjust as necessary
      gff=annotation,
      conditionsNames = c("SQL", "FLL"), # Adjust as necessary
      main = paste("Methylation Profile -", gene_name) # Add title to the plot
    )
  }, error = function(e) {
    print(paste("Error processing gene:", gene_name, "-", conditionMessage(e)))
  })
  # Close the PDF device
  dev.off()
  print(paste("Saved plot for gene:", gene_name, "to file:", pdf_filename))
}

for (i in 1:nrow(LL_pos_cocordant_2)) {
  # Extract the current gene information
  gene <- LL_pos_cocordant_1[i, ]
  # Ensure gene$gene is a character and handle missing values
  gene_name <- as.character(gene$gene)
  if (is.na(gene_name)) {
    print(paste("Skipping row", i, "due to missing gene name."))
    next
  }
  # Define the region of interest (10kb upstream and downstream)
  region_start <- max(1, gene$start - 10000)  # Ensure start is not below 1
  region_end <- gene$end + 10000
  # Create a GRanges object for the region
  region <- GRanges(
    seqnames = gene$seqnames,
    ranges = IRanges(start = region_start, end = region_end)
  )
  # Informative message about the current gene being processed
  print(paste("Processing gene:", gene_name))
  # Create a unique file name for each gene
  pdf_filename <- paste0("Methylation_ProfileLL_", gene_name, ".pdf")
  # Open a PDF device for the current gene
  pdf(pdf_filename, width = 9, height = 3.5, pointsize=12)
  # Try to plot the local methylation profile
  tryCatch({
    plotLocalMethylationProfile(
      SQL_all, FLL_all,
      region = region,
      DMRs = list("CG"=LL_DMR_CG), 
      context = "CG",  # Adjust as necessary
      gff=annotation,
      conditionsNames = c("SQL", "FLL"), # Adjust as necessary
      main = paste("Methylation Profile -", gene_name) # Add title to the plot
    )
  }, error = function(e) {
    print(paste("Error processing gene:", gene_name, "-", conditionMessage(e)))
  })
  # Close the PDF device
  dev.off()
  print(paste("Saved plot for gene:", gene_name, "to file:", pdf_filename))
}


AT1g66370<-GRanges(seqnames= Rle("Chr1"), ranges=IRanges(6755000, 6765000))
pdf("AT1g66370_localProfile_edit.pdf", width= 10, height = 4, pointsize = 12)
plotLocalMethylationProfile(
  SQL_all, FLL_all,
  region = AT1g66370,
  DMRs = list("CG"=LL_DMR_CG), 
  context = "CG", 
  gff=annotation,
  conditionsNames = c("SQL", "FLL"), 
  main = paste("Methylation Profile AT1G66370" ) 
)
dev.off()

At1g19520<-GRanges(seqnames = Rle("Chr1"), ranges=IRanges(24750000,24755000))
pdf("AT1g19520_localProfile_edit.pdf", width= 10, height = 4, pointsize = 12)
plotLocalMethylationProfile(
  SQL_all, FLL_all,
  region = At1g19520,
  DMRs = list("CG"=LL_DMR_CG), 
  context = "CG", 
  gff=annotation,
  conditionsNames = c("SQL", "FLL"), 
  main = paste("Methylation Profile AT1G19520" ) )
dev.off()


At3g27640<-GRanges(seqnames = Rle("Chr3"), ranges=IRanges(10230000,10240000))
pdf("AT3g27640_localProfile.pdf", width= 10, height = 4, pointsize = 12)
plotLocalMethylationProfile(
  SQL_all, FLL_all,
  region = At3g27640,
  DMRs = list("CG"=LL_DMR_CG), 
  context = "CG", 
  gff=annotation,
  conditionsNames = c("SQL", "FLL"), 
  main = paste("Methylation Profile AT3G27640" ) )
dev.off()

At4g09510<-GRanges(seqnames = Rle("Chr4"), ranges=IRanges(6020000,6025000))
pdf("AT4g09510_localProfile_edit.pdf", width= 10, height = 4, pointsize = 12)
plotLocalMethylationProfile(
  SQL_all, FLL_all,
  region = At4g09510,
  DMRs = list("CG"=LL_DMR_CG), 
  context = "CG", 
  gff=annotation,
  conditionsNames = c("SQL", "FLL"), 
  main = paste("Methylation Profile AT4G09510" ) )
dev.off()

