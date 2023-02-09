#!/usr/bin/env Rscript

# script: MGI_mergeStats4.R
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-01-23 v1.01

suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))

vectorBarcodeStat <- list.files(path = ".", 
                pattern = "BarcodeStat.txt" ,
                full.names = FALSE, 
                recursive = TRUE,
                ignore.case = FALSE,
                include.dirs = TRUE)

# merge all BarcodeStat.txt files to a single file
BarcodeStatMerge <- rbindlist(lapply(vectorBarcodeStat, fread))
# remove Total rows
BarcodeStatMerge <- BarcodeStatMerge[!grepl("Total", BarcodeStatMerge$`#SpeciesNO`),]
BarcodeStatMerge[,5] <- rep("0",nrow(BarcodeStatMerge))

# summarize and sum
BarcodeStat <- ddply(BarcodeStatMerge, ~ `#SpeciesNO`, numcolwise(sum))

# add back Pct column
BarcodeStat$Pct <- sprintf("%1.2f%%", 100*BarcodeStat$Total/sum(BarcodeStat$Total))

# save to csv file
write.csv(BarcodeStat, file = "BarcodeStat_merged.csv",row.names = FALSE)

# plot barcode bars
pdf(file="BarcodeStat_density.pdf", width = 10, height = 10, bg = "white")
ggplot(data=BarcodeStat, aes(x=`#SpeciesNO`, y=Total)) +
  geom_bar(width=0.7, stat="identity") +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y=element_text(angle = 0, hjust = 0, size=5)) +
  xlab("barcodes") +
  ylab("read (pair) count")
null <- dev.off()

##############
vectorTagStat <- list.files(path = ".",
                            pattern = "TagStat.txt" ,
                            full.names = FALSE, 
                            recursive = TRUE,
                            ignore.case = FALSE,
                            include.dirs = TRUE)

# merge all TagStat.txt files to a single file
TagStatMerge <- rbindlist(lapply(vectorTagStat, fread))
TagStatMerge[,4] <- rep("0",nrow(TagStatMerge))

# summarize and sum
TagStat <- ddply(TagStatMerge, .(`#Sequence`,SpeciesNO), numcolwise(sum))

# add back Pct column
TagStat$Pct <- sprintf("%1.2f%%", 100*TagStat$readCount/sum(TagStat$readCount))

# reorder descending
TagStat <- TagStat[order(-TagStat$readCount),]

# save to csv file
write.csv(TagStat, file = "TagStat_merged.csv",row.names = FALSE)


