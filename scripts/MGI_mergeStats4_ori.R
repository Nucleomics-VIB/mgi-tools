#!/usr/bin/env Rscript

# script: MGI_mergeStats4.R
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-10 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# compatible with the output of the original SplitDualbarcodes.pl
# merge multiple BarcodeStat.txt and TagStat.txt obtained after demultiplexing MGI data
# save resulting merged tables to file
# plot barcode frequency from the merged BarcodeStat

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggplot2"))

version <- "v1.00, 2023-02-10"

##############################
# merge BarcodeStat.txt files
##############################

vectorBarcodeStat <- list.files(path = '.', 
                pattern = "BarcodeStat.txt" ,
                full.names = FALSE, 
                recursive = TRUE,
                ignore.case = FALSE,
                include.dirs = TRUE)

# merge all BarcodeStat.txt files to a single file using data.table fread
BarcodeStatMerge <- data.table::rbindlist(lapply(vectorBarcodeStat, data.table::fread))
# remove Total rows and zero percentages
BarcodeStatMerge <- BarcodeStatMerge[!grepl("Total", BarcodeStatMerge$`#SpeciesNO`),]
BarcodeStatMerge$Pct <- rep("0",nrow(BarcodeStatMerge))

# summarize and sum using dplyr
BarcodeStat <- BarcodeStatMerge %>% 
  dplyr::group_by(`#SpeciesNO`) %>%
  dplyr::summarise(across(c(Correct, Corrected, Total), sum))

# add back Pct column
BarcodeStat$Pct <- sprintf("%1.2f%%", 100*BarcodeStat$Total/sum(BarcodeStat$Total))

# save to csv file
outfile <- "BarcodeStat_merged.csv"
write.csv(BarcodeStat, file = outfile, row.names = FALSE)

##############################
# plot barcode densities
##############################

plot.data <- BarcodeStat[,1:3] %>% reshape2::melt(id.vars = "#SpeciesNO")
colnames(plot.data) <- c("barcodes", "type", "count")

outfile <- "BarcodeStat_density.pdf"
pdf(file=outfile, width = 10, height = 10, bg = "white")

ggplot(plot.data, aes(x = barcodes, y= count, fill = forcats::fct_rev(type))) + 
  geom_bar(width=0.7,stat = "identity") +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y=element_text(angle = 0, hjust = 0, size=5)) +
  theme(legend.title=element_blank()) +
  theme(legend.position = "top") +
  theme(legend.key.size = unit(0.25, 'cm')) +
  guides(fill = guide_legend(reverse = TRUE)) +
  xlab("barcodes") +
  ylab("read (pair) count") +
  ggtitle("MGI400 Barcode frequency")

null <- dev.off()

##############################
# merge TagStat.txt files
##############################

vectorTagStat <- list.files(path = '.',
                            pattern = "TagStat.txt" ,
                            full.names = FALSE, 
                            recursive = TRUE,
                            ignore.case = FALSE,
                            include.dirs = TRUE)

# merge all TagStat.txt files to a single file
TagStatMerge <- data.table::rbindlist(lapply(vectorTagStat, data.table::fread))
TagStatMerge$Pct <- rep("0",nrow(TagStatMerge))

# summarize and sum using dplyr
TagStat <- TagStatMerge %>% 
  dplyr::group_by(`#Sequence`,SpeciesNO) %>%
  dplyr::summarise(across(c(readCount), sum), .groups = "drop") %>%
  arrange(desc(readCount))

# add back Pct column
TagStat$Pct <- sprintf("%1.2f%%", 100*TagStat$readCount/sum(TagStat$readCount))

# save to csv files
outfile <- "TagStat_merged.csv"
write.csv(TagStat, file = outfile, row.names = FALSE)

# top 100 unknown barcodes
topUnknown <- TagStat %>%
  filter(SpeciesNO == "unknown") %>%
  arrange(desc(readCount)) %>%
  head(100)
  
outfile <- "TagStat_top100_unknown.csv"
write.csv(topUnknown, file = outfile, row.names = FALSE)
