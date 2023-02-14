#!/usr/bin/env Rscript

# script: MGI_mergeStats4.R
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-10 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# merge multiple BarcodeStat.txt and SequenceStat.txt obtained after demultiplexing MGI data
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
BarcodeStatMerge <- BarcodeStatMerge[!grepl("Total", BarcodeStatMerge$`#Barcode`),]
BarcodeStatMerge$`Percentage(%)` <- rep("0",nrow(BarcodeStatMerge))

# summarize and sum using dplyr
BarcodeStat <- BarcodeStatMerge %>% 
  dplyr::group_by(`#Barcode`) %>%
  dplyr::summarise(across(c(Correct, Corrected, Total), sum))

# add back Percentage(%) column
BarcodeStat$`Percentage(%)` <- sprintf("%1.2f%%", 100*BarcodeStat$Total/sum(BarcodeStat$Total))

# save to csv file
outfile <- "BarcodeStat_merged.csv"
write.csv(BarcodeStat, file = outfile, row.names = FALSE)

##############################
# plot barcode densities
##############################

plot.data <- BarcodeStat[,1:3] %>% reshape2::melt(id.vars = "#Barcode")
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
# merge SequenceStat.txt files
##############################

vectorSequenceStat <- list.files(path = '.',
                            pattern = "SequenceStat.txt" ,
                            full.names = FALSE, 
                            recursive = TRUE,
                            ignore.case = FALSE,
                            include.dirs = TRUE)

# merge all SequenceStat.txt files to a single file
SequenceStatMerge <- data.table::rbindlist(lapply(vectorSequenceStat, data.table::fread))
SequenceStatMerge$`Percentage(%)` <- rep("0",nrow(SequenceStatMerge))

# summarize and sum using dplyr
SequenceStat <- SequenceStatMerge %>% 
  dplyr::group_by(`#Sequence`,Barcode) %>%
  dplyr::summarise(across(c(Count), sum), .groups = "drop") %>%
  arrange(desc(Count))

# add back Percentage(%) column
SequenceStat$`Percentage(%)` <- sprintf("%1.2f%%", 100*SequenceStat$Count/sum(SequenceStat$Count))

# save to csv files
outfile <- "SequenceStat_merged.csv"
write.csv(SequenceStat, file = outfile, row.names = FALSE)

# top 100 unknown barcodes
topUnknown <- SequenceStat %>%
  filter(Barcode == "undecoded") %>%
  arrange(desc(Count)) %>%
  head(100)
  
outfile <- "SequenceStat_top100_unknown.csv"
write.csv(topUnknown, file = outfile, row.names = FALSE)
