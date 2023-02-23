#!/usr/bin/env Rscript

# script: MGI_barcode_distance.R
# Stéphane Plaisance - VIB-Nucleomics Core - 2023-02-17 v1.00
# visit our Git: https://github.com/Nucleomics-VIB

# compute distance between all MGI barcode@1's, barcode#2's
# plot as heatmaps to identify barcode maximal edits for efficient demultiplexing

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringdist"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("cowplot"))

version <- "v1.00, 2023-02-17"

option_list <- list(
  make_option(c("-l", "--list"), type="character", default=NA,
              help="a space-separated list of barcodes")
)

# parse options
opt <- parse_args(OptionParser(option_list=option_list))

if ( is.na(opt$list) ) {
  stop("Usage: 'MGI_barcode_distance.R -l barcode.list'")
}

# compute pairwise edit distance
func <-function(x,y){z<-stringdist(x, y, method= "lv")
return(z)
}

# load log data in
data <- read.delim(opt$list, sep = "",
                   header=FALSE,
                   comment.char = "#",
                   stringsAsFactors=FALSE)

# data <- read.delim("/data/analyses/MGI_demux/new.barcode.list",
#                     sep="\t",
#                     header=FALSE,
#                     comment.char = "#",
#                     stringsAsFactors=FALSE)

#######################
# barcode1 with itself
#######################

ed1 <- crossing(data$V2, data$V2)
colnames(ed1) <- c("bc1", "bc1_")

# compute all pairwise
ed1$dist <- mapply(func, ed1$bc1, ed1$bc1_)

# find minimum >0
min.d1 <- min(ed1$dist[ed1$dist> 0])
plot.title.d1 <- paste0("Mismatch number between any two barcode-1 sequences (min:", min.d1, ")", sep="")

p1 <- ggplot(data = ed1, aes(bc1, bc1_, fill = dist)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "black",  mid="white", high = "white", 
                       limit = c(0,10), space = "Lab",
                       midpoint = 4,
                       name="edit distance") +
  geom_text(aes(bc1, bc1_, label = dist), color = "black", size = 1.5) +
  theme_minimal() + 
  theme(axis.text.x = element_text(color = "grey20", size = 3, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 3, angle = 0, hjust = 1, vjust = 0, face = "plain")) +
  coord_fixed() +
  scale_color_continuous(trans = 'reverse') +
  ggtitle(plot.title.d1) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# get max count to add some space above plot2 and show the top count
max.cnt1 <- max(hist(ed1$dist, breaks=10, plot=FALSE)$counts)

p1.hist <- ggplot(data = ed1, aes(x=dist, fill=factor(ifelse(dist>0,"different","self"),levels = c("self", "different")))) +
  geom_bar(stat="count") +
  scale_fill_manual(name = "type", values=c("red","grey50")) +
  scale_x_continuous(breaks= pretty_breaks()) +
  theme_minimal() +
  xlab('mutual edit distance') +
  ylab('barcode pair number') +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1, colour = "black") +
  theme(legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.position = "top") +
  guides(size=guide_legend(direction='horizontal')) +
  ylim(0, max.cnt1 * 1.25) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p1.tot <- plot_grid(p1, p1.hist, ncol = 1, nrow = 2, rel_heights = c(5,2))

suppressMessages(ggsave(plot=p1.tot, filename="bc1xself_distance.pdf"))

#######################
# barcode2 with itself
#######################

ed2 <- crossing(data$V3, data$V3)
colnames(ed2) <- c("bc2", "bc2_")

# compute all pairwise
ed2$dist <- mapply(func, ed2$bc2, ed2$bc2_)

# find minimum >0
min.d2 <- min(ed2$dist[ed2$dist> 0])
plot.title.d2 <- paste0("Mismatch number between any two barcode-2 sequences (min:", min.d2, ")", sep="")

p2 <- ggplot(data = ed2, aes(bc2, bc2_, fill = dist)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "black",  mid="white", high = "white", 
                       limit = c(0,10), space = "Lab",
                       midpoint = 4,
                       name="edit distance") +
  geom_text(aes(bc2, bc2_, label = dist), 
            color = "black", size = 1.5) +
  theme_minimal() + 
  theme(axis.text.x = element_text(color = "grey20", size = 3, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 3, angle = 0, hjust = 1, vjust = 0, face = "plain")) +
  coord_fixed() +
  scale_color_continuous(trans = 'reverse') +
  ggtitle(plot.title.d2)

# get max count to add some space above plot2 and show the top count
max.cnt2 <- max(hist(ed2$dist, breaks=10, plot=FALSE)$counts)

p2.hist <- ggplot(data = ed2, aes(x=dist, fill=factor(ifelse(dist>0,"different","self"),levels = c("self", "different")))) +
  geom_bar(stat="count") +
  scale_fill_manual(name = "type", values=c("red","grey50")) +
  scale_x_continuous(breaks= pretty_breaks()) +
  theme_minimal() +
  xlab('mutual edit distance') +
  ylab('barcode pair number') +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1, colour = "black") +
  theme(legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.position = "top") +
  guides(size=guide_legend(direction='horizontal')) +
  ylim(0, max.cnt2 * 1.25) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p2.tot <- plot_grid(p2, p2.hist, ncol = 1, nrow = 2, rel_heights = c(5,2))

suppressMessages(ggsave(plot=p2.tot, filename="bc2xself_distance.pdf"))

#######################
# barcode1 with barcode2
#######################

ed12 <- crossing(data$V2, data$V3)
colnames(ed12) <- c("bc1", "bc2")

# compute all pairwise
ed12$dist <- mapply(func, ed12$bc1, ed12$bc2)

# find minimum >0
min.d12 <- min(ed12$dist[ed12$dist> 0])
plot.title.d12 <- paste0("Mismatch number between barcode-1 and barcode-2 sequences (min:", min.d12, ")", sep="")

p3 <- ggplot(data = ed12, aes(bc1, bc2, fill = dist)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "black",  mid="white", high = "white", 
                       limit = c(0,10), space = "Lab",
                       midpoint = 4,
                       name="edit distance") +
  geom_text(aes(bc1, bc2, label = dist), 
            color = "black", size = 1.5) +
  theme_minimal() + 
  theme(axis.text.x = element_text(color = "grey20", size = 3, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 3, angle = 0, hjust = 1, vjust = 0, face = "plain")) +
  coord_fixed() +
  scale_color_continuous(trans = 'reverse') +
  ggtitle(plot.title.d12)


# get max count to add some space above plot2 and show the top count
max.cnt3 <- max(hist(ed12$dist, breaks=10, plot=FALSE)$counts)

p3.hist <- ggplot(data = ed12, aes(x=dist, fill=factor(ifelse(dist>0,"different","self")))) +
  geom_bar(stat="count") +
  scale_fill_manual(name = "type", values=c("grey50","red")) +
  scale_x_continuous(breaks= pretty_breaks()) +
  theme_minimal() +
  xlab('mutual edit distance') +
  ylab('barcode pair number') +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1, colour = "black") +
  theme(legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.position = "top") +
  guides(size=guide_legend(direction='horizontal')) +
  ylim(0, max.cnt3 * 1.25) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p3.tot <- plot_grid(p3, p3.hist, ncol = 1, nrow = 2, rel_heights = c(5,2))

suppressMessages(ggsave(plot=p3.tot, filename="bc1xbc2_distance.pdf"))

################################
# barcode1_barcode2 with itself
################################

data$merged <- paste0(data$V2, data$V3, sep="")

edmerged <- crossing(data$merged, data$merged) 
colnames(edmerged) <- c("dual_bc", "dual_bc_")

# compute all pairwise
edmerged$dist <- mapply(func, edmerged$dual_bc, edmerged$dual_bc_)

# find minimum >0
min.dual <- min(edmerged$dist[edmerged$dist> 0])
plot.title.dual <- paste0("Mismatch number between any two merged barcodes (min:", min.dual, ")", sep="")

p4 <- ggplot(data = edmerged, aes(dual_bc, dual_bc_, fill = dist)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "black",  mid="white", high = "white", 
                       limit = c(0,20), space = "Lab",
                       midpoint = 10,
                       name="edit distance") +
  geom_text(aes(dual_bc, dual_bc_, label = dist), 
            color = "black", size = 1.2) +
  theme_minimal() + 
  theme(axis.text.x = element_text(color = "grey20", size = 3, angle = 45, hjust = 1, vjust = 1, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 3, angle = 0, hjust = 1, vjust = 0, face = "plain")) +
  coord_fixed() +
  scale_color_continuous(trans = 'reverse') +
  ggtitle(plot.title.dual)


# get max count to add some space above plot2 and show the top count
max.cnt4 <- max(hist(edmerged$dist, breaks=20, plot=FALSE)$counts)

p4.hist <- ggplot(data = edmerged, aes(x=dist, fill=factor(ifelse(dist>0,"different","self")))) +
  geom_bar(stat="count") +
  scale_fill_manual(name = "type", values=c("grey50","red")) +
  scale_x_continuous(breaks= pretty_breaks()) +
  theme_minimal() +
  xlab('mutual edit distance') +
  ylab('merged-barcode pair number') +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1, colour = "black") +
  theme(legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.position = "top") +
  guides(size=guide_legend(direction='horizontal')) +
  ylim(0, max.cnt4 * 1.25) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

p4.tot <- plot_grid(p4, p4.hist, ncol = 1, nrow = 2, rel_heights = c(5,2))

suppressMessages(ggsave(plot=p4.tot, filename="bcdualxself-distance.pdf"))
