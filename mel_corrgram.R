## Script for making corrgrams of gene or isoform fpkm data from many samples
## Import data; quote option is needed because some genes include '
fpkm.data = read.table ("~/Desktop/normalized_out2/genes.fpkm_table.txt", header = TRUE, sep = "\t", quote = " ")

## Remove all zero expression values by converting to NA
zer.NA.fpkmdata  = fpkm.data[,2:13]
zer.NA.fpkmdata[zer.NA.fpkmdata==0] <- NA

## Make a corrgram showing all comparisons.
library(corrgram)
## Write correlogram into png file
png(filename = "~/Desktop/fpkmdata_all_zerNA_corrgraph_ordered.png", width = 1200, height = 1200, units = "px")
corrgram(zer.NA.fpkmdata, order = TRUE, lower.panel = panel.shade, upper.panel = panel.conf)
dev.off()

## Make a corrgram with fewer samples, and order the samples by similarity
png(filename = "~/Desktop/fpkmdata_pnsplus_zerNA_corrgraphs_UnAbl.png", width = 1200, height = 1200, units = "px")
corrgram(zer.NA.fpkmdata[,c(7:12)], order = TRUE, lower.panel = panel.shade, upper.panel = panel.conf)
dev.off()
