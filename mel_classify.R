library(pamr)

## Read in sample dataset: mel data, 12 samples  
mel.data <- pamr.from.excel("/Users/laurensaunders/Desktop/normalized_out2/genes.fpkm_table_pam.txt", 14, sample.labels=TRUE)
mel.data2 <- pamr.knnimpute(mel.data)
