library(cummeRbund)
library(pheatmap)

mel_cuff_data <- readCufflinks('/Users/laurensaunders/Desktop/diff_out')
csDensity(genes(mel_cuff_data), pseudocount = .1)
csScatter(genes(mel_cuff_data), 'UA', 'A')
csVolcano(genes(mel_cuff_data), 'UA', 'A')

# plot expression levels for genes of interest
rb1 <- getGene(mel_cuff_data, 'rb1')
expressionBarplot(rb1)

sox11a <- getGene(mel_cuff_data, 'sox11a')
expressionBarplot(sox11a)

# Plot volcano plots for all pairwise combinations of coditions
csVolcano(genes(mel_cuff_data),'A','UA',alpha=0.05, ...
          showSignificant=TRUE,xlimits=c(-10,10))

# Find and Plot expression of senescence biomarkers
mel_genes <- getGenes(mel_cuff_data,c("tyrp1b","pmela","mitfa","kita","gch2","dct", "mlphb"))
expressionBarplot(mel_genes,logMode=TRUE)

NC_genes <- getGenes(mel_cuff_data,c("sox10","foxd3","snai2"))
expressionBarplot(NC_genes, logMode=TRUE)

# record DEGS and put data in easy to view files
gene_diff_data <- diffData(genes(mel_cuff_data))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
nrow(sig_gene_data)
write.table(sig_gene_data, 'diff_genes_cd.txt', sep = '\t', row.names = F, col.names = T, quote = F)

# Record fpkm matrix of differentially expressed genes for heatmap
fpkm_data<-fpkmMatrix(genes(mel_cuff_data))
fpkm_diff<-(subset(fpkm_data,rownames(fpkm_data) %in% sig_gene_data$gene_id))

# Make a data.frame CuffGeneSet of differentially expressed genes for CummeRbund package
diff_expr=getGenes(mel_cuff_data,sig_gene_data$gene_id)

# Apply a aq-value, fold-change and FPKM cutoff to keep differentially expressed genes whose expression in at least one sample is >10 FPKM
diff_sig_cutoff<-getGenes(mel_cuff_data,sig_gene_data$gene_id[which(sig_gene_data$q_value<=0.05 & abs(sig_gene_data$log2_fold_change)>=2 & apply(cbind(sig_gene_data$value_1,sig_gene_data$value_2),1,max)>5)])
diff_Matrix_cutoff<-fpkmMatrix(diff_sig_cutoff)
csHeatmap(diff_sig_cutoff,clustering='row')