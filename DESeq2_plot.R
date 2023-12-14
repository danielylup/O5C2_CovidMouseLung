library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(goseq)
library(patchwork)

#Create a dataframe describing the four dataset
Treatment <- c("PBS", "PBS", "O5C2", "O5C2")
bulkCML_GG_sample_C2 <- as.data.frame(Treatment)
rownames(bulkCML_GG_sample_C2) <- colnames(bulkCML_GG_count_C2) 
bulkCML_GG_sample2_C2 <- bulkCML_GG_sample_C2
bulkCML_GG_sample2_C2$Sample <- rownames(bulkCML_GG_sample2_C2)

#Import raw read count dataset
bulkCML_GG_count_C2 <- read.table("./bulkCML_GG_count_C2", header = T, sep = "\t")

#Import the read count into DEseq2 format and normalize
bulkCML_GG_dds_C2 <- DESeqDataSetFromMatrix(countData = bulkCML_GG_count_C2, colData = bulkCML_GG_sample_C2, design = ~ Treatment)
bulkCML_GG_vsd_C2 <- vst(bulkCML_GG_dds_C2, blind = FALSE)

#Calculate distance matrix between samples
bulkCML_GG_assay_C2 <- assay(bulkCML_GG_vsd_C2)
colnames(bulkCML_GG_assay_C2) <- rownames(bulkCML_GG_sample_C2)
bulkCML_GG_Dists_C2 <- dist(t(bulkCML_GG_assay_C2))
bulkCML_GG_DistMatrix_C2 <- as.matrix(bulkCML_GG_Dists_C2)

#Performing DGE analysis using DESeq2
bulkCML_GG_dds_C2 <- DESeq(bulkCML_GG_dds_C2)
bulkCML_GG_res_C2 <- results(bulkCML_GG_dds_C2, contrast = c("Treatment", "O5C2", "PBS"))


#####Plot figure#####

#Figure 6f
png("./Publication_figure/Fig6f_volplot.png", units="in", width=5, height=6, res=600)
EnhancedVolcano(bulkCML_GG_res_C2, lab = rownames(bulkCML_GG_res_C2), x = 'log2FoldChange', y = 'pvalue', title = '', axisLabSize = 18, pCutoff = 0.01, subtitle = '', subtitleLabSize = 0.1)
dev.off()

#Figure 6g
bulkCML_GG_res_C2_DF <- as.data.frame(bulkCML_GG_res_C2)
bulkCML_GG_res_C2_DF$significant <- ifelse(bulkCML_GG_res_C2_DF$padj < .05, "Significant", NA)
bulkCML_GG_res_C2_DF$Gene <- rownames(bulkCML_GG_res_C2_DF)
sigGenes <- rownames(bulkCML_GG_res_C2_DF[bulkCML_GG_res_C2_DF$padj <= .05 & abs(bulkCML_GG_res_C2_DF$log2FoldChange) > 3,])
bulkCML_GG_deseq2_C2 <- assay(bulkCML_GG_vsd_C2)
bulkCML_GG_deseq2_C2 <- as.data.frame(bulkCML_GG_deseq2_C2)
bulkCML_GG_deseq2_C2$Gene <- rownames(bulkCML_GG_deseq2_C2)
bulkCML_GG_heatmap_C2 <- bulkCML_GG_deseq2_C2[bulkCML_GG_deseq2_C2$Gene %in% sigGenes,]
bulkCML_GG_heatmap_C2 <- melt(bulkCML_GG_heatmap_C2, id.vars = c("Gene"))
bulkCML_GG_heatmap_matrix_C2 <- dcast(bulkCML_GG_heatmap_C2, Gene ~ variable)
rownames(bulkCML_GG_heatmap_matrix_C2) <- bulkCML_GG_heatmap_matrix_C2$Gene
bulkCML_GG_heatmap_matrix_C2$Gene <- NULL

png("./Publication_figure/Fig6g_DEGheatmap_.png", units = "in", height = 15, width = 5, res = 600)
pheatmap(bulkCML_GG_heatmap_matrix_C2, cluster_rows = T, show_rownames = T, annotation = bulkCML_GG_sample_C2, fontsize = 7, cutree_rows = 2, scale = )
dev.off()

#Figure 6h
upreg_C2_genes <- bulkCML_GG_res_C2$padj < 0.05 & bulkCML_GG_res_C2$log2FoldChange > 0 & !is.na(bulkCML_GG_res_C2$padj)
names(upreg_C2_genes) <- rownames(bulkCML_GG_res_C2)
Mmus_bias_upC2 <- Mmus_bias[names(upreg_C2_genes)]
pwf_upC2 <- nullp(upreg_C2_genes, bias.data = as.numeric(Mmus_bias_upC2))
GOres_upC2 <- goseq(pwf_upC2, "mm9", "geneSymbol", test.cats = "GO:BP")

downreg_C2_genes <- bulkCML_GG_res_C2$padj < 0.05 & bulkCML_GG_res_C2$log2FoldChange < 0 & !is.na(bulkCML_GG_res_C2$padj)
names(downreg_C2_genes) <- rownames(bulkCML_GG_res_C2)
Mmus_bias_downC2 <- Mmus_bias[names(downreg_C2_genes)]
pwf_downC2 <- nullp(downreg_C2_genes, bias.data = as.numeric(Mmus_bias_downC2))
GOres_downC2 <- goseq(pwf_downC2, "mm9", "geneSymbol", test.cats = "GO:BP")

GOplot_upC2 <- GOres_upC2 %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term, colour=over_represented_pvalue, size=numDEInCat)) +
    ggtitle("Enriched GO term of up-regulated genes") +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Gene ratio", y="Biological process GO term", colour="p value", size="Count")

GOplot_downC2 <- GOres_downC2 %>% 
  top_n(20, wt=-over_represented_pvalue) %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>% 
  ggplot(aes(x=hitsPerc, y=term, colour=over_represented_pvalue, size=numDEInCat)) +
    ggtitle("Enriched GO term of down-regulated genes") +
    geom_point() +
    expand_limits(x=0) +
    labs(x="Gene ratio", y="Biological process GO term", colour="p value", size="Count")

png("./Publication_figure/Fig6h_GOplot_.png", units = "in", height = 8, width = 12, res = 600)
GOplot_upC2 | GOplot_downC2
dev.off()
