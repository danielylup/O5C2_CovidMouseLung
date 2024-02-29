library(DESeq2)
library(ggplot2)
library(pheatmap)
library(patchwork)

##Supplementary Fig. 10A
png("./Bulk_seq/Publication_figure/FigS10a_PCA_C2.png", units="in", width=4, height=3, res=600)
ggplot(bulkCML_GG_pca_C2, aes(PC1, PC2, color=Treatment)) +
    geom_point(size=3) +
    scale_color_manual(values = c("cyan", "salmon")) +
    xlim(-30, 60) +
    xlab(paste0("PC1: ",percentVar_C2[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar_C2[2],"% variance")) +
    guides(col = guide_legend(order = 1)) +
    geom_text(aes(label = rownames(bulkCML_GG_sample_C2)), position = position_nudge(y=2.5), color = "black", size = 2.5) +
    coord_fixed() +
    theme_bw() +
    theme(axis.title = element_text(size = 8), legend.key.height = unit(0.4, 'cm'), legend.title = element_text(size = 8), legend.text = element_text(size = 7))
dev.off())

##Supplementary Fig. 10B
PBSdotplot <- ggplot(bulkCML_GG_deseq2_C2, aes(x=`PBS-1`, y=`PBS-2`)) +
    geom_point(size = 0.5) +
    ggtitle("PBS") +
    scale_fill_viridis_c(trans = "log") +
    xlab("PBS-1") +
    ylab("PBS-2") +
    geom_text(x = 9, y = 19, size = 2.5, label = paste0("italic(R) == ", "0.991"), parse = TRUE) +
    geom_smooth(method = "lm", color = "red", linewidth = 0.5) +
    theme_bw() +
    theme(legend.position = "None", axis.text = element_text(size = 8), plot.title = element_text(size = 8))

O5C2dotplot <- ggplot(bulkCML_GG_deseq2_C2, aes(x=`O5C2-1`, y=`O5C2-2`)) +
    geom_point(size = 0.5) +
    ggtitle("O5C2") +
    scale_fill_viridis_c(trans = "log") +
    xlab("O5C2-1") +
    ylab("O5C2-2") +
    geom_text(x = 9, y = 19, size = 2.5, label = paste0("italic(R) == ", "0.955"), parse = TRUE) +
    geom_smooth(method = "lm", color = "red", linewidth = 0.5) +
    theme_bw() +
    theme(legend.position = "None", axis.text = element_text(size = 8), plot.title = element_text(size = 8))

png("./Publication_figure/FigS10b_scatter_C2.png", units="in", width=4, height=2.5, res=600)
(PBSdotplot | O5C2dotplot)
dev.off()

##Supplementary Fig. 10C
png("./Publication_figure/FigS10c_DistMatrix_C2.png", units="in", width=4, height=3, res=600)
pheatmap(bulkCML_GG_DistMatrix_C2, 
         clustering_distance_rows = bulkCML_GG_Dists_C2, 
         clustering_distance_cols = bulkCML_GG_Dists_C2, 
         treeheight_row = 20,
         treeheight_col = 20,
         col = colors, 
         show_colnames = FALSE,
         display_numbers = TRUE,
         fontsize_number = 8, 
         fontsize = 8)
dev.off()

##Supplementary Fig. 10D
png("./Publication_figure/FigS10d_clheatmap_C2.png", units="in", width=4, height=7, res=600)
pheatmap(assay(bulkCML_GG_vsd_C2)[select_C2,], cluster_row=FALSE, show_rownames = FALSE, annotation_col = bulkCML_GG_sample_C2, fontsize = 8)
dev.off()
