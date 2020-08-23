#load ac_mat and ms_tab
source("/cluster/home/chenhy/project/Brainspan/src/get_acmat_and_mstab.R")
#get sig_driver
sig_driver <- draw.volcanoPlot(dat=DE,label_col='ID',
                               logFC_col='logFC',
                               Pv_col='P.Value',
                               logFC_thre=0.5,
                               Pv_thre=1e-2,
                               main='Volcano Plot for ASD.Vs.TD',
                               show_label=F,
                               label_type = 'origin',
                               label_cex = 0.5,
                               show_plot = T)
sig_driver <- sig_driver[1:300,]


#Plot the heatmap of GSE113834
gse_phe <- pData(mstab_analysis.par$merge.ac.eset)
exp_mat <- Biobase::exprs(mstab_analysis.par$cal.eset)
draw.heatmap(mat=exp_mat,use_genes=sig_driver$ID,
             use_gene_label=sig_driver$ID,
             use_samples=colnames(exp_mat),
             #use_sample_label=phe_info[colnames(ac_mat),'donor_id'],
             phenotype_info=gse_phe,use_phe=c('group'),
             main='Activity for Top drivers',scale='row',
             show_column_names = FALSE,
             cluster_rows=TRUE,cluster_columns=TRUE,
             clustering_distance_rows='euclidean',
             clustering_distance_columns='euclidean',
             clustering_method_rows = "complete",
             row_names_gp = gpar(fontsize = 10),
             show_row_names = TRUE,
             #show_column_names = FALSE,
             col=circlize::colorRamp2(breaks=seq(-2,2,length.out=10),
                                      colors=rev(brewer.pal(10,'Spectral')))
             #pre_define=c('ASD'='blue','control'='red')
)

draw.heatmap(mat=ac_mat,use_genes=ms_tab[rownames(sig_driver[which(sig_driver$logFC>0),]),'external_gene_name'],
             use_gene_label=ms_tab[rownames(sig_driver[which(sig_driver$logFC>0),]),'external_gene_name'],
             use_samples=colnames(ac_mat),
             #use_sample_label=phe_info[colnames(ac_mat),'donor_id'],
             phenotype_info=pre_phe,use_phe=c('gender','class1'),
             main='Activity for Top drivers',scale='row',
             show_column_names = FALSE,
             cluster_rows=TRUE,cluster_columns=FALSE,
             clustering_distance_rows='euclidean',
             clustering_distance_columns='euclidean',
             clustering_method_rows = "complete",
             row_names_gp = gpar(fontsize = 10),
             show_row_names = TRUE,
             col=circlize::colorRamp2(breaks=seq(-2,2,length.out=10),
                                      colors=rev(brewer.pal(10,'Spectral')))
             #pre_define=c('ASD'='blue','control'='red')
             )


