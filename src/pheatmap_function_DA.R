plot_gse_pheatmap_DA <- function(gse,gse_DA,status,use_phe){
  
  analysis.par <- gse
  gse_phe <- pData(analysis.par$merge.ac.eset)
  ac_mat <- Biobase::exprs(analysis.par$merge.ac.eset)
  
  draw.heatmap(mat=ac_mat,use_genes=gse_DA$originalID_label,
               use_gene_label=gse_DA$gene_label,
               use_samples=colnames(ac_mat),
               #use_sample_label=phe_info[colnames(ac_mat),'donor_id'],
               phenotype_info=gse_phe,use_phe=use_phe,
               main=sprintf('Activity for %s overlap top drivers',status),scale='row',
               show_column_names = FALSE,
               cluster_rows=TRUE,cluster_columns=TRUE,
               clustering_distance_rows='euclidean',
               clustering_distance_columns='euclidean',
               clustering_method_rows = "complete",
               row_names_gp = gpar(fontsize = 10),
               show_row_names = FALSE,
               col=circlize::colorRamp2(breaks=seq(-2,2,length.out=10),
                                        colors=rev(brewer.pal(10,'Spectral')))
               #pre_define=c('ASD'='blue','control'='red')
  )
}