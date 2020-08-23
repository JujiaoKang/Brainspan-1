#load ac_mat and ms_tab
source("/cluster/home/chenhy/project/Brainspan/src/get_acmat_and_mstab.R")
#get sig_driver
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',
                               logFC_col='logFC.ASD.Vs.TD_DA',
                               Pv_col='P.Value.ASD.Vs.TD_DA',
                               logFC_thre=0.05,
                               Pv_thre=1e-3,
                               main='Volcano Plot for ASD.Vs.TD',
                               show_label=FALSE,
                               label_type = 'origin',
                               label_cex = 0.5,
                               show_plot = FALSE)
gs.preload(use_spe='Homo sapiens',update=FALSE)
res1 <- funcEnrich.Fisher(input_list=ms_tab[rownames(sig_driver),'geneSymbol'],
                          bg_list=ms_tab[,'geneSymbol'],
                          use_gs=c('H','C5'),
                          Pv_thre=0.1,Pv_adj = 'none')
res2 <- res1[str_detect(rownames(res1),c('CELL','NEU')),]
draw.funcEnrich.bar(funcEnrich_res=res1,top_number=20,
                    main='Function Enrichment for Top drivers',
                    gs_cex=0.7,gene_cex=0.7,
                    pdf_file="../result/plot/funcEnrich_bar.pdf")
draw.funcEnrich.bar(funcEnrich_res=res1,top_number=5,
                    main='Function Enrichment for Top drivers',
                    display_genes = TRUE,eg_num=3,
                    gs_cex=0.5,gene_cex=0.5)
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=30,gs_cex = 0.8,
                        gene_cex=0.8,pv_cex=0.8,pdf_file="../result/funcEnrich_cluster.pdf")
DA_Z <- z2col(ms_tab[rownames(sig_driver),'Z.ASD.Vs.TD_DA'],
              blue_col=brewer.pal(9,'Blues')[3],
              red_col=brewer.pal(9,'Reds')[3],col_max_thre=6)
names(DA_Z) <- ms_tab[rownames(sig_driver),'geneSymbol']
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=20,gs_cex = 0.5,
                        gene_cex=0.9,pv_cex=0.8,inner_color=DA_Z,pdf_file="../result/plot/funcEnrich_cluster.pdf")
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=20,gs_cex = 0.6,
                        gene_cex=1,pv_cex=1,
                        cluster_gs=TRUE,cluster_gene = TRUE,
                        pdf_file="../result/plot/funcEnrich_cluster.pdf")
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=15,gs_cex = 0.8,
                        gene_cex=0.9,pv_cex=0.8,
                        cluster_gs=TRUE,cluster_gene = FALSE,pdf_file="./result/funcEnrich_cluster.pdf")
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=20,gs_cex = 0.8,
                        gene_cex=0.9,pv_cex=0.8,
                        cluster_gs=FALSE,cluster_gene = TRUE
                        ,pdf_file="./result/funcEnrich_cluster.pdf")
draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=20,gs_cex = 1,
                        gene_cex=1,pv_cex=0.8,
                        cluster_gs=FALSE,cluster_gene = FALSE
                        ,pdf_file="../result/funcEnrich_cluster.pdf")
