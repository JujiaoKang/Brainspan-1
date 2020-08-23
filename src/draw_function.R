draw_boxplot <- function(gse_data,group_var,gene_symbol,gse_name){
  phe_info <- Biobase::pData(gse_data$cal.eset)
  group_var <- group_var[group_var %in% colnames(phe_info)]
  expr_dt <- exprs(gse_data$cal.eset) %>% .[which(rownames(.) == gene_symbol),] %>% as.matrix()
  expr_dt <- cbind(expr_dt,pData(gse_data$cal.eset)[,group_var]) %>% as.data.frame()
  colnames(expr_dt) <- c("y","group");expr_dt$y <- as.numeric(expr_dt$y)
  expr_p_value <- gse_data$final_ms_tab %>% 
    .[which(.$geneSymbol==gene_symbol),"P.Value.ASD.Vs.TD_DE"] %>% round(3)
  
  expr_box <- ggplot(expr_dt,aes(x=group,y=y))+ylab("Expression value")+
    geom_boxplot(aes(color=group))+geom_point(position = position_jitter(0.1))+ggtitle(str_c(gene_symbol,"-",gse_name))+
    geom_segment(aes(x=1,y=max(expr_dt$y)+mean(expr_dt$y)/2,xend=1,yend=max(expr_dt$y)+mean(expr_dt$y)))+
    geom_segment(aes(x=2,y=max(expr_dt$y)+mean(expr_dt$y)/2,xend=2,yend=max(expr_dt$y)+mean(expr_dt$y)))+
    geom_segment(aes(x=1,y=max(expr_dt$y)+mean(expr_dt$y),xend=2,yend=max(expr_dt$y)+mean(expr_dt$y)))+
    annotate("text",x=1.5,y=max(expr_dt$y)+mean(expr_dt$y)*1.3,label=sprintf("P = %s",expr_p_value),size=5)+
    theme_classic()
  
  ac_dt <- exprs(gse_data$merge.ac.eset) %>% 
    .[which(str_split(rownames(.),"_",simplify = T)[,1] == gene_symbol),] %>% as.matrix() %>% .[,1]
  ac_dt <- cbind(ac_dt,pData(gse_data$cal.eset)[,group_var]) %>% as.data.frame()
  colnames(ac_dt) <- c("y","group");ac_dt$y <- as.numeric(as.character(ac_dt$y))
  ac_p_value <- gse_data$final_ms_tab %>% 
    .[which(.$geneSymbol==gene_symbol),"P.Value.ASD.Vs.TD_DA"] %>% round(3)
  
  ac_box <- ggplot(ac_dt,aes(x=group,y=y))+ylab("Activity value")+
    geom_boxplot(aes(color=group))+geom_point(position = position_jitter(0.1))+ggtitle(str_c(gene_symbol,"-",gse_name))+
    geom_segment(aes(x=1,y=abs(max(ac_dt$y))*1.1,xend=1,yend=abs(max(ac_dt$y))*1.3))+
    geom_segment(aes(x=2,y=abs(max(ac_dt$y))*1.1,xend=2,yend=abs(max(ac_dt$y))*1.3))+
    geom_segment(aes(x=1,y=abs(max(ac_dt$y))*1.3,xend=2,yend=abs(max(ac_dt$y))*1.3))+
    annotate("text",x=1.5,y=abs(max(ac_dt$y))*1.6,label=sprintf("P = %s",ac_p_value),size=5)+
    theme_classic()
  box_map <- ggpubr::ggarrange(expr_box,ac_box,widths = c(1,1), align = "hv",hjust = 10,vjust = 10)
  
  return(box_map)
}

draw_network <- function(gse_data,gene_symbol){
  use_driver <- exprs(gse_data$merge.ac.eset) %>% rownames(.) %>% 
    .[which(str_split(.,"_",simplify = T)[,1] == gene_symbol)] %>%.[1]
  analysis.par <- gse_data
  comp_name <- "ASD.Vs.TD"
  # Define edges of the network
  edge_score <- analysis.par$merge.network$target_list[[use_driver]]$MI*sign(analysis.par$merge.network$target_list[[use_driver]]$spearman)
  names(edge_score) <- analysis.par$merge.network$target_list[[use_driver]]$target
  
  # Draw sub-network structure of selected driver
  draw.targetNet(source_label=analysis.par$final_ms_tab[use_driver,'gene_label'],
                 source_z=analysis.par$final_ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],
                 edge_score = edge_score,pdf_file=sprintf('/cluster/home/chenhy/project/autism/result/plot/%s/targetNet_out.pdf',gene_symbol),
                 label_cex = 1,n_layer=2, alphabetical_order=TRUE)
  
  draw.targetNet(source_label=analysis.par$final_ms_tab[use_driver,'gene_label'],
                 source_z=analysis.par$final_ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],
                 edge_score = edge_score,pdf_file=sprintf('/cluster/home/chenhy/project/autism/result/plot/%s/targetNet_in.pdf',gene_symbol),
                 label_cex = 1,arrow_direction = 'in',n_layer=2)
  
  target_gene <- analysis.par$merge.network$target_list[[use_driver]]$target
  res1 <- funcEnrich.Fisher(input_list=target_gene,
                            bg_list=rownames(exprs(analysis.par$cal.eset)),
                            use_gs=c('H','C5'),Pv_thre=0.1,Pv_adj = 'none')
  draw.funcEnrich.cluster(funcEnrich_res=res1,top_number=20,gs_cex = 0.5,
                          gene_cex=0.9,pv_cex=0.8,cluster_gs=TRUE,cluster_gene = TRUE,
                          pdf_file=sprintf('/cluster/home/chenhy/project/autism/result/plot/%s/targetGene_funEnrich.pdf',gene_symbol))
  
}

draw_pattern <- function(pattern_ac,pattern_expr,gene_symbol){
  ### The expression and activity pattern of driver in Brainspan
  #scatter plot for activity of gene 
  point_plot_ac <- ggplot(data=pattern_ac)+
    geom_point(aes(x=age_num,y=value,color = class1),size=1.5,alpha=.98)+
    scale_color_manual(values = c("red","blue","green","purple","orange", "brown"))+
    geom_smooth(aes(x=age_num,y=value),method='loess',se = T,color = "black",size=0.3,span=0.25)+
    guides(shape = guide_legend(override.aes = list(size = 2)),color = guide_legend(override.aes = list(size = 2)))+
    scale_x_continuous(limits=c(min(pattern_ac$age_num),max(pattern_ac$age_num)),breaks=unique(pattern_ac$age_num),labels=unique(pattern_ac$age))+
    scale_y_continuous(limits=c(min(pattern_ac[,"value"]),max(pattern_ac[,"value"])),breaks=round(seq(min(pattern_ac[,"value"]),max(pattern_ac[,"value"]),length.out=7),1))+
    labs(x="Age",y="Activity level",title = sprintf("Activity level for %s in brain with different period",gene_symbol))+
    theme(legend.title = element_blank(),legend.position = c(0.85, 0.8),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(panel.background=element_rect(fill='transparent',color='transparent'),
          legend.key=element_rect(fill='transparent', color='transparent'),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(linetype="solid",fill = "transparent"))
  ggsave(point_plot_ac, file=sprintf("/cluster/home/chenhy/project/autism/result/plot/%s/pattern_activity_plot.pdf",gene_symbol), width=7, height=5.5)
  
  
  #### Expression level
  #scatter plot for expression of gene 
  point_plot_expr <- ggplot(data=pattern_expr)+
    geom_point(aes(x=age_num,y=value,color = class1),size=1.5,alpha=.98)+
    scale_color_manual(values = c("red","blue","green","purple","orange", "brown"))+
    geom_smooth(aes(x=age_num,y=value),method='loess',se = T,color = "black",size=0.3,span=0.25)+
    guides(shape = guide_legend(override.aes = list(size = 2)),color = guide_legend(override.aes = list(size = 2)))+
    scale_x_continuous(limits=c(min(pattern_expr$age_num),max(pattern_expr$age_num)),breaks=unique(pattern_expr$age_num),labels=unique(pattern_expr$age))+
    scale_y_continuous(limits=c(min(pattern_expr[,"value"]),max(pattern_expr[,"value"])),breaks=round(seq(min(pattern_expr[,"value"]),max(pattern_expr[,"value"]),length.out=7),1))+
    labs(x="Age",y="Expression level",title = sprintf("Expression level for %s in brain with different period",gene_symbol))+
    theme(legend.title = element_blank(),legend.position = c(0.85, 0.8),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(panel.background=element_rect(fill='transparent',color='transparent'),
          legend.key=element_rect(fill='transparent', color='transparent'),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(linetype="solid",fill = "transparent"))
  ggsave(point_plot_expr, file=sprintf("/cluster/home/chenhy/project/autism/result/plot/%s/pattern_expr_plot.pdf",gene_symbol), width=7, height=5.5)
  
}



if(F){
  analysis.par <- gse_list[["GSE27919067"]]
  gene_symbol <- candidate_driver[1]
  group_var <- "ASD.CTL"
  ms_tab <- analysis.par$final_ms_tab
  use_driver <- exprs(analysis.par$merge.ac.eset) %>% rownames(.) %>% 
    .[which(str_split(.,"_",simplify = T)[,1] == gene_symbol)] %>% .[1]
  exp_mat <- Biobase::exprs(analysis.par$cal.eset)
  ## expression,the rownames could match originalID
  ac_mat <- Biobase::exprs(analysis.par$merge.ac.eset)
  ## activity,the rownames could match originalID_label
  phe_info <- Biobase::pData(analysis.par$cal.eset)
  use_obs_class <- get_obs_label(phe_info = phe_info,group_var)
  draw.categoryValue(ac_val=ac_mat[use_driver,],
                     exp_val=exp_mat[ms_tab[use_driver,'originalID'],],
                     use_obs_class=use_obs_class,
                     class_srt=30,
                     main_ac = ms_tab[use_driver,'gene_label'],
                     main_exp=ms_tab[use_driver,'geneSymbol'])
}
