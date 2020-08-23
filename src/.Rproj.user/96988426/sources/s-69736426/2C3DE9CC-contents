library(easyPubMed)
gene <- openxlsx::read.xlsx("/cluster/home/chenhy/project/autism/result/variant_ewas.xlsx",1)$Gene
article <- data.frame()
for (i in gene) {
  #i <- intersect_gene[1]
  final_df <- NULL
  my_query <- sprintf('autism %s',i)
  my_query <- get_pubmed_ids(my_query)
  my_batches <- if(my_query$Count>1)seq(from = 1, to = my_query$Count, by = 10)
  
  if(!is.null(my_batches)){
    my_abstracts_xml <- lapply(my_batches,  function(i) {
      fetch_pubmed_data(my_query, retmax = 1000, retstart = i)  
    })
    ## 储存为 list
    all_xml <- list()
    for(x in my_abstracts_xml) {
      xx <- articles_to_list(x)
      for(y in xx) {
        all_xml[[(1 + length(all_xml))]] <- y
      }  
    }
    ## max_chars = -1 即提取全部摘要
    final_df <- do.call(rbind, lapply(all_xml, article_to_df, 
                                      max_chars = -1, getAuthors = FALSE))
  }
  if(!is.null(final_df)){
    final_df$gene <- i
    final_df <- final_df[,c("gene","pmid","title","year","jabbrv")]
    article <- plyr::rbind.fill(article,as.data.frame(final_df))
  }
  Sys.sleep(2)
}
write.csv(article,file = "/cluster/home/chenhy/project/autism/result/variant_gene_article.csv",row.names = F)
