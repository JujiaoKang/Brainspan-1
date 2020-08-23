############### Step 0: Preparation ###############
# Load NetBID2 package
library(NetBID2)
library(tidyverse)
# Try parallel
# library(BiocParallel)
# library(BatchJobs)
# param <- MulticoreParam(workers = 8)
# Redefine the load and save function
load.list <- function(analysis.par,analysis.par.name,step){NetBID.loadRData(analysis.par=analysis.par,step=sprintf("%s-%s",step,analysis.par.name))}
save.list <- function(analysis.par,analysis.par.name,step){NetBID.saveRData(analysis.par=analysis.par,step=sprintf("%s-%s",step,analysis.par.name))}
#Create initial analysis.par.list
create_par.list <- function(){
  if(T){# Get network data
  network.dir <- "./network/brainspan"
  # Define main working directory and project name
  project_main_dir <- 'analysis/' 
  project_name <- 'brainspan' 
  sjar_dir <- dir(sprintf("%s/SJAR",network.dir))
  analysis.par.list <- list()
  for (i in 1:length(sjar_dir)) {
    network.project.name <- sjar_dir[i]
    analysis.par.list[[i]]  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                                network_dir=network.dir, network_project_name=network.project.name,
                                                tf.network.file=sprintf("%s/SJAR/%s/structure_acronym_TF/consensus_network_ncol_.txt",network.dir,network.project.name),
                                                sig.network.file=sprintf("%s/SJAR/%s/structure_acronym_SIG/consensus_network_ncol_.txt",network.dir,network.project.name))
  }
  names(analysis.par.list) <- sjar_dir
    }
  return(analysis.par.list)
}
############### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############
analysis.par.list <- create_par.list()
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) # RData saved after QC in the network construction step
phe <- pData(network.par$net.eset)
mat <- exprs(network.par$net.eset)
for (i in sjar_dir) {
  sub_phe <- phe[which(phe$structure_acronym==i),]
  sub_mat <- mat[,rownames(sub_phe)]
  analysis.par.list[[i]]$cal.eset <- generate.eset(exp_mat = sub_mat,phenotype_info = sub_phe)
  # Save Step 1 network.par as RData
  NetBID.saveRData(analysis.par=analysis.par.list[[i]],step=sprintf('exp-QC-%s',i))
}
############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

# Reload network.par RData from Step 1
analysis.par.list <- create_par.list()
for (i in sjar_dir) {load.list(analysis.par.list[[i]], i, step='exp-QC');analysis.par.list[[i]] <- analysis.par;rm(analysis.par)}

# Creat QC report for the network
# draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_',html_info_limit=FALSE)
# draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_',html_info_limit=TRUE)

# Merge network function 
merge_cal_eset <- function(analysis.par){
  #suppressMessages(library(NetBID2))
  # Get network information
  analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
  analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)
  analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)
  ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
  # Create eset using activity matrix
  analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                                  feature_info=NULL,annotation_info='activity in net-dataset')
}

analysis.par.list <- lapply(analysis.par.list, merge_cal_eset)
# # QC plot for activity eset
# draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=c("structure_acronym","class1"),do.logtransform=FALSE,prefix='AC_',
             # choose_plot=c('heatmap','pca','density'),pre_define=c('prenatal'='blue','adult'='red'),pca_plot_type='2D.interactive')
# Save Step 2 analysis.par as RData
map2(analysis.par.list, as.list(names(analysis.par.list)), save.list, step='act-get')
############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############

# Reload network.par RData from Step 2
analysis.par.list <- create_par.list()
for (i in sjar_dir) {load.list(analysis.par.list[[i]], i, step='act-get');analysis.par.list[[i]] <- analysis.par;rm(analysis.par)}

# Comparision function
de_da_comp <- function(analysis.par){
  # Create empty list to store comparison result
  analysis.par$DE <- list()
  analysis.par$DA <- list()
  # First comparison: adult vs. prenatal
  comp_name <- 'adult.Vs.prenatal' # Each comparison must has a name
  # Get sample names from each compared group
  phe_info <- pData(analysis.par$cal.eset)
  G1  <- rownames(phe_info)[which(phe_info$`class1`=='adult')] # Experiment group
  G0  <- rownames(phe_info)[which(phe_info$`class1`=='prenatal')] # Control group
  DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='adult',G0_name='prenatal')
  DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='adult',G0_name='prenatal')
  # Save comparison result to list element in analysis.par, with comparison name
  analysis.par$DE[[comp_name]] <- DE_gene_bid
  analysis.par$DA[[comp_name]] <- DA_driver_bid
  return(analysis.par)
}
for (i in sjar_dir) {analysis.par.list[[i]] <- de_da_comp(analysis.par.list[[i]])}
map2(analysis.par.list, as.list(names(analysis.par.list)), save.list, step='act-DA')
# Visualize top drivers
for (i in sjar_dir) {
  draw.NetBID(DA_list=analysis.par.list[[i]]$DA,DE_list=analysis.par.list[[i]]$DE,main_id='adult.Vs.prenatal',
              DA_display_col='P.Value',DE_display_col='logFC',z_col='Z-statistics',
              pdf_file=sprintf('%s/%s_NetBID_TOP_adult.Vs.prenatal.pdf',analysis.par.list[[i]]$out.dir.PLOT,i),text_cex=0.8) # Save as PDF
}

############### Step 4: Generate a master table for drivers (ms-tab) ###############
# Reload analysis.par RData from Step 3
analysis.par.list <- create_par.list()
for (i in sjar_dir) {load.list(analysis.par.list[[i]], i, step='act-DA');analysis.par.list[[i]] <- analysis.par;rm(analysis.par)}

db.preload(use_level='gene',use_spe='human',update=FALSE)
for (i in sjar_dir) {
  # Get all comparison names
  all_comp <- names(analysis.par.list[[i]]$DE) # Users can use index or name to get target ones
  # Prepare the conversion table (OPTIONAL)
  use_genes <- unique(c(analysis.par.list[[i]]$merge.network$network_dat$source.symbol,analysis.par.list[[i]]$merge.network$network_dat$target.symbol))
  transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes)
  analysis.par.list[[i]]$transfer_tab <- transfer_tab
  # Creat the final master table
  for (j in all_comp) {
    analysis.par.list[[i]]$final_ms_tab[[j]] <- generate.masterTable(use_comp=j,DE=analysis.par.list[[i]]$DE,DA=analysis.par.list[[i]]$DA,
                                                           target_list=analysis.par.list[[i]]$merge.network$target_list,
                                                           tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                                           main_id_type='external_gene_name')
  }
}

# Highlight marker genes that may be related to brain development
all_markgene <- c("E2F1", "SETDB1","MTSS1","SSRP1","EGR4","SETDB1","PURA","SLC4A1","FANCG","PCGF5","ADRA1D","NPAS4","PIN1","SOX12","TEF","EIF4A3","IGF2BP2")
mark_gene <- list(prenatal=all_markgene,
                  adult=all_markgene)
# Customize highlight color codes
mark_col <- get.class.color(names(mark_gene),
                            pre_define=c('prenatal'='blue','adult'='red'))
# Save the final master table as EXCEL file
for (i in sjar_dir) {
    # Path and file name of the output EXCEL file
    out_file <- sprintf('%s/ms_tab_%s.xlsx',analysis.par.list[[i]]$out.dir.DATA,i)
    out2excel(analysis.par.list[[i]]$final_ms_tab,out.xlsx = out_file,mark_gene = mark_gene,mark_col = mark_col)
}
# Save Step 4 analysis.par as RData, ESSENTIAL
map2(analysis.par.list, as.list(names(analysis.par.list)), save.list, step='ms-tab')


# Look the overlap driver
tmp <- analysis.par.list[[1]]$final_ms_tab[[1]][1:1000,"geneSymbol"]
for (i in 2:16) {
  tmp <- intersect(tmp,analysis.par.list[[i]]$final_ms_tab[[1]][1:1000,"geneSymbol"])
}
