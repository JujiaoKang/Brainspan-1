library(NetBID2)
library(tidyverse)
library(stringr)
project_main_dir <- '/cluster/home/chenhy/project/Brainspan/network' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
project_name <- "brainspan_prefrontal" # project name for the project folders under main directory.
# Create a hierarchcial working directory and return a list contains the hierarchcial working directory information
# This list object (network.par) is an ESSENTIAL variable in network construction pipeline
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)
NetBID.loadRData(network.par = network.par,step='exp-QC')
net_eset <- network.par$net.eset
phe <- pData(net_eset)

# Load database
db.preload(use_level='transcript',use_spe='human',update=FALSE)

# Converts gene ID into the corresponding TF/SIG list
use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)
#只要前额叶皮层的数据
phe <- phe[which(phe$structure_acronym %in% c("MFC","OFC","DFC","VFC")),]
phe$sample_name <- rownames(phe)
net_eset <- update_eset.phenotype(use_eset=net_eset,
                                  use_phenotype_info=phe,
                                  use_sample_col='sample_name')
network.par$net.eset <- net_eset
use.samples <- rownames(phe) 
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.1,
                 SJAR.project_name="prefrontal",SJAR.main_dir=network.par$out.dir.SJAR)

#===========================================analysis=================================================
######### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############
# Get the demo's constructed network data
network.dir <- sprintf("/cluster/home/chenhy/project/Brainspan/network/%s",project_name) # use demo network in the package
network.project.name <- 'prefrontal'
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)
analysis.par$cal.eset <- network.par$net.eset

# Save Step 1 network.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

# Reload network.par RData from Step 1
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')
# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

# Creat QC report for the network
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_',html_info_limit=FALSE)
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_',html_info_limit=TRUE)

# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')
# Create eset using activity matrix
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')
# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')
