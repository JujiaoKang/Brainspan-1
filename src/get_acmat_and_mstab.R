# get ac_mat
library(NetBID2)
library(tidyverse)
library(stringr)
ac_analysis.par <- list()
ac_project_main_dir <- '/cluster/home/chenhy/project/Brainspan/network/' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
ac_project_name <- "brainspan_prefrontal/DATA" # project name for the project folders under main directory.
ac_analysis.par$out.dir.DATA <- str_c(ac_project_main_dir,ac_project_name)
NetBID.loadRData(analysis.par = ac_analysis.par,step='act-get')
ac_analysis.par <- analysis.par
rm(analysis.par)
ac_eset <- ac_analysis.par$merge.ac.eset
phe_info <- Biobase::pData(ac_eset) ## phenotype information
ac_mat <- Biobase::exprs(ac_eset)


# get ms_tab
gse <- 'GSE113834'
gpl <- 'GPL15207'
mstab_analysis.par <- list()
mstab_project_name <- sprintf('Autism_%s_%s',gse,gpl) 
mstab_analysis.par$out.dir.DATA <- sprintf('/cluster/home/chenhy/project/neuron_drivers_landscape/data/%s/DATA/',mstab_project_name)
NetBID.loadRData(analysis.par=mstab_analysis.par,step='ms-tab')
mstab_analysis.par <- analysis.par
rm(analysis.par)
DE<- mstab_analysis.par$DE$ASD.Vs.TD
ms_tab <- mstab_analysis.par$final_ms_tab
