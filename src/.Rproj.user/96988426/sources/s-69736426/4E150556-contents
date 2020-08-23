load("RiboHistGeneCounts.rda")#load GeneCount file
Pheno<-read.table("RiboPheno.txt", header = T)#load demographic/technical info for subjects

#remove "," from numbers to use
Pheno$TotalNumMapped<-as.numeric(gsub(",", "", Pheno$TotalNumMapped)) 

#parse out gene counts and gene annotation info
Geneanno = GeneCounts[,1:6]
gCounts = GeneCounts[, 7:length(colnames(GeneCounts))]

#get RPKM values
genebgE= matrix(rep(Pheno$TotalNumMapped), nc = nrow(Pheno), nr = nrow(GeneCounts), byrow = T)
genewidE = matrix(rep(Geneanno$Length), nr = nrow(GeneCounts), 
                  nc = nrow(Pheno),  byrow=F)
geneRpkm<- GeneCounts/(genewidE/1000)/(genebgE/1e6)
geneRpkm_withanno <- data.frame(Geneanno, geneRpkm)

#get log2(Rpkm values +1)
geneRpkm2 <- data.frame(Geneanno, log2(geneRpkm + 1))

save(Pheno, geneRpkm2, geneRpkm, geneRpkm_withanno, file = "geneRpkm2.rda")
