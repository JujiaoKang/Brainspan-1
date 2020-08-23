#library(NetBID2)
load("/cluster/home/chenhy/project/Brainspan/db/brainspan.RData")
pre_phe <- pData(pre_eset)
pre_phe$week <- pre_phe$age
w1 <- grep('pcw',pre_phe$age);pre_phe$week[w1] <- as.numeric(gsub("(.*) pcw","\\1",pre_phe$age[w1]))
w1 <- grep('mos',pre_phe$age);pre_phe$week[w1] <- 42+4*as.numeric(gsub("(.*) mos","\\1",pre_phe$age[w1]))
w1 <- grep('yrs',pre_phe$age);pre_phe$week[w1] <- 42+52*as.numeric(gsub("(.*) yrs","\\1",pre_phe$age[w1]))
pre_phe$week_numeric <- as.numeric(pre_phe$week)
pre_phe$cortex <- ifelse(grepl("cortex",pre_phe$structure_name),'cortex',pre_phe$structure_name)
pre_phe$period <- ifelse(pre_phe$week_numeric<=9,'Embryo',
                         ifelse(pre_phe$week_numeric<=14,'Prenatal_1st-Trimester',
                                ifelse(pre_phe$week_numeric<=28,'Prenatal_2nd-Trimester',
                                       ifelse(pre_phe$week_numeric<=42,'Prenatal_3rd-Trimester',
                                              ifelse(pre_phe$week_numeric<=42+52*13,'Childhood','Adulthood')))))
