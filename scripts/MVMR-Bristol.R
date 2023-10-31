
library(remotes)
library(devtools)
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
install_github(c("MRCIEU/TwoSampleMR","MRCIEU/MRInstruments"))
install_github("WSpiller/MRPracticals",build_opts = c("--no-resave-data", "--no-manual"),build_vignettes = TRUE)
install.packages("MendelianRandomization")
#library(MRPracticals)
library(plyr)

library(MendelianRandomization)
library(TwoSampleMR)
library(MVMR)


###### details on MVMR 
vignette("MVMR Tutorial")
#https://www.youtube.com/watch?v=xNsXNzUAjSk
# youtube video is a tuorial on MVMR r pakcage

########################

#getting genome wide sig SNPs
dat_X <- extract_instruments(outcomes='ukb-b-5192',clump=F)
ieu_M <- extract_instruments(outcomes='ieu-b-25',clump=F)

#clumping instruments
snps<-rbind(dat_X,ieu_M)
snps<-clump_data(snps) 
#this gives us a list of SNPs which are indepednet and predict either the exposure or meidator

#getting SNP data for the exposure, outcome, and meidator
exp<- extract_outcome_data(snps$SNP,  "ukb-b-5192",  proxies = TRUE,   rsq = 0.8, align_alleles = 1,  palindromes = 1,  maf_threshold = 0.3,  access_token = ieugwasr::check_access_token(),   splitsize = 10000, proxy_splitsize = 500)
exp<- convert_outcome_to_exposure(exp)

med<- extract_outcome_data(snps$SNP,  "ieu-b-25",  proxies = TRUE,   rsq = 0.8, align_alleles = 1,  palindromes = 1,  maf_threshold = 0.3,  access_token = ieugwasr::check_access_token(),   splitsize = 10000, proxy_splitsize = 500)
med<-harmonise_data(exp,med, action = 2 )
med<-rename(med, replace = c("beta.outcome"="b_adj", "se.outcome"="se_adj"))
med<-med[c("b_adj", "SNP", "se_adj")] 

out<- extract_outcome_data(snps$SNP,  "ieu-a-966",  proxies = TRUE,   rsq = 0.8, align_alleles = 1,  palindromes = 1,  maf_threshold = 0.3,  access_token = ieugwasr::check_access_token(),   splitsize = 10000, proxy_splitsize = 500)
out<-harmonise_data(exp,out, action = 2)
out<-out[c("SNP",  "beta.outcome", "se.outcome")]

exp<-rename(exp, replace = c("beta.exposure"="b_exp", "se.exposure"="se_exp"))
exp<-exp[c("b_exp", "SNP", "se_exp")] 


# merging data
mvmr_raw1<-merge(exp, med, by=c("SNP"), all=TRUE) # in order to do this cleanly I have created 3 data frames, each with the SNP beta and exposure varaibesl for one pheotype, and then merged them (got issues when trying to replace existing data)
mvmr_raw1<-merge(out, mvmr_raw1, by=c("SNP"), all=TRUE) # if doing yourself  you want this list to be ranked based of the lowest p value each snp has with any exposure

nrow(mvmr_raw1)

#remvoing any missing data (latter functions don't work with non missing data)
mvmr_raw1<-subset(mvmr_raw1, !is.na(beta.outcome))
mvmr_raw1<-subset(mvmr_raw1, !is.na(b_adj))
mvmr_raw1<-subset(mvmr_raw1, !is.na(b_exp))
nrow(mvmr_raw1)

# formating data for mvmr #n.b. you can't have any missing data in this.
F.data<-format_mvmr( BXGs = mvmr_raw1[,c("b_exp","b_adj")], BYG = mvmr_raw1$beta.outcome,
                     seBXGs = mvmr_raw1[,c("se_exp","se_adj")],   seBYG = mvmr_raw1$se.outcome,  RSID = mvmr_raw1$SNP)


#pheno cov
cov<-0#phenocov_mvmr(matrix(c(?,?,?,?),2), F.data[,c("sebetaX1", "sebetaX2")])
# The F stat makes assumptions about how the variables are correlated. See documentation for more detail. 
# IF estimteated in spearte samples you can assume it to be zero. 
# cov needs the phenoypic correlation between the exposures. these should be ordered x1,x2...
# can assume to be zero when exposures and outcomes are created in their own data sets
# where possible I have found it is better to know the covaraince matrix. 

#weak instrument test
sres<-strength_mvmr(r_input=F.data,gencov=cov) #false postive rate for politropy will be increased by using weak instruments


#testing for plytipy 
pres<-pleiotropy_mvmr(r_input=F.data,gencov=cov) #probably some plitropy given low p value

#MVMR IVW
res<-ivw_mvmr(r_input=F.data, gencov = cov)

# ajust for weak instrument bias where present
#qhet_mvmr(r_input = F.data, pcor=matrix(c(1,? ,? ,1), ncol = 2, nrow = 2), CI=TRUE, iterations=100)
#pcor = correlation between the exposures, should be ordered x1, x2 etc 
# n.b. any file which changes f.data has to be used after the mvmr r packge ones. 

#MVMR egger
mr_mvegger(mr_mvinput(bx = cbind(F.data$betaX1, F.data$betaX2), bxse = cbind(F.data$sebetaX1, F.data$sebetaX2),
                      by = F.data$betaYG, byse = F.data$sebetaYG), orientate = 1)
#MVMR median
mr_mvmedian(mr_mvinput(bx = cbind(F.data$betaX1, F.data$betaX2), bxse = cbind(F.data$sebetaX1, F.data$sebetaX2),
                       by = F.data$betaYG, byse = F.data$sebetaYG), iterations = 1000)
