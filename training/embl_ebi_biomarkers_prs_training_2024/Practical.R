###################################################
###################################################
####        EMBL-EBI Dementia Course		   ####
#### 			Biomarkers and PRS			   ####
#### Emily Simmonds & Valentina Escott-Price   ####
#### 		UKDRI at Cardiff University		   ####
#### 			10th April 2024				   ####
###################################################
###################################################

###################################################
## Set working directory and install libraries   ##
###################################################

setwd("/home/training/VEP_ES_PRS/")

library(ROCR)
library(pROC)
library(data.table)

###################################################

###################################################
## 			Compute PRS, using PRS(C+T)  		 ##
###################################################

### Clump genotype data


command <- sprintf("cd /home/training/VEP_ES_PRS/")
system(command)

command <- sprintf("plink --bfile AD_GWAS_data_final --clump Kunkle_summary_stats.txt --clump-r2 0.1 --clump-kb 1000 --out AD_GWAS_data_final_clumped")
system(command)

command <- sprintf("plink --bfile AD_GWAS_data_final --extract AD_GWAS_data_final_clumped.clumped --make-bed --out AD_GWAS_data_final_clumped")
system(command)


### Subset summary statistics with pT=0.1 and 0.0001

library(data.table)
kun <- fread("Kunkle_summary_stats.txt", header=T, data.table=F)
write.table(kun[which(kun$P <= 0.1), c(2,4,6)], "Kunkle_summary_stats_p0.1.txt", col.names=F, row.names=F, quote=F, sep="\t") 
write.table(kun[which(kun$P <= 0.00001), c(2,4,6)], "Kunkle_summary_stats_p0.00001.txt", col.names=F, row.names=F, quote=F, sep="\t") 


### Compute PRS

command <- sprintf("plink --bfile AD_GWAS_data_final_clumped --score Kunkle_summary_stats_p0.1.txt --out AD_GWAS_data_final_clumped_p0.1")
system(command)

command <- sprintf("plink --bfile AD_GWAS_data_final_clumped --score Kunkle_summary_stats_p0.00001.txt --out AD_GWAS_data_final_clumped_p0.00001")
system(command)

###################################################

###################################################
## 			Compute PRS, using PRS-CS   		 ##
###################################################

### Output summary stats in correct format

kun <- fread("Kunkle_summary_stats.txt", header=T, data.table=F)
write.table(kun[,c(2,4,5,6,7)], "Kunkle_summary_stats_PRScs.txt", col.names=T, row.names=F, quote=F, sep="\t") 


### Use PRS-CS with default parameters to compute adjusted effect sizes

command <- sprintf("python /home/training/PRScs/PRScs.py --ref_dir /home/training/PRScs/ldblk_ukbb_eur --bim_prefix=AD_GWAS_data_final --sst_file=Kunkle_summary_stats_PRScs.txt --n_gwas=63926 --out_dir=Kunkle_PRScs_weights")
system(command)



### Compute PRS score

summ <- data.frame(matrix(nrow=0, ncol=6))
for (chr in 1:22){
	data <- fread(paste("Kunkle_PRScs_weights_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt", sep=""), header=F, data.table=F);print(nrow(data))
	summ <- rbind(summ, data)
}
write.table(summ[, c(2,4,6)], "Kunkle_PRScs_weights.txt", col.names=F, row.names=F, quote=F, sep="\t")

command <- sprintf("plink --bfile AD_GWAS_data_final --score Kunkle_PRScs_weights.txt --out AD_GWAS_data_final_PRScs")
system(command)


###################################################

###################################################
## 	    	Merge PRS data and covariates   		 ##
###################################################

### Merge PRS files

prs0.1 <- fread("AD_GWAS_data_final_clumped_p0.1.profile", header=T, data.table=F)
prs0.1 <- prs0.1[,c(1,2,6)]
colnames(prs0.1)[3] <- "PRSp0.1"

prs0.00001 <- fread("AD_GWAS_data_final_clumped_p0.00001.profile", header=T, data.table=F)
prs0.00001 <- prs0.00001[,c(1,2,6)]
colnames(prs0.00001)[3] <- "PRSp0.00001"

prscs <- fread("AD_GWAS_data_final_PRScs.profile", header=T, data.table=F)
prscs <- prscs[,c(1,2,6)]
colnames(prscs)[3] <- "PRScs"

prs <- merge(prs0.1, prs0.00001, by=c("FID", "IID"), sort=F, all=F)
prs <- merge(prs, prscs, by=c("FID", "IID"), sort=F, all=F)


### Adjust PRS for principal components (PCs)

cov <- fread("AD_GWAS_covars.txt", header=T, data.table=F)

m <- merge(cov, prs, by=c("FID", "IID"), sort=F, all=F)

fit <- glm(PRSp0.1 ~ PC1 + PC2 + PC3, family="gaussian", data=m)
m$PRSp0.1adj <- fit$residuals
m$PRSp0.1adjNorm <- (m$PRSp0.1adj- mean(m$PRSp0.1adj))/sd(m$PRSp0.1adj)

fit <- glm(PRSp0.00001 ~ PC1 + PC2 + PC3, family="gaussian", data=m)
m$PRSp0.00001adj <- fit$residuals
m$PRSp0.00001adjNorm <- (m$PRSp0.00001adj- mean(m$PRSp0.00001adj))/sd(m$PRSp0.00001adj)

fit <- glm(PRScs ~ PC1 + PC2 + PC3, family="gaussian", data=m)
m$PRScsadj <- fit$residuals
m$PRScsadjNorm <- (m$PRScsadj- mean(m$PRScsadj))/sd(m$PRScsadj)

###################################################

###################################################
## Look at descriptive statistics for covariates ##
###################################################

summary(m[,c(5,10:15)])

###################################################

###################################################
## 			Prediction Accuracy of PRS           ##
###################################################


### Compute AUC for three PRSes

m$AD <- m$AD - 1

model <- glm(AD ~ PRSp0.1adjNorm, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ PRSp0.00001adjNorm, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ PRScsadjNorm, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")


### Keep most predictive PRS

m <- m[, c(1:6,10:15,24)]

###################################################

###################################################
## 		Prediction Accuracy of Biomarkers        ##
###################################################

### Compute AUC for biomarkers

model <- glm(AD ~ AB40, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ AB42, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ GFAP, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ NfL, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ Ptau181, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ AB40_42, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

###################################################

###################################################
## 				Stepwise Regression    		     ##
###################################################

library(MASS)

full.model <- glm(AD ~ SEX + AGE + APOEe4 + AB40 + AB42 + GFAP +NfL + Ptau181 + AB40_42 + PRScsadjNorm , family="binomial", data=m)  
step.model <- stepAIC(full.model, direction="both", trace=F)
summary(step.model)

auc(step.model$y, step.model$fitted.values, levels = c(0, 1), direction = "<")



###################################################

###################################################
## 	Remove APOE region, compute PRS.no.APOE      ##
###################################################

kun <- fread("Kunkle_summary_stats.txt", header=T, data.table=F)
a <- which(kun$CHR==19 & kun$BP >= 44400000 & kun$BP <= 46500000)
write.table(kun[-a ,c(2,4,5,6,7)], "Kunkle_noAPOE_summary_stats_PRScs.txt", col.names=T, row.names=F, quote=F, sep="\t") 


command <- sprintf("python /home/training/PRScs/PRScs.py --ref_dir /home/training/PRScs/ldblk_ukbb_eur --bim_prefix=AD_GWAS_data_final --sst_file=Kunkle_noAPOE_summary_stats_PRScs.txt --n_gwas=63926 --out_dir=Kunkle_noAPOE_PRScs_weights")
system(command)


### Compute PRS score

summ <- data.frame(matrix(nrow=0, ncol=6))
for (chr in 1:22){
	data <- fread(paste("Kunkle_noAPOE_PRScs_weights_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt", sep=""), header=F, data.table=F);print(nrow(data))
	summ <- rbind(summ, data)
}
write.table(summ[, c(2,4,6)], "Kunkle_noAPOE_PRScs_weights.txt", col.names=F, row.names=F, quote=F, sep="\t")

command <- sprintf("plink --bfile AD_GWAS_data_final --score Kunkle_noAPOE_PRScs_weights.txt --out AD_GWAS_noAPOE_data_final_PRScs")
system(command)


prscs <- fread("AD_GWAS_noAPOE_data_final_PRScs.profile", header=T, data.table=F)
prscs <- prscs[,c(1,2,6)]
colnames(prscs)[3] <- "PRScs.noAPOE"

### Adjust PRS for principal components (PCs)

cov <- fread("AD_GWAS_covars.txt", header=T, data.table=F)

m <- merge(cov, prscs, by=c("FID", "IID"), sort=F, all=F)

fit <- glm(PRScs.noAPOE ~ PC1 + PC2 + PC3, family="gaussian", data=m)
m$PRScs.noAPOEadj <- fit$residuals
m$PRScs.noAPOEadjNorm <- (m$PRScs.noAPOEadj- mean(m$PRScs.noAPOEadj))/sd(m$PRScs.noAPOEadj)


m$AD <- m$AD - 1

model <- glm(AD ~ PRScs.noAPOEadjNorm, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")

model <- glm(AD ~ APOEe4 + PRScs.noAPOEadjNorm, family="binomial", data=m)
auc(model$y, model$fitted.values, levels = c(0, 1), direction = "<")



