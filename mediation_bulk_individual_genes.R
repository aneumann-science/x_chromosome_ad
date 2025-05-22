library(data.table)
library(parallel)
library(lavaan)
library(glmnet)
library(psych)

###### Phenotype prep
pheno.data <- fread("phenotype/amp_ad_pheno.csv")

# Center gender and code female as 1, male as 0
pheno.data$gender <- as.numeric(as.factor(pheno.data$Sex)) - 1
pheno.data$gender <- 1-pheno.data$gender

# Specify continuous age variable 
pheno.data$age <- pheno.data$AgeAtDeath
pheno.data$age[pheno.data$age == "90_or_above"] <- 90
pheno.data$age <- as.numeric(pheno.data$age)

# Define ADvsCon as ordered for lavaan
pheno_pathological.data <- pheno.data[pheno.data$FK_final_CaseControl != "Other", ]
pheno_pathological.data$FK_final_CaseControl <- ordered(pheno_pathological.data$FK_final_CaseControl, levels = c("Control","AD"))

# Load imputed genetic components
load("phenotype/gPC_imputed.Rdata")
pheno_pathological.data <- merge(pheno_pathological.data, gPC_imputed.data, by = "specimenID")

# Stratify by cohort
pheno_rosmap.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "ROSMAP_DLPFC", ]
pheno_mayo.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "MayoRNASeq_TCX", ]
pheno_msbb_ba10.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "MSBB_BA10", ]
pheno_msbb_ba22.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "MSBB_BA22", ]
pheno_msbb_ba36.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "MSBB_BA36", ]
pheno_msbb_ba44.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "MSBB_BA44", ]

### Load transcriptomes
rosmap_transcriptome.data <- fread("transcriptome/ROSMAP/rosmap_dlpfc_chrX.txt")
mayo_transcriptome.data <- fread("transcriptome/Mayo/mayo_tcx_chrX.txt")
msbb_ba10_transcriptome.data <- fread("transcriptome/MSBB_BA10/msbb_ba10_chrX.txt")
msbb_ba22_transcriptome.data <- fread("transcriptome/MSBB_BA22/msbb_ba22_chrX.txt")
msbb_ba36_transcriptome.data <- fread("transcriptome/MSBB_BA36/msbb_ba36_chrX.txt")
msbb_ba44_transcriptome.data <- fread("transcriptome/MSBB_BA44/msbb_ba44_chrX.txt")

## Load peer factors
rosmap_peer.data <- fread("expression_matrices_and_PEER/ePEER_60_ROSMAP_DLPFC.autosomes.allSamples_n636.14April2025.PEER_covariates.txt", data.table = F)
mayo_peer.data <- fread("expression_matrices_and_PEER/ePEER_45_MayoRNASeq_TCX.autosomes.allSamples_n264.14April2025.PEER_covariates.txt", data.table = F)
msbb_ba10_peer.data <- fread("expression_matrices_and_PEER/ePEER_30_MSBB_BA10.autosomes.allSamples_n227.14April2025.PEER_covariates.txt", data.table = F)
msbb_ba22_peer.data <- fread("expression_matrices_and_PEER/ePEER_30_MSBB_BA22.autosomes.allSamples_n202.14April2025.PEER_covariates.txt", data.table = F)
msbb_ba36_peer.data <- fread("expression_matrices_and_PEER/ePEER_30_MSBB_BA36.autosomes.allSamples_n180.14April2025.PEER_covariates.txt", data.table = F)
msbb_ba44_peer.data <- fread("expression_matrices_and_PEER/ePEER_30_MSBB_BA44.autosomes.allSamples_n193.14April2025.PEER_covariates.txt", data.table = F)

# Prepare for merging
peer_prep <- function(peer.data) {
  rownames(peer.data) <- peer.data$ID
  peer.data <- peer.data[,2:ncol(peer.data)]
  peer.data <- as.data.frame(t(peer.data))
  peer.data$specimenID <- rownames(peer.data)
  return(peer.data)
}

rosmap_peer.data <- peer_prep(rosmap_peer.data)
mayo_peer.data <- peer_prep(mayo_peer.data)
msbb_ba10_peer.data <- peer_prep(msbb_ba10_peer.data)
msbb_ba22_peer.data <- peer_prep(msbb_ba22_peer.data)
msbb_ba36_peer.data <- peer_prep(msbb_ba36_peer.data)
msbb_ba44_peer.data <- peer_prep(msbb_ba44_peer.data)

# Count number on unique transcripts
length(unique(c(names(rosmap_transcriptome.data),names(mayo_transcriptome.data),names(msbb_ba10_transcriptome.data),names(msbb_ba22_transcriptome.data),names(msbb_ba36_transcriptome.data),names(msbb_ba44_transcriptome.data))))

# Count number of cases and controls across all samples
dim(unique(pheno_pathological.data[pheno_pathological.data$FK_final_CaseControl == "AD", "individualID"]))
dim(unique(pheno_pathological.data[pheno_pathological.data$FK_final_CaseControl == "Control", "individualID"]))

### Rosmap
pheno_transcriptome_rosmap.data <- merge(pheno_rosmap.data, rosmap_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_rosmap.data <- merge(pheno_transcriptome_rosmap.data, rosmap_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
rosmap_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_rosmap.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:60))]), pheno_transcriptome_rosmap.data$FK_final_CaseControl, family = "binomial")
plot(rosmap_lasso.fit)
coef(rosmap_lasso.fit, s = "lambda.1se")

mediation_rosmap.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov2 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov8 + InferredCov10 + InferredCov18 + InferredCov20 + InferredCov24 + InferredCov25 + InferredCov27 + InferredCov29 + InferredCov34 + InferredCov35 + InferredCov36 + InferredCov41 + InferredCov47 + InferredCov49 + InferredCov51 + InferredCov56 + InferredCov60
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov2 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov8 + InferredCov10 + InferredCov18 + InferredCov20 + InferredCov24 + InferredCov25 + InferredCov27 + InferredCov29 + InferredCov34 + InferredCov35 + InferredCov36 + InferredCov41 + InferredCov47 + InferredCov49 + InferredCov51 + InferredCov56 + InferredCov60

  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

rosmap_transcripts <- names(rosmap_transcriptome.data)[2:750]
mediation_rosmap.list <- mclapply(pheno_transcriptome_rosmap.data[,..rosmap_transcripts], function(transcript){
  pheno_transcriptome_rosmap.data$transcript <- transcript
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = pheno_transcriptome_rosmap.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 8)
mediation_rosmap.data <- do.call(rbind, mediation_rosmap.list)

names(mediation_rosmap.data) <- paste(names(mediation_rosmap.data),"rosmap",sep="_")

### Mayo
pheno_transcriptome_mayo.data <- merge(pheno_mayo.data, mayo_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_mayo.data <- merge(pheno_transcriptome_mayo.data, mayo_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
mayo_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_mayo.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:45))]), pheno_transcriptome_mayo.data$FK_final_CaseControl, family = "binomial")
plot(mayo_lasso.fit)
coef(mayo_lasso.fit, s = "lambda.1se")


mediation_mayo.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov2 + InferredCov3 + InferredCov4 + InferredCov5 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov11 + InferredCov14 + InferredCov17 + InferredCov19 + InferredCov21 + InferredCov25
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov2 + InferredCov3 + InferredCov4 + InferredCov5 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov11 + InferredCov14 + InferredCov17 + InferredCov19 + InferredCov21 + InferredCov25
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

mayo_transcripts <- names(mayo_transcriptome.data)[2:900]
mediation_mayo.list <- mclapply(pheno_transcriptome_mayo.data[,..mayo_transcripts], function(transcript){
  pheno_transcriptome_mayo.data$transcript <- transcript
  mediation_mayo.fit <- sem(mediation_mayo.model, data = pheno_transcriptome_mayo.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_mayo.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 8)
mediation_mayo.data <- do.call(rbind, mediation_mayo.list)

names(mediation_mayo.data) <- paste(names(mediation_mayo.data),"mayo",sep="_")

### MSBB BA10
pheno_transcriptome_msbb_ba10.data <- merge(pheno_msbb_ba10.data, msbb_ba10_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba10.data <- merge(pheno_transcriptome_msbb_ba10.data, msbb_ba10_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba10_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba10.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba10.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba10_lasso.fit)
coef(msbb_ba10_lasso.fit, s = "lambda.1se")

mediation_msbb_ba10.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov8 + InferredCov11 + InferredCov12 + InferredCov19 + InferredCov20 + InferredCov21 + InferredCov22 + InferredCov26
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov4 + InferredCov5 + InferredCov6 + InferredCov11 + InferredCov12 + InferredCov19 + InferredCov20 + InferredCov21 + InferredCov22 + InferredCov26
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba10_transcripts <- names(msbb_ba10_transcriptome.data)[2:799]
mediation_msbb_ba10.list <- mclapply(pheno_transcriptome_msbb_ba10.data[,..msbb_ba10_transcripts], function(transcript){
  pheno_transcriptome_msbb_ba10.data$transcript <- transcript
  mediation_msbb_ba10.fit <- sem(mediation_msbb_ba10.model, data = pheno_transcriptome_msbb_ba10.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba10.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 8)
mediation_msbb_ba10.data <- do.call(rbind, mediation_msbb_ba10.list)

names(mediation_msbb_ba10.data) <- paste(names(mediation_msbb_ba10.data),"msbb_ba10",sep="_")

### MSBB BA22
pheno_transcriptome_msbb_ba22.data <- merge(pheno_msbb_ba22.data, msbb_ba22_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba22.data <- merge(pheno_transcriptome_msbb_ba22.data, msbb_ba22_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba22_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba22.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba22.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba22_lasso.fit)
coef(msbb_ba22_lasso.fit, s = "lambda.1se")

mediation_msbb_ba22.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov3 + InferredCov4 + InferredCov11 + InferredCov12 + InferredCov13 + InferredCov14 + InferredCov15 + InferredCov16 + InferredCov19 + InferredCov20 + InferredCov22 + InferredCov23 + InferredCov25 + InferredCov28
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov3 + InferredCov4 + InferredCov11 + InferredCov12 + InferredCov13 + InferredCov14 + InferredCov15 + InferredCov16 + InferredCov19 + InferredCov20 + InferredCov22 + InferredCov23 + InferredCov25 + InferredCov28
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba22_transcripts <- names(msbb_ba22_transcriptome.data)[2:762]
mediation_msbb_ba22.list <- mclapply(pheno_transcriptome_msbb_ba22.data[,..msbb_ba22_transcripts], function(transcript){
  pheno_transcriptome_msbb_ba22.data$transcript <- transcript
  mediation_msbb_ba22.fit <- sem(mediation_msbb_ba22.model, data = pheno_transcriptome_msbb_ba22.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba22.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 8)
mediation_msbb_ba22.data <- do.call(rbind, mediation_msbb_ba22.list)

names(mediation_msbb_ba22.data) <- paste(names(mediation_msbb_ba22.data),"msbb_ba22",sep="_")

### MSBB BA36
pheno_transcriptome_msbb_ba36.data <- merge(pheno_msbb_ba36.data, msbb_ba36_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba36.data <- merge(pheno_transcriptome_msbb_ba36.data, msbb_ba36_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba36_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba36.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba36.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba36_lasso.fit)
coef(msbb_ba36_lasso.fit, s = "lambda.1se")

mediation_msbb_ba36.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov3 + InferredCov5 + InferredCov10 + InferredCov14 + InferredCov15 + InferredCov19 + InferredCov20 + InferredCov22 + InferredCov24 + InferredCov28
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov3 + InferredCov5 + InferredCov10 + InferredCov14 + InferredCov15 + InferredCov19 + InferredCov20 + InferredCov22 + InferredCov24 + InferredCov28
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba36_transcripts <- names(msbb_ba36_transcriptome.data)[2:784]
mediation_msbb_ba36.list <- mclapply(pheno_transcriptome_msbb_ba36.data[,..msbb_ba36_transcripts], function(transcript){ 
  pheno_transcriptome_msbb_ba36.data$transcript <- transcript
  mediation_msbb_ba36.fit <- sem(mediation_msbb_ba36.model, data = pheno_transcriptome_msbb_ba36.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba36.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 8)
mediation_msbb_ba36.data <- do.call(rbind, mediation_msbb_ba36.list)

names(mediation_msbb_ba36.data) <- paste(names(mediation_msbb_ba36.data),"msbb_ba36",sep="_")

### MSBB BA44
pheno_transcriptome_msbb_ba44.data <- merge(pheno_msbb_ba44.data, msbb_ba44_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba44.data <- merge(pheno_transcriptome_msbb_ba44.data, msbb_ba44_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba44_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba44.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba44.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba44_lasso.fit)
coef(msbb_ba44_lasso.fit, s = "lambda.1se")

mediation_msbb_ba44.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov2 + InferredCov4 + InferredCov6 + InferredCov7 + InferredCov11 + InferredCov12 + InferredCov13 + InferredCov15 + InferredCov18 + InferredCov23 + InferredCov25 + InferredCov26 + InferredCov29
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov2 + InferredCov4 + InferredCov6 + InferredCov7 + InferredCov11 + InferredCov12 + InferredCov13 + InferredCov15 + InferredCov18 + InferredCov23 + InferredCov25 + InferredCov26 + InferredCov29
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba44_transcripts <- names(msbb_ba44_transcriptome.data)[2:753]
mediation_msbb_ba44.list <- mclapply(pheno_transcriptome_msbb_ba44.data[,..msbb_ba44_transcripts], function(transcript){
  pheno_transcriptome_msbb_ba44.data$transcript <- transcript
  mediation_msbb_ba44.fit <- sem(mediation_msbb_ba44.model, data = pheno_transcriptome_msbb_ba44.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba44.fit)
  estimates$pvalue <- 2*pnorm(q=abs(estimates$z), lower.tail=FALSE)
  results_gender_transcript <- estimates[estimates$label == "gender_transcript",c("est","se","pvalue")]
  names(results_gender_transcript) <- c("est_gender_transcript","se_gender_transcript","pvalue_gender_transcript")
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  results_mediation <- estimates[estimates$label == "mediation",c("est","se","pvalue")]
  names(results_mediation) <- c("est_mediation","se_mediation","pvalue_mediation")
  results_direct <- estimates[estimates$label == "gender_ad",c("est","se","pvalue")]
  names(results_direct) <- c("est_direct","se_direct","pvalue_direct")
  results_total <- estimates[estimates$label == "total",c("est","se","pvalue")]
  names(results_total) <- c("est_total","se_total","pvalue_total")
  cbind(results_gender_transcript,results_transcript_ad,results_mediation,results_direct,results_total)
}, mc.cores = 8)
mediation_msbb_ba44.data <- do.call(rbind, mediation_msbb_ba44.list)

names(mediation_msbb_ba44.data) <- paste(names(mediation_msbb_ba44.data),"msbb_ba44",sep="_")

# Make column for row.names
# Mediation estimates
mediation.list <- list(mediation_rosmap.data,mediation_mayo.data,mediation_msbb_ba10.data,mediation_msbb_ba22.data,mediation_msbb_ba36.data,mediation_msbb_ba44.data)
mediation.list <- lapply(mediation.list, function(x) data.frame(x, transcript = row.names(x)))

results_all_mediation.data <- Reduce(function(x, y) merge(x, y, by="transcript", all=TRUE), mediation.list)

# Write file for meta analysis
fwrite(results_all_mediation.data[,c("transcript","est_mediation_rosmap","se_mediation_rosmap","est_mediation_mayo","se_mediation_mayo")], file = "meta_input/mediation_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_gender_transcript_rosmap","se_gender_transcript_rosmap","est_gender_transcript_mayo","se_gender_transcript_mayo")], file = "meta_input/gender_transcript_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_transcript_ad_rosmap","se_transcript_ad_rosmap","est_transcript_ad_mayo","se_transcript_ad_mayo")], file = "meta_input/transcript_ad_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_direct_rosmap","se_direct_rosmap","est_direct_mayo","se_direct_mayo")], file = "meta_input/direct_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_total_rosmap","se_total_rosmap","est_total_mayo","se_total_mayo")], file = "meta_input/total_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)

# Perform meta-analysis
java -jar /home/aneumann/Metasoft/Metasoft.jar -pvalue_table /home/aneumann/Metasoft/HanEskinPvalueTable.txt -input meta_input/mediation_rosmap_mayo_pathological.txt -output meta_output/mediation_rosmap_mayo_meta_pathological.txt
java -jar /home/aneumann/Metasoft/Metasoft.jar -pvalue_table /home/aneumann/Metasoft/HanEskinPvalueTable.txt -input meta_input/gender_transcript_rosmap_mayo_pathological.txt -output meta_output/gender_transcript_rosmap_mayo_meta_pathological.txt
java -jar /home/aneumann/Metasoft/Metasoft.jar -pvalue_table /home/aneumann/Metasoft/HanEskinPvalueTable.txt -input meta_input/transcript_ad_rosmap_mayo_pathological.txt -output meta_output/transcript_ad_rosmap_mayo_meta_pathological.txt
java -jar /home/aneumann/Metasoft/Metasoft.jar -pvalue_table /home/aneumann/Metasoft/HanEskinPvalueTable.txt -input meta_input/direct_rosmap_mayo_pathological.txt -output meta_output/direct_rosmap_mayo_meta_pathological.txt
java -jar /home/aneumann/Metasoft/Metasoft.jar -pvalue_table /home/aneumann/Metasoft/HanEskinPvalueTable.txt -input meta_input/total_rosmap_mayo_pathological.txt -output meta_output/total_rosmap_mayo_meta_pathological.txt

# Read in meta-analysis results
mediation_rosmap_mayo_meta.data <- fread("meta_output/mediation_rosmap_mayo_meta_pathological.txt")

# Number of total genes tested
dim(mediation_rosmap_mayo_meta.data)
# Number of genes overlapping
table(!is.na(mediation_rosmap_mayo_meta.data$PVALUE_FE))
# Bonferroni corrected p-value threshold
bonferroni <- 0.05/731

# Format discovery results
hits <- mediation_rosmap_mayo_meta.data[mediation_rosmap_mayo_meta.data$PVALUE_FE < bonferroni, RSID]
hits.data <- results_all_mediation.data[results_all_mediation.data$transcript %in% hits, ]
hits_formatted.data <- hits.data[,c("transcript","est_mediation_rosmap","se_mediation_rosmap","pvalue_mediation_rosmap","est_mediation_mayo","se_mediation_mayo","pvalue_mediation_mayo","est_mediation_msbb_ba10","se_mediation_msbb_ba10","pvalue_mediation_msbb_ba10","est_mediation_msbb_ba22","se_mediation_msbb_ba22","pvalue_mediation_msbb_ba22","est_mediation_msbb_ba36","se_mediation_msbb_ba36","pvalue_mediation_msbb_ba36","est_mediation_msbb_ba44","se_mediation_msbb_ba44","pvalue_mediation_msbb_ba44")]
write.csv(hits_formatted.data, file = "results/hits_sex_bias.csv", quote = F, row.names = F)
write.csv(mediation_rosmap_mayo_meta.data[mediation_rosmap_mayo_meta.data$PVALUE_FE < bonferroni, ], file = "results/hits_sex_bias_meta.csv", quote = F, row.names = F)

# gender_transcript
gender_transcript_rosmap_mayo_meta.data <- fread("meta_output/gender_transcript_rosmap_mayo_meta_pathological.txt")
hits_formatted_gender_transcript.data <- hits.data[,c("transcript","est_gender_transcript_rosmap","se_gender_transcript_rosmap","pvalue_gender_transcript_rosmap","est_gender_transcript_mayo","se_gender_transcript_mayo","pvalue_gender_transcript_mayo","est_gender_transcript_msbb_ba10","se_gender_transcript_msbb_ba10","pvalue_gender_transcript_msbb_ba10","est_gender_transcript_msbb_ba22","se_gender_transcript_msbb_ba22","pvalue_gender_transcript_msbb_ba22","est_gender_transcript_msbb_ba36","se_gender_transcript_msbb_ba36","pvalue_gender_transcript_msbb_ba36","est_gender_transcript_msbb_ba44","se_gender_transcript_msbb_ba44","pvalue_gender_transcript_msbb_ba44")]
write.csv(hits_formatted_gender_transcript.data, file = "results/hits_sex_bias_gender_transcript.csv", quote = F, row.names = F)
write.csv(gender_transcript_rosmap_mayo_meta.data[gender_transcript_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_gender_transcript.csv", quote = F, row.names = F)

# transcript_ad
transcript_ad_rosmap_mayo_meta.data <- fread("meta_output/transcript_ad_rosmap_mayo_meta_pathological.txt")
hits_formatted_transcript_ad.data <- hits.data[,c("transcript","est_transcript_ad_rosmap","se_transcript_ad_rosmap","pvalue_transcript_ad_rosmap","est_transcript_ad_mayo","se_transcript_ad_mayo","pvalue_transcript_ad_mayo","est_transcript_ad_msbb_ba10","se_transcript_ad_msbb_ba10","pvalue_transcript_ad_msbb_ba10","est_transcript_ad_msbb_ba22","se_transcript_ad_msbb_ba22","pvalue_transcript_ad_msbb_ba22","est_transcript_ad_msbb_ba36","se_transcript_ad_msbb_ba36","pvalue_transcript_ad_msbb_ba36","est_transcript_ad_msbb_ba44","se_transcript_ad_msbb_ba44","pvalue_transcript_ad_msbb_ba44")]
write.csv(hits_formatted_transcript_ad.data, file = "results/hits_sex_bias_transcript_ad.csv", quote = F, row.names = F)
write.csv(transcript_ad_rosmap_mayo_meta.data[transcript_ad_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_transcript_ad.csv", quote = F, row.names = F)

# direct
direct_rosmap_mayo_meta.data <- fread("meta_output/direct_rosmap_mayo_meta_pathological.txt")
hits_formatted_direct.data <- hits.data[,c("transcript","est_direct_rosmap","se_direct_rosmap","pvalue_direct_rosmap","est_direct_mayo","se_direct_mayo","pvalue_direct_mayo","est_direct_msbb_ba10","se_direct_msbb_ba10","pvalue_direct_msbb_ba10","est_direct_msbb_ba22","se_direct_msbb_ba22","pvalue_direct_msbb_ba22","est_direct_msbb_ba36","se_direct_msbb_ba36","pvalue_direct_msbb_ba36","est_direct_msbb_ba44","se_direct_msbb_ba44","pvalue_direct_msbb_ba44")]
write.csv(hits_formatted_direct.data, file = "results/hits_sex_bias_direct.csv", quote = F, row.names = F)
write.csv(direct_rosmap_mayo_meta.data[direct_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_direct.csv", quote = F, row.names = F)

# total
total_rosmap_mayo_meta.data <- fread("meta_output/total_rosmap_mayo_meta_pathological.txt")
hits_formatted_total.data <- hits.data[,c("transcript","est_total_rosmap","se_total_rosmap","pvalue_total_rosmap","est_total_mayo","se_total_mayo","pvalue_total_mayo","est_total_msbb_ba10","se_total_msbb_ba10","pvalue_total_msbb_ba10","est_total_msbb_ba22","se_total_msbb_ba22","pvalue_total_msbb_ba22","est_total_msbb_ba36","se_total_msbb_ba36","pvalue_total_msbb_ba36","est_total_msbb_ba44","se_total_msbb_ba44","pvalue_total_msbb_ba44")]
write.csv(hits_formatted_total.data, file = "results/hits_sex_bias_total.csv", quote = F, row.names = F)
write.csv(total_rosmap_mayo_meta.data[total_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_total.csv", quote = F, row.names = F)

### Look-up of XWAS suggestive genes
xwas <- c("ENSG00000146938.14","ENSG00000169933.12","ENSG00000198947.15","ENSG00000125675.17","ENSG00000196632.10","ENSG00000181191.11","ENSG00000126733.20","ENSG00000065923.9","ENSG00000169306.9","ENSG00000131171.12","ENSG00000122122.9","ENSG00000196664.4","ENSG00000173198.5")
xwas.data <- results_all_mediation.data[results_all_mediation.data$transcript %in% xwas, ]
xwas_formatted.data <- xwas.data[,c("transcript","est_mediation_rosmap","se_mediation_rosmap","pvalue_mediation_rosmap","est_mediation_mayo","se_mediation_mayo","pvalue_mediation_mayo","est_mediation_msbb_ba10","se_mediation_msbb_ba10","pvalue_mediation_msbb_ba10","est_mediation_msbb_ba22","se_mediation_msbb_ba22","pvalue_mediation_msbb_ba22","est_mediation_msbb_ba36","se_mediation_msbb_ba36","pvalue_mediation_msbb_ba36","est_mediation_msbb_ba44","se_mediation_msbb_ba44","pvalue_mediation_msbb_ba44")]
write.csv(xwas_formatted.data, file = "results/xwas_sex_bias.csv", quote = F, row.names = F)

# gender_transcript
xwas_formatted_gender_transcript.data <- xwas.data[,c("transcript","est_gender_transcript_rosmap","se_gender_transcript_rosmap","pvalue_gender_transcript_rosmap","est_gender_transcript_mayo","se_gender_transcript_mayo","pvalue_gender_transcript_mayo","est_gender_transcript_msbb_ba10","se_gender_transcript_msbb_ba10","pvalue_gender_transcript_msbb_ba10","est_gender_transcript_msbb_ba22","se_gender_transcript_msbb_ba22","pvalue_gender_transcript_msbb_ba22","est_gender_transcript_msbb_ba36","se_gender_transcript_msbb_ba36","pvalue_gender_transcript_msbb_ba36","est_gender_transcript_msbb_ba44","se_gender_transcript_msbb_ba44","pvalue_gender_transcript_msbb_ba44")]
write.csv(xwas_formatted_gender_transcript.data, file = "results/xwas_sex_bias_gender_transcript.csv", quote = F, row.names = F)

# transcript_ad
xwas_formatted_transcript_ad.data <- xwas.data[,c("transcript","est_transcript_ad_rosmap","se_transcript_ad_rosmap","pvalue_transcript_ad_rosmap","est_transcript_ad_mayo","se_transcript_ad_mayo","pvalue_transcript_ad_mayo","est_transcript_ad_msbb_ba10","se_transcript_ad_msbb_ba10","pvalue_transcript_ad_msbb_ba10","est_transcript_ad_msbb_ba22","se_transcript_ad_msbb_ba22","pvalue_transcript_ad_msbb_ba22","est_transcript_ad_msbb_ba36","se_transcript_ad_msbb_ba36","pvalue_transcript_ad_msbb_ba36","est_transcript_ad_msbb_ba44","se_transcript_ad_msbb_ba44","pvalue_transcript_ad_msbb_ba44")]
write.csv(xwas_formatted_transcript_ad.data, file = "results/xwas_sex_bias_transcript_ad.csv", quote = F, row.names = F)

# Descriptives
cohort.list <- list(pheno_transcriptome_rosmap.data,pheno_transcriptome_mayo.data,pheno_transcriptome_msbb_ba10.data,pheno_transcriptome_msbb_ba22.data,pheno_transcriptome_msbb_ba36.data,pheno_transcriptome_msbb_ba44.data)

descriptives_females <- function(data) {
  data <- data[data$Sex == "Female", ]
  name <- data$Cohort[1] 
  n_cases <- table(data$FK_final_CaseControl)["AD"]
  n_control <- table(data$FK_final_CaseControl)["Control"]
  age <- mean(data$age)
  eur_ancestry <- prop.table(table(data$FK_Ancestry))["EUR"]
  descriptives.data <- data.frame(cohort = name, n_cases = n_cases, n_control=n_control, age=age, ancestry=eur_ancestry)
}

descriptives_males <- function(data) {
  data <- data[data$Sex == "Male", ]
  name <- data$Cohort[1] 
  n_cases <- table(data$FK_final_CaseControl)["AD"]
  n_control <- table(data$FK_final_CaseControl)["Control"]
  age <- mean(data$age)
  eur_ancestry <- prop.table(table(data$FK_Ancestry))["EUR"]
  descriptives.data <- data.frame(cohort = name, n_cases = n_cases, n_control=n_control, age=age, ancestry=eur_ancestry)
}

descriptives_females.data <- do.call(rbind, lapply(cohort.list, descriptives_females))
descriptives_males.data <- do.call(rbind, lapply(cohort.list, descriptives_males))
write.csv(descriptives_females.data, file = "results/descriptives_females.csv", row.names = F)
write.csv(descriptives_males.data, file = "results/descriptives_males.csv", row.names = F)




