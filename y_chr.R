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
rosmap_transcriptome.data <- fread("transcriptome/ROSMAP/rosmap_dlpfc_chrY.txt")
mayo_transcriptome.data <- fread("transcriptome/Mayo/mayo_tcx_chrY.txt")
msbb_ba10_transcriptome.data <- fread("transcriptome/MSBB_BA10/msbb_ba10_chrY.txt")
msbb_ba22_transcriptome.data <- fread("transcriptome/MSBB_BA22/msbb_ba22_chrY.txt")
msbb_ba36_transcriptome.data <- fread("transcriptome/MSBB_BA36/msbb_ba36_chrY.txt")
msbb_ba44_transcriptome.data <- fread("transcriptome/MSBB_BA44/msbb_ba44_chrY.txt")

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

### Rosmap
pheno_transcriptome_rosmap.data <- merge(pheno_rosmap.data, rosmap_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_rosmap.data <- merge(pheno_transcriptome_rosmap.data, rosmap_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
rosmap_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_rosmap.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:60))]), pheno_transcriptome_rosmap.data$FK_final_CaseControl, family = "binomial")
plot(rosmap_lasso.fit)
coef(rosmap_lasso.fit, s = "lambda.1se")

y_chr_rosmap.model <- '
  FK_final_CaseControl ~ transcript_ad*transcript + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov4 + InferredCov8 + InferredCov35 + InferredCov36
'

rosmap_transcripts <- names(rosmap_transcriptome.data)[2:27]
y_chr_rosmap.list <- mclapply(pheno_transcriptome_rosmap.data[,..rosmap_transcripts], function(transcript){
  pheno_transcriptome_rosmap.data$transcript <- transcript
  y_chr_rosmap.fit <- sem(y_chr_rosmap.model, data = pheno_transcriptome_rosmap.data, estimator = "WLSMV")
  estimates <- parameterestimates(y_chr_rosmap.fit)
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  return(results_transcript_ad)
}, mc.cores = 8)
y_chr_rosmap.data <- do.call(rbind, y_chr_rosmap.list)

names(y_chr_rosmap.data) <- paste(names(y_chr_rosmap.data),"rosmap",sep="_")

### Mayo
pheno_transcriptome_mayo.data <- merge(pheno_mayo.data, mayo_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_mayo.data <- merge(pheno_transcriptome_mayo.data, mayo_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
mayo_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_mayo.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:45))]), pheno_transcriptome_mayo.data$FK_final_CaseControl, family = "binomial")
plot(mayo_lasso.fit)
coef(mayo_lasso.fit, s = "lambda.1se")

y_chr_mayo.model <- '
  FK_final_CaseControl ~ transcript_ad*transcript + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov2 + InferredCov3 + InferredCov4 + InferredCov5 + InferredCov7
'

mayo_transcripts <- names(mayo_transcriptome.data)[2:46]
y_chr_mayo.list <- mclapply(pheno_transcriptome_mayo.data[,..mayo_transcripts], function(transcript){
  pheno_transcriptome_mayo.data$transcript <- transcript
  y_chr_mayo.fit <- sem(y_chr_mayo.model, data = pheno_transcriptome_mayo.data, estimator = "WLSMV")
  estimates <- parameterestimates(y_chr_mayo.fit)
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  return(results_transcript_ad)
}, mc.cores = 8)
y_chr_mayo.data <- do.call(rbind, y_chr_mayo.list)

names(y_chr_mayo.data) <- paste(names(y_chr_mayo.data),"mayo",sep="_")

### MSBB BA10
pheno_transcriptome_msbb_ba10.data <- merge(pheno_msbb_ba10.data, msbb_ba10_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba10.data <- merge(pheno_transcriptome_msbb_ba10.data, msbb_ba10_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba10_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba10.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba10.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba10_lasso.fit)
coef(msbb_ba10_lasso.fit, s = "lambda.1se")

y_chr_msbb_ba10.model <- '
  FK_final_CaseControl ~ transcript_ad*transcript + age + gPC1 + gPC2 + gPC3 + gPC4
'

msbb_ba10_transcripts <- names(msbb_ba10_transcriptome.data)[2:37]
y_chr_msbb_ba10.list <- mclapply(pheno_transcriptome_msbb_ba10.data[,..msbb_ba10_transcripts], function(transcript){
  pheno_transcriptome_msbb_ba10.data$transcript <- transcript
  y_chr_msbb_ba10.fit <- sem(y_chr_msbb_ba10.model, data = pheno_transcriptome_msbb_ba10.data, estimator = "WLSMV")
  estimates <- parameterestimates(y_chr_msbb_ba10.fit)
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  return(results_transcript_ad)
}, mc.cores = 4)
y_chr_msbb_ba10.data <- do.call(rbind, y_chr_msbb_ba10.list)

names(y_chr_msbb_ba10.data) <- paste(names(y_chr_msbb_ba10.data),"msbb_ba10",sep="_")

### MSBB BA22
pheno_transcriptome_msbb_ba22.data <- merge(pheno_msbb_ba22.data, msbb_ba22_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba22.data <- merge(pheno_transcriptome_msbb_ba22.data, msbb_ba22_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba22_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba22.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba22.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba22_lasso.fit)
coef(msbb_ba22_lasso.fit, s = "lambda.1se")

y_chr_msbb_ba22.model <- '
  FK_final_CaseControl ~ transcript_ad*transcript + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov3 + InferredCov4
'

msbb_ba22_transcripts <- names(msbb_ba22_transcriptome.data)[2:37]
y_chr_msbb_ba22.list <- mclapply(pheno_transcriptome_msbb_ba22.data[,..msbb_ba22_transcripts], function(transcript){
  pheno_transcriptome_msbb_ba22.data$transcript <- transcript
  y_chr_msbb_ba22.fit <- sem(y_chr_msbb_ba22.model, data = pheno_transcriptome_msbb_ba22.data, estimator = "WLSMV")
  estimates <- parameterestimates(y_chr_msbb_ba22.fit)
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  return(results_transcript_ad)
}, mc.cores = 8)
y_chr_msbb_ba22.data <- do.call(rbind, y_chr_msbb_ba22.list)

names(y_chr_msbb_ba22.data) <- paste(names(y_chr_msbb_ba22.data),"msbb_ba22",sep="_")

### MSBB BA36
pheno_transcriptome_msbb_ba36.data <- merge(pheno_msbb_ba36.data, msbb_ba36_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba36.data <- merge(pheno_transcriptome_msbb_ba36.data, msbb_ba36_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba36_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba36.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba36.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba36_lasso.fit)
coef(msbb_ba36_lasso.fit, s = "lambda.1se")

y_chr_msbb_ba36.model <- '
  FK_final_CaseControl ~ transcript_ad*transcript + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov3
'

msbb_ba36_transcripts <- names(msbb_ba36_transcriptome.data)[2:38]
y_chr_msbb_ba36.list <- mclapply(pheno_transcriptome_msbb_ba36.data[,..msbb_ba36_transcripts], function(transcript){
  pheno_transcriptome_msbb_ba36.data$transcript <- transcript
  y_chr_msbb_ba36.fit <- sem(y_chr_msbb_ba36.model, data = pheno_transcriptome_msbb_ba36.data, estimator = "WLSMV")
  estimates <- parameterestimates(y_chr_msbb_ba36.fit)
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  return(results_transcript_ad)
}, mc.cores = 8)
y_chr_msbb_ba36.data <- do.call(rbind, y_chr_msbb_ba36.list)

names(y_chr_msbb_ba36.data) <- paste(names(y_chr_msbb_ba36.data),"msbb_ba36",sep="_")

### MSBB BA44
pheno_transcriptome_msbb_ba44.data <- merge(pheno_msbb_ba44.data, msbb_ba44_transcriptome.data, by.x = "specimenID", by.y = "IID")
pheno_transcriptome_msbb_ba44.data <- merge(pheno_transcriptome_msbb_ba44.data, msbb_ba44_peer.data, by = "specimenID")

# Lasso
set.seed(20250105)
msbb_ba44_lasso.fit <- cv.glmnet(as.matrix(pheno_transcriptome_msbb_ba44.data[,c(c("gender","age","gPC1","gPC2","gPC3","gPC4"),paste0("InferredCov",1:30))]), pheno_transcriptome_msbb_ba44.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba44_lasso.fit)
coef(msbb_ba44_lasso.fit, s = "lambda.1se")

y_chr_msbb_ba44.model <- '
  FK_final_CaseControl ~ transcript_ad*transcript + age + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1
'

msbb_ba44_transcripts <- names(msbb_ba44_transcriptome.data)[2:35]
y_chr_msbb_ba44.list <- mclapply(pheno_transcriptome_msbb_ba44.data[,..msbb_ba44_transcripts], function(transcript){
  pheno_transcriptome_msbb_ba44.data$transcript <- transcript
  y_chr_msbb_ba44.fit <- sem(y_chr_msbb_ba44.model, data = pheno_transcriptome_msbb_ba44.data, estimator = "WLSMV")
  estimates <- parameterestimates(y_chr_msbb_ba44.fit)
  results_transcript_ad <- estimates[estimates$label == "transcript_ad",c("est","se","pvalue")]
  names(results_transcript_ad) <- c("est_transcript_ad","se_transcript_ad","pvalue_transcript_ad")
  return(results_transcript_ad)
}, mc.cores = 4)
y_chr_msbb_ba44.data <- do.call(rbind, y_chr_msbb_ba44.list)

names(y_chr_msbb_ba44.data) <- paste(names(y_chr_msbb_ba44.data),"msbb_ba44",sep="_")

# Make column for row.names
# Mediation estimates
y_chr.list <- list(y_chr_rosmap.data,y_chr_mayo.data,y_chr_msbb_ba10.data,y_chr_msbb_ba22.data,y_chr_msbb_ba36.data,y_chr_msbb_ba44.data)
y_chr.list <- lapply(y_chr.list, function(x) data.frame(x, transcript = row.names(x)))

results_all_y_chr.data <- Reduce(function(x, y) merge(x, y, by="transcript", all=TRUE), y_chr.list)
write.csv(results_all_y_chr.data[results_all_y_chr.data$transcript %in% c("ENSG00000012817.15","ENSG00000183878.15"),], file = "results/y_chr.csv", quote = F, row.names = F)

# Meta
fwrite(results_all_y_chr.data[,c("transcript","est_transcript_ad_rosmap","se_transcript_ad_rosmap","est_transcript_ad_mayo","se_transcript_ad_mayo")], file = "meta_input/transcript_ad_rosmap_mayo_y_chr.txt", sep = "\t", na = "NA", col.names = F)

java -jar /home/aneumann/Metasoft/Metasoft.jar -pvalue_table /home/aneumann/Metasoft/HanEskinPvalueTable.txt -input meta_input/transcript_ad_rosmap_mayo_y_chr.txt -output meta_output/transcript_ad_rosmap_mayo_meta_y_chr.txt

