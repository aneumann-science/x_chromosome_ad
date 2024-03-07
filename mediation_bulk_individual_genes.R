library(data.table)
library(parallel)
library(lavaan)
library(glmnet)
library(psych)

###### Phenotype prep
wgs_info.data <- fread("covariates/AMPAD_and_EADBBelgian_FullCovariatesAndSampleInfo_Correlation_19May2022.csv")

# Center gender and code female as 1, male as 0
wgs_info.data$gender <- as.numeric(as.factor(wgs_info.data$Sex)) - 1
wgs_info.data$gender <- 1-wgs_info.data$gender

# Specify continous age variable 
wgs_info.data$age <- wgs_info.data$AgeAtDeath_or_Sampling
wgs_info.data$age[wgs_info.data$age == "90_or_above"] <- 90
wgs_info.data$age <- as.numeric(wgs_info.data$age)

# Define ADvsCon as ordered for lavaan
wgs_info_pathological.data <- wgs_info.data[wgs_info.data$FK_final_CaseControl != "Other", ]
wgs_info_pathological.data$FK_final_CaseControl <- ordered(wgs_info_pathological.data$FK_final_CaseControl, levels = c("Control","AD"))

# Stratify by cohort
wgs_info_rosmap.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "ROSMAP_DLPFC", ]
wgs_info_mayo.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "MayoRNASeq_TCX", ]
wgs_info_msbb_ba10.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "MSBB_BA10", ]
wgs_info_msbb_ba22.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "MSBB_BA22", ]
wgs_info_msbb_ba36.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "MSBB_BA36", ]
wgs_info_msbb_ba44.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "MSBB_BA44", ]

### Load transcriptoms
rosmap_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/ROSMAP/rosmap_dlpfc_chrX.txt")
mayo_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/Mayo/mayo_tcx_chrX.txt")
msbb_ba10_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/MSBB_BA10/msbb_ba10_chrX.txt")
msbb_ba22_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/MSBB_BA22/msbb_ba22_chrX.txt")
msbb_ba36_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/MSBB_BA36/msbb_ba36_chrX.txt")
msbb_ba44_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/MSBB_BA44/msbb_ba44_chrX.txt")

# Count number on unique transcripts for p-value calculations
length(unique(c(names(rosmap_transcriptome.data),names(mayo_transcriptome.data),names(msbb_ba10_transcriptome.data),names(msbb_ba22_transcriptome.data),names(msbb_ba36_transcriptome.data),names(msbb_ba44_transcriptome.data))))

# 1076 not counting
p_bonferroni <- 0.05/994

### Rosmap
wgs_info_transcriptome_rosmap.data <- merge(wgs_info_rosmap.data, rosmap_transcriptome.data, by.x = "ID", by.y = "IID")

# Lasso
set.seed(20200516)
rosmap_lasso.fit <- cv.glmnet(as.matrix(wgs_info_transcriptome_rosmap.data[,c(c("gender","age"),paste0("gPC",1:4),paste0("ePEER",1:60))]), wgs_info_transcriptome_rosmap.data$FK_final_CaseControl, family = "binomial")
plot(rosmap_lasso.fit)
coef(rosmap_lasso.fit, s = "lambda.1se")
# wgs_info_transcriptome_rosmap.data$peer_predicted <- predict(rosmap_lasso.fit, as.matrix(wgs_info_transcriptome_rosmap.data[,paste0("ePEER",1:60)]))
# peer.fit <- glm(FK_final_CaseControl ~ peer_predicted, data = wgs_info_transcriptome_rosmap.data, family = "binomial")
# summary(peer.fit)
# exp(coef(peer.fit))

mediation_rosmap.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER4 + ePEER6 + ePEER7 + ePEER10 + ePEER18 + ePEER24 + ePEER40 + ePEER41 + ePEER42 + ePEER52 + ePEER60
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER4 + ePEER6 + ePEER7 + ePEER10 + ePEER18 + ePEER24 + ePEER40 + ePEER41 + ePEER42 + ePEER52 + ePEER60
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

rosmap_transcripts <- names(rosmap_transcriptome.data)[2:749]
mediation_rosmap.list <- mclapply(wgs_info_transcriptome_rosmap.data[,..rosmap_transcripts], function(transcript){
  wgs_info_transcriptome_rosmap.data$transcript <- transcript
  mediation_rosmap.fit <- sem(mediation_rosmap.model, data = wgs_info_transcriptome_rosmap.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_rosmap.fit)
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
}, mc.cores = 4)
mediation_rosmap.data <- do.call(rbind, mediation_rosmap.list)

names(mediation_rosmap.data) <- paste(names(mediation_rosmap.data),"rosmap",sep="_")

# # Average transcription
# wgs_info_transcriptome_rosmap.data$transcript <- rowMeans(wgs_info_transcriptome_rosmap.data[,..rosmap_transcripts])
# summary(lm(transcript~ gender + age + peer_predicted, data = wgs_info_transcriptome_rosmap.data))
# mediation_rosmap_avg.fit <- sem(mediation.model, data = wgs_info_transcriptome_rosmap.data, estimator = "WLSMV")
# summary(mediation_rosmap_avg.fit)

### Mayo
wgs_info_mayo.data$ID <- as.numeric(wgs_info_mayo.data$ID) 
wgs_info_transcriptome_mayo.data <- merge(wgs_info_mayo.data, mayo_transcriptome.data, by.x = "ID", by.y = "IID")

# Lasso
set.seed(20200516)
mayo_lasso.fit <- cv.glmnet(as.matrix(wgs_info_transcriptome_mayo.data[,c(c("gender","age"),paste0("gPC",1:4),paste0("ePEER",1:45))]), wgs_info_transcriptome_mayo.data$FK_final_CaseControl, family = "binomial")
plot(mayo_lasso.fit)
coef(mayo_lasso.fit, s = "lambda.1se")
# wgs_info_transcriptome_mayo.data$peer_predicted <- predict(mayo_lasso.fit, as.matrix(wgs_info_transcriptome_mayo.data[,paste0("ePEER",1:45)]))
# peer.fit <- glm(FK_final_CaseControl ~ peer_predicted, data = wgs_info_transcriptome_mayo.data, family = "binomial")
# plot(wgs_info_transcriptome_mayo.data$peer_predicted, as.numeric(wgs_info_transcriptome_mayo.data$FK_final_CaseControl))
# summary(peer.fit)
# exp(coef(peer.fit))

mediation_mayo.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER2 + ePEER3 + ePEER4 + ePEER7 + ePEER9 + ePEER11 + ePEER13 + ePEER14 + ePEER17 + ePEER21 + ePEER28
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER2 + ePEER3 + ePEER4 + ePEER7 + ePEER9 + ePEER11 + ePEER13 + ePEER14 + ePEER17 + ePEER21 + ePEER28
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

mayo_transcripts <- names(mayo_transcriptome.data)[2:900]
mediation_mayo.list <- mclapply(wgs_info_transcriptome_mayo.data[,..mayo_transcripts], function(transcript){
  wgs_info_transcriptome_mayo.data$transcript <- transcript
  mediation_mayo.fit <- sem(mediation_mayo.model, data = wgs_info_transcriptome_mayo.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_mayo.fit)
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
}, mc.cores = 4)
mediation_mayo.data <- do.call(rbind, mediation_mayo.list)

names(mediation_mayo.data) <- paste(names(mediation_mayo.data),"mayo",sep="_")

# # Average transcription
# wgs_info_transcriptome_mayo.data$transcript <- rowMeans(wgs_info_transcriptome_mayo.data[,..mayo_transcripts])
# summary(lm(transcript~ gender + age + ePEER1 + ePEER2 + ePEER3 + ePEER4 + ePEER5 + ePEER6 + ePEER7 + ePEER8 + ePEER9 + ePEER10 + ePEER11 + ePEER12 + ePEER13 + ePEER14 + ePEER15, data = wgs_info_transcriptome_mayo.data))
# mediation_mayo_avg.fit <- sem(mediation.model, data = wgs_info_transcriptome_mayo.data, estimator = "WLSMV")
# summary(mediation_mayo_avg.fit)

### MSBB BA10
wgs_info_msbb_ba10.data$ID <- as.numeric(wgs_info_msbb_ba10.data$ID)
wgs_info_transcriptome_msbb_ba10.data <- merge(wgs_info_msbb_ba10.data, msbb_ba10_transcriptome.data, by.x = "ID", by.y = "IID")

# Lasso
set.seed(20200516)
msbb_ba10_lasso.fit <- cv.glmnet(as.matrix(wgs_info_transcriptome_msbb_ba10.data[,c(c("gender","age"),paste0("gPC",1:4),paste0("ePEER",1:30))]), wgs_info_transcriptome_msbb_ba10.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba10_lasso.fit)
coef(msbb_ba10_lasso.fit, s = "lambda.1se")
# wgs_info_transcriptome_msbb_ba10.data$peer_predicted <- predict(msbb_ba10_lasso.fit, as.matrix(wgs_info_transcriptome_msbb_ba10.data[,paste0("ePEER",1:30)]))
# peer.fit <- glm(FK_final_CaseControl ~ peer_predicted, data = wgs_info_transcriptome_msbb_ba10.data, family = "binomial")
# summary(peer.fit)
# exp(coef(peer.fit))

mediation_msbb_ba10.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER4 + ePEER5 + ePEER6 + ePEER7 + ePEER8 + ePEER11 + ePEER12 + ePEER13 + ePEER19 + ePEER20 + ePEER21 + ePEER23 + ePEER27
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER4 + ePEER5 + ePEER6 + ePEER7 + ePEER8 + ePEER11 + ePEER12 + ePEER13 + ePEER19 + ePEER20 + ePEER21 + ePEER23 + ePEER27
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba10_transcripts <- names(msbb_ba10_transcriptome.data)[2:796]
mediation_msbb_ba10.list <- mclapply(wgs_info_transcriptome_msbb_ba10.data[,..msbb_ba10_transcripts], function(transcript){
  wgs_info_transcriptome_msbb_ba10.data$transcript <- transcript
  mediation_msbb_ba10.fit <- sem(mediation_msbb_ba10.model, data = wgs_info_transcriptome_msbb_ba10.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba10.fit)
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
}, mc.cores = 4)
mediation_msbb_ba10.data <- do.call(rbind, mediation_msbb_ba10.list)

names(mediation_msbb_ba10.data) <- paste(names(mediation_msbb_ba10.data),"msbb_ba10",sep="_")

# # Average transcription
# wgs_info_transcriptome_msbb_ba10.data$transcript <- rowMeans(wgs_info_transcriptome_msbb_ba10.data[,..msbb_ba10_transcripts])
# summary(lm(transcript~ gender + age + ePEER1 + ePEER2 + ePEER3 + ePEER4 + ePEER5 + ePEER6 + ePEER7 + ePEER8 + ePEER9 + ePEER10 + ePEER11 + ePEER12 + ePEER13 + ePEER14 + ePEER15, data = wgs_info_transcriptome_msbb_ba10.data))
# mediation_msbb_ba10_avg.fit <- sem(mediation.model, data = wgs_info_transcriptome_msbb_ba10.data, estimator = "WLSMV")
# summary(mediation_msbb_ba10_avg.fit)

### MSBB BA22
wgs_info_msbb_ba22.data$ID <- as.numeric(wgs_info_msbb_ba22.data$ID)
wgs_info_transcriptome_msbb_ba22.data <- merge(wgs_info_msbb_ba22.data, msbb_ba22_transcriptome.data, by.x = "ID", by.y = "IID")

# Lasso
set.seed(20200516)
msbb_ba22_lasso.fit <- cv.glmnet(as.matrix(wgs_info_transcriptome_msbb_ba22.data[,c(c("gender","age"),paste0("gPC",1:4),paste0("ePEER",1:30))]), wgs_info_transcriptome_msbb_ba22.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba22_lasso.fit)
coef(msbb_ba22_lasso.fit, s = "lambda.1se")
# wgs_info_transcriptome_msbb_ba22.data$peer_predicted <- predict(msbb_ba22_lasso.fit, as.matrix(wgs_info_transcriptome_msbb_ba22.data[,paste0("ePEER",1:30)]))
# peer.fit <- glm(FK_final_CaseControl ~ peer_predicted, data = wgs_info_transcriptome_msbb_ba22.data, family = "binomial")
# summary(peer.fit)
# exp(coef(peer.fit))

mediation_msbb_ba22.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER3 + ePEER4 + ePEER15 + ePEER16 + ePEER19 + ePEER22 + ePEER23 + ePEER24 + ePEER27 + ePEER29 + ePEER30
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER3 + ePEER4 + ePEER15 + ePEER16 + ePEER19 + ePEER22 + ePEER23 + ePEER24 + ePEER27 + ePEER29 + ePEER30
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba22_transcripts <- names(msbb_ba22_transcriptome.data)[2:761]
mediation_msbb_ba22.list <- mclapply(wgs_info_transcriptome_msbb_ba22.data[,..msbb_ba22_transcripts], function(transcript){
  wgs_info_transcriptome_msbb_ba22.data$transcript <- transcript
  mediation_msbb_ba22.fit <- sem(mediation_msbb_ba22.model, data = wgs_info_transcriptome_msbb_ba22.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba22.fit)
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
}, mc.cores = 4)
mediation_msbb_ba22.data <- do.call(rbind, mediation_msbb_ba22.list)

names(mediation_msbb_ba22.data) <- paste(names(mediation_msbb_ba22.data),"msbb_ba22",sep="_")

# # Average transcription
# wgs_info_transcriptome_msbb_ba22.data$transcript <- rowMeans(wgs_info_transcriptome_msbb_ba22.data[,..msbb_ba22_transcripts])
# summary(lm(transcript~ gender + age + ePEER1 + ePEER2 + ePEER3 + ePEER4 + ePEER5 + ePEER6 + ePEER7 + ePEER8 + ePEER9 + ePEER10 + ePEER11 + ePEER12 + ePEER13 + ePEER14 + ePEER15, data = wgs_info_transcriptome_msbb_ba22.data))
# mediation_msbb_ba22_avg.fit <- sem(mediation.model, data = wgs_info_transcriptome_msbb_ba22.data, estimator = "WLSMV")
# summary(mediation_msbb_ba22_avg.fit)

### MSBB BA36
wgs_info_msbb_ba36.data$ID <- as.numeric(wgs_info_msbb_ba36.data$ID)
wgs_info_transcriptome_msbb_ba36.data <- merge(wgs_info_msbb_ba36.data, msbb_ba36_transcriptome.data, by.x = "ID", by.y = "IID")

# Lasso
set.seed(20200516)
msbb_ba36_lasso.fit <- cv.glmnet(as.matrix(wgs_info_transcriptome_msbb_ba36.data[,c(c("gender","age"),paste0("gPC",1:4),paste0("ePEER",1:30))]), wgs_info_transcriptome_msbb_ba36.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba36_lasso.fit)
coef(msbb_ba36_lasso.fit, s = "lambda.1se")
# wgs_info_transcriptome_msbb_ba36.data$peer_predicted <- predict(msbb_ba36_lasso.fit, as.matrix(wgs_info_transcriptome_msbb_ba36.data[,paste0("ePEER",1:30)]))
# peer.fit <- glm(FK_final_CaseControl ~ peer_predicted, data = wgs_info_transcriptome_msbb_ba36.data, family = "binomial")
# summary(peer.fit)
# exp(coef(peer.fit))

mediation_msbb_ba36.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER3 + ePEER5 + ePEER9 + ePEER10 + ePEER14 + ePEER19 + ePEER21 + ePEER23 + ePEER26 + ePEER30
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER3 + ePEER5 + ePEER9 + ePEER10 + ePEER14 + ePEER19 + ePEER21 + ePEER23 + ePEER26 + ePEER30
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba36_transcripts <- names(msbb_ba36_transcriptome.data)[2:783]
mediation_msbb_ba36.list <- mclapply(wgs_info_transcriptome_msbb_ba36.data[,..msbb_ba36_transcripts], function(transcript){ 
  wgs_info_transcriptome_msbb_ba36.data$transcript <- transcript
  mediation_msbb_ba36.fit <- sem(mediation_msbb_ba36.model, data = wgs_info_transcriptome_msbb_ba36.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba36.fit)
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
}, mc.cores = 4)
mediation_msbb_ba36.data <- do.call(rbind, mediation_msbb_ba36.list)

names(mediation_msbb_ba36.data) <- paste(names(mediation_msbb_ba36.data),"msbb_ba36",sep="_")

# # Average transcription
# wgs_info_transcriptome_msbb_ba36.data$transcript <- rowMeans(wgs_info_transcriptome_msbb_ba36.data[,..msbb_ba36_transcripts])
# summary(lm(transcript~ gender + age + ePEER1 + ePEER2 + ePEER3 + ePEER4 + ePEER5 + ePEER6 + ePEER7 + ePEER8 + ePEER9 + ePEER10 + ePEER11 + ePEER12 + ePEER13 + ePEER14 + ePEER15, data = wgs_info_transcriptome_msbb_ba36.data))
# mediation_msbb_ba36_avg.fit <- sem(mediation.model, data = wgs_info_transcriptome_msbb_ba36.data, estimator = "WLSMV")
# summary(mediation_msbb_ba36_avg.fit)

### MSBB BA44
wgs_info_msbb_ba44.data$ID <- as.numeric(wgs_info_msbb_ba44.data$ID)
wgs_info_transcriptome_msbb_ba44.data <- merge(wgs_info_msbb_ba44.data, msbb_ba44_transcriptome.data, by.x = "ID", by.y = "IID")

# Lasso
set.seed(20200516)
msbb_ba44_lasso.fit <- cv.glmnet(as.matrix(wgs_info_transcriptome_msbb_ba44.data[,c(c("gender","age"),paste0("gPC",1:4),paste0("ePEER",1:30))]), wgs_info_transcriptome_msbb_ba44.data$FK_final_CaseControl, family = "binomial")
plot(msbb_ba44_lasso.fit)
coef(msbb_ba44_lasso.fit, s = "lambda.1se")
# wgs_info_transcriptome_msbb_ba44.data$peer_predicted <- predict(msbb_ba44_lasso.fit, as.matrix(wgs_info_transcriptome_msbb_ba44.data[,paste0("ePEER",1:30)]))
# peer.fit <- glm(FK_final_CaseControl ~ peer_predicted, data = wgs_info_transcriptome_msbb_ba44.data, family = "binomial")
# summary(peer.fit)
# exp(coef(peer.fit))

mediation_msbb_ba44.model <- '
  transcript ~ gender_transcript*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER2 + ePEER4 + ePEER6 + ePEER7 + ePEER11 + ePEER12 + ePEER13 + ePEER17 + ePEER25
  FK_final_CaseControl ~ transcript_ad*transcript + gender_ad*gender + age + gPC1 + gPC2 + gPC3 + gPC4 + ePEER1 + ePEER2 + ePEER4 + ePEER6 + ePEER7 + ePEER11 + ePEER12 + ePEER13 + ePEER17 + ePEER25
  
  mediation := gender_transcript*transcript_ad
  total := mediation+gender_ad
'

msbb_ba44_transcripts <- names(msbb_ba44_transcriptome.data)[2:754]
mediation_msbb_ba44.list <- mclapply(wgs_info_transcriptome_msbb_ba44.data[,..msbb_ba44_transcripts], function(transcript){
  wgs_info_transcriptome_msbb_ba44.data$transcript <- transcript
  mediation_msbb_ba44.fit <- sem(mediation_msbb_ba44.model, data = wgs_info_transcriptome_msbb_ba44.data, estimator = "WLSMV")
  estimates <- parameterestimates(mediation_msbb_ba44.fit)
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
}, mc.cores = 4)
mediation_msbb_ba44.data <- do.call(rbind, mediation_msbb_ba44.list)

names(mediation_msbb_ba44.data) <- paste(names(mediation_msbb_ba44.data),"msbb_ba44",sep="_")

# # Average transcription
# wgs_info_transcriptome_msbb_ba44.data$transcript <- rowMeans(wgs_info_transcriptome_msbb_ba44.data[,..msbb_ba44_transcripts])
# summary(lm(transcript~ gender + age + ePEER1 + ePEER2 + ePEER3 + ePEER4 + ePEER5 + ePEER6 + ePEER7 + ePEER8 + ePEER9 + ePEER10 + ePEER11 + ePEER12 + ePEER13 + ePEER14 + ePEER15, data = wgs_info_transcriptome_msbb_ba44.data))
# mediation_msbb_ba44_avg.fit <- sem(mediation.model, data = wgs_info_transcriptome_msbb_ba44.data, estimator = "WLSMV")
# summary(mediation_msbb_ba44_avg.fit)

# Make column for row.names
# Mediation estimates
mediation.list <- list(mediation_rosmap.data,mediation_mayo.data,mediation_msbb_ba10.data,mediation_msbb_ba22.data,mediation_msbb_ba36.data,mediation_msbb_ba44.data)
mediation.list <- lapply(mediation.list, function(x) data.frame(x, transcript = row.names(x)))

results_all_mediation.data <- Reduce(function(x, y) merge(x, y, by="transcript", all=TRUE), mediation.list)

# Write file for meta analysis
fwrite(results_all_mediation.data[,c("transcript","est_mediation_rosmap","se_mediation_rosmap","est_mediation_mayo","se_mediation_mayo")], file = "meta/mediation_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_gender_transcript_rosmap","se_gender_transcript_rosmap","est_gender_transcript_mayo","se_gender_transcript_mayo")], file = "meta/gender_transcript_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_transcript_ad_rosmap","se_transcript_ad_rosmap","est_transcript_ad_mayo","se_transcript_ad_mayo")], file = "meta/transcript_ad_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_direct_rosmap","se_direct_rosmap","est_direct_mayo","se_direct_mayo")], file = "meta/direct_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)
fwrite(results_all_mediation.data[,c("transcript","est_total_rosmap","se_total_rosmap","est_total_mayo","se_total_mayo")], file = "meta/total_rosmap_mayo_pathological.txt", sep = "\t", na = "NA", col.names = F)

# 
java -jar /home/alexander/Metasoft/Metasoft.jar -pvalue_table /home/alexander/Metasoft/HanEskinPvalueTable.txt -input meta/mediation_rosmap_mayo_pathological.txt -output meta/mediation_rosmap_mayo_meta_pathological.txt
java -jar /home/alexander/Metasoft/Metasoft.jar -pvalue_table /home/alexander/Metasoft/HanEskinPvalueTable.txt -input meta/gender_transcript_rosmap_mayo_pathological.txt -output meta/gender_transcript_rosmap_mayo_meta_pathological.txt
java -jar /home/alexander/Metasoft/Metasoft.jar -pvalue_table /home/alexander/Metasoft/HanEskinPvalueTable.txt -input meta/transcript_ad_rosmap_mayo_pathological.txt -output meta/transcript_ad_rosmap_mayo_meta_pathological.txt
java -jar /home/alexander/Metasoft/Metasoft.jar -pvalue_table /home/alexander/Metasoft/HanEskinPvalueTable.txt -input meta/direct_rosmap_mayo_pathological.txt -output meta/direct_rosmap_mayo_meta_pathological.txt
java -jar /home/alexander/Metasoft/Metasoft.jar -pvalue_table /home/alexander/Metasoft/HanEskinPvalueTable.txt -input meta/total_rosmap_mayo_pathological.txt -output meta/total_rosmap_mayo_meta_pathological.txt

mediation_rosmap_mayo_meta.data <- fread("meta/mediation_rosmap_mayo_meta_pathological.txt")

dim(mediation_rosmap_mayo_meta.data)
table(!is.na(mediation_rosmap_mayo_meta.data$PVALUE_FE))
bonferroni <- 0.05/730

hits <- mediation_rosmap_mayo_meta.data[mediation_rosmap_mayo_meta.data$PVALUE_FE < bonferroni, RSID]
hits.data <- results_all_mediation.data[results_all_mediation.data$transcript %in% hits, ]
hits_formatted.data <- hits.data[,c("transcript","est_mediation_rosmap","se_mediation_rosmap","pvalue_mediation_rosmap","est_mediation_mayo","se_mediation_mayo","pvalue_mediation_mayo","est_mediation_msbb_ba10","se_mediation_msbb_ba10","pvalue_mediation_msbb_ba10","est_mediation_msbb_ba22","se_mediation_msbb_ba22","pvalue_mediation_msbb_ba22","est_mediation_msbb_ba36","se_mediation_msbb_ba36","pvalue_mediation_msbb_ba36","est_mediation_msbb_ba44","se_mediation_msbb_ba44","pvalue_mediation_msbb_ba44")]
write.csv(hits_formatted.data, file = "results/hits_sex_bias.csv", quote = F, row.names = F)
write.csv(mediation_rosmap_mayo_meta.data[mediation_rosmap_mayo_meta.data$PVALUE_FE < bonferroni, ], file = "results/hits_sex_bias_meta.csv", quote = F, row.names = F)

# gender_transcript
gender_transcript_rosmap_mayo_meta.data <- fread("meta/gender_transcript_rosmap_mayo_meta_pathological.txt")
hits_formatted_gender_transcript.data <- hits.data[,c("transcript","est_gender_transcript_rosmap","se_gender_transcript_rosmap","pvalue_gender_transcript_rosmap","est_gender_transcript_mayo","se_gender_transcript_mayo","pvalue_gender_transcript_mayo","est_gender_transcript_msbb_ba10","se_gender_transcript_msbb_ba10","pvalue_gender_transcript_msbb_ba10","est_gender_transcript_msbb_ba22","se_gender_transcript_msbb_ba22","pvalue_gender_transcript_msbb_ba22","est_gender_transcript_msbb_ba36","se_gender_transcript_msbb_ba36","pvalue_gender_transcript_msbb_ba36","est_gender_transcript_msbb_ba44","se_gender_transcript_msbb_ba44","pvalue_gender_transcript_msbb_ba44")]
write.csv(hits_formatted_gender_transcript.data, file = "results/hits_sex_bias_gender_transcript.csv", quote = F, row.names = F)
write.csv(gender_transcript_rosmap_mayo_meta.data[gender_transcript_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_gender_transcript.csv", quote = F, row.names = F)

# transcript_ad
transcript_ad_rosmap_mayo_meta.data <- fread("meta/transcript_ad_rosmap_mayo_meta_pathological.txt")
hits_formatted_transcript_ad.data <- hits.data[,c("transcript","est_transcript_ad_rosmap","se_transcript_ad_rosmap","pvalue_transcript_ad_rosmap","est_transcript_ad_mayo","se_transcript_ad_mayo","pvalue_transcript_ad_mayo","est_transcript_ad_msbb_ba10","se_transcript_ad_msbb_ba10","pvalue_transcript_ad_msbb_ba10","est_transcript_ad_msbb_ba22","se_transcript_ad_msbb_ba22","pvalue_transcript_ad_msbb_ba22","est_transcript_ad_msbb_ba36","se_transcript_ad_msbb_ba36","pvalue_transcript_ad_msbb_ba36","est_transcript_ad_msbb_ba44","se_transcript_ad_msbb_ba44","pvalue_transcript_ad_msbb_ba44")]
write.csv(hits_formatted_transcript_ad.data, file = "results/hits_sex_bias_transcript_ad.csv", quote = F, row.names = F)
write.csv(transcript_ad_rosmap_mayo_meta.data[transcript_ad_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_transcript_ad.csv", quote = F, row.names = F)

# direct
direct_rosmap_mayo_meta.data <- fread("meta/direct_rosmap_mayo_meta_pathological.txt")
hits_formatted_direct.data <- hits.data[,c("transcript","est_direct_rosmap","se_direct_rosmap","pvalue_direct_rosmap","est_direct_mayo","se_direct_mayo","pvalue_direct_mayo","est_direct_msbb_ba10","se_direct_msbb_ba10","pvalue_direct_msbb_ba10","est_direct_msbb_ba22","se_direct_msbb_ba22","pvalue_direct_msbb_ba22","est_direct_msbb_ba36","se_direct_msbb_ba36","pvalue_direct_msbb_ba36","est_direct_msbb_ba44","se_direct_msbb_ba44","pvalue_direct_msbb_ba44")]
write.csv(hits_formatted_direct.data, file = "results/hits_sex_bias_direct.csv", quote = F, row.names = F)
write.csv(direct_rosmap_mayo_meta.data[direct_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_direct.csv", quote = F, row.names = F)

# total
total_rosmap_mayo_meta.data <- fread("meta/total_rosmap_mayo_meta_pathological.txt")
hits_formatted_total.data <- hits.data[,c("transcript","est_total_rosmap","se_total_rosmap","pvalue_total_rosmap","est_total_mayo","se_total_mayo","pvalue_total_mayo","est_total_msbb_ba10","se_total_msbb_ba10","pvalue_total_msbb_ba10","est_total_msbb_ba22","se_total_msbb_ba22","pvalue_total_msbb_ba22","est_total_msbb_ba36","se_total_msbb_ba36","pvalue_total_msbb_ba36","est_total_msbb_ba44","se_total_msbb_ba44","pvalue_total_msbb_ba44")]
write.csv(hits_formatted_total.data, file = "results/hits_sex_bias_total.csv", quote = F, row.names = F)
write.csv(total_rosmap_mayo_meta.data[total_rosmap_mayo_meta.data$RSID %in% hits, ], file = "results/hits_sex_bias_meta_total.csv", quote = F, row.names = F)

# Descriptives
cohort.list <- list(wgs_info_transcriptome_rosmap.data,wgs_info_transcriptome_mayo.data,wgs_info_transcriptome_msbb_ba10.data,wgs_info_transcriptome_msbb_ba22.data,wgs_info_transcriptome_msbb_ba36.data,wgs_info_transcriptome_msbb_ba44.data)

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

# Total gender effects
# Gender effect
cohort_ad_gender.list <- list(wgs_info_pathological.data, wgs_info_transcriptome_rosmap.data,wgs_info_transcriptome_mayo.data,wgs_info_transcriptome_msbb_ba10.data,wgs_info_transcriptome_msbb_ba22.data,wgs_info_transcriptome_msbb_ba36.data,wgs_info_transcriptome_msbb_ba44.data)

ad_gender_model <- "FK_final_CaseControl ~ gender_ad*gender + age"

fit_ad_gender <- function(data) {
  ad_gender.fit <- sem(ad_gender_model, data = data, estimator = "WLSMV")
  parameterestimates(ad_gender.fit)[1,]
}

apply(cohort_ad_gender.list, fit_ad_gender)
ad_gender.data <- do.call(rbind, lapply(cohort_ad_gender.list, fit_ad_gender))

