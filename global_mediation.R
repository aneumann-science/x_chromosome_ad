library(hdmed)
library(data.table)
library(dplyr)
library(lavaan)
library(parallel)
library(glmnet)
library(boot)

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

# Create dummy variables for Mayo and MSBB BA10
pheno_pathological.data$mayo <- 0
pheno_pathological.data$mayo[pheno_pathological.data$Cohort == "MayoRNASeq_TCX"] <- 1
pheno_pathological.data$msbb <- 0
pheno_pathological.data$msbb[pheno_pathological.data$Cohort == "MSBB_BA10"] <- 1

# Reduce to ROSMAP, Mayo and MSBB BA10
pheno_pathological.data <- pheno_pathological.data[pheno_pathological.data$Cohort == "ROSMAP_DLPFC" | 
                                                     pheno_pathological.data$Cohort == "MayoRNASeq_TCX" |
                                                     pheno_pathological.data$Cohort == "MSBB_BA10"]

# Determine optimal PEER factor combination
peer.data <- fread("expression_matrices_and_PEER/ePEER_60_rosmapDLPFC_mayoTCX_msbbBA10_innerJoint_normalizedSeparately_allSamples_14April2025.PEER_covariates.txt", data.table = F)
names(peer.data) <- gsub("rosmapDLPFC_","",names(peer.data))
names(peer.data) <- gsub("mayoTCX_","",names(peer.data))
names(peer.data) <- gsub("msbbBA10_","",names(peer.data))
rownames(peer.data) <- peer.data$ID
peer.data <- peer.data[,2:ncol(peer.data)]
peer.data <- as.data.frame(t(peer.data))
peer.data$specimenID <- rownames(peer.data)
pheno_pathological_peer.data <- merge(pheno_pathological.data, peer.data, by = "specimenID")

set.seed(20250501)
rosmap_lasso.fit <- cv.glmnet(as.matrix(pheno_pathological_peer.data[,c(c("gender","age","mayo","msbb"),paste0("gPC",1:4),paste0("InferredCov",1:60))]), pheno_pathological_peer.data$FK_final_CaseControl, family = "binomial")
plot(rosmap_lasso.fit)
coef(rosmap_lasso.fit, s = "lambda.1se")
selected_covariates <- rownames(coef(rosmap_lasso.fit, s = 'lambda.1se'))[coef(rosmap_lasso.fit, s = 'lambda.1se')[,1]!= 0]
selected_peer <- selected_covariates[grep("InferredCov",selected_covariates)]
covariates <- c("age","mayo","msbb",paste0("gPC",1:4),selected_peer)

# Test total effect of gender, no peer
total_no_peer.fit <- glm(as.numeric(FK_final_CaseControl) - 1 ~ gender + age + mayo + msbb + gPC1 + gPC2 + gPC3 + gPC4, data = pheno_pathological_peer.data, family = "binomial")
summary(total_no_peer.fit)

total_coef_no_peer <- coef(total_no_peer.fit)["gender"]
exp(total_coef_no_peer)
exp(confint(total_no_peer.fit))

# Test total effect of gender, peer
total.fit <- glm(as.numeric(FK_final_CaseControl) - 1 ~ gender + age + mayo + msbb + gPC1 + gPC2 + gPC3 + gPC4 + InferredCov1 + InferredCov2 + InferredCov4 + InferredCov8 + InferredCov10 + InferredCov11 + InferredCov13 + InferredCov15 + InferredCov17 + InferredCov19 + InferredCov20 + InferredCov23 + InferredCov24 + InferredCov25 + InferredCov26 + InferredCov27 + InferredCov29 + InferredCov31 + InferredCov33 + InferredCov34 + InferredCov35 + InferredCov39 + InferredCov42 + InferredCov43 + InferredCov48 + InferredCov57 + InferredCov60, data = pheno_pathological_peer.data, family = "binomial")
summary(total.fit)
total_coef <- coef(total.fit)["gender"]

### Load and merge transcriptome
transcriptome.data <- fread("transcriptome/amp_ad/amp_ad_chrX.txt")
transcriptome.data$IID <- gsub("rosmapDLPFC_","",transcriptome.data$IID)
transcriptome.data$IID <- gsub("mayoTCX_","",transcriptome.data$IID)
transcriptome.data$IID <- gsub("msbbBA10_","",transcriptome.data$IID)
pheno_transcriptome.data <- merge(pheno_pathological_peer.data, transcriptome.data, by.x = "specimenID", by.y = "IID")

### Mediation model
# Exposure
A <- pheno_transcriptome.data$gender

# Mediator
transcripts_discovery <- names(transcriptome.data)[2:dim(transcriptome.data)[2]]
M <- pheno_transcriptome.data[,..transcripts_discovery]

# Outcome 
Y <- as.numeric(pheno_transcriptome.data$FK_final_CaseControl) - 1

# Covariates
C <- as.matrix(pheno_transcriptome.data[,..covariates])

### HDMA analysis with full covariates and without bootstrap
set.seed(20250501)
out_all <- mediate_hdma(A, M, Y, C, binary_y = T, n_include=659)
save(out_all, file = "results/out_all.RData")

### HDMA analysis with full covariates and with bootstrap
mediate_hdma_boot <- function(data, indices) {
  tryCatch(
    expr = {
      A_boot <- data[indices,gender]
      M_boot <- data[indices,..transcripts_discovery]
      Y_boot <- as.numeric(data[indices,]$FK_final_CaseControl) - 1
      C_boot <- as.matrix(data[indices,..covariates])
      # If no mediator is selected, set direct effect to total effect
      out_all_boot <- mediate_hdma(A_boot, M_boot, Y_boot, C_boot, binary_y = T, n_include=659)
      if(!is.null(out_all_boot)) return(out_all_boot$effects$estimate) else return(c(0,total_coef,total_coef))
    },
    error = function(e){
      return(c(0,total_coef,total_coef))
    })
}

set.seed(20250501)
mediate_hdma_boot_results <- boot(data=pheno_transcriptome.data, statistic=mediate_hdma_boot,
                                  R=1000, parallel = "multicore", ncpus = 70)

save(mediate_hdma_boot_results, file = "results/mediate_hdma_boot_results.Rdata")

boot.ci(mediate_hdma_boot_results, index = 1)
boot.ci(mediate_hdma_boot_results, index = 2)
boot.ci(mediate_hdma_boot_results, index = 3)

png(file = "figures/mediate_hdma_boot_results.png")
plot(mediate_hdma_boot_results)
dev.off()

### HDMA analysis with basic covariates and without bootstrap
set.seed(20250501)
out_all_no_peer <- mediate_hdma(A, M, Y, C[,c("age","mayo","msbb","gPC1","gPC2","gPC3","gPC4")], binary_y = T, n_include=659)
save(out_all_no_peer, file = "results/out_all_no_peer.RData")

### HDMA analysis with no peer factors and bootstrap
covariates_no_peer <- c("age","mayo","msbb","gPC1","gPC2","gPC3","gPC4")

mediate_hdma_boot_no_peer <- function(data, indices) {
  tryCatch(
    expr = {
      A_boot <- data[indices,gender]
      M_boot <- data[indices,..transcripts_discovery]
      Y_boot <- as.numeric(data[indices,]$FK_final_CaseControl) - 1
      C_boot <- as.matrix(data[indices,..covariates_no_peer])
      out_all_boot <- mediate_hdma(A_boot, M_boot, Y_boot, C_boot, binary_y = T, n_include=659)
      # If no mediator is selected, set direct effect to total effect
      if(!is.null(out_all_boot)) return(out_all_boot$effects$estimate) else return(c(0,total_coef_no_peer,total_coef_no_peer))
    },
    error = function(e){
      return(c(0,total_coef_no_peer,total_coef_no_peer))
    })
}

set.seed(20250501)
mediate_hdma_boot_results_no_peer <- boot(data=pheno_transcriptome.data, statistic=mediate_hdma_boot_no_peer,
                                          R=1000, parallel = "multicore", ncpus = 70)
save(mediate_hdma_boot_results_no_peer, file = "results/mediate_hdma_boot_results_no_peer.Rdata")

boot.ci(mediate_hdma_boot_results_no_peer, index = 1)
boot.ci(mediate_hdma_boot_results_no_peer, index = 2)
boot.ci(mediate_hdma_boot_results_no_peer, index = 3)

### Export individual genes
# Gene ID - ENSG dictionary
load("ENSG_relatedFiles/ensg_gene_names.Rdata")
ensg_gene_names.data <- unique(ensg_gene_names.data)

mediators_no_peer.data <- out_all_no_peer$contributions
mediators_no_peer.data <- merge(mediators_no_peer.data, ensg_gene_names.data, by.x = "mediator", by.y = "gene_id")
xci_no_peer <- mediators_no_peer.data[mediators_no_peer.data$alpha > 0 & mediators_no_peer.data$beta > 0, ]
dim(xci_no_peer)
write.csv(mediators_no_peer.data, file = "results/mediators_no_peer.csv")

mediators.data <- out_all$contributions
mediators.data <- merge(mediators.data, ensg_gene_names.data, by.x = "mediator", by.y = "gene_id")

# Genes with XCIe and positive mediation 
xci <- mediators.data[mediators.data$alpha > 0 & mediators.data$beta > 0, ]
dim(xci)

# Higher female-specific expression
mean(mediators.data[mediators.data$alpha > 0, "alpha"])
length(mediators.data[mediators.data$alpha > 0, "alpha"])

# Lower female-specific expression
mean(mediators.data[mediators.data$alpha < 0, "alpha"])
length(mediators.data[mediators.data$alpha < 0, "alpha"])

# Mean effect gene -> AD
exp(mean(abs(mediators.data$beta)))

write.csv(mediators.data, file = "results/mediators.csv")

# Associations between covariates and transcripts
covariates_reg.list <- mclapply(pheno_transcriptome.data[,..transcripts_discovery], function(transcript){
  pheno_transcriptome.data$transcript <- transcript
  covariates.fit <- glm(transcript ~ FK_final_CaseControl + gender + age + mayo + msbb + scale(gPC1) + scale(gPC2) + scale(gPC3) + scale(gPC4) + scale(InferredCov1) + scale(InferredCov2) + scale(InferredCov3) + scale(InferredCov4) + scale(InferredCov7) + scale(InferredCov8) + scale(InferredCov10) + scale(InferredCov11) + scale(InferredCov13) + scale(InferredCov14) + scale(InferredCov15) + scale(InferredCov17) + scale(InferredCov19) + scale(InferredCov20) + scale(InferredCov21) + scale(InferredCov24) + scale(InferredCov25) + scale(InferredCov26) + scale(InferredCov28) + scale(InferredCov29) + scale(InferredCov30) + scale(InferredCov31) + scale(InferredCov32) + scale(InferredCov33) + scale(InferredCov34) + scale(InferredCov36) + scale(InferredCov37) + scale(InferredCov39) + scale(InferredCov40) + scale(InferredCov41) + scale(InferredCov42) + scale(InferredCov44) + scale(InferredCov47) + scale(InferredCov50) + scale(InferredCov55) + scale(InferredCov56) + scale(InferredCov57) + scale(InferredCov60), data = pheno_transcriptome.data)
  coef(summary(covariates.fit))[,"Estimate"]
}, mc.cores = 8)
covariates_reg.data <- do.call(cbind, covariates_reg.list)

write.csv(covariates_reg.data, file = "results/covariates_reg.csv")
