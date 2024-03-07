library(hdmed)
library(data.table)
library(dplyr)
library(lavaan)
library(parallel)
library(glmnet)
library(boot)

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
wgs_info_pathological.data$FK_final_CaseControl

# Create dummy variable for Mayo and MSBB
wgs_info_pathological.data$mayo <- 0
wgs_info_pathological.data$mayo[wgs_info_pathological.data$Cohort == "MayoRNASeq_TCX"] <- 1
wgs_info_pathological.data$msbb <- 0
wgs_info_pathological.data$msbb[wgs_info_pathological.data$Cohort == "MSBB_BA10"] <- 1

# Reduce to ROSMAP, Mayo and MSBB BA10
wgs_info_pathological.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "ROSMAP_DLPFC" | 
                                                           wgs_info_pathological.data$Cohort == "MayoRNASeq_TCX" |
                                                           wgs_info_pathological.data$Cohort == "MSBB_BA10"]


# Determine optimal PEER factor combination
peer.data <- fread("../../fahri/Alex_xChr/PEER_rosmapDLPFC_mayoTCX_msbbBA10_innerJoint_OLDnormalizedExp_27June2023.PEER_covariates.txt", data.table = F)
names(peer.data) <- gsub("rosmapDLPFC_","",names(peer.data))
names(peer.data) <- gsub("mayoTCX_","",names(peer.data))
names(peer.data) <- gsub("msbbBA10_","",names(peer.data))
rownames(peer.data) <- peer.data$ID
peer.data <- peer.data[,-1]
peer.data <- as.data.frame(t(peer.data))
peer.data$ID <- row.names(peer.data)
wgs_info_pathological.data <- merge(wgs_info_pathological.data, peer.data, by = "ID")

set.seed(20230711)
rosmap_lasso.fit <- cv.glmnet(as.matrix(wgs_info_pathological.data[,c(c("gender","age","mayo","msbb"),paste0("gPC",1:4),paste0("InferredCov",1:60))]), wgs_info_pathological.data$FK_final_CaseControl, family = "binomial")
plot(rosmap_lasso.fit)
coef(rosmap_lasso.fit, s = "lambda.1se")
selected_covariates <- rownames(coef(rosmap_lasso.fit, s = 'lambda.1se'))[coef(rosmap_lasso.fit, s = 'lambda.1se')[,1]!= 0]
selected_peer <- selected_covariates[grep("InferredCov",selected_covariates)]
covariates <- c("age","mayo","msbb",paste0("gPC",1:4),selected_peer)

# Test total effect of gender
glm(FK_final_CaseControl ~ gender + age + mayo + msbb, data = wgs_info_pathological.data)

# Stratify by cohort
wgs_info_rosmap.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "ROSMAP_DLPFC", ]
wgs_info_mayo.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "MayoRNASeq_TCX", ]
wgs_info_msbb_ba10.data <- wgs_info_pathological.data[wgs_info_pathological.data$Cohort == "MSBB_BA10", ]

### Load transcriptoms
rosmap_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/ROSMAP/rosmap_dlpfc_chrX.txt")
mayo_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/Mayo/mayo_tcx_chrX.txt")
msbb_ba10_transcriptome.data <- fread("../amp_ad_genotypes/transcriptome/MSBB_BA10/msbb_ba10_chrX.txt")

# Rosmap
wgs_info_transcriptome_rosmap.data <- merge(wgs_info_rosmap.data, rosmap_transcriptome.data, by.x = "ID", by.y = "IID")

# Mayo
wgs_info_mayo.data$ID <- as.numeric(wgs_info_mayo.data$ID) 
wgs_info_transcriptome_mayo.data <- merge(wgs_info_mayo.data, mayo_transcriptome.data, by.x = "ID", by.y = "IID")
wgs_info_transcriptome_mayo.data$ID <- as.character(wgs_info_transcriptome_mayo.data$ID)

# Merge ROSMAP and MAYO
wgs_info_transcriptome_rosmap_mayo.data <- bind_rows(wgs_info_transcriptome_rosmap.data, wgs_info_transcriptome_mayo.data)

### MSBB BA10
wgs_info_msbb_ba10.data$ID <- as.numeric(wgs_info_msbb_ba10.data$ID)
wgs_info_transcriptome_msbb_ba10.data <- merge(wgs_info_msbb_ba10.data, msbb_ba10_transcriptome.data, by.x = "ID", by.y = "IID")
wgs_info_transcriptome_msbb_ba10.data$ID <- as.character(wgs_info_transcriptome_msbb_ba10.data$ID)

wgs_info_transcriptome_rosmap_mayo_msbb_ba10.data <- bind_rows(wgs_info_transcriptome_rosmap.data, wgs_info_transcriptome_mayo.data, wgs_info_transcriptome_msbb_ba10.data)

# Exposure
A <- wgs_info_transcriptome_rosmap_mayo_msbb_ba10.data$gender

# Mediator
transcripts <- intersect(names(rosmap_transcriptome.data), names(mayo_transcriptome.data))
transcripts_discovery <- intersect(transcripts, names(msbb_ba10_transcriptome.data))[-1]
M <- wgs_info_transcriptome_rosmap_mayo_msbb_ba10.data[,..transcripts_discovery]

# Outcome 
Y <- as.numeric(wgs_info_transcriptome_rosmap_mayo_msbb_ba10.data$FK_final_CaseControl) - 1

# Covariates
C <- as.matrix(wgs_info_transcriptome_rosmap_mayo_msbb_ba10.data[,..covariates])

### HDMA analysis with full covariates and without bootstrap
set.seed(20230711)
out_all <- mediate_hdma(A, M, Y, C, binary_y = T, n_include=656)
save(out_all, file = "results/out_all.RData")
mediators <- out_all$contributions$mediator

### HDMA analysis with full covariates and with bootstrap
mediate_hdma_boot <- function(data, indices) {
  tryCatch(
    expr = {
      A_boot <- data[indices,gender]
      M_boot <- data[indices,..mediators]
      Y_boot <- as.numeric(data[indices,]$FK_final_CaseControl) - 1
      C_boot <- as.matrix(data[indices,..covariates])
      out_all_boot <- mediate_hdma(A_boot, M_boot, Y_boot, C_boot, binary_y = T, n_include=21)
      if(!is.null(out_all_boot)) return(out_all_boot$effects$estimate) else return(c(0,0,0.09278296))
    },
    error = function(e){
      return(c(0,0,0.09278296))
    })
}

set.seed(20230730)
mediate_hdma_boot_results <- boot(data=wgs_info_transcriptome_rosmap_mayo_msbb_ba10.data, statistic=mediate_hdma_boot,
                                   R=50, parallel = "multicore", ncpus = 60)
boot.ci(mediate_hdma_boot_results)

save(mediate_hdma_boot_results, file = "results/mediate_hdma_boot_results.Rdata")
boot.ci(mediate_hdma_boot_results_no_peer)

png(file = "figures/mediate_hdma_boot_results.png")
plot(mediate_hdma_boot_results)
dev.off()

### HDMA analysis with basic covariates and without bootstrap
set.seed(20230711)
out_all_no_peer <- mediate_hdma(A, M, Y, C[,c("age","mayo","msbb","gPC1","gPC2","gPC3","gPC4")], binary_y = T, n_include=656)
save(out_all_no_peer, file = "results/out_all_no_peer.RData")
mediators_no_peer <- out_all_no_peer$contributions$mediator

### HDMA analysis with full covariates and bootstrap
covariates_no_peer <- c("age","mayo","msbb","gPC1","gPC2","gPC3","gPC4")
                                 
mediate_hdma_boot_no_peer <- function(data, indices) {
  tryCatch(
    expr = {
      A_boot <- data[indices,gender]
      M_boot <- data[indices,..mediators_no_peer]
      Y_boot <- as.numeric(data[indices,]$FK_final_CaseControl) - 1
      C_boot <- as.matrix(data[indices,..covariates_no_peer])
      out_all_boot <- mediate_hdma(A_boot, M_boot, Y_boot, C_boot, binary_y = T, n_include=21)
      if(!is.null(out_all_boot)) return(out_all_boot$effects$estimate) else return(c(0,0,0.09278296))
    },
    error = function(e){
      return(c(0,0,0.09278296))
    })
}

set.seed(20230730)
mediate_hdma_boot_results_no_peer <- boot(data=wgs_info_transcriptome_rosmap_mayo_msbb_ba10.data, statistic=mediate_hdma_boot_no_peer,
                                          R=1000, parallel = "multicore", ncpus = 40)
save(mediate_hdma_boot_results_no_peer, file = "results/mediate_hdma_boot_results_no_peer.Rdata")
boot.ci(mediate_hdma_boot_results_no_peer)

### Export individual genes
# Gene ID - ENSG dictionary
load("../amp_ad_genotypes/ENSG_relatedFiles/ensg_gene_names.Rdata")
ensg_gene_names.data <- unique(ensg_gene_names.data)

mediators_no_peer.data <- out_all_no_peer$contributions
mediators_no_peer.data <- merge(mediators_no_peer.data, ensg_gene_names.data, by.x = "mediator", by.y = "gene_id")
xci_no_peer <- mediators_no_peer.data[mediators_no_peer.data$alpha > 0 & mediators_no_peer.data$beta > 0, ]
dim(xci_no_peer)
write.csv(mediators_no_peer.data, file = "results/mediators_no_peer.csv")

mediators.data <- out_all$contributions
mediators.data <- merge(mediators.data, ensg_gene_names.data, by.x = "mediator", by.y = "gene_id")
xci <- mediators.data[mediators.data$alpha > 0 & mediators.data$beta > 0, ]
dim(xci)
write.csv(mediators.data, file = "results/mediators.csv")

