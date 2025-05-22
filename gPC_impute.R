library(missMDA)
library(data.table)

###### Phenotype prep
pheno.data <- fread("phenotype/amp_ad_pheno.csv")

# Center gender and code female as 1, male as 0
pheno.data$gender <- as.numeric(as.factor(pheno.data$Sex)) - 1
pheno.data$gender <- 1-pheno.data$gender

# Specify continuous age variable 
pheno.data$age <- pheno.data$AgeAtDeath
pheno.data$age[pheno.data$age == "90_or_above"] <- 90
pheno.data$age <- as.numeric(pheno.data$age)

# Define ADvsCon
pheno_pathological.data <- pheno.data[pheno.data$FK_final_CaseControl != "Other", ]
pheno_pathological.data$FK_final_CaseControl <- ordered(pheno_pathological.data$FK_final_CaseControl, levels = c("Control","AD"))

# Individuals with missing WGS data
dim(unique(pheno_pathological.data[is.na(pheno_pathological.data$WGS_Availability), "individualID"]))

# Define race variables
pheno_pathological.data$RaceAncestry[pheno_pathological.data$RaceAncestry == "Black_or_AfricanAmerican"] <- "Black"
pheno_pathological.data$RaceAncestry <- as.factor(pheno_pathological.data$RaceAncestry)

old_phenotype.data <- fread("phenotype/AMPAD_and_EADBBelgian_FullCovariatesAndSampleInfo_Correlation_19May2022.csv")
old_phenotype.data <- old_phenotype.data[,c("Sample","gPC1","gPC2","gPC3","gPC4")]

pheno_pathological_wgs.data <- merge(pheno_pathological.data, old_phenotype.data, by.x = "specimenID", by.y = "Sample", all.x = T)

pre_imputation.data <- pheno_pathological_wgs.data[,c("FK_final_CaseControl","gender","age","RaceAncestry","gPC1","gPC2","gPC3","gPC4")]
imputed <- imputeFAMD(pre_imputation.data, ncp = 6)
gPC_imputed.data <- imputed$completeObs[c("gPC1","gPC2","gPC3","gPC4")]
gPC_imputed.data$specimenID <- pheno_pathological_wgs.data$specimenID

save(gPC_imputed.data, file = "phenotype/gPC_imputed.Rdata")