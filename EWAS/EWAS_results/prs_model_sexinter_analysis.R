library(limma)
# Load cell type adjusted beta values
load("adjusted_betas.Rdata")

# Perform analysis of full model:
# Methylation (beta values adjusted for cell type) ~ 1 + age + sex + PRS + 3 Genomic PC's + Smoking Score + Medication History + Missing PRS genotypes + PRS x Sex
sample.annot <- sample.annot[match(colnames(adj.betas),sample.annot$patientID),]

temp.status <- factor(sample.annot$status, levels=c("control","case"))
temp.sex <- factor(sample.annot$sex, levels=c("Female","Male"))
temp.design <- model.matrix(~sample.annot$PRS+temp.sex+sample.annot$Age+sample.annot$PC1+sample.annot$PC2+sample.annot$PC3+sample.annot$smoking.score+sample.annot$Medication+sample.annot$Missing.genos+sample.annot$PRS:temp.sex)
results.fit <- lmFit(adj.betas, temp.design)
results.fit <- eBayes(results.fit)

# Save the results object
save(results.fit,file="prs_model_fit_object.Rdata")
