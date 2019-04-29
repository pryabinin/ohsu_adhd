# This script desribes the process of cleaning genotype data. This script is based off of the GWASTools vignette.
# Peter Ryabinin (ryabinip@ohsu.edu)

library(GWASTools)

# Load genotype data
gds.fn <- "genotypes.gds"
(gds <- GdsGenotypeReader(gds.fn))
snp.annotation <- read.delim("/data/Peter/projects/ADHD/data/snpAnnotation.txt",stringsAsFactors = F)
scan.annotation <- read.delim("/data/Peter/projects/ADHD/data/sampleAnnotation.txt",stringsAsFactors = F)
names(scan.annotation)[1]<-"subjectID"

snpAnnot <- SnpAnnotationDataFrame(snp.annotation)
scanAnnot <- ScanAnnotationDataFrame(scan.annotation)
# The data was run in two batches
scanAnnot$batch <- ifelse(scanAnnot$plate %in% c("ADHD_Nigg_Psych_RW_01","ADHD_Nigg_Psych_RW_02"),"rerun","original")
genoData <- GenotypeData(gds,snpAnnot=snpAnnot,scanAnnot = scanAnnot)

# Remove probes that are missing from the latest version of the Illumina manifest
illum <- read.delim("InfiniumPsychArray-24v1-1_A1_StrandReport_FDT.txt",skip=5,stringsAsFactors = F)
removed.snps <- snpAnnot$snpName[!(snpAnnot$snpName %in% illum$SNP_Name)]
removed.snp.ids <- snpAnnot$snpID[!(snpAnnot$snpName %in% illum$SNP_Name)]
filteredsnps <- data.frame(snpID = removed.snp.ids, snpName=removed.snps,
                           reason="Removed by Illumina")

# Calculate missing call rate of each SNP and how many SNPs would be removed at various thresholds
miss <- missingGenotypeBySnpSex(genoData,verbose = F)
allequal(snpAnnot$snpID, as.numeric(names(miss$missing.fraction)))
snpAnnot$missing.n1 <- miss$missing.fraction
varMetadata(snpAnnot)["missing.n1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all samples",
  "except that females are excluded for Y chr SNPs")
num_snps_lost <- c(sum(snpAnnot$missing.n1>0.02),sum(snpAnnot$missing.n1>0.03),sum(snpAnnot$missing.n1>0.04),sum(snpAnnot$missing.n1>0.05))
perc_snps_lost <- c(sum(snpAnnot$missing.n1>0.02)/nsnp(genoData)*100,sum(snpAnnot$missing.n1>0.03)/nsnp(genoData)*100,sum(snpAnnot$missing.n1>0.04)/nsnp(genoData)*100,sum(snpAnnot$missing.n1>0.05)/nsnp(genoData)*100)

lostdf <- data.frame("SNPs_lost" = num_snps_lost, "Perc_SNPs_Lost" = round(perc_snps_lost,2))
graphics.off()
theme <- ttheme_default(core=list(fg_params=list(cex=2.0)))
grid.table(lostdf,rows=c("98% Call Rate","97% Call Rate","96% Call Rate","95% Call Rate"),cols=c("# of SNPs lost","% of SNPs lost"),theme=theme)

hist((1-snpAnnot$missing.n1[!(snpAnnot$snpID %in% filteredsnps[,1])])*100, ylim=c(0,10000),
     xlab="SNP call rate (%)",
     main="Histogram of SNP Call Rate")

# Remove 0% call rate SNPs and calculate call rate of each sample

sum(snpAnnot$missing.n1 == 1)
# Create a variable that contains the IDs of these SNPs to exclude
snpexcl <- snpAnnot$snpID[snpAnnot$missing.n1 == 1 | snpAnnot$snpID %in% filteredsnps[,1]]
length(snpexcl)
# Calculate the missing call rate per sample
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
names(miss)
head(miss$missing.counts)
head(miss$snps.per.chr)
# Check to make sure that the correct number of SNPs were excluded
sum(miss$snps.per.chr)
nrow(snpAnnot) - sum(miss$snps.per.chr)

head(miss$missing.fraction)
# Check the ordering matches the sample annotation file
allequal(names(miss$missing.fraction), scanAnnot$scanID)
# Add the missing call rates vector to the sample annotation file
scanAnnot$missing.e1 <- miss$missing.fraction
varMetadata(scanAnnot)["missing.e1", "labelDescription"] <- paste(
  "fraction of genotype calls missing over all snps with missing.n1<1",
  "except that Y chr SNPs are excluded for females")
summary(scanAnnot$missing.e1)

hist(100*(1-scanAnnot$missing.e1),
     xlab="Sample call rate (%)",
     main="Call Rate for Samples")

snpexcl <- snpAnnot$snpID[snpAnnot$missing.n1 >= 0.03 | snpAnnot$snpID %in% filteredsnps[,1]]
miss <- missingGenotypeByScanChrom(genoData, snp.exclude=snpexcl)
hist(100*(1-miss$missing.fraction),
     xlab="Sample call rate (%)",
     main="Call Rate for Samples using SNPs with >97% call rate")

# Since no samples had low call rate in this experiment, remove low call rate probes (call rate < 97%) according to first calculation

temp <- snpAnnot$snpID[snpAnnot$missing.n1 >=0.03]
temp <- temp[!(temp %in% filteredsnps$snpID)]

filteredsnps <- rbind(filteredsnps,data.frame(snpID=temp,snpName=snpAnnot$snpName[match(temp,snpAnnot$snpID)],
                                              reason="Call rate < 97%"))

# Calculate call rates for autosomes and X chromosome for future anaylsis.
# Also, find duplicate samples and make all but the one with the highest call rate as "duplicate" so they can be found and excluded from some future analysis.
auto <- colnames(miss$missing.counts) %in% 1:22
missa <- rowSums(miss$missing.counts[,auto]) / sum(miss$snps.per.chr[auto])
summary(missa)
missx <- miss$missing.counts[,"X"] / miss$snps.per.chr["X"]
summary(missx)
# save the order of scanAnnot
scanAnnot$order <- 1:length(getScanID(scanAnnot))
# check they match sample annotation file
allequal(names(missa), scanAnnot$scanID)
allequal(names(missx), scanAnnot$scanID)
# Add these separate sample missing call rates to the sample annotation
scanAnnot$miss.e1.auto <- missa
scanAnnot$miss.e1.xchr <- missx
# Order scanAnnot by missing.e1 so duplicate subjectIDs with a higher missing rate are marked as duplicates
scanAnnot <- scanAnnot[order(scanAnnot$subjectID, scanAnnot$missing.e1),]
scanAnnot$duplicated <- duplicated(scanAnnot$subjectID)
table(scanAnnot$duplicated, useNA="ifany")
# Put scanAnnot back in original order
scanAnnot <- scanAnnot[order(scanAnnot$order),]
scanAnnot$order<-NULL
allequal(scanAnnot$scanID, getScanID(genoData))
varMetadata(scanAnnot)["duplicated", "labelDescription"] <- "TRUE for duplicate scan with higher missing.e1"

# Calculate missing call rate per plate
batches <- unique(scanAnnot$plate)
bmiss <- rep(NA,length(batches)); names(bmiss) <- batches
bn <- rep(NA,length(batches)); names(bn) <- batches
for(i in 1:length(batches)) {
  x <- scanAnnot$missing.e1[is.element(scanAnnot$plate, batches[i])]
  bmiss[i] <- mean(x)
  bn[i] <- length(x)
}
barplot(100*(1-bmiss),names.arg = 1:11,main="Sample call rates across plates",xlab="Plate number",ylab="Call rate",ylim=c(0,100))
bmiss.plate <- bmiss
bn.plate <- bn

# Calculate missing call rate per batch
batches <- unique(scanAnnot$batch)
bmiss <- rep(NA,length(batches)); names(bmiss) <- batches
bn <- rep(NA,length(batches)); names(bn) <- batches
for(i in 1:length(batches)) {
  x <- scanAnnot$missing.e1[is.element(scanAnnot$batch, batches[i])]
  bmiss[i] <- mean(x)
  bn[i] <- length(x)
}
barplot(100*(1-bmiss),main="Sample call rates across batches",xlab="Batch",ylab="Call rate",ylim=c(0,100))
100*(1-bmiss)
bn
# All plates and batches have about the same call rate, the reruns are not significantly different.

# Calculate regression between mean call rate and size of batch.
y <- lm(100*(1-bmiss.plate) ~ bn.plate)
summary(y)
plot(bn.plate, 100*(1-bmiss.plate),
     xlab="Number of samples per plate", ylab="Mean call rate",
     main="Mean Call Rate vs\nSamples per Plate")
abline(y$coefficients)

# Calculate chi-squared test for allelic frequency between plates by looking at how each plate differs from all other plates.

res <- batchChisqTest(genoData, batchVar="plate", return.by.snp=TRUE,
                      snp.include = snpAnnot$snpID[!(snpAnnot$snpID %in% filteredsnps[,1]) &
                                                     snpAnnot$chromosome %in% c(1:22,24,26,27)])
res$pval <- pchisq(res$chisq,1,lower.tail=F)
res$mean.pval <- apply(res$pval,2,mean,na.rm=T)
sum(res$pval[!is.na(res$pval)] < 0.05)
res$adj.p <- apply(res$pval,2,p.adjust,method="BH",n=(length(res$pval)))
res$mean.adj.p <- apply(res$adj.p,2,mean,na.rm=T)
sum(res$adj.p[!is.na(res$adj.p)] < 0.05)

# There are many instances of SNPs being over or underrepresented in each batch but after correcting for multiple testing that number becomes very small

# Calculate chi-squared statistic to find relationship between plate and self reported race.
x <- table(scanAnnot$race1, useNA="ifany")
x
x[1] / sum(x)
x[2] / sum(x)
x <- table(scanAnnot$race1, scanAnnot$plate)
x
# Run an approximate chi-square test to see if there are ethnic effects
chisq <- chisq.test(x)
chisq$p.value
# Calculate the fraction of samples in each batch that are CEU
batches <- unique(scanAnnot$plate)
eth <- rep(NA,length(batches)); names(eth) <- batches 
for(i in 1:length(batches)){
  x <- scanAnnot$race1[is.element(scanAnnot$plate, batches[i])]
  xl <- length(x[x == "White/Middle Eastern"])
  eth[i] <- xl / length(x)
}
eth <- eth[order(names(eth))] 
allequal(names(eth), names(res$mean.chisq))

# Plot ethnicity against average chi-square statistic for each plates, looking for plates with similar ethnicity representation but different genotyping to indiciate a possible batch effect.
# Plot the average Chi-Square test statistic against the
#     fraction of samples that are CEU
l.mod <- lm(res$mean.chisq ~ eth)
summary(l.mod)
plot(eth, res$mean.chisq, xlab="Fraction of White/Middle Eastern Samples per Plate",
     ylab="Average Chi-square Test Statistic",
     main="Fraction of White/Middle Eastern Samples per Plate
     vs Average Chi-square Test Statistic\n P-value: .53")
text(eth, res$mean.chisq, 1:11, cex=1, pos=4, col="blue")
text(0.925,0.775,"plate number",col="blue")
abline(l.mod$coefficients)

# Check how many SNPs remain over/under represented per plate after multiple testing correction since plate with high average chi square statstic may not necessarily have a higher number of significanly different SNPs

res$sig.per.plate <- apply(res$adj.p,2, function(x) { sum(x[!is.na(x)]<0.05)})
res$sig.per.plate
graphics.off()
grid.table(as.data.frame(res$sig.per.plate),cols="Number of SNPs with differing allele frequencies",
           rows=paste("plate",1:11,sep=" "))

# Check if minor allele counts are very low, if so, consider using Fisher exact test since the chi-squared test be inaccurate.
resFish <- batchFisherTest(genoData, batchVar="plate", return.by.snp=TRUE,
                           snp.include = snpAnnot$snpID[!(snpAnnot$snpID %in% filteredsnps[,1]) &
                                                          snpAnnot$chromosome %in% c(1:22,24,26,27)])

minallele <- apply(resFish$allele.counts,1,min)
hist(minallele,xlab= "Minor allele counts",main="Histogram of minor allele counts")

hist(minallele[minallele<50 & minallele>0],xlab= "Minor allele counts",main="Histogram of minor allele counts > 0 and < 50")

sig.snps <- which(apply(res$adj.p,1,min)<0.05)
hist(minallele[sig.snps],main="Histogram of minor allele counts for SNPs that differ across plates",xlab="Minor allele counts")

barplot(table(minallele[1:length(minallele) %in% sig.snps & minallele<10]),main="Histogram of minor allele counts for SNPs that differ\n across plates with minor allele counts < 10",xlab="Minor allele counts")

# Many SNPs have very low minor allele counts, Fisher exact test may be appropriate instead of the chi-squared test
resFish$mean.pval <- apply(resFish$pval,2,mean,na.rm=T)
resFish$adj.p <- apply(resFish$pval,2,p.adjust,method="BH",n=(length(res$pval)))
resFish$mean.adj.p <- apply(resFish$adj.p,2,mean,na.rm=T)
resFish$sig.per.plate <- apply(resFish$adj.p,2, function(x) { sum(x[!is.na(x)]<0.05)})
resFish$sig.per.plate
graphics.off()
grid.table(as.data.frame(resFish$sig.per.plate),cols="Number of SNPs with differing allele frequencies",
           rows=paste("plate",1:11,sep=" "))

sig.snps <- which(apply(resFish$adj.p,1,min)<0.05)
hist(minallele[sig.snps],xlab= "Minor allele counts",main="Histogram of minor allele counts for SNPs that differ across plates by Fisher exact test")
barplot(table(minallele[1:length(minallele) %in% sig.snps & minallele<20]),xlab= "Minor allele counts",main="Histogram of minor allele counts < 20 by Fisher exact test")

resFish <- batchFisherTest(genoData, batchVar="batch", return.by.snp=TRUE,
                           snp.include = snpAnnot$snpID[!(snpAnnot$snpID %in% filteredsnps[,1]) &
                                                          snpAnnot$chromosome %in% c(1:22,24,26,27)])
resFish$mean.pval <- apply(resFish$pval,2,mean,na.rm=T)
resFish$adj.p <- apply(resFish$pval,2,p.adjust,method="BH",n=(length(res$pval)))
resFish$mean.adj.p <- apply(resFish$adj.p,2,mean,na.rm=T)
resFish$sig.per.plate <- apply(resFish$adj.p,2, function(x) { sum(x[!is.na(x)]<0.05)})
resFish$sig.per.plate
resFish$sig.per.plate.ids <- apply(resFish$adj.p,2, function(x) { names(x[!is.na(x) & x<0.05])})

sig.snps <- which(apply(resFish$adj.p,1,min)<0.05)
barplot(table(minallele[1:length(minallele) %in% sig.snps]),xlab= "Minor allele counts",main="Histogram of minor allele counts using Fisher exact test between batches")

# Now look at X SNPs using only female samples (due to differences in sexes)
resFishX <- batchFisherTest(genoData, batchVar="batch", return.by.snp=TRUE,
                            snp.include = snpAnnot$snpID[!(snpAnnot$snpID %in% filteredsnps[,1]) &
                                                           snpAnnot$chromosome == 23],sex.include = "F")

resFishX$mean.pval <- apply(resFishX$pval,2,mean,na.rm=T)
resFishX$adj.p <- apply(resFishX$pval,2,p.adjust,method="BH",n=(length(resFishX$pval)))
resFishX$mean.adj.p <- apply(resFishX$adj.p,2,mean,na.rm=T)
resFishX$sig.per.plate <- apply(resFishX$adj.p,2, function(x) { sum(x[!is.na(x)]<0.05)})
resFishX$sig.per.plate
resFishX$sig.per.plate.ids <- apply(resFishX$adj.p,2, function(x) { names(x[!is.na(x) & x<0.05])})

# Now look at SNPs on the Y chromosome using only male samples
resFishY <- batchFisherTest(genoData, batchVar="batch", return.by.snp=TRUE,
                            snp.include = snpAnnot$snpID[!(snpAnnot$snpID %in% filteredsnps[,1]) &
                                                           snpAnnot$chromosome == 25],sex.include = "M")

resFishY$mean.pval <- apply(resFishY$pval,2,mean,na.rm=T)
resFishY$adj.p <- apply(resFishY$pval,2,p.adjust,method="BH",n=(length(resFishY$pval)))
resFishY$mean.adj.p <- apply(resFishY$adj.p,2,mean,na.rm=T)
resFishY$sig.per.plate <- apply(resFishY$adj.p,2, function(x) { sum(x[!is.na(x)]<0.05)})
resFishY$sig.per.plate
resFishY$sig.per.plate.ids <- apply(resFishY$adj.p,2, function(x) { as.integer(names(x[!is.na(x) & x<0.05]))})

# Create list of SNPs that differ between batches and plot their SNP clusters to see if they should be removed (no SNPs on the Y chromosome need to be removed)
batcheffectsnpids <- unique(as.integer(c(resFish$sig.per.plate.ids[,1], resFishX$sig.per.plate.ids[1])))

qxyfile <- "intensity.gds"
qualGDS <- GdsIntensityReader(qxyfile)
qualData <- IntensityData(qualGDS, scanAnnot=scanAnnot,snpAnnot = snpAnnot)

named.chrom <- sapply(snpAnnot$chromosome[batcheffectsnpids], function(x) {
  switch(as.character(x),"23"="X","24"="XY","25"="Y","26"="MT","27"="UNK",as.character(x))
})
main.txt <- paste(snpAnnot$snpName[batcheffectsnpids],paste("chrom",named.chrom,sep = "_"),
                  paste("pos",snpAnnot$position[batcheffectsnpids],sep="_"),sep=" ")
pdf("batch_effected_snps.pdf")
genoClusterPlot(qualData,genoData,snpID=batcheffectsnpids,main.txt=main.txt,scan.hilite = scanAnnot$scanID[scanAnnot$batch == "rerun"],start.axis.at.0 = T)
dev.off()

# Add batch effected SNPs to filtered SNP list
temp <- batcheffectsnpids
temp <- batcheffectsnpids[!(batcheffectsnpids %in% filteredsnps$snpID)]
filteredsnps <- rbind(filteredsnps,data.frame(snpID=temp,snpName=snpAnnot$snpName[match(temp,as.character(snpAnnot$snpID))],reason="Batch effected"))

# Calculate missing genotypes and look at missingness across chromosomes.
miss <- missingGenotypeByScanChrom(genoData,snp.exclude = filteredsnps[,1])
miss.rate <- t(apply(miss$missing.counts, 1, function(x) {
  x / miss$snps.per.chr}))
miss.rate <- as.data.frame(miss.rate)
miss.rate$Y[scanAnnot$sex=="F" | is.na(scanAnnot$sex)] <- NA
cols <- names(miss.rate) %in% c(1:22, "X", "XY","Y","M","U")
boxplot(100*(1-miss.rate[,cols]), main="Call Rate by Chromosome",
        ylab="Call Rate", xlab="Chromosome")

# Calculate missing call rate per sex on the X chromosome
boxplot(100*(1-miss.rate$X) ~ scanAnnot$sex,
        main="X Chromosome Call Rate by Sex",
        ylab="Call Rate")

# Calculate heterozygosity
# Calculate heterozygosity by scan by chromosome
het.results <- hetByScanChrom(genoData,snp.exclude = filteredsnps[,1])
#close(genoData)
# Ensure heterozygosity results are ordered correctly
allequal(scanAnnot$scanID, rownames(het.results))
# Write autosomal and X chr heterozygosity to sample annot
scanAnnot$het.A <- het.results[,"A"]
scanAnnot$het.X <- het.results[,"X"]
varMetadata(scanAnnot)["het.A", "labelDescription"] <-
  "fraction of heterozygotes for autosomal SNPs"
varMetadata(scanAnnot)["het.X", "labelDescription"] <-
  "fraction of heterozygotes for X chromosome SNPs"

# Plot heterozygosity by race
boxplot(scanAnnot$het.A ~ scanAnnot$race1,
        main="Autosomal Heterozygosity")

# Plot heterozygosity on the X chromosome of female samples
female <- scanAnnot$sex == "F"
boxplot(scanAnnot$het.X[female] ~ scanAnnot$race1[female],
        main="X Chromosome Heterozygosity in Females")

# Mis-annotated sex check
inten.by.chrom <- meanIntensityByScanChrom(qualGDS,snp.exclude = filteredsnps[,1])
names(inten.by.chrom)

mninten <- inten.by.chrom[[1]]  # mean intensities
dim(mninten)
# Check to be sure sample ordering is consistent
allequal(scanAnnot$scanID, rownames(mninten))
# Assign each sex a color
xcol <- rep(NA, nrow(scanAnnot))
xcol[scanAnnot$sex == "M"] <- "blue"
xcol[scanAnnot$sex == "F"] <- "red"
nx <- sum(snpAnnot$chromosome == 23)
ny <- sum(snpAnnot$chromosome == 25)

# Create plots to check for misannotated sexes
x1 <-mninten[,"X"]; y1 <- mninten[,"Y"]
main1 <- "Mean X vs \nMean Y Chromosome Intensity"
#Het on X vs X intensity
x2 <- mninten[,"X"]; y2 <- scanAnnot$het.X
main2 <- "Mean X Chromosome Intensity vs
Mean X Chromosome Heterozygosity"
# Het on X vs Y intensity
y3 <- mninten[,"Y"]; x3 <- scanAnnot$het.X
main3 <- "Mean X Chromosome Heterozygosity vs
Mean Y Chromosome Intensity"
# X vs A het
x4 <- scanAnnot$het.A[scanAnnot$sex == "F"]
y4 <- scanAnnot$het.X[scanAnnot$sex == "F"]
main4 <- "Mean Autosomal Heterozygosity vs
Mean X Chromosome Heterozygosity\nFemales only"
cols <- c("blue","red")
mf <- c("male", "female")
xintenlab <- paste("X intensity (n=", nx, ")", sep="")
yintenlab <- paste("Y intensity (n=", ny, ")", sep="")
pdf("/data/Peter/projects/ADHD/data/DataCleaning-sex.pdf")
par(mfrow=c(2,2))
plot(x1, y1, xlab=xintenlab, ylab=yintenlab,
     main=main1, col=xcol, cex.main=0.8)
legend("topright",mf,col=cols,pch=c(1,1))
plot(x2, y2, col=xcol, xlab=xintenlab,
     ylab="X heterozygosity", main=main2, cex.main=0.8)
plot(x3, y3, col=xcol, ylab=yintenlab,
     xlab="X heterozygosity", main=main3, cex.main=0.8)
col4 <- ifelse(scanAnnot$race1=="White/Middle Eastern","red","green")
plot(x4,y4, col=col4, xlab="Autosomal heterozygosity",
     ylab="X heterozygosity", main=main4, cex.main=0.8)
legend("topleft",c("White","Other"),col=c("red","green"),pch=c(1,1))
dev.off()
dev.off()

# Summarize the SNPs removed so far
lostdf <- data.frame("Reason for removal" = names(table(filteredsnps$reason)),
                     "Number of SNPs removed" = as.integer(table(filteredsnps$reason)),
                     "Percentage of SNPs Removed" = as.numeric(round(table(filteredsnps$reason)/nrow(manifest)*100,2)))
graphics.off()
theme <- ttheme_default(core=list(fg_params=list(cex=2.0)))
grid.table(lostdf,rows=c("","",""),cols=c("Reason for removal","Number of SNPs removed","Percentage of all SNPs"),theme=theme)

# Calculate Population structure and relatedness. Initial estimate is from the KING method in SNPRelate, and then use PC-AiR and PC-Relate to calculate relatedness
library(GENESIS)
library(SNPRelate)
library(data.table)

close(genoData)
gdsfile <- "genotypes.gds"
gdsobj <- snpgdsOpen(gdsfile)
snpset <- snpgdsLDpruning(gdsobj,
                          snp.id=snp.sel[!(snp.sel %in% filteredsnps[,1])],
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(0.1))
snp.pruned <- unlist(snpset, use.names=FALSE)
ibdobj <- snpgdsIBDKING(gdsobj,num.thread=16,autosome.only = T,
                        snp.id=snp.pruned,type="KING-robust")
snpgdsClose(gdsobj)

# Setup pedigree, note that the mother's ID ends with 4 and father's ID ends with 5
scanAnnot$mother <- paste(substr(as.character(scanAnnot$subjectID),1,nchar(as.character(scanAnnot$subjectID))-1),"4",sep="")
scanAnnot$father <- paste(substr(as.character(scanAnnot$subjectID),1,nchar(as.character(scanAnnot$subjectID))-1),"5",sep="")
ped <- pData(scanAnnot)[,c("family", "subjectID", "father", "mother", "sex")]
for (i in 1:ncol(ped))
  ped[,i] <- as.character(ped[,i])
names(ped) <- c("family", "individ", "father", "mother", "sex")
(chk <- pedigreeCheck(ped))
dups <- chk$duplicates
uni.ped <- pedigreeDeleteDuplicates(ped, dups)
(chk <- pedigreeCheck(uni.ped))

ni <- chk$parent.no.individ.entry
temp.ped <- lapply(unique(ni$family),function(x) {
  data.frame(family=c(x[1],x[1]), individ=c(paste(x[1],"4",sep=""),paste(x[1],"5",sep="")), father=c(0,0),mother=c(0,0),sex=c("F","M"),stringsAsFactors = F)
})
temp.ped <- as.data.frame(rbindlist(temp.ped))
ped.complete <- rbind(uni.ped, temp.ped)
(chk <- pedigreeCheck(ped.complete))

# Calculate expected relatedness
rels <- pedigreePairwiseRelatedness(ped.complete)
length(rels$inbred.fam)
relprs <- rels$relativeprs
table(relprs$relation)

# Calculate relatedness using PC-AiR and PC-Relate. Need to determine how many PC's to use with PC-Relate, and how many iterations of PC-AiR and PC-Relate to use
gds.fn <- "genotypes.gds"
(gds <- GdsGenotypeReader(gds.fn))
genoData <- GenotypeData(gds,snpAnnot=snpAnnot,scanAnnot = scanAnnot)

kinmat <- ibdobj$kinship
row.names(kinmat) <- ibdobj$sample.id
colnames(kinmat) <- ibdobj$sample.id
mypcair <- pcair(genoData = genoData,snp.include = snp.pruned,
                 kinMat = kinmat, divMat = kinmat)

mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:2], 
                       snp.include = snp.pruned, training.set = mypcair$unrels)
pc.iter2 <- as.list(rep(NA,4))
pc.iter2[[1]] <- list(mypcair,mypcrelate)

for (i in 2:length(pc.iter2)) {
  mypcair <- pcair(genoData = genoData,snp.include = snp.pruned,
                   kinMat = pc.iter2[[i-1]][[2]]$kinship, divMat= pc.iter2[[i-1]][[2]]$kinship)
  
  mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:2], 
                         snp.include = snp.pruned,
                         training.set = mypcair$unrels)
  
  pc.iter2[[i]] <- list(mypcair,mypcrelate)
}

mypcair <- pcair(genoData = genoData,snp.include = snp.pruned,v = NULL,
                 kinMat = kinmat, divMat = kinmat)

mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:3], 
                       snp.include = snp.pruned, training.set = mypcair$unrels)
pc.iter3 <- as.list(rep(NA,4))
pc.iter3[[1]] <- list(mypcair,mypcrelate)

for (i in 2:length(pc.iter3)) {
  mypcair <- pcair(genoData = genoData,snp.include = snp.pruned,
                   kinMat = pc.iter3[[i-1]][[2]]$kinship, divMat= pc.iter3[[i-1]][[2]]$kinship)
  
  mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:3], 
                         snp.include = snp.pruned,
                         training.set = mypcair$unrels)
  
  pc.iter3[[i]] <- list(mypcair,mypcrelate)
}

mypcair <- pcair(genoData = genoData,snp.include = snp.pruned,
                 kinMat = kinmat, divMat = kinmat)

mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:4], 
                       snp.include = snp.pruned, training.set = mypcair$unrels)
pc.iter4 <- as.list(rep(NA,4))
pc.iter4[[1]] <- list(mypcair,mypcrelate)

for (i in 2:length(pc.iter4)) {
  mypcair <- pcair(genoData = genoData,snp.include = snp.pruned,
                   kinMat = pc.iter4[[i-1]][[2]]$kinship, divMat= pc.iter4[[i-1]][[2]]$kinship)
  
  mypcrelate <- pcrelate(genoData = genoData, pcMat = mypcair$vectors[,1:4], 
                         snp.include = snp.pruned,
                         training.set = mypcair$unrels)
  
  pc.iter4[[i]] <- list(mypcair,mypcrelate)
}

# Analyze relatedness when using 2 PCs
pcrelatekin <- pcrelateReadKinship(pcrelObj = pc.iter2[[3]][[1]])
ibd <- pcrelatekin[,c("ID1","ID2","k0","kin")]
samp <- pData(scanAnnot)[,c("scanID", "subjectID")]
samp <- samp[match(pc.iter2[[3]][[2]]$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")

examine.relatedness <- function() {
ibd <- merge(ibd, samp, by.x="ID1", by.y="scanID")
ibd <- merge(ibd, samp, by.x="ID2", by.y="scanID", suffixes=c("1","2"))
ibd$ii <- pasteSorted(ibd$Individ1, ibd$Individ2)
relprs$ii <- pasteSorted(relprs$Individ1, relprs$Individ2)
ibd <- merge(ibd, relprs[,c("ii","relation")], all.x=TRUE)
names(ibd)[names(ibd) == "relation"] <- "exp.rel"
ibd$exp.rel[ibd$Individ1 == ibd$Individ2] <- "Dup"
ibd$exp.rel[is.na(ibd$exp.rel)] <- "U"
ibd$obs.rel <- ibdAssignRelatednessKing(ibd$k0, ibd$kin)
return(ibd)
}
ibd <- examine.relatedness()
table(ibd$exp.rel, useNA="ifany")
table(ibd$obs.rel, useNA="ifany")

# Analyze relatedness when using 3 PCs
pcrelatekin <- pcrelateReadKinship(pcrelObj = pc.iter3[[3]][[1]])
ibd <- pcrelatekin[,c("ID1","ID2","k0","kin")]
samp <- pData(scanAnnot)[,c("scanID", "subjectID")]
samp <- samp[match(pc.iter3[[3]][[2]]$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")

ibd <- examine.relatedness()
table(ibd$exp.rel, useNA="ifany")
table(ibd$obs.rel, useNA="ifany")

# Analyze relatedness when using 4 PCs
pcrelatekin <- pcrelateReadKinship(pcrelObj = pc.iter4[[3]][[1]])
ibd <- pcrelatekin[,c("ID1","ID2","k0","kin")]
samp <- pData(scanAnnot)[,c("scanID", "subjectID")]
samp <- samp[match(pc.iter4[[3]][[2]]$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")

ibd <- examine.relatedness()
table(ibd$exp.rel, useNA="ifany")
table(ibd$obs.rel, useNA="ifany")

# Relatedness does not change from 3 to 4 PCs, therefore 3 PCs will be used.

# Determine how many iterations of PC-AiR and PC-Relate to use:
races <- names(table(scanAnnot$race1))
cols <- c("American Indian/Alaska Native/Eskimo"="blue",
          "Asian/East Indian"="red",
          "Black"="green",
          "Native Hawaiian/Pacific Islander"="magenta",
          "White/Middle Eastern"="black")
colorbyrace <- cols[scanAnnot$race1]

# First iteration:
plot(pc.iter3[[1]][[1]],main="PC 1 vs PC 2 of PC-AiR, second iteration",col=colorbyrace)
legend("topleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(pc.iter3[[1]][[1]],main="PC 3 vs PC 4 of PC-AiR, second iteration",col=colorbyrace, vx = 3, vy = 4)
legend("bottomleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(pc.iter3[[1]][[1]],main="PC 5 vs PC 6 of PC-AiR, second iteration",col=colorbyrace, vx = 5, vy = 6)
legend("topleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(x=1:length(pc.iter3[[1]][[1]]$values),y=pc.iter3[[1]][[1]]$values,xaxp=c(1,length(pc.iter3[[1]][[1]]$values),length(pc.iter3[[1]][[1]]$values)-1),
     main="Scree plot of Eigenvalues from PC-AiR, second iteration",
     xlab="Principal Component Number",
     ylab="Eigenvalue")

pcrelatekin <- pcrelateReadKinship(pcrelObj = pc.iter3[[1]][[2]])
ibd <- pcrelatekin[,c("ID1","ID2","k0","kin")]
samp <- pData(scanAnnot)[,c("scanID", "subjectID")]
samp <- samp[match(mypcrelate$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")

ibd <- examine.relatedness()
table(ibd$exp.rel, useNA="ifany")
table(ibd$obs.rel, useNA="ifany")

plot(ibd$k0, ibd$kin, col=cols[ibd$exp.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient",main="Sample Relatedness by kinship and IBS=0\ncolored by expected relatedness, plotting all comparisons")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=c("Duplicate","Unrelated","Full Siblings"), col=cols, pch=1)

cols <- c(Dup="magenta", U="black", FS="blue",Deg2="red",Deg3="orange")
plot(ibd$k0, ibd$kin, col=cols[ibd$obs.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient",main="Sample Relatedness using top 2 PCs from PC-AiR\n followed by PC-Relate, first iteration\n colored by observed relatedness, plotting all comparisons")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=c("Duplicate","Unrelated","Full Siblings","Degree 2","Degree 3"), col=cols, pch=1)

# Second iteration:
plot(pc.iter3[[2]][[1]],main="PC 1 vs PC 2 of PC-AiR, second iteration",col=colorbyrace)
legend("topleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(pc.iter3[[2]][[1]],main="PC 3 vs PC 4 of PC-AiR, second iteration",col=colorbyrace, vx = 3, vy = 4)
legend("bottomleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(pc.iter3[[2]][[1]],main="PC 5 vs PC 6 of PC-AiR, second iteration",col=colorbyrace, vx = 5, vy = 6)
legend("topleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(x=1:length(pc.iter3[[2]][[1]]$values),y=pc.iter3[[2]][[1]]$values,xaxp=c(1,length(pc.iter3[[2]][[1]]$values),length(pc.iter3[[2]][[1]]$values)-1),
     main="Scree plot of Eigenvalues from PC-AiR, second iteration",
     xlab="Principal Component Number",
     ylab="Eigenvalue")

pcrelatekin <- pcrelateReadKinship(pcrelObj = pc.iter[[2]][[2]])
ibd <- pcrelatekin[,c("ID1","ID2","k0","kin")]
samp <- pData(scanAnnot)[,c("scanID", "subjectID")]
samp <- samp[match(mypcrelate$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")

ibd <- examine.relatedness()
table(ibd$exp.rel, useNA="ifany")
table(ibd$obs.rel, useNA="ifany")

plot(ibd$k0, ibd$kin, col=cols[ibd$exp.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient",main="Sample Relatedness by kinship and IBS=0\ncolored by expected relatedness, plotting all comparisons")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=c("Duplicate","Unrelated","Full Siblings"), col=cols, pch=1)

cols <- c(Dup="magenta", U="black", FS="blue",Deg2="red",Deg3="orange")
plot(ibd$k0, ibd$kin, col=cols[ibd$obs.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient",main="Sample Relatedness using top 2 PCs from PC-AiR\n followed by PC-Relate, second iteration\n colored by observed relatedness, plotting all comparisons")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=c("Duplicate","Unrelated","Full Siblings","Degree 2","Degree 3"), col=cols, pch=1)

# Third iteration:
plot(pc.iter3[[3]][[1]],main="PC 1 vs PC 2 of PC-AiR, second iteration",col=colorbyrace)
legend("topleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(pc.iter3[[3]][[1]],main="PC 3 vs PC 4 of PC-AiR, second iteration",col=colorbyrace, vx = 3, vy = 4)
legend("bottomleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(pc.iter3[[3]][[1]],main="PC 5 vs PC 6 of PC-AiR, second iteration",col=colorbyrace, vx = 5, vy = 6)
legend("topleft",text.col=cols,legend=races,cex=.75,title="Self reported race",title.col = "black")
legend("topright",legend=c("Unrelated subset","Related subset"),pch=c(19,3),cex=.75,title="PC-AiR subsets")

plot(x=1:length(pc.iter3[[3]][[1]]$values),y=pc.iter3[[3]][[1]]$values,xaxp=c(1,length(pc.iter3[[3]][[1]]$values),length(pc.iter3[[3]][[1]]$values)-1),
     main="Scree plot of Eigenvalues from PC-AiR, second iteration",
     xlab="Principal Component Number",
     ylab="Eigenvalue")

pcrelatekin <- pcrelateReadKinship(pcrelObj = pc.iter3[[3]][[2]])
ibd <- pcrelatekin[,c("ID1","ID2","k0","kin")]
samp <- pData(scanAnnot)[,c("scanID", "subjectID")]
samp <- samp[match(mypcrelate$sample.id, samp$scanID),]
names(samp) <- c("scanID", "Individ")

ibd <- examine.relatedness()
table(ibd$exp.rel, useNA="ifany")
table(ibd$obs.rel, useNA="ifany")

plot(ibd$k0, ibd$kin, col=cols[ibd$exp.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient",main="Sample Relatedness by kinship and IBS=0\ncolored by expected relatedness, plotting all comparisons")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=c("Duplicate","Unrelated","Full Siblings"), col=cols, pch=1)

cols <- c(Dup="magenta", U="black", FS="blue",Deg2="red",Deg3="orange")
plot(ibd$k0, ibd$kin, col=cols[ibd$obs.rel],
     xlab="Fraction of IBS=0", ylab="Kinship coefficient",main="Sample Relatedness using top 2 PCs from PC-AiR\n followed by PC-Relate, first iteration\n colored by observed relatedness, plotting all comparisons")
abline(h=c(cut.deg1, cut.deg2, cut.deg3, cut.dup), lty=2, col="gray")
legend("topright", legend=c("Duplicate","Unrelated","Full Siblings","Degree 2","Degree 3"), col=cols, pch=1)

# Relatedness does not change after the second iteration, therefore use the results of the second iteration

# Combine the data with HapMap reference samples and calculate genomic principal components. HapMap genotypes have been filtered to match the alleles of the experiment's data set and only unambiguous SNPs:

snpgdsCombineGeno(gds.fn=c("genotypes.gds",
                           "HapMap_SNPs_to_combine.gds"),
                  out.fn="combined_hapmap_gwas_flip.gds")

combined.gds.fn <- "combined_hapmap_gwas_flip.gds"
(combined.gds <- GdsGenotypeReader(combined.gds.fn))
combined.genoData <- GenotypeData(combined.gds)

combined.scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=getScanID(combined.genoData),
                                                         sex=c(scanAnnot$sex,rep("F",nscan(combined.genoData)-nrow(scanAnnot)))))
close(combined.genoData)
(combined.gds <- GdsGenotypeReader(combined.gds.fn))
combined.genoData <- GenotypeData(combined.gds,scanAnnot = combined.scanAnnot)

miss <- missingGenotypeByScanChrom(combined.genoData)
hist(100*(1-miss$missing.fraction),xlab="Call Rate",main="Missing call rates of combined GWAS and HapMap data set")

# Some HapMap samples have a call rate that is too low to use, identify those samples and remove them from subsequence analysis.
good.call.rate.samps <- getScanID(combined.genoData)[miss$missing.fraction<0.4]

filt <- get(data(pcaSnpFilters.hg19))
chrom <- getChromosome(combined.genoData)
pos <- getPosition(combined.genoData)
snpID <- getSnpID(combined.genoData)
snp.filt <- rep(TRUE, length(snpID))
for (f in 1:nrow(filt)) {
  snp.filt[chrom == filt$chrom[f] & filt$start.base[f] < pos
           & pos < filt$end.base[f]] <- FALSE
}
snp.sel <- snpID[snp.filt]
length(snp.sel)

sample.sel <- good.call.rate.samps[good.call.rate.samps %in% getScanID(combined.genoData)[!duplicated.samps]]
length(sample.sel)

close(combined.genoData)
gdsfile <- combined.gds.fn
gdsobj <- snpgdsOpen(gdsfile)
snpset <- snpgdsLDpruning(gdsobj, sample.id=sample.sel,
                          snp.id=snp.sel[!(snp.sel %in% filteredsnps[,1])],
                          autosome.only=TRUE, maf=0.05, missing.rate=0.05,
                          method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(0.1))

snp.pruned <- unlist(snpset, use.names=FALSE)

ibdobj <- snpgdsIBDKING(gdsobj,num.thread=16,autosome.only = T,
                        snp.id=snp.pruned,sample.id = good.call.rate.samps,
                        type="KING-robust")
kinmat <- ibdobj$kinship
row.names(kinmat) <- ibdobj$sample.id
colnames(kinmat) <- ibdobj$sample.id

snpgdsClose(gdsobj)
(combined.gds <- GdsGenotypeReader(combined.gds.fn))
combined.genoData <- GenotypeData(combined.gds,scanAnnot = combined.scanAnnot)

# load HapMap ethnicity information
hapmap_race <- read.delim("hapmap_relationships.txt",
                          stringsAsFactors = F)
combined.scanAnnot$race=c(scanAnnot$race1,
                          hapmap_race$population[match(
                            combined.scanAnnot$scanID[790:nscan(combined.genoData)],hapmap_race$IID)])


mypcair.combined <- pcair(genoData = combined.genoData,v=NULL,scan.include = good.call.rate.samps,
                 snp.include = snp.pruned,kinMat = kinmat, divMat = kinmat)
cols <- c("American Indian/Alaska Native/Eskimo"="blue",
          "Asian/East Indian"="red",
          "Black"="green",
          "Native Hawaiian/Pacific Islander"="magenta",
          "White/Middle Eastern"="black",
          "ASW"="green4",
          "CEU"="gray48",
          "CHB"="pink",
          "CHD"="palevioletred",
          "GIH"="paleturquoise1",
          "JPT"="orange",
          "LWK"="mediumspringgreen",
          "MEX"="seashell3",
          "MKK"="palegreen",
          "TSI"="slategray1",
          "YRI"="turquoise")
races <- names(cols)
races.used.in.order <- sapply(good.call.rate.samps,function(x) {
  combined.scanAnnot$race[as.character(combined.scanAnnot$scanID)==as.character(x)]
})
colorbyrace <- cols[races.used.in.order]
# Make alpha values of the colors and point symbols be different between the experiment and the HapMap samples
colorbyrace <- c(alpha(colorbyrace[1:789],1),alpha(colorbyrace[789:length(colorbyrace)],.1))
pchval <- c(rep(9,789),rep(19,length(good.call.rate.samps)-789))

plot(mypcair.combined,main="PC 1 vs PC 2 of PC-AiR, using HapMap and ADHD GWAS samples",col=colorbyrace,pch=pchval)
legend("topright",text.col=cols,legend=races,cex=.6,title.col = "black",bty = "n")
legend("topleft",legend=c("ADHD GWAS samples","HapMap samples"),pch=c(9,19),cex=.75,bty = "n")

colorbyrace <- cols[races.used.in.order]
colorbyrace <- c(alpha(colorbyrace[1:789],1),alpha(colorbyrace[789:length(colorbyrace)],.2))
plot(mypcair.combined,main="PC 3 vs PC 4 of PC-AiR, using HapMap and ADHD GWAS samples",vx=3,vy=4 ,col=colorbyrace,pch=pchval)
legend("bottomright",text.col=cols,legend=races,cex=.6,title.col = "black",bty = "n")
legend("bottomleft",legend=c("ADHD GWAS samples","HapMap samples"),pch=c(9,19),cex=.75,bty = "n")

# Compare the PC's between the combined experiment and GWAS data set with the PC's using only the experiment data
plot(mypcair.combined$vectors[1:789,1],pc.iter3[[3]][[2]]$vectors[,1],main="Using PC-AiR: PC 1 of GWAS+HapMap data set vs PC 1 of GWAS only",xlab="PC 1 of GWAS+HapMap",ylab="PC 1 of GWAS only")
plot(mypcair.combined$vectors[1:789,2],pc.iter3[[3]][[2]]$vectors[,2],main="Using PC-AiR: PC 2 of GWAS+HapMap data set vs PC 2 of GWAS only",xlab="PC 2 of GWAS+HapMap",ylab="PC 2 of GWAS only")
plot(mypcair.combined$vectors[1:789,3],pc.iter3[[3]][[2]]$vectors[,3],main="Using PC-AiR: PC 3 of GWAS+HapMap data set vs PC 2 of GWAS only",xlab="PC 3 of GWAS+HapMap",ylab="PC 3 of GWAS only")

# Graph parallel coordinates plots for the experiment PC's
gds.fn <- "genotypes.gds"
(gds <- GdsGenotypeReader(gds.fn))
genoData <- GenotypeData(gds,snpAnnot=snpAnnot,scanAnnot = scanAnnot)

par.coord <- pc.iter3[[3]][[2]]$vectors
rangel <- apply(par.coord, 2, function(x) range(x)[1])
rangeh <- apply(par.coord, 2, function(x) range(x)[2])
std.coord <- par.coord
num.pcs.plot <- 10
for (i in 1:num.pcs.plot)
  std.coord[,i] <- (par.coord[,i] - rangel[i])/(rangeh[i]-rangel[i])
plot(c(0,num.pcs.plot+1), c(0,1), type = 'n', axes = FALSE, ylab = "", xlab = "",
     main = "Parallel Coordinates Plot of PC-AiR PC's on GWAS data set")
for (j in 1:(num.pcs.plot-1))
  for (i in sample(1:nrow(std.coord)) )
    lines(c(j,j+1), std.coord[i,c(j,j+1)], col=colorbyrace[i], lwd=0.25)
axis(1, at = 1:num.pcs.plot, labels = paste("PC",1:num.pcs.plot, sep = "."))
legend("topright",text.col=cols,legend=races,cex=.75,title.col = "black",bty = "n",y.intersp=.5)

# Calculate the correlations between SNP genotypes and the principal components

rv <- list(sample.id = row.names(pc.iter3[[3]][[2]]$vectors), snp.id = snp.pruned, 
           eigenval = pc.iter3[[3]][[2]]$values, eigenvect = pc.iter3[[3]][[2]]$vectors, varprop = pc.iter3[[3]][[2]]$values/pc.iter3[[3]][[2]]$sum.values, 
           TraceXTX = NULL, Bayesian = NULL, genmat = NULL)
class(rv) <- "snpgdsPCAClass"
close(genoData)
gdsobj <- snpgdsOpen(gds.fn)

num.eig <- 50
corr <- snpgdsPCACorr(rv, gdsobj, eig.which=1:num.eig,verbose = T,
                      snp.id = snpAnnot$snpID[!(snpAnnot$snpID %in% filteredsnps[,1])])
snp <- snpAnnot[match(corr$snp.id, snpID),]
chrom <- getChromosome(snp, char=TRUE)
pdf("DataCleaning-corr-cleanedsnps.pdf")
for (i in 1:num.eig) {
  snpCorrelationPlot(abs(corr$snpcorr[i,]), chrom,
                     main=paste("Eigenvector",i), ylim=c(0,1))
}
dev.off()
snpgdsClose(gdsobj)

# Find if there is a relationship between principal components and status
pca <-pc.iter3[[3]][[2]]

princomp <- as.data.frame(pca$vectors)
princomp$scanID <- as.factor(scanAnnot$scanID)
princomp$case.ctrl.status <- as.factor(scanAnnot$status)
princomp$race <- as.factor(scanAnnot$race1)

pc.percent <- 100 * pca$values/pca$sum.values
pc.percent
lbls <- paste("PC", 1:3, "\n", format(pc.percent[1:3], digits=2), "%", " variability", sep="")
cols <- rep(NA, nrow(scanAnnot))
cols[scanAnnot$status == 1] <- "green"
cols[scanAnnot$status == 0] <- "magenta"

pairs(pca$vectors[,1:3], col=cols, labels=lbls,
      main = "First Three PCs by PC-AiR Case = Green, Control = Purple")
boxplot(princomp[, 1] ~ princomp$case.ctrl.status,
        ylab = "PC1", main = "PC-AiR PC1 vs. Case-control Status")
boxplot(princomp[, 2] ~ princomp$case.ctrl.status,
        ylab = "PC2", main = "PC-AiR PC2 vs. Case-control Status")
boxplot(princomp[, 3] ~ princomp$case.ctrl.status,
        ylab = "PC3", main = "PC-AiR PC3 vs. Case-control Status")

lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", " variability", sep="")
cols <- c("American Indian/Alaska Native/Eskimo"="blue",
          "Asian/East Indian"="red",
          "Black"="green",
          "Native Hawaiian/Pacific Islander"="magenta",
          "White/Middle Eastern"="black",
          "ASW"="green4",
          "CEU"="gray48",
          "CHB"="pink",
          "CHD"="palevioletred",
          "GIH"="paleturquoise1",
          "JPT"="orange",
          "LWK"="mediumspringgreen",
          "MEX"="seashell3",
          "MKK"="palegreen",
          "TSI"="slategray1",
          "YRI"="turquoise")
races <- names(cols)
colorbyrace <- alpha(cols[scanAnnot$race1],1)
pairs(pca$vectors[,1:4], col=colorbyrace, labels=lbls,
      main = "First Four PCs by PC-AiR, colored by race")
boxplot(princomp[, 1] ~ princomp$case.ctrl.status,
        ylab = "PC1", main = "PC-AiR PC1 vs. Case-control Status")
boxplot(princomp[, 2] ~ princomp$case.ctrl.status,
        ylab = "PC2", main = "PC-AiR PC2 vs. Case-control Status")
boxplot(princomp[, 3] ~ princomp$case.ctrl.status,
        ylab = "PC3", main = "PC-AiR PC3 vs. Case-control Status")
boxplot(princomp[, 4] ~ princomp$case.ctrl.status,
        ylab = "PC4", main = "PC-AiR PC4 vs. Case-control Status")

aov.p1 <- aov(princomp[,1] ~ princomp$race *
                princomp$case.ctrl.status, princomp)
summary(aov.p1)
aov.p2 <- aov(princomp[,2] ~ princomp$race *
                princomp$case.ctrl.status, princomp)
summary(aov.p2)
aov.p3 <- aov(princomp[,3] ~ princomp$race *
                princomp$case.ctrl.status, princomp)
summary(aov.p3)
aov.p4 <- aov(princomp[,4] ~ princomp$race *
                princomp$case.ctrl.status, princomp)
summary(aov.p4)

aov.all <- pbapply(princomp[,1:length(pca$values)],2, function(x) {
  summary(aov(x~princomp$race * princomp$case.ctrl.status,princomp))
})
pvals<-sapply(aov.all,function(x){x[[1]][2,5]})
plot(1:length(pca$values),pvals,xlab="Principal Component Number",
     ylab="P-value of ANOVA model: PC~status",main="P-values of the ANOVA model: PC~status for all
     PCs of PC-AiR on GWAS data set")

# Calculate missing call rate for cases and controls separately:

controls <- scanAnnot$scanID[scanAnnot$status == 0 | is.na(scanAnnot$status)]
miss<-missingGenotypeBySnpSex(genoData,scan.exclude = controls)
snpAnnot$missing.case <- miss$missing.fraction
cases <- scanAnnot$scanID[scanAnnot$status == 1 | is.na(scanAnnot$status)]
miss<-missingGenotypeBySnpSex(genoData,scan.exclude = cases)
snpAnnot$missing.controls <- miss$missing.fraction

differences <- 100*(1-snpAnnot$missing.case[!(snpAnnot$snpID %in%  filteredsnps[,1])]) -100*(1-snpAnnot$missing.controls[!(snpAnnot$snpID %in% filteredsnps[,1])])
sum(differences==0)
hist(differences[differences!=0],main="Case sample call rate minus Control samples call rate over cleaned SNPs\nzero difference SNPs (396,593) excluded",xlab="Call rate difference")

# remove all SNPs which differ in call rate between cases and controls by at least 2% (absolute difference in call rate)

temp <- snpAnnot$snpID[abs(snpAnnot$missing.case - snpAnnot$missing.controls)> 0.02 & !(snpAnnot$snpID %in% filteredsnps[,1])]
tempnames <- sapply(temp,function(x) {snpAnnot$snpName[snpAnnot$snpID==x]})

filteredsnps <- rbind(filteredsnps,data.frame(snpID=temp,snpName=tempnames,
                                              reason="Call Rate Diff in Status 0.02"))

# Perform chromosomal anomaly detection using BAF and circular binary segmentation
# Note, the final parameter set was determined by parameter exploration

chrom <- getChromosome(snpAnnot, char=TRUE)
pos <- getPosition(snpAnnot)
hla.df <- get(data(HLA.hg19))
hla <- chrom == "6" & pos >= hla.df$start.base & pos <= hla.df$end.base
xtr.df <- get(data(pseudoautosomal.hg19))
xtr <- chrom == "X" & pos >= xtr.df["X.XTR", "start.base"] &
  pos <= xtr.df["X.XTR", "end.base"]
centromeres <- get(data(centromeres.hg19))
gap <- rep(FALSE, nrow(snpAnnot))
for (i in 1:nrow(centromeres)) {
  ingap <- chrom == centromeres$chrom[i] & pos > centromeres$left.base[i] &
    pos < centromeres$right.base[i]
  gap <- gap | ingap
}
ignore <- snpAnnot$missing.n1 == 1 #ignore includes intensity-only and failed snps
snp.exclude <- ignore | hla | xtr | gap
snp.ok <- snpAnnot$snpID[!snp.exclude]
snp.ok <- snp.ok[!(snp.ok %in% filteredsnps[,1])]
chromsnps <- snpAnnot$snpID[snpAnnot$chromosome %in% 1:23]
snp.ok <- snp.ok[snp.ok %in% chromsnps]
scan.ids <- scanAnnot$scanID

baf.seg <- anomSegmentBAF(qualData, genoData, scan.ids=scan.ids,chrom.ids=1:24, 
                          snp.ids=snp.ok,alpha=0.01,min.width = 2, verbose=TRUE)

dim(baf.seg)
baf.anom <- anomFilterBAF(qualData, genoData, segments=baf.seg,num.mark.thresh = 10,
                          snp.ids=snp.ok, centromere=centromeres, low.qual.ids=low.qual.ids,
                          verbose=T)
names(baf.anom)
baf.filt <- baf.anom$filtered
dim(baf.filt)
baf.filt$method <- "BAF"
anoms <- baf.filt
anoms$anom.id <- 1:nrow(anoms)
snp.not.ok <- snpAnnot$snpID[!(snpAnnot$snpID %in% snp.ok)]

stats <- anomSegStats(qualData, genoData, snp.ids=snp.ok, anom=anoms,
                      centromere=centromeres)

a<-table(table(baf.filt$scanID))
a<-c(nscan(genoData)-length(unique(baf.filt$scanID)),a)
names(a)[1]<-"0"
barplot(a,main="Number of Samples with a certain number of anomalies",xlab="Number of Anomalies in a Sample",ylab="Number of Samples")
barplot(table(stats$chromosome),main="Number of Anomalies identified by BAF per Chromosome",xlab="Chromosome",ylab="Number of Anomalies")


# Perform chromosomal anomaly detection using circular binary segmentation and LRR
# Parameter exploration is not shown

loh.anom <- anomDetectLOH(qualData, genoData, scan.ids=scan.ids,
                          chrom.ids=1:23, snp.ids=snp.ok,known.anoms=baf.filt,
                          small.num=20,verbose=T)

loh.filt <- loh.anom$filtered
loh.filt$method <- "LOH"
anoms <- loh.filt
anoms$anom.id <- 1:nrow(anoms)
snp.not.ok <- snpAnnot$snpID[!(snpAnnot$snpID %in% snp.ok)]

stats <- anomSegStats(qualData, genoData, snp.ids=snp.ok, anom=anoms,
                      centromere=centromeres)

# The LRR method of anomaly detection fails to remove to centromere spanning segments. Remove these segments then combine these anomalies with the ones found using BAF.

not.centromere <- apply(loh.filt, 1, function(x) {
  left.base <- as.integer(x[17])
  right.base <- as.integer(x[18])
  chr <- trimws(x[2])
  
  cent.left <- centromeres$left.base[centromeres$chrom==chr]
  cent.right <- centromeres$right.base[centromeres$chrom==chr]
  
  ifelse(abs(left.base-cent.left) < abs(left.base-cent.right) & abs(right.base-cent.right) < abs(right.base-cent.left),F,T)
})

loh.filt <- loh.filt[not.centromere,]

if (!is.null(loh.filt)) {
  loh.filt$method <- "LOH"
  cols <- intersect(names(baf.filt), names(loh.filt))
  anoms <- rbind(baf.filt[,cols], loh.filt[,cols])
} else {
  anoms <- baf.filt
}
anoms$anom.id <- 1:nrow(anoms)
anoms$color <- ifelse(anoms$method=="BAF","red","blue")

stats <- anomSegStats(qualData, genoData, snp.ids=snp.ok, anom=anoms,
                      centromere=centromeres)
anom.filt <- stats[,c("scanID", "chromosome", "left.base", "right.base")]
anom.filt$whole.chrom <- FALSE

# Create data set with the anomalous regions zeroed out
# select subjects
subj <- scanAnnot$scanID
subj.filt.file <- "genotypes_anomfilt.gds"
genofile <- "genotypes.gds"
close(genoData)
setMissingGenotypes(genofile, subj.filt.file, anom.filt,
                    file.type="gds", sample.include=subj, verbose=T)

(gds <- GdsGenotypeReader(subj.filt.file))
genoData <- GenotypeData(gds,snpAnnot=snpAnnot,scanAnnot = scanAnnot)

# Find SNPs that were discordantly called between duplicate samples

snp.excl <- filteredsnps[,1]
length(snp.excl)
dupdisc <- duplicateDiscordance(genoData, subjName.col="subjectID",
                                one.pair.per.subj=F,snp.exclude=snp.excl)
names(dupdisc)
head(dupdisc$discordance.by.snp)
length(dupdisc$discordance.by.subject)

disc.snps <- dupdisc$discordance.by.snp
names(disc.snps)[1] <- "snpID"
temp <- data.frame(snpID=getSnpID(snpAnnot),snpName=snpAnnot$snpName,
                   chromosome=getChromosome(snpAnnot),plus_allele_A=getAlleleA(snpAnnot),
                   plus_allele_B=getAlleleB(snpAnnot))
disc.snps <- merge(disc.snps,temp,by="snpID",all.x=T)
theme <- ttheme_default(core=list(fg_params=list(cex=1.0)))
graphics.off()
grid.table(disc.snps,theme=theme,rows=rep("",nrow(disc.snps)))

# Add discordant SNPs in duplicate samples to list of SNPs to remove

temp <- data.frame(snpID=disc.snps$snpID,snpName=disc.snps$snpName,reason="discordant")
filteredsnps <- rbind(filteredsnps,temp)

# Check the HWE for each SNP on autosomes and the X chromosome:
# Only use unrelated samples of the same self-reported race
scan.excl <- scanAnnot$scanID[scanAnnot$race1 != "White/Middle Eastern" | 
                                scanAnnot$duplicated | scanAnnot$sibling]

chr <- getChromosome(genoData)
auto <- range(which(chr %in% 1:22))
X <- range(which(chr == 23))

hwe <- exactHWE(genoData, scan.exclude=scan.excl, snpStart=auto[1], snpEnd=auto[2],verbose = T)
hweX <- exactHWE(genoData, scan.exclude=scan.excl, snpStart=X[1], snpEnd=X[2],verbose = T)
hwe <- rbind(hwe, hweX)
hwe <- hwe[!(hwe$snpID %in% filteredsnps[,1]),]
hwe.0 <- hwe[hwe$MAF > 0,]

# Add SNPs with HWE p-value < 1e-06 to filtered list
temp <- data.frame(snpID=hwe.0$snpID[hwe.0$pval<1e-06],
                   snpName=snpAnnot$snpName[match(hwe.0$snpID[hwe.0$pval<1e-06],snpAnnot$snpID)],
                   reason="HWE p-val < 1e-06")
filteredsnps <- rbind(filteredsnps,temp)

# Create data set with only cleaned SNPs
close(genoData)
gdsSubset(subj.filt.file,"/data/Peter/projects/ADHD/data/ADHD_GWAS_genotypes_new_manifest_cleaned.gds",
          snp.include=snpAnnot$snpID[!(snpAnnot$snpID %in% filteredsnps[,1])])

snpAnnot.c<-snpAnnot[!(snpAnnot$snpID %in% filteredsnps[,1])]
write.table(getAnnotation(snpAnnot.c),"/data/Peter/projects/ADHD/data/snpAnnotation_new_manifest_cleaned.txt",quote = F,sep="\t",row.names = F)
(gds <- GdsGenotypeReader("/data/Peter/projects/ADHD/data/ADHD_GWAS_genotypes_new_manifest_cleaned.gds"))
genoData.c <- GenotypeData(gds,snpAnnot=snpAnnot.c,scanAnnot = scanAnnot)

qxyfile <- "/data/Peter/projects/ADHD/data/ADHD_GWAS_intensity_new_manifest_cleaned_dat.gds"
qualGDS <- GdsIntensityReader(qxyfile)
qualData.c <- IntensityData(qualGDS, scanAnnot=scanAnnot,snpAnnot = snpAnnot.c)






