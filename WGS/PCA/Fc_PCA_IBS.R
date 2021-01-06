# This version is using subsets of T0 as replicates
#
# We are looking for outliers, mislabelled samples etc.
# Do PCA, SNP correlation and Identity-by-state analysis.

setwd("D://BCF_VAR_CALLS/noMisMatchFilteringAllSamples/vcfFile/FC_subT0")


# Load the R packages: gdsfmt and SNPRelate
library("gdsfmt")
library("SNPRelate")
library("scales")
library("dendextend")
if (!require("RColorBrewer")){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}


# Set up colour ranges
blueScale <- brewer.pal(9,"Blues")
redScale <- brewer.pal(9,"Reds")

cols <- c(rep(blueScale[2],3),
          rep(blueScale[4],3),
          rep(redScale[4],3),
          blueScale[8],
          redScale[8])
  

repNumbers <- c("1","2","3",
                "1","3","5",
                "1","3","5",
                "1",
                "1")


# Read in vcf file
#vcf.fn <- "FC_mergedFiltered_NoMismatchFilter_allBatches_SNPsOnly.vcf"
vcf.fn <- "FC_mergedFiltered_noMismatchFilter_allBatches_subT0_SNPsOnly.vcf"
snpgdsVCF2GDS(vcf.fn, "vcfNoMismatchFilterAllBatches.gds", method="bi")

# Summary
snpgdsSummary("vcfNoMismatchFilterAllBatches.gds")
#The file name: E:\BCF_VAR_CALLS\noMisMatchFilteringAllSamples\vcfFile\FC_subT0\vcfNoMismatchFilterAllBatches.gds 
#The total number of samples: 11 
#The total number of SNPs: 77303 
#SNP genotypes are stored in SNP-major mode (Sample X SNP).

# Open the GDS file
genofile <- snpgdsOpen("vcfNoMismatchFilterAllBatches.gds")

# Prune SNPs by missing rate and MAF (just for PCA)
snpset <- snpgdsSelectSNP(genofile,
                          autosome.only=FALSE,
                          missing.rate=0.3,
                          maf=0.05)

#Excluding 2,708 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: 0.3)

# Get all selected snp id
snpset.id <- unlist(snpset)


pcaFiltered <- snpgdsPCA(genofile, 
                         snp.id=snpset.id,
                         missing.rate=0.3,
                         maf=0.05,
                         autosome.only=FALSE)

# Plot filtered PCA
pdf("PCA/PCA.noOutliers.noLD.MAF005_MISS03.pdf")
  plot(pcaFiltered,
       col=alpha(cols,0.5),
       pch=16,
       cex=2.5)
  text(pcaFiltered$eigenvect[,1],pcaFiltered$eigenvect[,2], 
       repNumbers,
       pos=3,cex=0.6,
       col="black")
dev.off()



#--- SNP correlation

# Replace character based chromosome names with numbers 1:27
# to use as colours
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

newvals <- c(1:125)
foo <- newvals[as.factor(chr)]

# Get correlations between SNPs and Eigenvectors. Which SNPs correlate with which axis?
# Use the set of snp.ids used in the PCA
CORR <- snpgdsPCACorr(pcaFiltered, 
                      genofile,
                      snp.id=pcaFiltered$snp.id,
                      eig.which=1:3)

pdf("PCA/PCA.correlation.noOutlier.noLD.MAF005.MISS03.pdf",width=12,height=8)
par(mfrow=c(3,1), mai=c(0.3, 0.55, 0.1, 0.25))
for (i in 1:3)
{
  plot(abs(CORR$snpcorr[i,]), 
       ylim=c(0,1), 
       xlab="", 
       ylab=paste("PC", i),
       col=foo, 
       pch=1)
}
dev.off()



#--- Get details of SNPs correlated with Principal components

# Get list of SNP details for all SNPs
SNPlist <- snpgdsSNPList(genofile)

# Get SNP loadings - Remember: some SNPs not used (210177 out of 242884)
snpLoads <- snpgdsPCASNPLoading(pcaFiltered,genofile)

# Filter SNPlist to include just those used in PCA
PCA_SNPS <- SNPlist[SNPlist$snp.id %in% snpLoads$snp.id,]


# Add correlations for PCs 1 & 2
PCA_SNPS$corr1 <- CORR$snpcorr[1,]

# Add SNP loadings for PCs 1 & 2
PCA_SNPS$pc1Load <- snpLoads$snploading[1,]

# Get SNPs with correlation >= 0.9
snpCorr1 <- na.omit(PCA_SNPS[PCA_SNPS$corr1>=0.9,])

write.table(snpCorr1, "PCA/PC1correlatedSNPs.txt",
            row.names=FALSE,quote=FALSE)

snpNegCorr1 <- na.omit(PCA_SNPS[PCA_SNPS$corr1<=-0.9,])

write.table(snpNegCorr1, "PCA/PC1negcorrelatedSNPs.txt",
            row.names=FALSE,quote=FALSE)



#--- Identity by state analysis

set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile,
                                   autosome.only=FALSE,
                                   maf=0.05,
                                   missing.rate=0.3,
                                   snp.id=snpset.id))

rv <- snpgdsCutTree(ibs.hc)

# Relabel and recolor leaves
tree <- rv$dendrogram

# Get current sample order
tree %>% labels

# Using subsets of T0
#
#[1] "T23-11-17_8C_FC1" "T28-11-17_4C_FC1" "T47_8C_FC3"       "T0_3pt5C_FC1_S2"  "T0_3pt5C_FC1_S1"  "T0_3pt5C_FC1_S3" 
#[7] "T47_8C_FC1"       "T47_8C_FC5"       "T66_3pt5C_FC3"    "T66_3pt5C_FC1"    "T66_3pt5C_FC5"   

newLabels <- c("T23-11-17_8C_FC1","T28-11-17_4C_FC1","T47_8C_FC3",
               "T0_3pt5C_FC1_S2","T0_3pt5C_FC1_S1","T0_3pt5C_FC1_S3",
               "T47_8C_FC1","T47_8C_FC5","T66_3.5C_FC3",
               "T66_3.5C_FC1","T66_3.5C_FC5")


newCols <- c(redScale[8],blueScale[8],redScale[4],
             blueScale[2],blueScale[2],blueScale[2],
             redScale[4],redScale[4],blueScale[4],
             blueScale[4],blueScale[4])


pdf("PCA/IBS.dendrogram.noOutlier.noLD.MAF005_MISS03.pdf")
tree %>% set("labels", newLabels) %>% 
  set("labels_cex",0.7) %>% 
  set("leaves_cex", 2) %>%
  set("leaves_col", newCols) %>% plot(main="Identity by state")
dev.off()



# Close the GDS file
snpgdsClose(genofile)
