# This version is using the complete T. pseudonana dataset. The latest batch
# included two more replicates for T0.
#
# We are looking for outliers, mislabelled samples etc.
# Do PCA, SNP correlation and Identity-by-state analysis.

setwd("D://BCF_VAR_CALLS/noMisMatchFilteringAllSamples/vcfFile/TP")


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
# Create color palettes
# Green for 22C, blue for 9C and red for 32C
greenScale <- brewer.pal(9,"Greens")
blueScale <- brewer.pal(9,"Blues")
redScale <- brewer.pal(9,"Reds")

# Commented out all 5 T450 replicates, 2 appear to be
# outliers, likely mislabelled. Just use replicates
# 2,3 and 5 for this group
cols <- c(rep(greenScale[2],3),
          redScale[6],
          rep(greenScale[4],5),
          rep(redScale[2],5),
          rep(blueScale[2],5),
          rep(greenScale[6],4),
          rep(redScale[4],5),
          rep(blueScale[4],5),
          rep(blueScale[6],5),
          #rep(redScale[8],5))
          rep(redScale[8],3))


sampleNamesMinusWeird <- c("1568","PRO2022_S1","PRO2022_S2",
                           "1568_2",
                           "SAM37194","SAM37195","SAM37196","SAM37197","SAM37198",
                           "SAM37199","SAM37200","SAM37201","SAM37202","SAM37203",
                           "SAM37204","SAM37205","SAM37206","SAM37207","SAM37208",
                           "SAM37209","SAM37210","SAM37211","SAM37212",
                           "SAM37213","SAM37214","SAM37215","SAM37216","SAM37217",
                           "SAM37218","SAM37219","SAM37220","SAM37221","SAM37222",
                           "SAM37223","SAM37224","SAM37225","SAM37226","SAM37227",
                           #"SAM37228","SAM37229","SAM37230","SAM37231","SAM37232")
                           "SAM37229","SAM37230","SAM37232")


repNumbers <- c("1","2","3",
                "1",
                "1","2","3","4","5",
                "1","2","3","4","5",
                "1","2","3","4","5",
                "2","3","4","5",
                "1","2","3","4","5",
                "1","2","3","4","5",
                "1","2","3","4","5",
                #"1","2","3","4","5")
                "2","3","5")


# Read in vcf file
vcf.fn <- "mergedFiltered_NoMismatchFilter_allBatches_SNPsOnly.vcf"
snpgdsVCF2GDS(vcf.fn, "vcfNoMismatchFilterAllBatches.gds", method="bi")

# Summary
snpgdsSummary("vcfNoMismatchFilterAllBatches.gds")
# The file name: \\ueahome3\stfsci3\bva09npu\data\Documents\BCF_VAR_CALLS\noMisMatchFilteringAllSamples\vcfFile\vcfNoMismatchFilterAllBatches.gds 
# The total number of samples: 43 
# The total number of SNPs: 183473 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).

# Open the GDS file
genofile <- snpgdsOpen("vcfNoMismatchFilterAllBatches.gds")

# Prune SNPs by missing rate and MAF (just for PCA)
snpset <- snpgdsSelectSNP(genofile,
                          autosome.only=FALSE,
                          sample.id=sampleNamesMinusWeird,
                          missing.rate=0.3,
                          maf=0.05)

# Get all selected snp id
snpset.id <- unlist(snpset)


pcaFiltered <- snpgdsPCA(genofile, 
                         snp.id=snpset.id,
                         sample.id=sampleNamesMinusWeird,
                         missing.rate=0.3,
                         maf=0.05,
                         autosome.only=FALSE)

# Plot filtered PCA
pdf("PCA/PCA.noOutliers.noLD.MAF005_MISS03.pdf")
  plot(pcaFiltered,
       col=alpha(cols,0.5),
       pch=16,
       cex=3.5)
  text(pcaFiltered$eigenvect[,1],pcaFiltered$eigenvect[,2], 
       repNumbers,
       pos=3,cex=0.6,
       col="black")
dev.off()



#--- SNP correlation

# Replace character based chromosome names with numbers 1:27
# to use as colours
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

newvals <- c(1:27)
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
                                   maf=0.1,
                                   missing.rate=0.3,
                                   snp.id=snpset.id,
                                   sample.id=sampleNamesMinusWeird))
rv <- snpgdsCutTree(ibs.hc)

# Relabel and recolor leaves
tree <- rv$dendrogram

# Get current sample order
tree %>% labels

#[1] "PRO2022_S1" "SAM37208"   "SAM37204"   "SAM37198"   "SAM37195"   "SAM37197"   "SAM37194"   "PRO2022_S2"
#[9] "SAM37196"   "SAM37218"   "SAM37221"   "1568"       "SAM37211"   "SAM37222"   "SAM37210"   "SAM37220"  
#[17] "SAM37212"   "SAM37223"   "SAM37209"   "SAM37226"   "SAM37227"   "SAM37224"   "SAM37219"   "SAM37225"  
#[25] "SAM37207"   "SAM37205"   "SAM37206"   "SAM37229"   "SAM37217"   "SAM37232"   "SAM37199"   "SAM37202"  
#[33] "SAM37203"   "SAM37201"   "1568_2"     "SAM37200"   "SAM37230"   "SAM37213"   "SAM37216"   "SAM37214"  
#[41] "SAM37215"


newLabels <- c("T0_22C_TP2","T32_9C_TP5","T32_9C_TP1","T70_22C_TP5","T70_22C_TP2",
               "T70_22C_TP4","T70_22C_TP1","T0_22C_TP3","T70_22C_TP3","T144_9C_TP1",
               "T144_9C_TP4","T0_22C_TP1","T270_22C_TP4","T144_9C_TP5","T270_22C_TP3",
               "T144_9C_TP3","T270_22C_TP5","T250_9C_TP1","T270_22C_TP2","T250_9C_TP4",
               "T250_9C_TP5","T250_9C_TP2","T144_9C_TP2","T250_9C_TP3","T32_9C_TP4",
               "T32_9C_TP2","T32_9C_TP3",
                 
               "T450_32C_TP2","T210_32C_TP5","T450_32C_TP5","T50_32C_TP1","T50_32C_TP4",
               "T50_32C_TP5","T50_32C_TP3","T300_32C_TP1","T50_32C_TP2","T450_32C_TP3",
               "T210_32C_TP1","T210_32C_TP4","T210_32C_TP2","T210_32C_TP3")


newCols <- c(greenScale[2],blueScale[2],blueScale[2],greenScale[4],greenScale[4],
             greenScale[4],greenScale[4],greenScale[2],greenScale[4],blueScale[4],
             blueScale[4],greenScale[2],greenScale[6],blueScale[4],greenScale[6],
             blueScale[4],greenScale[6],blueScale[6],greenScale[6],blueScale[6],
             blueScale[6],blueScale[6],blueScale[4],blueScale[6],blueScale[2],
             blueScale[2],blueScale[2],
             
             redScale[8],redScale[4],redScale[8],redScale[2],redScale[2],
             redScale[2],redScale[2],redScale[6],redScale[2],redScale[8],
             redScale[4],redScale[4],redScale[4],redScale[4])


pdf("PCA/IBS.dendrogram.noOutlier.noLD.MAF005_MISS03.pdf")
tree %>% set("labels", newLabels) %>% 
  set("labels_cex",0.7) %>% 
  set("leaves_cex", 2) %>%
  set("leaves_col", newCols) %>% plot(main="Identity by state")
dev.off()



# Close the GDS file
snpgdsClose(genofile)
