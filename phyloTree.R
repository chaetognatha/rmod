#Phylogenetic trees with SNPRelate



#setting up the workspace
setwd("C:/Users/V/Documents/LundUni/BINP28/Project")
getwd()
#loading packages
#library("vcfR")
library("devtools")
library("gdsfmt")
library("SNPRelate")
library("ggplot2")
library("ggtree")
library("ape")


#reading in the VCF file, and transferring the data to the gds format
#vcf_file <- read.vcfR( "MysteryTaxa.recode.vcf", verbose = FALSE ) 
#had that worked, this would have been what the vcfR package was for
#file type conversion: 
#source: https://www.rdocumentation.org/packages/SNPRelate/versions/1.6.2/topics/snpgdsVCF2GDS
#vcf.fn <- system.file("extdata", "C:/Users/V/Documents/LundUni/BINP28/Project/MysteryTaxa.recode.vcf", package="SNPRelate")
#couldn't get the system.file() to read the vcf, so I just used the file in teh next line directly, instead of putting it in a vecor first
snpgdsVCF2GDS("ProjTaxa.vcf", "ProjTaxa.gds", method= "biallelic.only")
#Returns: 
# Start file conversion from VCF to SNP GDS ...
# Method: exacting biallelic SNPs
# Number of samples: 16
# Parsing "ProjTaxa.vcf" ...
# import 3241886 variants.
# + genotype   { Bit2 16x3241886, 12.4M } *
#     Optimize the access efficiency ...
# Clean up the fragments of GDS file:
#     open the file 'ProjTaxa.gds' (26.7M)
# # of fragments: 409
# save to 'ProjTaxa.gds.tmp'
# rename 'ProjTaxa.gds.tmp' (26.6M, reduced: 4.6K)
# # of fragments: 20


#following tutorial from: https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html#preparing-data
genofile <- snpgdsOpen("ProjTaxa.gds")
snpgdsSummary(genofile)
#Returns: 
# The file name: C:\Users\V\Documents\LundUni\BINP28\Project\ProjTaxa.gds 
# The total number of samples: 16 
# The total number of SNPs: 3241886 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).
pop_code <- scan("pop.txt", what=character())
#Returns: 
# Read 16 items
table(pop_code)
#Returns: 
# pop_code
# Group8N     GroupK0 GroupLesina  GroupNaxos 
# 5           5           5           1 


#LD-based SNP pruning
set.seed(1000)
# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
#Returns: 
# SNP pruning based on LD:
#     Excluding 1,370,870 SNPs on non-autosomes
# Excluding 12,059 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# # of samples: 16
# # of SNPs: 1,858,957
# using 1 thread
# sliding window: 500,000 basepairs, Inf SNPs
# |LD| threshold: 0.2
# method: composite
# Chromosome 5: 0.10%, 1,875/1,871,016
# 1,875 markers are selected in total.
names(snpset)
#Returns: 
# [1] "chr5"
snpset.id <- unlist(snpset)
snpset.id

#PCA 
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
#Returns: 
# Principal Component Analysis (PCA) on genotypes:
#     Excluding 3,240,011 SNPs (non-autosomes or non-selection)
# Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# # of samples: 16
# # of SNPs: 1,875
# using 2 threads
# # of principal components: 32
# PCA:    the sum of all selected genotypes (0,1,2) = 56016
# CPU capabilities: Double-Precision SSE2
# Mon Feb 15 15:14:36 2021    (internal increment: 28416)
# [==================================================] 100%, completed, 0s  
# Mon Feb 15 15:14:36 2021    Begin (eigenvalues and eigenvectors)
# Mon Feb 15 15:14:36 2021    Done.
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)
# Draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#if there is population information
# get sample id
sample.id <- scan("individuals.txt", what=character())
# assume the order of sample ids is the same as the population codes
head(cbind(sample.id, pop_code))
# Make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
#this works!


#create a phylogeny
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2))
#Returns: 
# Identity-By-State (IBS) analysis on genotypes:
#     Excluding 1,370,870 SNPs on non-autosomes
# Excluding 12,059 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
# # of samples: 16
# # of SNPs: 1,858,957
# using 2 threads
# IBS:    the sum of all selected genotypes (0,1,2) = 47698233
# Mon Feb 15 15:16:26 2021    (internal increment: 65536)
# [==================================================] 100%, completed, 0s  
# Mon Feb 15 15:16:26 2021    Done.

# Determine groups of individuals automatically
rv <- snpgdsCutTree(ibs.hc)
#Returns: 
# Determine groups by permutation (Z threshold: 15, outlier threshold: 5):
#     Create 1 groups.

plot(rv$dendrogram, leaflab="none", main="Phylogeny", ylim=c(-0.1, 0.4))

table(rv$samp.group)
#Returns: 
# G001 
# 16 

# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))
#Returns: 
# Create 4 groups.

sampleGroup <- as.factor(pop_code)

plot(rv2$dendrogram, ylab = "Degrees of Relatedness", xlab = "Species", leaflab="perpendicular", main="Phylogeny", ylim=c(-0.25, 0.4))
legend("bottomleft", legend=levels(sampleGroup), col=1:nlevels(sampleGroup), pch=19, ncol=3)
#this works! 


#creating a Newick tree
#source: https://rpubs.com/adel922/560260
#some cool graphics here, come back to it later
#required packages: ape, ggtree
treePlot2 <- rv$dendrogram

ggtree(as.phylo(as.hclust(treePlot2)), layout="circular",color='darkgreen', branch.length="branch.length") + 
    geom_tiplab(size=2.5, aes(angle=angle)) + 
    ggtitle("Phylogeny")

# converting dendrograms to class hclust and to new variables
hcProject <- as.hclust(rv$dendrogram)
# Making the hclust object into a phylo object in ape
thisProject <- as.phylo(hcProject) 
# Writing to a newick tree file
write.tree(phy=thisProject, file="ProjTaxa.newick") 
