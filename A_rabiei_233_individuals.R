## Import the filtered vcf dataset including 233 individuals (including the replicated samples)

library(vcfR)
library(poppr)

#for MAC
Ar_maf1_vcf <- read.vcfR("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics/Data/filtered_233/snp.ann.all_filtered.vcf")

#Windows
#Ar_maf1_vcf <- read.vcfR("G:/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics/Data/filtered_233/snp.ann.all_filtered.vcf")

Ar_genind_maf1 <- vcfR2genind(Ar_maf1_vcf)
Ar_genind_maf1
indNames(Ar_genind_maf1)

X <- genind2df(Ar_genind_maf1)
Arsnps_maf1 <- df2genind(X, ploidy = 1, ncode=1)
indNames(Arsnps_maf1)

Ar_maf1 <- poppr::as.genclone(Arsnps_maf1)
Ar_maf1
indNames(Ar_maf1)

## Check the genetic distance between replicated samples
#  this will calculate the pairwise euclidean distances between samples

Ar_maf1.dist_all <- bitwise.dist(Ar_maf1)
Ar_maf1.dist_all
library(MASS)
write.matrix (Ar_maf1.dist_all, file = "Ar233_distance.matrix.bitwise.txt")

# In the exported file you will see distance between three pairs of replicates is zero
# this shows minimum genotyping error and no need for further filtering
