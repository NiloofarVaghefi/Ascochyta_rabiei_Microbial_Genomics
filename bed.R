Ar230_maf1_vcf

install.packages("bedr")
library("bedr")
vcf2bed(Ar230_maf1_vcf, filename = NULL, header = FALSE, other = NULL, verbose = TRUE)



vcf.file <- system.file("G:/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/SNPeff/snp.ann.vcf", package = "bedr")
install.packages("pegas")
library("pegas")
x <- read.vcf(vcf.file)
x <- read.vcf("G:/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/SNPeff/snp.ann.vcf")
Ar230_maf1.bed <- vcf2bed(x, header = TRUE, other = NULL)

Ar230_maf1.bed

write.table(Ar230_maf1.bed,file="Ar230_maf1.bed",sep="\t", row.names=FALSE)


 