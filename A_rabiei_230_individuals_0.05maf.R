## IMPORT DATA WITH REPLICATES REMOVED

library(vcfR)
library(poppr)

#Windows
#Ar230_maf1_vcf <- read.vcfR("G:/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/SNPeff/snp.ann.vcf")
#MAC
Ar230_maf5_vcf <- read.vcfR("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/SNPeff/snp.ann5.vcf")
Ar230_genind_maf5 <- vcfR2genind(Ar230_maf5_vcf)
Ar230_genind_maf5
indNames(Ar230_genind_maf5)

X <- genind2df(Ar230_genind_maf5)
Ar230snps_maf5 <- df2genind(X, ploidy = 1, ncode=1)
indNames(Ar230snps_maf5)

Ar230_maf5 <- poppr::as.genclone(Ar230snps_maf5)
Ar230_maf5
# 659 loci
indNames(Ar230_maf5)


## SET STRATA

#MAC
strata.df <- read.csv("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/230_strata.csv", head = FALSE, sep = ",")
strata(Ar230_maf1) <- strata.df
splitStrata(Ar230_maf1) <- ~Year/State/Location/Farm/Host/Zone/ICC3996_rating/Genesis090_rating/HatTrick_rating/Seamer_rating/Path1/Path2
Ar230_maf1
# This is a genclone object
# 230 original multilocus genotypes 
# 230 haploid individuals
# 2759 codominant loci
# 12 strata - Year, State, Location, ..., Seamer_rating, Path1, Path2
# 0 populations defined. 

# To view the table showing all strata
library(dplyr)
Allstrata <- strata(Ar230_maf1) %>% group_by(Year,State, Location, Farm, Host, Zone, ICC3996_rating, Genesis090_rating, HatTrick_rating, Seamer_rating, Path1, Path2) %>% summarize(Count = n())
View(Allstrata)
nameStrata(Ar230_maf1)


## DISTANCE BETWEEN INDIVIDUALS (PAIRWISE)

library("poppr")
library("pegas")
library("ape")
library("adegenet")
library("ade4")

# This is applied on allele frequency within individuals (genind object). 
# Function dist() from adegenet provides different options. 
# using the euclidean distance among vector of allele frequencies.
distgenEUCL <- dist(Ar230_genind_maf1, method = "euclidean", diag = FALSE, upper = FALSE, p = 1)
hist(distgenEUCL)

X <- genind2loci(Ar230_genind_maf1)
pairwise_genetic_distance <- dist.gene(X, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(pairwise_genetic_distance)
# The option pairwise.deletion = FALSE removes all loci with one missing value
# can see on the histogram that we get a maximum distance of around 400 loci
# but the majority of individuals are different at around 200-250 loci

#gives density on y axis
ratio_of_allelic_differences <- diss.dist(Ar230_genind_maf1, percent = TRUE, mat = FALSE)
hist(ratio_of_allelic_differences, freq=FALSE)
#gives density on y axis
Number_of_allelic_differences <- diss.dist(Ar230_genind_maf1, percent = FALSE, mat = FALSE)
hist(Number_of_allelic_differences, freq=FALSE)

#gives frequency on y axis
ratio_of_allelic_differences <- diss.dist(Ar230_genind_maf1, percent = TRUE, mat = FALSE)
hist(ratio_of_allelic_differences)
#gives frequency on y axis
Number_of_allelic_differences <- diss.dist(Ar230_genind_maf1, percent = FALSE, mat = FALSE)
hist(Number_of_allelic_differences)


## CHECK MISSING DATA

setPop(Ar230_maf1) <- ~Zone
Missing_data <- info_table(Ar230_maf1, type = "missing")
sum(Missing_data["Total", 1:2759] > 0)
barplot(Missing_data["Total", 1:2759], xlab = "Locus", ylab = "Proportion Missing", las=2)
# no locus has more than 0.1 missing data
summary(Ar230_maf1)
# your will see that heterozygosity observed is zero and missing data is only 0.35%


## IMPORT DATA WITH REPLICATES REMOVED - LOW AND MODIFIER LOCI FOR POP GEN ANALYSES

Ar230_LM_maf1_vcf <- read.vcfR("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/SNPeff/lowmodifierimpact.vcf")
Ar230_LM_genind_maf1 <- vcfR2genind(Ar230_LM_maf1_vcf)
Ar230_LM_genind_maf1
indNames(Ar230_LM_genind_maf1)

X <- genind2df(Ar230_LM_genind_maf1)
Ar230snps_LM_maf1 <- df2genind(X, ploidy = 1, ncode=1)
indNames(Ar230snps_LM_maf1)

Ar230_LM_maf1 <- poppr::as.genclone(Ar230snps_LM_maf1)
Ar230_LM_maf1
indNames(Ar230_LM_maf1)


##   SET STRATA
strata.df <- read.csv("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/230_strata.csv", head = FALSE, sep = ",")
strata(Ar230_LM_maf1) <- strata.df
splitStrata(Ar230_LM_maf1) <- ~Year/State/Location/Farm/Host/Zone/ICC3996_rating/Genesis090_rating/HatTrick_rating/Seamer_rating/Path1/Path2
Ar230_LM_maf1

#To view the table showing all strata
library(dplyr)
Allstrata <- strata(Ar230_LM_maf1) %>% group_by(Year,State, Location, Farm, Host, Zone, ICC3996_rating, Genesis090_rating, HatTrick_rating, Seamer_rating, Path1, Path2) %>% summarize(Count = n())
View(Allstrata)
nameStrata(Ar230_LM_maf1)

#define populations based on strata
popNames(Ar230_LM_maf1)
setPop(Ar230_LM_maf1) <- ~Year
Ar230_LM_maf1


##  AMOVA 
setPop(Ar230_LM_maf1) <- ~Year
popNames(Ar230_LM_maf1)
# make a subpob without 2013 and 2014 as these are too small
Ar230_LM_maf1_year <- popsub (Ar230_LM_maf1, c(2,3,4,5))
indNames(Ar230_LM_maf1_year)

(Ar230_LM_maf1_amova <- poppr.amova(Ar230_LM_maf1_year, ~Year/State/Location, filter = TRUE, threshold = 0.00))
(Ar230_LM_maf1_amova_Pval   <- randtest(Ar230_LM_maf1_amova, nrepet = 999))
plot(Ar230_LM_maf1_amova_Pval)

(Ar230_LM_maf1_amova <- poppr.amova(Ar230_LM_maf1_year, ~Year/Zone/Location, filter = TRUE, threshold = 0.00))
(Ar230_LM_maf1_amova_Pval   <- randtest(Ar230_LM_maf1_amova, nrepet = 999))
plot(Ar230_LM_maf1_amova_Pval)

#without removing 2013 and 2014
(Ar230_LM_maf1_amova <- poppr.amova(Ar230_LM_maf1, ~Year/State/Location, filter = TRUE, threshold = 0.00))
(Ar230_LM_maf1_amova_Pval   <- randtest(Ar230_LM_maf1_amova, nrepet = 999))
plot(Ar230_LM_maf1_amova_Pval)

(Ar230_LM_maf1_amova <- poppr.amova(Ar230_LM_maf1, ~Year/Zone/Location, filter = TRUE, threshold = 0.00))
(Ar230_LM_maf1_amova_Pval   <- randtest(Ar230_LM_maf1_amova, nrepet = 999))
plot(Ar230_LM_maf1_amova_Pval)


##  Discriminant Analysis of Principal components
# need to define the no. of PCs retained for DAPC, which can have a substantial impact 
# on the results. Cross-validation (xvalDapc) is an objective way to identify no. of PCs 
# to retain.

# DAPC with populations pre-defined based on the year of collection
setPop(Ar230_LM_maf1) <- ~Year
popNames(Ar230_LM_maf1)
cols3 <- c("darkorange1", "purple", "light blue","cyan", "blue", "Yellow2")

set.seed(999)
xvalAr2021 <- xvalDapc(tab(Ar230_LM_maf1, NA.method = "mean"), pop(Ar230_LM_maf1), training.set = 0.7)
set.seed(999)
xvalAr2021 <- xvalDapc(tab(Ar230_LM_maf1, NA.method = "mean"), pop(Ar230_LM_maf1), parallel = "multicore", ncpus = 4L, training.set = 0.9, n.rep = 1000, n.pca = 70:90)
xvalAr2021 [2:6]
#80
popNames(Ar230_LM_maf1)

dapc_year <- dapc(Ar230_LM_maf1, var.contrib = TRUE, scale = FALSE, n.pca = 80, n.da = nPop(Ar230_LM_maf1) - 1)
#final:
scatter(dapc_year, col = cols3, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1,cell=0, cstar=0, legend = TRUE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
contrib_Year_axis1 <- loadingplot(dapc_year$var.contr, axis = 1, thres = 0.005)
contrib_Year_axis1

contrib_Year_axis2 <- loadingplot(dapc_year$var.contr, axis = 2, thres = 0.005)
contrib_Year_axis2

# DAPC with populations pre-defined based on AE Zone
setPop(Ar230_LM_maf1) <- ~Zone
popNames(Ar230_LM_maf1)

#removing isolates with unknown Zone
Ar230_LM_maf1_zone <- popsub(Ar230_LM_maf1, c(1,2,3,5,6,7,8,9,10,11))

set.seed(999)
xvalAr2021 <- xvalDapc(tab(Ar230_LM_maf1_zone, NA.method = "mean"), pop(Ar230_LM_maf1_zone), training.set = 0.7)
set.seed(999)
xvalAr2021 <- xvalDapc(tab(Ar230_LM_maf1_zone, NA.method = "mean"), pop(Ar230_LM_maf1_zone), parallel = "multicore", ncpus = 4L, training.set = 0.9, n.rep = 1000, n.pca = 30:60)
xvalAr2021 [2:6]
#54
popNames(Ar230_LM_maf1_zone)
cols <- c("Yellow2", "dark green", "green", "turquoise", "sky blue", "magenta", "orange", "purple", "blue", "red")
dapc.zone <- dapc(Ar230_LM_maf1_zone, var.contrib = TRUE, scale = FALSE, n.pca = 54, n.da = nPop(Ar230_LM_maf1_zone) - 1)
#final:
scatter(dapc.zone, col = cols, clabel=FALSE, xax = 1, yax = 2,inset.solid = 1,cell=0, cstar=0, legend = TRUE, cleg = 0.75, posi.leg = "bottomright", cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.03, inset.pca = 0.03, posi.da = "bottomleft")
contrib_zone_axis1 <- loadingplot(dapc.zone$var.contr, axis = 1, thres = 0.005)
contrib_zone_axis1

contrib_zone_axis2 <- loadingplot(dapc.zone$var.contr, axis = 2, thres = 0.003)
contrib_zone_axis2


## POPULATION DIFFERENTIATION

setPop(Ar230_LM_maf1) <- ~Year
popNames(Ar230_LM_maf1)

library(mmod)
#gobal GST
Gst_Hedrick(Ar230_LM_maf1_year)
#pairwise
(D_Ar2021 <- pairwise_D(Ar230_LM_maf1_year))
(GST_Ar2021 <- pairwise_Gst_Nei(Ar230_LM_maf1_year))
(GST_Ar2021 <- pairwise_Gst_Hedrick(Ar230_LM_maf1_year))


##  DISTANCE TREE GENLIGHT OBJECT

#MAC
Ar230_LM_maf1_vcf <- read.vcfR("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/SNPeff/lowmodifierimpact.vcf")
Ar230_LM_maf1_vcf
Ar230_LM_maf1_genlight <- vcfR2genlight(Ar230_LM_maf1_vcf)
Ar230_LM_maf1_genlight
indNames(Ar230_LM_maf1_genlight)
ploidy(Ar230_LM_maf1_genlight) <- 1

strata.df <- read.csv("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/Data/230_strata.csv", head = FALSE, sep = ",")
strata(Ar230_LM_maf1_genlight) <- strata.df
splitStrata(Ar230_LM_maf1_genlight) <- ~Year/State/Location/Farm/Host/Zone/ICC3996_rating/Genesis090_rating/HatTrick_rating/Seamer_rating/Path1/Path2
Ar230_LM_maf1_genlight

setPop(Ar230_LM_maf1_genlight) <- ~Year
popNames(Ar230_LM_maf1_genlight)

ind_dist <- bitwise.dist(Ar230_LM_maf1_genlight)
Ar2021_genlight_tree <- aboot(Ar230_LM_maf1_genlight, tree = "nj", distance = bitwise.dist, sample = 100, showtree = T, cutoff = 50, quiet = T)

# Visualize - colour by year of collection
cols <- c("Red3", "darkorange1", "cyan", "yellow2", "light blue","purple")
popNames(Ar230_LM_maf1_genlight)
library(ape)
plot.phylo(Ar2021_genlight_tree, cex = 0.9, font = 1, adj = 0, tip.color = cols[Ar230_LM_maf1_genlight$pop], label.offset = 0.0001)
plot.phylo(Ar2021_genlight_tree, type="fan", cex = 0.9, font = 2, adj = 0, tip.color = cols[Ar230_LM_maf1_genlight$pop], label.offset = 0.0001)
nodelabels(Ar2021_genlight_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.6, font = 0.5, xpd = TRUE)
axisPhylo(3)
legend('bottomright', legend = c("2014", "2015", "2016", "2017", "2020", "2013"), fill = cols, border = FALSE, bty = "n", cex = 0.8)


##  INDICES OF GENETIC DIVERSITY

setPop(Ar230_LM_maf1) <- ~Zone
popNames(Ar230_LM_maf1)

(Ar2021_diversity <- poppr(Ar230_LM_maf1, sample = 100))
(N <- Ar2021_diversity$N)
(lambda <- Ar2021_diversity$lambda) 
#corrected simpson's index for sample size
(N/(N-1))*lambda
(MLG <- Ar2021_diversity$MLG)
(Clonal_Fraction <- (1-(MLG/N)))

setPop(Ar230_LM_maf1) <- ~State
popNames(Ar230_LM_maf1)

(Ar2021_diversity <- poppr(Ar230_LM_maf1, sample = 100))
(N <- Ar2021_diversity$N)
(lambda <- Ar2021_diversity$lambda) 
#corrected simpson's index for sample size
(N/(N-1))*lambda
(MLG <- Ar2021_diversity$MLG)
(Clonal_Fraction <- (1-(MLG/N)))


## VISUALIZE OUTPUTS FROM VCFTOOLS

library(tidyverse)

#Variant mean depth
#the mean depth for each variant. 
#This is essentially the number of reads that have mapped to this position. 
#the mean of the read depth across all inds - for both alleles at a position and is not partitioned between the reference and the alternative
#10x is a good rule of thumb as a minimum cutoff for read depth
#although if we wanted to be conservative, we could go with 15x
# for max depth, a good rule of thumb is something the mean depth x 2

var_depth <- read_delim("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics/vcfTools_on_final_LM_maf1_230/out.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0,250)
summary(var_depth)

#Variant missingness
#proportion of missingness at each variant. 
#This is a measure of how many individuals lack a genotype at a call site

var_miss <- read_delim("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics/vcfTools_on_final_LM_maf1_230/out.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_miss$fmiss)

#Distribution of allele frequencies will help inform  minor-allele frequency (MAF) thresholds
var_freq <- read_delim("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics/vcfTools_on_final_LM_maf1_230/out.frq", delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
# It is an important measure because low MAF alleles may only occur in one or two individuals. 
# It is possible that some of these low frequency alleles are in fact unreliable base calls - i.e. a source of error.
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
summary(var_freq$maf)

#Distribution of mean depth among individuals.
ind_depth <- read_delim("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics/vcfTools_on_final_LM_maf1_230/out.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, binwidth = 1)
a + theme_light()
view(ind_depth)
summary(ind_depth)

#Proportion of missing data per individual
ind_miss  <- read_delim("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics/vcfTools_on_final_LM_maf1_230/out.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3, bins = 233)
a + theme_light()
view (ind_miss)
summary(ind_miss)


##  FIND CLUSTERS - underlying population structure

library(adegenet)
library(vcfR)
Ar230_LM_maf1_vcf <- read.vcfR("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/data/SNPeff/lowmodifierimpact.vcf")
Ar230_LM_maf1_vcf

Ar230_maf1_LM_genlight <- vcfR2genlight(Ar230_LM_maf1_vcf)
Ar230_maf1_LM_genlight
indNames(Ar230_maf1_LM_genlight)

Ar230_maf1_LM_genlight.clusters <- find.clusters(Ar230_maf1_LM_genlight, max.n.clust = 20, n.start =1000)
#250
#3
dapcAr2021_LM_maf1 <- dapc(Ar230_maf1_LM_genlight, Ar230_maf1_LM_genlight.clusters$grp, var.contrib = TRUE)
# Choose the number PCs to retain (>=1): 
#75
#2
scatter(dapcAr2021_LM_maf1, legend = TRUE, cleg = 1, col = rainbow(6), clabel = 0, cex = 1.5, scree.da = TRUE, scree.pca = TRUE, posi.pca = "topleft", posi.da = "bottomleft", ratio.da = 0.15, ratio.pca = 0.15, inset.da = 0.01, inset.pca = 0.01, posi.leg = "bottomright")

# check which isolates belong to each of the three clusters
dapcAr2021_LM_maf1
dapcAr2021_LM_maf1$grp

###   DAPC again for the figure
library(adegenet)
library(vcfR)
Ar230_LM_maf1_vcf <- read.vcfR("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/data/SNPeff/lowmodifierimpact.vcf")
Ar230_LM_maf1_vcf

pop.data <- read.table("~/Library/CloudStorage/GoogleDrive-vaghefi.n@gmail.com/My Drive/GRDC-Ascochyta/Paper_Microbial_Genomics2/data/230_samples_data_table.txt", sep = "\t", header = TRUE)

all(colnames(Ar230_LM_maf1_vcf@gt)[-1] == pop.data$ISOLATE_ID)
# [1] TRUE
Ar230_maf1_LM_genlight <- vcfR2genlight(Ar230_LM_maf1_vcf)
Ar230_maf1_LM_genlight
indNames(Ar230_maf1_LM_genlight)

maxK <- 10
myMat <- matrix(nrow=10, ncol=maxK)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(Ar230_maf1_LM_genlight, n.pca = 75, choose.n.clust = FALSE,  max.n.clust = maxK, n.start=1000)
  myMat[i,] <- grp$Kstat
}

library(ggplot2)
library(reshape2)
my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1

my_k <- 3:3

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  #set.seed(9)
  grp_l[[i]] <- find.clusters(Ar230_maf1_LM_genlight, n.pca = 75, n.clust = my_k[i], n.start=1000)
  dapc_l[[i]] <- dapc(Ar230_maf1_LM_genlight, pop = grp_l[[i]]$grp, n.pca = 75, n.da = my_k[i])
}

my_df <- as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group <- dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

my_pal <- RColorBrewer::brewer.pal(n=5, name = "Dark2")

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 <- p2 + geom_point(size = 4, shape = 21)
p2 <- p2 + theme_bw()
p2 <- p2 + scale_color_manual(values=c(my_pal))
p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Zone <- pop.data$Zone
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Zone <- pop.data$Zone
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Zone, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_color_brewer(palette="Dark2")
p3 <- p3 + scale_fill_manual(values=c(my_pal))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3

library("ggpubr")
ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(1, 2)
)

dapc_l$grp


## MSN

cols <- c("darkgreen", "lightgreen", "darkorange3", "yellow1", "blue","magenta","green", "cyan", "orange", "yellow3")
setPop(Ar230_maf1) <- ~Zone
popNames(Ar230_maf1)
#removing indivuduals with unknown zone
Ar230_zones <- popsub(Ar230_maf1, c(1,2,3,5,6,7,8,9,10,11))

library(igraph)
adist <- bitwise.dist(Ar230_zones, mat=TRUE, euclidean = FALSE, percent = FALSE)
amsn <- poppr.msn(Ar230_zones, adist, showplot = TRUE, threshold = 0, palette=cols, vertex.label = NA, margin=rep(-0.1,4), wscale = FALSE)
plot_poppr_msn(Ar230_zones, amsn, gadj = 6, layfun = layout_with_kk)
set.seed(500)
plot_poppr_msn(Ar230_zones, amsn, inds="NONE",  palette = cols, nodescale = 15, margin=rep(0,4), wscale = FALSE)
#( inds = "059")
#write.matrix (adist, file = "adist.matrix.bitwise.txt")
#amsn$weight
