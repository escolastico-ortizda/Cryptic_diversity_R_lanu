### --------------------- 0. Clonal analyses in poppR of R. lanuginosum --------------------- ###

# install.packages("poppr")
library("poppr")
packageVersion("poppr")
library("adegenet")
packageVersion("adegenet")

# install.packages("vcfR")
library(vcfR)
packageVersion("vcfR")

# Import vcf file from populations-Stacks for analyses.
setwd("DIRECTORY_NAME_REPLACE")

### --------------------- 1. Load R. lanuginosum Species dataset  (230 samples) --------------------- ###

pop_230 <- read.vcfR("populations_fineRADstructure/R20_230_samples_t125/populations.snps.vcf")
class(pop_230)
pop_230_snps <- vcfR2genlight(pop_230, n.cores = 10)
pop_230_snps@ind.names
pop_230_map <- read.delim("populations_fineRADstructure/R20_230_samples_t125/popmap_230_data.txt",header=T)

# Insert population and other information to the genlight object
pop_230_snps@pop <- as.factor(pop_230_map[,2])
pop_230_snps@strata <- pop_230_map


### --------------------- 2. Principal components analysis of the Species dataset (230 samples) --------------------- ###

# PCA
library(rlang)
library(vcfR)
library(ade4)
library(adegenet)

pop_230_snps
PCA_230 <- glPca(pop_230_snps,nf=2)
barplot(PCA_230$eig, main="eigenvalues",col=heat.colors(length(PCA_230$eig)))
scatter(PCA_230)
sum(PCA_230$eig)
PCA_230$scores

library(ggplot2)

# Transform the scores into a data frame to used it in ggplot.
PCA_230$scores<- as.data.frame(PCA_230$scores)

PCA_groups_230 <- ggplot(PCA_230$scores [c(1,2)],aes(x = PCA_230$scores[,1],y = PCA_230$scores[,2],col = pop_230_snps@strata[,2]))+
  geom_point()+
  geom_vline(xintercept=0,linetype=2,colour="gray")+ 
  geom_hline(yintercept=0,linetype=2,colour="gray")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_ellipse(linetype = 1)+
  scale_color_manual(values = c("dodgerblue1","red","orange", "forestgreen","black","darkgray","lightgray"))


PCA_groups_230 + labs(x="PCA 1 (5.87.82%)",y="PCA 2 (2.94%)",colour="Genetic groups")


### --------------------- 3. Load R. lanuginosum Tundra dataset (125 samples) --------------------- ###

pop_125 <- read.vcfR("populations_fineRADstructure/all_r80_125samples/populations.snps.vcf")
class(pop_125)
pop_125_samples_snps <- vcfR2genlight(pop_125, n.cores = 10)
ploidy(pop_125_samples_snps) <- 1
pop_125_samples_snps@ind.names
# write.table(pop_125_snpclone@ind.names,"Samples_order.txt", sep= "\t")

pop_125_map <- read.delim("populations_fineRADstructure/all_r80_125samples/popmap_125_samples.txt",header=T)
# YOU HAVE TO VERIFY THAT THE ORDER OF YOUR SAMPLES IN THE VCF FILE MATCH THAT OF THE POPMAP !!!

# Insert population and other information to the genlight object
pop_125_samples_snps@pop <- as.factor(pop_125_map[,2])
pop_125_samples_snps@strata <- pop_125_map

# setPop(pop_125_samples_snps) <- ~Group
# pop_125_samples_snps<- popsub(pop_125_samples_snps, exclude=c("AB")) # To separate genetic groups
# pop_AB
# setPop(pop_AB) <- ~Habitat


### --------------------- 4. Principal component analysis of the Tundra dataset (125 samples) --------------------- ###

library(rlang)
library(vcfR)
library(ade4)
library(adegenet)

pop_125_samples_snps
PCA_125 <- glPca(pop_125_samples_snps,nf=2)
barplot(PCA_125$eig, main="eigenvalues",col=heat.colors(length(PCA_125$eig)))
scatter(PCA_125)

library(ggplot2)

# Transform the scores into a data frame to used it in ggplot.
PCA_125$scores<- as.data.frame(PCA_125$scores)

PCA_125_df<- as.data.frame(PCA_125$scores)
PCA_125_df<- cbind(PCA_125_df,pop_125_samples_snps@strata[,13])
names(PCA_125_df)[3] ="Group_Habitat"

# PCA Optional

PCA_125_groups_opt <- ggplot(PCA_125_df,aes(x = PC1,y = PC2, col = Group_Habitat))+
  geom_point(aes(shape=Group_Habitat), size = 2) +
  scale_shape_manual(values=c(17,16,17,16)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_ellipse(linetype = 2) +
  scale_color_manual(values=c("#666666","#666666","#cccccc","#cccccc"))


PCA_125_groups_opt + labs(x="PCA 1 (49.82%)",y="PCA 2 (9.41%)")

ggsave("Figure_S3_PCA_125_samples_groups_habitat.tiff", units="mm", width=140, height=100, dpi=300, compression = "lzw")
ggsave("Figure_S3_PCA_125_samples_groups_habitats.svg",units="mm", width=140, height=100, dpi=300)

# PCA Normal

PCA_125_groups <- ggplot(PCA_125$scores [c(1,2)],aes(x = PCA_125$scores[,1],y = PCA_125$scores[,2],col = pop_125_samples_snps@strata[,3]))+
  geom_point()+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_ellipse(linetype = 2) +
  scale_color_manual(values=c("#666666","#cccccc"))+
  #scale_shape_manual(values=c(17,16))
 
PCA_125_groups + labs(x="PCA 1 (49.82%)",y="PCA 2 (9.41%)",colour="Genetic groups")

#ggsave("Figure_S3_PCA_125_samples_genetic_groups.tiff", units="mm", width=140, height=100, dpi=300, compression = 'lzw')
#ggsave("Figure_S3_PCA_125_samples_genetic_groups.svg",units="mm", width=140, height=100, dpi=300)

PCA_125_habitat <- ggplot(PCA_125$scores [c(1,2)],aes(x = PCA_125$scores[,1],y = PCA_125$scores[,2],col = pop_125_samples_snps@strata[,4]))+
  geom_point()+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_ellipse(linetype = 2)

PCA_125_habitat + labs(x="PCA 1 (49.82%)",y="PCA 2 (9.41%)",colour="Habitat")

#ggsave("Figure_S4_PCA_125_samples_habitat.tiff", units="mm", width=140, height=100, dpi=300, compression = 'lzw')
#ggsave("Figure_S4_PCA_125_samples_habitat.svg",units="mm", width=140, height=100, dpi=300)

# Both figures
library("ggpubr")
ggarrange(PCA_125_groups + 
            labs(x="PCA 1 (49.82%)",y="PCA 2 (9.41%)",colour="Genetic groups"),
          PCA_125_habitat + labs(x="PCA 1 (49.82%)",y="PCA 2 (9.41%)",colour="Habitat") +
            labs(x="PCA 1 (49.82%)",y="PCA 2 (9.41%)"), 
          labels = c("a", "b"),
          ncol = 1, nrow = 2)

#ggsave("Figure_S3_PCA_125_samples_genetic_groups_habitat.tiff", units="mm", width=140, height=200, dpi=300, compression = 'lzw')
#ggsave("Figure_S3_PCA_125_samples_genetic_groups_habitat.svg",units="mm", width=140, height=200, dpi=300)


### --------------------- 5. Genetic diversity comparison between groups AB vs CD (population scale) --------------------- ###
library(tidyr)
library(MASS)
library(dplyr)

# Set the working directory before loading the dataset
habitat_diversity <- read.delim("./Genetic_diversity_comparison/Summary_population_scale_r80_66_samples_even.txt", sep = "\t", header = TRUE)
AB_diversity <- read.delim("./Genetic_diversity_comparison/Summary_population_scale_20_samples_even_AB_habitat.txt", sep = "\t", header = TRUE)
CD_diversity <- read.delim("./Genetic_diversity_comparison/Summary_population_scale_46_samples_even_CD_habitat.txt", sep = "\t", header = TRUE)
  
  
#### Separate diversity per HABITAT
diversity_F <- habitat_diversity[which(habitat_diversity$Pop_ID=="Forest_tundra"),]
diversity_S <- habitat_diversity[which(habitat_diversity$Pop_ID=="Shrub_tundra"),]

# Haplotype diversity per habitat
truehist(diversity_F[,3],prob=F)
truehist(diversity_S[,3],prob=F)

# Normality tests
shapiro.test(diversity_F[,3]) # Not normal distributed
shapiro.test(diversity_S[,3]) # Not normal distributed

# Comparison using Mann-Whitney test
wilcox.test(diversity_F[,3],diversity_S[,3]) # Comparison using Mann-Whitney
# The haplotype diversity do not differ between habitats W = 98754, p-value = 0.6158

# Nlucleotide diversity
truehist(diversity_F[,4],prob=F)
truehist(diversity_S[,4],prob=F)

# Normality tests
shapiro.test(diversity_F[,4]) # Not normal distributed
shapiro.test(diversity_S[,4]) # Not normal distributed

# Comparison using Mann-Whitney test
wilcox.test(diversity_F[,4],diversity_S[,4]) # Comparison using Mann-Whitney
# The nucleotide diversity do not differ between habitats W = 96199, p-value = 0.2464


######## Group AB ######## 

#### Separate diversity per HABITAT
AB_diversity_F <- AB_diversity[which(AB_diversity$Pop_ID=="Forest_tundra"),]
AB_diversity_S <- AB_diversity[which(AB_diversity$Pop_ID=="Shrub_tundra"),]

# Haplotype diversity per habitat
truehist(AB_diversity_F[,3],prob=F)
truehist(AB_diversity_S[,3],prob=F)

# Normality tests
shapiro.test(AB_diversity_F[,3]) # Not normal distributed
shapiro.test(AB_diversity_S[,3]) # Not normal distributed

# Comparison using Mann-Whitney test
wilcox.test(AB_diversity_F[,3],AB_diversity_S[,3]) # Comparison using Mann-Whitney
# The haplotype diversity differs between habitats W = 75817, p-value = 1.809e-05

# Nlucleotide diversity
truehist(AB_diversity_F[,4],prob=F)
truehist(AB_diversity_S[,4],prob=F)

# Normality tests
shapiro.test(AB_diversity_F[,4]) # Not normal distributed
shapiro.test(AB_diversity_S[,4]) # Not normal distributed

# Comparison using Mann-Whitney test
wilcox.test(AB_diversity_F[,4],AB_diversity_S[,4]) # Comparison using Mann-Whitney
# The nucleotide diversity do not differ between habitats W = 75744, p-value = 1.652e-05

######## Group CD ######## 

#### Separate diversity per HABITAT
CD_diversity_F <- CD_diversity[which(CD_diversity$Pop_ID=="Forest_tundra"),]
CD_diversity_S <- CD_diversity[which(CD_diversity$Pop_ID=="Shrub_tundra"),]

# Haplotype diversity per habitat
truehist(CD_diversity_F[,3],prob=F)
truehist(CD_diversity_S[,3],prob=F)

# Normality tests
shapiro.test(CD_diversity_F[,3]) # Not normal distributed
shapiro.test(CD_diversity_S[,3]) # Not normal distributed

# Comparison using Mann-Whitney test
wilcox.test(CD_diversity_F[,3],CD_diversity_S[,3]) # Comparison using Mann-Whitney
# The haplotype diversity differs between habitats W = 153572, p-value = 2.2e-16

# Nucleotide diversity
truehist(CD_diversity_F[,4],prob=F)
truehist(CD_diversity_S[,4],prob=F)

# Normality tests
shapiro.test(CD_diversity_F[,4]) # Not normal distributed
shapiro.test(CD_diversity_S[,4]) # Not normal distributed

# Comparison using Mann-Whitney test
wilcox.test(CD_diversity_F[,4],CD_diversity_S[,4]) # Comparison using Mann-Whitney
# The nucleotide diversity do not differ between habitats W = 152993, p-value = < 2.2e-16


### --------------------- 6. Clonal assignment (MLGs) of the Tundra dataset (125 samples) --------------------- ###

# Create a snpclone object
pop_125_snpclone <- as.snpclone(pop_125_samples_snps)
haploid <- seq(1,1, length.out= 125)
pop_125_snpclone@ploidy <- as.integer(haploid)
pop_125_snpclone # 611 SNPs

# Genetic distance among individuals or samples
#nei.dist(pop_all_snpclone)
#(nei_all_snp<- nei.dist(pop_all_snpclone))


# Use genetic distance to assign genotypes
pop_125_snp_filtered <- filter_stats(pop_125_snpclone, plot=T)

# Save images that are not ggplot objects

#tiff(file="Figure_S1_Genetic_distance_cutoff_inference_125_samples_600_700_300dpi.tiff",
#     width=6, height=4, units="in", res=300)
#pop_125_snp_filtered <- filter_stats(pop_125_snpclone, plot=T)
#dev.off()

#svg("Figure_S1_Genetic_distance_cutoff_inference_125_samples_600_700_300dpi.svg")
#pop_125_snp_filtered <- filter_stats(pop_125_snpclone, plot=T)
#dev.off()

print(farthest_thresh <-cutoff_predictor(pop_125_snp_filtered$farthest$THRESHOLDS)) # 0.02864157
print(average_thresh <-cutoff_predictor(pop_125_snp_filtered$average$THRESHOLDS)) # 0.02741408
print(nearest_thresh <-cutoff_predictor(pop_125_snp_filtered$nearest$THRESHOLDS)) # 0.0008183306

mlg.filter(pop_125_snpclone, threads = 1L) <- 0.0286
pop_125_snpclone
pop_tab_snp<-mlg.table(pop_125_snpclone, strata = ~Habitat)


#tiff(file="Fig_S6_New_multilocus_genotype_counts_125_samples_300dpi.tiff",
#   width=6, height=4, units="in", res=300)
#pop_tab_snp<-mlg.table(pop_125_snpclone, strata = ~Habitat)
#dev.off()

#svg("Fig_S6_New_multilocus_genotype_counts_125_samples_300dpi.svg")
#pop_tab_snp<-mlg.table(pop_125_snpclone, strata = ~Habitat)
#dev.off()

# Using MLG information plot the clones at fine scale. 
mlg.vector(pop_125_snpclone)
pop_125_snpclone@strata
mlg.id(pop_125_snpclone)
#write.table(pop_125_snpclone@mlg@mlg[["contracted"]],"MLG__125_samples_snps.txt", sep= "\t")

new_MLGs <- read.delim("poppR/new_MLGs_125_samples.txt",header=F)
new_MLGs <- as.vector(new_MLGs[,1])
class(new_MLGs)

mll(pop_125_snpclone)
mll.custom(pop_125_snpclone, set=TRUE) <- new_MLGs

### --------------------- 7. Population scale (125 samples) - Clonal divesity and structure  --------------------- ###


# Filter by genetic group and habitat with an equal sample size. The fine-scale samples must be filtered out.
# Corrected sample size
setPop(pop_125_snpclone) <- ~Correction_Habitat_Group
pop_125_snp_corrected <- popsub(pop_125_snpclone, sublist = c("TRUE"))
pop_125_snp_corrected@pop
pop_125_snp_corrected@ind.names
mll(pop_125_snp_corrected)
pop_125_snp_corrected@strata
pop_tab_125_corrected <- mlg.table(pop_125_snp_corrected, strata = ~Group/Habitat)

#write.table(pop_tab_125_corrected,"MLGs_per_habitat_125_samples.txt", sep= "\t")


#tiff(file="Fig_S7_New_multilocus_genotype_counts_genetic_group_125_samples_300dpi.tiff",
#     width=6, height=4, units="in", res=300)
#pop_tab_125_corrected <- mlg.table(pop_125_snp_corrected, strata = ~Group/Habitat)
#dev.off()

#svg("Fig_S7_New_multilocus_genotype_counts_genetic_group_125_samples_300dpi.svg")
#pop_tab_125_corrected <- mlg.table(pop_125_snp_corrected, strata = ~Group/Habitat)
#dev.off()

(pop_stats_125_corrected <- diversity_stats(pop_tab_125_corrected))
diversity_ci(pop_tab_125_corrected, n = 100L, rarefy = TRUE, raw = FALSE)

# AMOVA
# DO not apply the new mlls otherwise the anova will not work
table(strata(pop_125_snp_corrected, ~Group/Habitat, combine = FALSE))  # Subpopulations
pop_125_snp_corrected 

# Split the dataset per genetic group

# Group AB
setPop(pop_125_snp_corrected) <- ~Group
pop_AB_snp_corrected <- popsub(pop_125_snp_corrected, sublist = c("AB"))

result_amova_AB <- poppr.amova(pop_AB_snp_corrected,~Habitat/Quadrat)
amova.test_AB <- randtest(result_amova_AB, nrepet = 999)
plot(amova.test_AB )

# Clone corrected group AB
result_amova_AB_clone <- poppr.amova(pop_AB_snp_corrected,~Habitat/Quadrat,clonecorrect = TRUE)
amova.test_AB_clone <- randtest(result_amova_AB_clone, nrepet = 999)
plot(amova.test_AB_clone)

# Group CD
pop_CD_snp_corrected <- popsub(pop_125_snp_corrected, sublist = c("CD"))

result_amova_CD <- poppr.amova(pop_CD_snp_corrected,~Habitat/Quadrat)
amova.test_CD <- randtest(result_amova_CD, nrepet = 999)
plot(amova.test_CD)

# Clone corrected group CD
result_amova_CD_clone <- poppr.amova(pop_CD_snp_corrected,~Habitat/Quadrat,clonecorrect = TRUE)
amova.test_CD_clone <- randtest(result_amova_CD_clone, nrepet = 999)
plot(amova.test_CD_clone)

### --------------------- 8. Isolation by distance analyses --------------------- ###

# install.packages("dartR")
# BiocManager::install("SNPRelate")
# install.packages("rgdal")
# install.packages("dismo")

library(dartR)
gl.install.vanilla.dartR()

# Work with the genlight object, not with the snpclone.

pop_125_samples_snps@pop
setPop(pop_125_samples_snps) <- ~Quadrat
pop_scale_125_snps <- popsub(pop_125_samples_snps,exclude = c("393","414"))
pop_scale_125_snps@pop
setPop(pop_scale_125_snps) <- ~Quadrat
ploidy(pop_scale_125_snps) <- 1


# Population scale using 125 samples

pop_scale_125_snps@strata
pop_scale_125_snps@other$latlon <- data.frame(pop_125_map[51:125,7:8])
pop_scale_125_snps@other$pop <- data.frame(pop_125_map[51:125,10])

#If any of the optional content slots indicated above are missing, consider running
gl <- gl.compliance.check(pop_scale_125_snps)


# Isolation by distance analyses

#ibd_all <- gl.ibd(pop_all_snps)
#str(ibd_all)
#plot(ibd_all$Dgeo, ibd_all$Dgen)

ibd_pop <- gl.ibd(pop_scale_125_snps) # (Options, distance = 'D'  paircols='pop')
# Mantel statistic r: 0.1138
# Significance: 0.012

# Upper quantiles of permutations (null model):
#  90%    95%  97.5%    99% 
#  0.0550 0.0745 0.0898 0.0984
# Permutation: free
# Number of permutations: 999

# Save images that are not ggplot objects

#tiff(file="Figure_XX_IBD_Genetic_Geographic_distances_A_group_700_600_300dpi.tiff",
#     width=6, height=4, units="in", res=300)
#ibd_pop <- gl.ibd(pop_scale_125_snps)
#dev.off()

#svg("Figure_XX_IBD_Genetic_Geographic_distances_A_group_700_600_300dpi.svg")
#ibd_pop <- gl.ibd(pop_scale_125_snps)
#dev.off()


str(ibd_pop)
plot(ibd_pop$Dgeo, ibd_pop$Dgen)

# Formatting IBD results
ibd.matrix <- as.matrix(ibd_pop$Dgen)
ind <- which(upper.tri(ibd.matrix),arr.ind = T)
ibd.df <- data.frame(Site1 = dimnames(ibd.matrix)[[2]][ind[,2]],
                    Site2 = dimnames(ibd.matrix)[[1]][ind[,1]],
                    ibd.scores=ibd.matrix[ ind ] %>% round(digits = 4))
# Distance is based on Fst
ibd.df

# Geographic distance
dgeo.matrix <- as.matrix(ibd_pop$Dgeo)
ind = which( upper.tri(dgeo.matrix), arr.ind = TRUE)
dgeo.df <- data.frame(Site1 = dimnames(dgeo.matrix)[[2]][ind[,2]],
                      Site2 = dimnames(dgeo.matrix)[[1]][ind[,1]],
                      dgeo.scores=dgeo.matrix[ ind ] %>% round(digits = 4))

dgeo.df
ibd.df$geodist <- dgeo.df$dgeo.scores

# Plotting IBD results
library(ggplot2)
library(ggpubr)


plot.ibd <- ggplot(ibd.df, aes(x=geodist, y=ibd.scores)) +
  geom_point() +
  labs(x="log (Geographic distance)", y ="Fst/1-Fst") +
  geom_smooth(method=lm, fill="lightgrey", se=T, level=0.95)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot.ibd

#ggsave("Figure_4_IBD_Genetic_Geographic_distances.tiff", units="mm", width=140, height=140, dpi=300, compression = 'lzw')
#ggsave("Figure_4_IBD_Genetic_Geographic_distances.svg",units="mm", width=140, height=140, dpi=300)


# IBD per genetic group

pop_scale_125_snps@pop
setPop(pop_scale_125_snps) <- ~Group
pop_scale_125_snps@strata
pop_scale_snps_AB <- popsub(pop_scale_125_snps,exclude = "AB")
pop_scale_snps_CD <- popsub(pop_scale_125_snps,exclude = "CD")

pop_scale_snps_AB
setPop(pop_scale_snps_AB) <- ~Quadrat

pop_scale_snps_CD
setPop(pop_scale_snps_CD) <- ~Quadrat

# Group AB
ibd_pop_AB <- gl.ibd(pop_scale_snps_AB,Dgeo_trans='log(Dgeo)' ,
                     Dgen_trans='Dgen/(1-Dgen)',paircols ='Group') # Mantel test statistic= 0.1191, p=0.019
#Mantel statistic r: 0.1987 
#Significance: 0.001 
str(ibd_pop_AB)
plot(ibd_pop_AB$Dgeo, ibd_pop_AB$Dgen)


# Group AB - Formatting IBD results
ibd.matrix_AB <- as.matrix(ibd_pop_AB $Dgen)
ind <- which(upper.tri(ibd.matrix_AB),arr.ind = T)
ibd.df_AB <- data.frame(Site1 = dimnames(ibd.matrix_AB)[[2]][ind[,2]],
                     Site2 = dimnames(ibd.matrix_AB)[[1]][ind[,1]],
                     ibd.scores=ibd.matrix_AB[ ind ] %>% round(digits = 4))

# Group AB - Geographic distance
dgeo.matrix_AB <- as.matrix(ibd_pop_AB$Dgeo)
ind = which( upper.tri(dgeo.matrix_AB), arr.ind = TRUE)
dgeo.df_AB <- data.frame(Site1 = dimnames(dgeo.matrix_AB)[[2]][ind[,2]],
                      Site2 = dimnames(dgeo.matrix_AB)[[1]][ind[,1]],
                      dgeo.scores=dgeo.matrix_AB[ ind ] %>% round(digits = 4))

dgeo.df
ibd.df_AB$geodist <- dgeo.df_AB$dgeo.scores

plot.ibd_AB <- ggplot(ibd.df_AB, aes(x=geodist, y=ibd.scores)) +
  geom_point() +
  labs(x="log (Geographic distance)", y ="Fst/1-Fst") +
  geom_smooth(method=lm, fill="lightgrey", se=T, level=0.95)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot.ibd_AB

# ggsave("Figure_S8_IBD_Genetic_Geographic_distances_group_AB.tiff", units="mm", width=140, height=140, dpi=300, compression = 'lzw')
# ggsave("Figure_S8_IBD_Genetic_Geographic_distances_group_AB.svg",units="mm", width=140, height=140, dpi=300)


# Group CD
ibd_pop_CD <- gl.ibd(pop_scale_snps_CD) 
# Mantel statistic r= 0.1325 ; p=0.066 
pop_scale_snps_CD@other
str(ibd_pop_CD)
plot(ibd_pop_CD$Dgeo, ibd_pop_CD$Dgen)

ibd_pop_CD$Dgen

# Group CD - Formatting IBD results
ibd.matrix_CD <- as.matrix(ibd_pop_CD $Dgen)
ind <- which(upper.tri(ibd.matrix_CD),arr.ind = T)
ibd.df_CD <- data.frame(Site1 = dimnames(ibd.matrix_CD)[[2]][ind[,2]],
                       Site2 = dimnames(ibd.matrix_CD)[[1]][ind[,1]],
                       ibd.scores=ibd.matrix_CD[ ind ] %>% round(digits = 4))

# Group CD - Geographic distance
dgeo.matrix_CD <- as.matrix(ibd_pop_CD$Dgeo)
ind = which( upper.tri(dgeo.matrix_CD), arr.ind = TRUE)
dgeo.df_CD <- data.frame(Site1 = dimnames(dgeo.matrix_CD)[[2]][ind[,2]],
                        Site2 = dimnames(dgeo.matrix_CD)[[1]][ind[,1]],
                        dgeo.scores=dgeo.matrix_CD[ ind ] %>% round(digits = 4))

dgeo.df
ibd.df_CD$geodist <- dgeo.df_CD$dgeo.scores

plot.ibd_CD <- ggplot(ibd.df_CD, aes(x=geodist, y=ibd.scores)) +
  geom_point() +
  labs(x="log (Geographic distance)", y ="Fst/1-Fst") +
  geom_smooth(method=lm, fill="lightgrey", se=T, level=0.95)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot.ibd_CD

# ggsave("Figure_S9_IBD_Genetic_Geographic_distances_group_CD.tiff", units="mm", width=140, height=140, dpi=300, compression = 'lzw')
# ggsave("Figure_S9_IBD_Genetic_Geographic_distances_group_CD.svg",units="mm", width=140, height=140, dpi=300)

