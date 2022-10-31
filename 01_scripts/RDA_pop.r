# RUN RDA FROM 012 IMP'ed GENOTYPE MATRIX

# Install libraries
#install.packages("vegan")
library(data.table)
library(vegan)

argv <- commandArgs(T)

# Specify paths
GEN_MAT_NAME <- argv[1]
GEN_MAT_PATH <- paste0(GEN_MAT_NAME, '.012.imp')
INDV_FILE <- paste0(GEN_MAT_NAME, '.012.indv')

OUT_DIR <- "rda"
ID_SEX_POP <- argv[2]
CORES <- argv[3]

# 1. Import 012 matrix 
gen_mat <- fread(GEN_MAT_PATH, header = TRUE, showProgress = TRUE) 
gen_mat <- as.data.frame(gen_mat)
print('done importing genotype matrix')

# Check missing data
sum(gen_mat == '-1')

# 2. Import 012 matrix samples' IDs
gen_mat_ID <- fread(INDV_FILE, header = FALSE, col.names = 'ID')


# 3. Import/build predictors dataframe
ID_sex_pop <- read.table(ID_SEX_POP, header = FALSE,
                         col.names = c("ID", "sex", "pop"))

#merge(gen_mat_ID, ID_sex_pop, by = 'ID')
                         
# Reorder predictors IDs according to gen_mat_ID
ID_sex_pop <-ID_sex_pop[match(gen_mat_ID$ID, ID_sex_pop$ID), ]

# Add IDs as row names for genotype matrix
rownames(gen_mat) <- gen_mat_ID$ID


# 4. Run rda
print('starting RDA')
pop.rda <- rda(gen_mat ~ ID_sex_pop$pop, scale = TRUE)

print(paste('completed RDA. saving to ', OUT_DIR))
saveRDS(pop.rda, file = paste0(OUT_DIR, '/pop_rda.rds'))

###pop.rda <- readRDS('rda/pop_rda.rds')

# Calculate R^2
print('Calculating R^2')
RsquareAdj(pop.rda) # r.squared=0.0181 # adj.r.squared=0.00118
summary(pop.rda)$concont

## export table
write.table(RsquareAdj(pop.rda), paste0(OUT_DIR, "/adjR2_pop.txt"), row.names = FALSE, quote = FALSE)
#pdf(file=paste0(OUT_DIR, "/screeplot.pdf"))

#dev.off()

# 6. Check RDA signifiance using ANOVA
print('checking sgnifiance using ANOVA')

#signif.full <- anova.cca(pop.rda, parallel = getOption("mc.cores")) # default is permutation=999
signif.full <- anova.cca(pop.rda, parallel = CORES)
signif.full
saveRDS(signif.full, file = paste0(OUT_DIR, '/signif.full_pop.rds'))
#signif.axis <- anova.cca(pop.rda, by = "axis", parallel = getOption("mc.cores"))
signif.axis <- anova.cca(pop.rda, by = "axis", parallel = CORES)
signif.axis
saveRDS(signif.axis, file = paste0(OUT_DIR, '/signif.axis_pop.rds'))



# 7. Identify candidate SNPs involved in local adaptation
# 7.1 Extract loadings
load.rda <- summary(pop.rda)$species[, 1:3]
#load.rda <- scores(pop.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
print('done extracting species')

pdf(file = paste0(OUT_DIR, "/hist_loadings_pop.pdf"))
par(mfrow = c(1, 3))
hist(load.rda[, 1], main = "Loadings on RDA1")
hist(load.rda[, 2], main = "Loadings on RDA2")
hist(load.rda[, 3], main = "Loadings on RDA3")
dev.off()
## export table
write.table (load.rda[, 1:3], paste0(OUT_DIR, "/rda_loading_pop.txt"), quote=FALSE)

# Outlier detection function
print('starting outlier detection')
outliers <- function(x, z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# 7.2 Detect outliers for RDA1 and 2
# 3 sd is standard, 2 sd is less conservative
cand1 <- outliers(load.rda[, 1], 3) # candidates SNPs for RDA1
cand2 <- outliers(load.rda[, 2], 3) # candidates SNPs for RDA2
#cand3 <- outliers(load.rda[, 3], 3) 

# 7.3 Number of candidate outlier SNPs
ncand <- length(cand1) + length(cand2) #+ length(cand3)
cat(paste('number of candidate SNPs : ', ncand, "\n"))

lim <- 3
if (ncand == 0){
  cand1 <- outliers(load.rda[, 1], 2.5) 
  cand2 <- outliers(load.rda[, 2], 2.5) 
  #cand3 <- outliers(load.rda[, 3], 2) 
  ncand<-length(cand1) + length (cand2) #+ length ( cand3) # total nb of snps
  lim <- 2.5
  print("no SNP outlier with 3 sd, outlier detection re-done with 2.5 sd")
}

# 7.4 Store candidate SNPs with axis, SNP name, loading, & correlation with each predictor
## For each axis/RDA
cand1 <- cbind.data.frame(axis = rep(1, times = length(cand1)), snp = names(cand1), loading = unname(cand1))
cand2 <- cbind.data.frame(axis = rep(2, times = length(cand2)), snp = names(cand2), loading = unname(cand2))
#cand3 <- cbind.data.frame(rep(3,times = length(cand3)), names(cand3), unname(cand3))
#colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis", "snp", "loading")

## Bind all together
cand <- rbind(cand1, cand2)#, cand3)
cand$snp <- as.character(cand$snp)
cand$CHROM <- sapply(X = cand$snp, FUN = function(x) unlist(strsplit(x, split = '_'))[1])
cand$POS <- sapply(X = cand$snp, FUN = function(x) as.numeric(unlist(strsplit(x, split = '_'))[2]))
cand$END <- cand$POS + 1

#saveRDS(cand, file = paste0(OUT_DIR, '/cand.rds')
write.table(cand[, c('CHROM', 'POS', 'END', 'loading')], file = paste0(OUT_DIR, '/cand_pop_SNPs.txt'), sep = "\t", quote = FALSE, row.names = FALSE)

## Check for duplicates
length(cand$snp[duplicated(cand$snp)])
cand <- cand[!duplicated(cand$snp), ]



# 6. Plot RDA SNPs
# SNPs at middle (red), samples are black, vectors are predictors
print('plotting')
jpeg(file = paste0(OUT_DIR, "/rda_1-2_1-3_pop.jpg"))
par(mfrow = c(1, 2)) 
plot(pop.rda, scaling = 3)
text(pop.rda, display = "sites", col = 1, scaling = 3)
plot(pop.rda, choices = c(1,3), scaling = 3) # axis 1 and 3
text(pop.rda, display = "sites", col = 1, scaling = 3)
dev.off()

print('done')

## Add correlation with each predictor : NOT RELEVANT FOR A SINGLE CATEGORICAL PREDICTOR
# n <- dim(ID_sex_pop)[2]
#foo <- matrix(nrow = (ncand), ncol = 1)  # ncol = number of predictors, nrow = nrow(cand)
#colnames(foo) <- c("sex")

#for (i in 1:length(cand$snp)) {
#  nam <- cand[i, 2] # SNP's name
#  snp.gen <- gen_mat[, ..nam] # # extract all GT for a given SNP
#  #foo[i, ] <- apply(ID_sex_pop, 2, function(x) cor(x, snp.gen)) # calculate correlation for each ID_SEX_POP column
#  foo[i, ] <- cor(ID_sex_pop$sex, snp.gen)
#}
#table of candidate snp with loading on each axis and correlation with env predictors
#cand <- cbind.data.frame(cand,foo)  
#head(cand)

#cand <- cbind.data.frame(cand,foo)  
#head(cand)

# 8. Investigate candidate SNPs
# Any SNPs associated with several axis? if yes remove them
#n_dupli <- length(cand$snp[duplicated(cand$snp)])
#n_dupli
#if (n_dupli >= 1){cand <- cand[!duplicated(cand$snp)]}

# Find the most stronly correlated predictor
#n < -dim(cand)[2]
#for (i in 1:length(cand$snp)) {
#  bar <- cand[i,]
#  cand[i, n+1] <- names(which.max(abs(bar[4:n]))) # gives the variable
#  cand[i, n+2] <- max(abs(bar[4:n]))              # gives the correlation
#}

#colnames(cand)[n+1] <- "predictor"
#colnames(cand)[n+2] <- "correlation"
#head(cand)
#table(cand$predictor) 
#write.table(cand, paste0(OUTPUT_FOLDER,ENV,"_candidate_SNP_",lim,"_sd.txt"), quote=FALSE, sep=" ", row.names=FALSE)


