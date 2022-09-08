# Run PCA with genotype matrix "$MAT_file"
# 

# Libraries
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(reshape2)
library(data.table)

# Access to files provided in command line arguments
argv <- commandArgs(T)
FILE <- argv[1] # "$MAT_file"
ID_SEX_POP <- argv[2]
#pca_path <- argv[2] # $pca_path

N_SAMPLES <- nrow(FILE)

# 1. Load SNPs info
## Create a new column for SNP info name (CHR + position)
geno.012.pos <- read.table(paste0(FILE, ".pos")) %>% #read genotype_matrix.pos; 1st column = CHR, 2nd column = position
  mutate(., locus = paste(V1, V2, sep = '_')) #paste CHR (V1) and position (V2) together and add this new "locus" column to genotype_matrix.pos
  
  
# 2. Load individuals info
geno.012.indv <- read.table(paste0(FILE, ".indv"), col.names = c("ID")) #read genotype_matrix.indv
#geno.012.indv <- read.table(paste0(FILE, ".pos")) 


# 3. Load genotype matrix 
geno.012 <- fread(FILE)[, -1] #read genotype_matrix and exclude 1st column from .012 matrix file 

## Set rownames and colnames for the genotype matrix
colnames(geno.012) <- geno.012.pos$locus #each column is a SNP
#rownames(geno.012) <- geno.012.indv #each row is a different indv

## Inspect matrix
#geno.012[1:12, 1:20] #show 12 first indv and 20 first SNPs 


# 4. Imputation of missing genotypes : replace missing geno (-1) with the most frequent genotype (0, 1, 2) for given SNP
geno.012.imp <- apply(
  geno.012, 2, function(x){
    replace(
      x, 
      which(x == "-1"), 
      max(
        0, 
          as.numeric(
            names(
              which.max(
              table(x)     # genotypes frequencies
              )
            )
          )
        )
      )
    }
  )


## Save as table to .imp file
write.table(geno.012.imp, paste0(FILE, ".imp"), sep = "\t", row.names = F, quote = F)


# 5. Run PCA
pca.all <- prcomp(geno.012.imp)

## Save screeplot
jpeg(paste0(FILE ,"screeplot.jpg"))
screeplot(pca.all)
dev.off()


# 6. Get PCA stats info
sum.pca <- summary(pca.all)

## Print stats info
#sum.pca$importance[ ,1:5]


# 7. Extract PCA % (proportion of variance * 100)
var1 <- round(sum.pca$importance[2, 1]*100, 1) #prop of variance of PC 1 * 100
var2 <- round(sum.pca$importance[2, 2]*100, 1) #prop of variance of PC 2 * 100
var3 <- round(sum.pca$importance[2, 3]*100, 1) #prop of variance of PC 3 * 100
var4 <- round(sum.pca$importance[2, 4]*100, 1) #prop of variance of PC 4 * 100


# 8. Save PCA results and loadings 
#replace Z by number of indv or samples
#write.table(pca.all$x[, 1:Z], paste0(FILE, ".pca"), sep = "\t", quote = F)

write.table(pca.all$x[, 1:N_SAMPLES], paste0(FILE, ".pca"), sep = "\t", quote = F)
write.table(pca.all$rotation[, 1:N_SAMPLES], paste0(FILE, ".loadings"), sep = "\t", quote = F)


# 9. Plot PCA

## Save as dataframe
SNP_df <- as.data.frame(pca.all$x[ , 1:N_SAMPLES])

## Add column for pop data
grouped <- fread(ID_SEX_POP, header = FALSE, col.names = c("ID", "sex", "pop")) # read ID - Pop correspondance file

ID_Pop <- merge(x = geno.012.indv, y = grouped, by = "ID") 

SNP_df <- cbind(SNP_df, ID_Pop)
saveRDS(SNP_df, file = paste0(FILE, "_SNP_df.rds"))


## PC 1 & 2
jpeg(paste0(FILE,".pc1-2.jpg"))

PC_1_2 <- ggplot(data = SNP_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Pop, shape = Pop)) +  
  scale_color_manual(values = c("red", "blue")) +
  guides(size = FALSE) +
  labs(title = "SNP merged", x = paste0("PC1  (", var1, ")"), y = paste0("PC2  (", var2, ")"))
PC_1_2
dev.off()
saveRDS(PC_1_2, file = paste0(FILE, "_PC1_2.plot.rds"))

## PC 1 & 3
jpeg(paste0(FILE,".pc1-3.jpg"))

PC_1_3 <- ggplot(data = SNP_df, aes(x = PC1, y = PC3)) +
  geom_point(aes(color = Pop, shape = Pop)) +  
  scale_color_manual(values = c("red", "blue")) +
  guides(size = FALSE) +
  labs(title = "SNP merged", x = paste0("PC1  (", var1, ")"), y = paste0("PC3  (", var3, ")"))
PC_1_3
dev.off()
saveRDS(PC_1_3, file = paste0(FILE, "_PC1_3.plot.rds"))

## PC 2 & 3
jpeg(paste0(FILE,".pc2-3.jpg"))

PC_2_3 <- ggplot(data = SNP_df, aes(x = PC2, y = PC3)) +
  geom_point(aes(color = Pop, shape = Pop)) +  
  scale_color_manual(values = c("red", "blue")) +
  guides(size = FALSE) +
  labs(title = "SNP merged", x = paste0("PC2  (", var2, ")"), y = paste0("PC3  (", var3, ")"))
PC_2_3
dev.off()
saveRDS(PC_2_3, file = paste0(FILE, "_PC2_3.plot.rds"))

## PC 3 & 4
jpeg(paste0(FILE,".pc3-4.jpg"))

PC_3_4 <- ggplot(data = SNP_df, aes(x = PC3, y = PC4)) +
  geom_point(aes(color = Pop, shape = Pop)) +  
  scale_color_manual(values = c("red", "blue")) +
  guides(size = FALSE) +
  labs(title = "SNP merged", x = paste0("PC3  (", var3, ")"), y = paste0("PC4  (", var4, ")"))
PC_3_4
dev.off()
saveRDS(PC_3_4, file = paste0(FILE, "_PC3_4.plot.rds"))


# 10. Plot loadings

## Generate loadings matrix
loading_mat <- cbind(geno.012.pos, pca.all$rotation[, 1:N_SAMPLES])
colnames(loading_mat)[1:2] <- c("CHR", "pos")

## PC1 loadings
jpeg(paste0(FILE, ".loadings_pc1.jpg"), width = 600, height = 600, quality = 90)
ggplot(data = loading_mat, aes(x = pos/1000000, y = PC1)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~CHR, scales = "free_x") + 
  theme_classic()
dev.off()

## PC2 loadings
jpeg(paste0(FILE, ".loadings_pc2.jpg"), width = 600, height = 600, quality = 90)
ggplot(data = loading_mat, aes(x = pos/1000000, y = PC2)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~CHR, scales = "free_x") + 
  theme_classic()
dev.off()

## PC3 loadings
jpeg(paste0(FILE, ".loadings_pc3.jpg"))
ggplot(data = loading_mat, aes(x = pos/1000000, y = PC3)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~CHR, scales = "free_x") + theme_classic()
dev.off()

## PC4 loadings
jpeg(paste0(FILE, ".loadings_pc4.jpg"))
ggplot(data = loading_mat, aes(x = pos/1000000, y = PC4)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~CHR, scales = "free_x") + 
  theme_classic()
dev.off()