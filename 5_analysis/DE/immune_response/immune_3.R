library(Seurat)
library(limma)
#library(devtools)
library(ggbiplot)
library(dplyr)
#library(sva)
library(SingleCellExperiment)
require(devtools)
#install_version("Seurat", version = "4.0.5", repos = "http://cran.us.r-project.org")


#setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/DE/child_vs_adult/")

load("sub_RNA_5.RData")


sub_5 <- RenameIdents(sub_5, 
                      
                      '0'='Suprabasal',
                      '1_0'='Secretory-1',
                      '1_1'='Secretory-1',
                      '1_2'='Secretory-1',
                      '1_3'='Secretory-1',
                      '1_4'='Secretory-1',
                      '1_5'='Goblet',
                      '1_6'='Secretory-1',
                      '1_7'='Goblet/Ionocyte',
                      '2'	="Basal-1",
                      '3'	="Ciliated-1",
                      '4_0'	="Goblet",
                      '4_1'	="Secretory-2",
                      '4_2'	="Goblet",
                      '4_3'	="Secretory-2",
                      '4_4'	="Secretory-2",
                      '4_5'	="Secretory-2",
                      '4_6'	="Secretory-2",
                      '4_7'	="Secretory-2",
                      '4_8'	="Secretory-2",
                      '5'	="Ciliated-3",
                      '6'	="Cycling-basal",
                      '7'	="Ciliated-2",
                      '8'	="Basal-2",
                      '9_0'	="Ciliated-SC2",
                      '9_1'	="Secretory-Ciliated-SC2",
                      '10_0' ="Goblet-IFN-stim",
                      '10_1' ="Secretory-3",
                      '10_2' ="Secretory-3",
                      '10_3' ="Secretory-3",
                      '10_4' ="Goblet-IFN-stim",
                      '10_5' ="Secretory-3",
                      '10_6' ="Secretory-3",
                      '10_7' ="Secretory-3",
                      '11'	="Deuterosomal",
                      '12_0' = "Ionocyte",
                      '12_1' = "Brush/Tuft",
                      '13' = "Secretory-Ciliated"
)

sub_5@meta.data$cluster <- Idents(sub_5)                                       

#------------pseudobulk approach (https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html)

##detach packages
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))


library(scater)
#library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
#library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
#library(DESeq2)
library(RColorBrewer)
#library(edgeR)
library(limma)
library(Glimma)
library(naniar)
library(ggfortify)
#library(limmaDE2)


# Extract raw counts and metadata to create SingleCellExperiment object
counts <- sub_5@assays$RNA@counts 

metadata <- sub_5@meta.data
metadata_2 <- sub_5@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(sub_5@active.ident)
metadata_2$cluster_id <- factor(sub_5@active.ident)
metadata$id <- paste(metadata$group,metadata$Infection,metadata$Sample,sep = "_")
metadata_2$id <- paste(metadata$cluster_id,metadata$group,metadata$Infection,metadata$Sample,sep = "_")
metadata_2$id <- factor(metadata_2$id)
metadata$group_id <- paste(metadata$group,metadata$Infection,sep = "_")
metadata_2$group_id <- paste(metadata$group,metadata$Infection,sep = "_")
metadata$id <- factor(metadata$id)
metadata$group_id <- factor(metadata$group_id)
metadata_2$group_id <- factor(metadata_2$group_id)


#check_rows:
#print(nrow(metadata_2[metadata_2$id == "Secretory_Ciliated_1_UK_72hpi_adult_Uninfected_Sample1_Sample1",]))


# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

sce_2 <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata_2)




# Explore the raw counts for the dataset

## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

counts(sce)[1:6, 1:6]

## Explore the cellular metadata for the dataset
dim(colData(sce))

head(colData(sce))



## new section filtering


# Perform QC if not already performed
dim(sce)
dim(sce_2)


## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
sce_2 <- sce_2[rowSums(counts(sce_2) > 1) >= 10, ]

dim(sce)
dim(sce_2)
# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$id))
sids_2 <- purrr::set_names(levels(sce_2$id))
# Total number of samples 
ns <- length(sids)
ns_2 <- length(sids_2)
ns
ns_2

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$id)
table(sce_2$id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$id))
n_cells_2 <- as.numeric(table(sce_2$id))
## Determine how to reorder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$id)
m_2 <- match(sids_2, sce_2$id)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")

ei_2 <- data.frame(colData(sce_2)[m_2, ], 
                 n_cells_2, row.names = NULL) %>% 
  select(-"cluster_id")


ei
ei_2


##get list of samples with greater than or equal to 15 counts
ei_2$n_cells_3 <- ifelse(ei_2$n_cells_2 >=15,print("=>15"),print("<15"))

new_count_list <- ei_2[ei_2$n_cells_2 >=15,]
new_count_list_2 <- new_count_list$id


##subset metadata_2 based on these

## PCA 
##PCA subset numeric columns

nums <- unlist(lapply(ei_2, is.numeric), use.names = FALSE)
subset <- ei_2[ , nums]

num2 <- unlist(lapply(metadata_2, is.numeric), use.names = FALSE)
subset2 <- metadata_2[ , num2]

##or PCA only nCount

#subset <- ei_2$nCount_RNA

pca_res <- prcomp(subset, scale. = TRUE)
summary(pca_res)


pca_res2 <- prcomp(subset2, scale. = TRUE)
summary(pca_res2)



autoplot(pca_res2,data=metadata_2,colour = 'Sample')


#---------------------------

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id","id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

pb[1:6, 1:6]



# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                   pattern = "_",  
                                 n = 2), 
                 `[`, 1)
# ##Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 #stringr::str_extract(rownames(u), "(?<=_)[:print:]+")))
                 stringr::str_extract(rownames(u), "[:print:]+")))


class(pb)
# Explore the different components of list
str(pb)



# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$id)

# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

# Create a data frame with the sample IDs, cluster IDs and condition

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    id = de_samples)

gg_df <- left_join(gg_df, ei_2[, c("id", "group_id", "Sample")]) 


gg_df_2 <- left_join(gg_df, new_count_list[, c("id", "group_id", "Sample")]) 
gg_df_2 <- gg_df_2 %>% drop_na()

#-----------------------------------
metadata <- gg_df_2 %>%
  dplyr::select(cluster_id, id, group_id, Sample) 
metadata

metadata$cluster_id <- as.factor(metadata$cluster_id)
# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
clusters



#-----batch correction

##no subset aggregate counts

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

count_matrix <- as.matrix(t(pb))


##get list of samples with greater than or equal to 15 counts
ei_2$n_cells_3 <- ifelse(ei_2$n_cells_2 >=15,print("=>15"),print("<15"))
new_count_list <- ei_2[ei_2$n_cells_2 >=15,]
new_count_list_2 <- new_count_list$id

count_matrix_2 <- count_matrix[, new_count_list_2] #######15 cells!!!!!


#---------------------------Try limiting to >= 15 before adjusting
batch_d2 <- data.frame(colnames(count_matrix_2))

##----include data about donor
batch_d2$batch <- ifelse(grepl("Sample1",batch_d2$colnames.count_matrix_2.)=="TRUE", print("1"),batch_d2$colnames.count_matrix_2.) 
batch_d2$batch <- ifelse(grepl("Sample2",batch_d2$colnames.count_matrix_2.)=="TRUE", print("2"),batch_d2$batch) 
batch_d2$batch <- ifelse(grepl("Sample3",batch_d2$colnames.count_matrix_2.)=="TRUE", print("3"),batch_d2$batch) 
batch_d2$batch <- ifelse(grepl("Sample4",batch_d2$colnames.count_matrix_2.)=="TRUE", print("4"),batch_d2$batch) 
batch_d2$batch <- ifelse(grepl("Sample5",batch_d2$colnames.count_matrix_2.)=="TRUE", print("5"),batch_d2$batch) 
batch_d2$batch <- ifelse(grepl("Sample6",batch_d2$colnames.count_matrix_2.)=="TRUE", print("6"),batch_d2$batch) 


##----include data about sex
batch_d2$sex <- ifelse(grepl("Sample1|Sample3|Sample6",batch_d2$colnames.count_matrix_2.)=="TRUE", print("1"), print("2")) 

##------include data about treatment
batch_d2$treatment <- ifelse(grepl("uninfected",batch_d2$colnames.count_matrix_2.)=="TRUE", print("1"),batch_d2$colnames.count_matrix_2.) 
batch_d2$treatment <- ifelse(grepl("UK_72hpi",batch_d2$colnames.count_matrix_2.)=="TRUE", print("2"),batch_d2$treatment) 
batch_d2$treatment <- ifelse(grepl("VIC01_72hpi",batch_d2$colnames.count_matrix_2.)=="TRUE", print("3"),batch_d2$treatment) 
batch_d2$treatment <- ifelse(grepl("VIC01_48hpi",batch_d2$colnames.count_matrix_2.)=="TRUE", print("4"),batch_d2$treatment) 

##-----include data about age

batch_d2$age <- ifelse(grepl("adult",batch_d2$colnames.count_matrix_2.)=="TRUE", print("1"),print("2")) 
#batch_d2$age <- ifelse(grepl("child",batch_d2$colnames.count_matrix_2.)=="TRUE", print("2"),batch_d2$inf) 



##---include data about infection status

batch_d2$inf <- ifelse(grepl("Infected",batch_d2$colnames.count_matrix_2.)=="TRUE", print("1"),batch_d2$colnames.count_matrix_2.) 
batch_d2$inf <- ifelse(grepl("Uninfected",batch_d2$colnames.count_matrix_2.)=="TRUE", print("2"),batch_d2$inf) 

##---include data about cell-type

#clusters
batch_d2$cluster <- ifelse(grepl("^Basal-1",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Basal-1"),batch_d2$colnames.count_matrix_2.) 
batch_d2$cluster <- ifelse(grepl("^Basal-2",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Basal-2"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("^Suprabasal",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Suprabasal"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("^Brush/Tuft",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Brush/Tuft"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("^Ciliated-1",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Ciliated-1"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("^Ciliated-2",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Ciliated-2"),batch_d2$cluster)
batch_d2$cluster <- ifelse(grepl("^Ciliated-3",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Ciliated-3"),batch_d2$cluster)
batch_d2$cluster <- ifelse(grepl("^Cycling-basal",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Cycling-basal"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Deuterosomal",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Deuterosomal"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Goblet",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Goblet"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Ionocyte",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Ionocyte"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Secretory-1",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Secretory-1"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Secretory-2",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Secretory-2"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Secretory-3",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Secretory-3"),batch_d2$cluster)
batch_d2$cluster <- ifelse(grepl("Secretory-Ciliated",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Secretory-Ciliated"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Ciliated-SC2",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Ciliated-SC2"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Secretory-Ciliated-SC2",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Secretory-Ciliated-SC2"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Goblet/Ionocyte",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Goblet/Ionocyte"),batch_d2$cluster) 
batch_d2$cluster <- ifelse(grepl("Goblet-IFN-stim",batch_d2$colnames.count_matrix_2.)=="TRUE", print("Goblet-IFN-stim"),batch_d2$cluster) 

#--------------------
##check with >=15 cells
list <- batch_d2$colnames.count_matrix_2.
list <- as.data.frame(list)

df = data.frame()
for (i in list$list) {
  output <-  sub("_Sample.*", "", i)
  # Using rbind() to append the output of one iteration to the dataframe
  df = rbind(df, output)
}



##DE unadjusted Differential expression: limma-voom
col <- as.data.frame(colnames(count_matrix_2))
##check lib size per sample
count_matrix_2 <- count_matrix[, new_count_list_2] #######15 cells!!!!!
batch_d2$cluster <- factor(batch_d2$cluster)
batch_d2$o <- order(batch_d2$cluster)
#count_matrix_2$o <- batch_d2
#col.concent <- c("gold" , "orange", "red", "blue")
barplot((colSums(count_matrix_2)*1e-6), names = batch_d2$colnames.count_matrix_2., ylab="Library size (millions)",
                                                         col = batch_d2$o, las = 2, cex.names = 0.3)
#legend("topleft", legend = c("1.5","2.5","3.5","4.5","5.5","6.5","7.5"), col = col.concent, pch = 15, cex = 0.7)

group <- factor(paste0(df$X.Basal.1_UK_72hpi_adult_Infected.))
group <- gsub("-", "_", group)
group <- gsub("Brush/Tuft", "Brush_Tuft",group)
group <- gsub("Goblet/Ionocyte", "Goblet_Ionocyte",group)
batch_d2$group <- group

#DE-------------------------------------
dge <- DGEList(counts=count_matrix_2,group=group)
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)


##----include data about donor
dge[["samples"]]$batch <- ifelse(grepl("Sample1",rownames(dge[["samples"]]))=="TRUE", print("Adult1"),rownames(dge[["samples"]])) 
dge[["samples"]]$batch <- ifelse(grepl("Sample2",rownames(dge[["samples"]]))=="TRUE", print("Adult2"),dge[["samples"]]$batch) 
dge[["samples"]]$batch <- ifelse(grepl("Sample3",rownames(dge[["samples"]]))=="TRUE", print("Adult3"),dge[["samples"]]$batch) 
dge[["samples"]]$batch <- ifelse(grepl("Sample4",rownames(dge[["samples"]]))=="TRUE", print("Child1"),dge[["samples"]]$batch) 
dge[["samples"]]$batch <- ifelse(grepl("Sample5",rownames(dge[["samples"]]))=="TRUE", print("Child2"),dge[["samples"]]$batch) 
dge[["samples"]]$batch <- ifelse(grepl("Sample6",rownames(dge[["samples"]]))=="TRUE", print("Child3"),dge[["samples"]]$batch) 


##----include data about sex
dge[["samples"]]$sex <- ifelse(grepl("Sample1|Sample3|Sample6",rownames(dge[["samples"]]))=="TRUE", print("Male"), print("Female")) 

##------include data about treatment
dge[["samples"]]$treatment <- ifelse(grepl("uninfected",rownames(dge[["samples"]]))=="TRUE", print("uninfected"),rownames(dge[["samples"]])) 
dge[["samples"]]$treatment <- ifelse(grepl("UK_72hpi",rownames(dge[["samples"]]))=="TRUE", print("2"),dge[["samples"]]$treatment) 
dge[["samples"]]$treatment <- ifelse(grepl("VIC01_72hpi",rownames(dge[["samples"]]))=="TRUE", print("3"),dge[["samples"]]$treatment) 
dge[["samples"]]$treatment <- ifelse(grepl("VIC01_48hpi",rownames(dge[["samples"]]))=="TRUE", print("4"),dge[["samples"]]$treatment) 

##---include data about infection status

dge[["samples"]]$inf <- ifelse(grepl("Infected",rownames(dge[["samples"]]))=="TRUE", print("1"),rownames(dge[["samples"]])) 
dge[["samples"]]$inf <- ifelse(grepl("Uninfected",rownames(dge[["samples"]]))=="TRUE", print("uninfected"),dge[["samples"]]$inf) 


##----------------include group

dge[["samples"]]$g <- ifelse(grepl("Sample6",rownames(dge[["samples"]]))=="TRUE", print("Donor6"), print("Donor1_5")) 


dge[["samples"]]$finf <- paste(dge[["samples"]]$g,dge[["samples"]]$treatment,dge[["samples"]]$inf, sep="_")


finf <- dge[["samples"]]$finf

#sex <- batch_d2$sex
sex <- dge[["samples"]]$sex
batch <- dge[["samples"]]$batch
design2 <- model.matrix(~0 + finf + sex)

v <- voom(dge, design2, plot=TRUE)
#plotMDS(v, labels=1:300, col=as.numeric(group))

##random effect
corfit <- duplicateCorrelation(v,design2,block=as.factor(batch))
corfit$consensus


#dge <- estimateDisp(dge, design2, robust=TRUE)
#plotBCV(dge)




#Likelihood ratio test

my.contrasts <- makeContrasts(
"donor_6_vs_all" =finfDonor6_uninfected_uninfected-finfDonor1_5_uninfected_uninfected,

levels = design2)


contrastlist <- c("donor_6_vs_all")



##Likelihood ratio + Withsex
fit <- lmFit(v,design2,block=batch,correlation=corfit$consensus)
#fit <- lmFit(v, design)
#fit <- lmFit(v, design2)
#fit <- eBayes(fit)


#topTable(fit, coef=ncol(design2))

fun <- function(i){
  
  fit2 <- contrasts.fit(fit,contrast=my.contrasts[,i])
  fit2 <- eBayes(fit2)
  contr1_withsex_table <- topTable(fit2, adjust.method="BH", sort.by="p",n = Inf)
  length(which(contr1_withsex_table$adj.P.Val < 0.05))
  write.csv(contr1_withsex_table, file = paste0(i,".csv"), row.names = TRUE)
} 


lapply(contrastlist,fun)







##-----------------run 2


dge[["samples"]]$g2 <- ifelse(grepl("Sample6",rownames(dge[["samples"]]))=="TRUE", print("Donor6"), print("Child")) 
dge[["samples"]]$g2 <- ifelse(grepl("Sample1|Sample2|Sample3",rownames(dge[["samples"]]))=="TRUE", print("Adult"), dge[["samples"]]$g2) 

dge[["samples"]]$fin2 <- paste(dge[["samples"]]$g2,dge[["samples"]]$treatment,dge[["samples"]]$inf, sep="_")

fin2 <- dge[["samples"]]$fin2 

#sex <- batch_d2$sex
sex <- dge[["samples"]]$sex
design2 <- model.matrix(~0 + fin2 + sex)

v <- voom(dge, design2, plot=TRUE)
#plotMDS(v, labels=1:300, col=as.numeric(group))

##random effect
corfit <- duplicateCorrelation(v,design2,block=as.factor(batch))
corfit$consensus
#Likelihood ratio test


my.contrasts <- makeContrasts(
  "donor_6_vs_4_5" =fin2Donor6_uninfected_uninfected-fin2Child_uninfected_uninfected,
  
  levels = design2)


contrastlist <- c("donor_6_vs_4_5")



##Likelihood ratio + Withsex
fit <- lmFit(v,design2,block=batch,correlation=corfit$consensus)
#fit <- lmFit(v, design)
#fit <- lmFit(v, design2)
#fit <- eBayes(fit)


#topTable(fit, coef=ncol(design2))

fun <- function(i){
  
  fit2 <- contrasts.fit(fit,contrast=my.contrasts[,i])
  fit2 <- eBayes(fit2)
  contr1_withsex_table <- topTable(fit2, adjust.method="BH", sort.by="p",n = Inf)
  length(which(contr1_withsex_table$adj.P.Val < 0.05))
  write.csv(contr1_withsex_table, file = paste0(i,".csv"), row.names = TRUE)
} 


lapply(contrastlist,fun)





