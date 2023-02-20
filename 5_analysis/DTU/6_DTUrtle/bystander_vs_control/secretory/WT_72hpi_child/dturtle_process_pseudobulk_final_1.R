#setwd("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt/Inf_vs_uninfected/secretory-Basal_Alpha_vs_uninf_child_save/")
#if(!requireNamespace("remotes", quietly = TRUE)){
#  install.packages("remotes")
#}
##remotes::install_github("TobiTekath/DTUrtle")
library("DTUrtle")
library("BiocParallel")
library("GenomicRanges") 
library("Gviz")
library("rtracklayer") 
library("stageR")
library("tximport")
#library("DESeq2")
library("stringr")
library("tidyr")
library("data.table")   
library("Matrix")
library("plyr")
library("dplyr")
library(stringr)
#library(Seurat)
library("remotes")
#require(devtools)
#BiocManager::install(c("BiocParallel", "GenomicRanges", "Gviz", "rtracklayer", "stageR", "tximport"))
#setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/DTU_15/inf_vs_uninf/Basal_Alpha_uninf_child/")
library("SingleCellExperiment")
library(SummarizedExperiment)
library(methods) 
library(Matrix.utils)
biocpar <- BiocParallel::MulticoreParam(1)
library("magrittr")
library("purrr")
library("ggplot2")
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
#remotes::install_github("traversc/trqwe")
#BiocManager::valid()

all_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/all_gtf.tsv", header=FALSE)
#write.table(all_gtf, file="all_gtf.tsv", sep="\t", row.names = FALSE, col.names = FALSE)
all_gtf$transcript_id <- gsub("transcript_id transcript:", "",all_gtf$V9)
all_gtf$transcript_id <- gsub(";.*", "",all_gtf$transcript_id)
all_gtf$gene_id <- gsub(".*gene_id ","",all_gtf$V9)
all_gtf$gene_id <- gsub(";","",all_gtf$gene_id)

txInfo <- all_gtf[,10:11]
txInfo <- unique(txInfo)
tx2gene <- txInfo

##Fix sampledata
sampledata <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/scp_meta_2_newcell_2.tsv")
colnames(sampledata)[1] <- "id"
sampledata$group <- paste(sampledata$Cohort, sampledata$Cell_State,sampledata$Patient ,sep="_")
sampledata$Cell_Type3 <- sub("/", "-", sampledata$Cell_Type3)
sampledata$Cell_Type3 <- as.factor(sampledata$Cell_Type3)
sampledata$group <- as.factor(sampledata$group)
sampledata$full_id <- paste(sampledata$Cell_Type3,sampledata$Cohort, sampledata$Cell_State,sampledata$Patient ,sep="_")
sampledata$full_id <- as.factor(sampledata$full_id)
sampledata$full_group <- paste(sampledata$Cell_Type3,sampledata$Cohort, sampledata$Cell_State,sep="_")
sampledata$full_group <- as.factor(sampledata$full_group)
load("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/DTU_15/mergelib_0.5_ENST.Rdata")

# Create single cell experiment object

sce <- SingleCellExperiment(assays=list(counts = mergelib2), 
                           colData = sampledata)

##TRY THIS LATER (REMOVE ANY TRANSCRIPTS WITH LESS THAN 10 CELLS WITH ANY COUNTS)
#sce <- sce[rowSums(counts(sce) > 1) >= 15, ]
#sce <- SummarizedExperiment(counts = Tasic_counts_vignette, 
 #                           colData = sampledata)


# Named vector of cluster names
kids <- purrr::set_names(levels(sce$Cell_Type3))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$full_id))
# Total number of samples 
ns <- length(sids)

ns

## Keep data with only 15 cells
# Generate sample level metadata

## Determine the number of cells per sample
ohh <- table(sce$full_id)
write.table(ohh, file="full_id_cells.tsv", sep="\t", col.names=TRUE, row.names=FALSE)
fg <- as.data.frame(table(sce$full_group))
write.table(fg, file="full_group_cells.tsv", sep="\t", col.names=TRUE, row.names=FALSE)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$full_id))
## Determine how to reorder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$full_id)
## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"id")
ei


save(ei, file="ei.RData")

#Keep samples with greater than equal to 15 cells
new_count_list <- ei[ei$n_cells >=15,]
new_count_list_2 <- new_count_list$full_id


# Identify groups for aggregation of counts
groups <- colData(sce)[, c("full_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

#save(pb, file="pb.Rdata")

count_matrix <- as.matrix(t(pb))


##subset with samples with at least 5 cells

#count_matrix_2 <- count_matrix[, new_count_list_2] 
save(count_matrix, file="count_matrix.Rdata")

colnames(new_count_list)[9] <- "id"
colnames(new_count_list)[7] <- "cluster_id"
colnames(new_count_list)[4] <- "Sample"
new_count_list$lol <- paste(new_count_list$cluster_id, new_count_list$Cell_State, sep="_")


write.table(new_count_list, "new_count_list.tsv", sep="\t", col.names=TRUE, row.names=TRUE)
#gg_df <- left_join(gg_df, ei[, c("id", "cluster_id", "Sample")]) 

p <-ggplot(data=new_count_list, aes(x=id, y=n_cells, fill=lol, pattern = Cell_State)) +
  geom_bar(stat="identity",color="black") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, face="bold"),
        axis.title = element_text(face="bold")) +
  geom_col_pattern(
    aes(id, n_cells, pattern_fill = Cell_State), 
    pattern_fill    = 'black',
    pattern_colour  = 'black',
    colour          = 'black', 
    pattern_density = 0.35
    )

pdf(file="Inf_uninft.pdf", width=20, height=15)
p
dev.off()


#==========================================
#create a sample data sheet, specifying which sample / cell belongs to which group
newsampledata <- new_count_list[,c(9:10)]
colnames(newsampledata) <- c("id", "group")


#use DRIMSeq for fitting a Dirichlet-multinomial model --- appropriate parameters are chosen based on the selected filtering_strategy ('bulk' or 'sc')
dturtle <- run_drimseq(counts = count_matrix, tx2gene = tx2gene, pd=newsampledata, id_col = "id",
                       cond_col = "group", filtering_strategy = "bulk", 
                       #BPPARAM = biocpar, cond_levels=cond_col[144]-cond_col[159])
                       #BPPARAM = biocpar, 
                       cond_levels=c("Secretory_VIC01_72hpi_child_Uninfected",
                                     "Secretory_uninfected_child_Uninfected"))
#run posthoc filtering and two-staged statistical correction with stageR
dturtle <- posthoc_and_stager(dturtle = dturtle, ofdr = 0.05)


##Flag priming

#Attention: calculation for all available genes might be time consuming.
#all_gtf <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/all_gtf.tsv", header=FALSE)

#all_gtf <- import_gtf(gtf_file = "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/all_gtf.gtf", feature_type = NULL)

#all_gtf$V4 <- as.numeric(all_gtf$V4)
#all_gtf$V5 <- as.numeric(all_gtf$V5)
#colnames(all_gtf) <- c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#colnames(all_gtf)[4:5] <- c("start", "end")
#priming_bias_df <- priming_bias_detection_probability(counts = count_matrix, gtf = all_gtf, tx2gene = tx2gene, one_to_one = NULL,
 #                                                 priming_enrichment = "3", BPPARAM = biocpar)



#highly flexible function to create a results data frame
dturtle <- create_dtu_table(dturtle = dturtle)

## View results data frame
#View(dturtle$dtu_table)

dtu_table <- dturtle$dtu_table
save(dtu_table, file="dtu_table.Rdata")
save(dturtle, file="dturtle.Rdata")
write.table(dtu_table, file="dtu_table.csv", sep=",",col.names=TRUE, row.names=FALSE)

write.table(dturtle[["FDR_table"]], file="FDR_table.tsv", sep="\t", col.names=TRUE)
dtu_table2 <- dtu_table[dtu_table$minimal_tx_qvalue<0.05,]
dtu_table2 <- dtu_table2[order(dtu_table2$minimal_tx_qvalue),]
write.table(dtu_table2, file="dtu_table2.tsv", sep="\t", col.names=TRUE, row.names=TRUE)
ID <- dtu_table2$gene_ID
top10 <- ID[1:10]

#change to results folder
#setwd("my_results_folder")    

#create plots, save them to disk and link them in the `dtu_table`.
dturtle <- plot_proportion_barplot(dturtle = dturtle, 
                                   savepath = "images", 
                                   add_to_table = "barplot",
                                   width=15,
                                   height=10,
                                   fit_line_color = "black",
                                   genes=top10)

dturtle <- plot_proportion_pheatmap(dturtle = dturtle, 
                                    savepath = "images", 
                                    include_expression = TRUE,
                                    add_to_table = "pheatmap",
                                    width=15,
                                   height=10,
                                    genes=top10)

