#setwd("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/")


##Test
dtu_table <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Inf_vs_uninfected/secretory-ciliated_Alpha_vs_uninf_adult_save/dtu_table.tsv")
FDR_table <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Inf_vs_uninfected/secretory-ciliated_Alpha_vs_uninf_adult_save/FDR_table.tsv")

genes_dtu <- dtu_table[dtu_table$gene_qvalue<0.5,]
genes_FDR <- FDR_table[FDR_table$gene<0.5,]
genes_FDR <- unique(genes_FDR$geneID)
genes_FDR <- as.data.frame(genes_FDR)

setequal(genes_dtu$gene_ID, genes_FDR$genes_FDR) ##should be the same


#==================== INF VS UNINFECTED 
#Load data
FDR_table1 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Inf_vs_uninfected/secretory-ciliated_Alpha_vs_uninf_adult_save/FDR_table.tsv")
FDR_table2 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Inf_vs_uninfected/secretory-ciliated_Alpha_vs_uninf_child_save/FDR_table.tsv")
FDR_table3 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Inf_vs_uninfected/secretory-ciliated_Alpha_vs_WT_child_save/FDR_table.tsv")
FDR_table4 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Inf_vs_uninfected/secretory-ciliated_WT_vs_uninf_child_save/FDR_table.tsv")


merge <- append(FDR_table$geneID, FDR_table2$geneID)
merge <- append(merge, FDR_table3$geneID)
merge <- append(merge, FDR_table4$geneID)

unimerge <- unique(merge)


write.table(unimerge, file="inf_vs_uninfected_0.5_firststage_background_ENST.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)


## ================ UNINFECTED

#Load data

FDR_table5 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Uninfected/Basal_save/FDR_table.tsv")
FDR_table6 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Uninfected/Basal_secretory_save/FDR_table.tsv")
FDR_table7 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Uninfected/Cycling_basal_save/FDR_table.tsv")
FDR_table8 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Uninfected/Deuterosomal_save/FDR_table.tsv")
FDR_table9 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Uninfected/Secretory_Ciliated_save/FDR_table.tsv")
FDR_table10 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Uninfected/Secretory_save/FDR_table.tsv")


merge <- append(FDR_table5$geneID, FDR_table6$geneID)
merge <- append(merge, FDR_table7$geneID)
merge <- append(merge, FDR_table8$geneID)
merge <- append(merge, FDR_table9$geneID)
merge <- append(merge, FDR_table10$geneID)

unimerge <- unique(merge)


write.table(unimerge, file="uninfected_0.5_firststage_background_ENST.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

##===================== Bystander vs control

FDR_table11 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/ciliated/ciliated_Alpha_adult_save/FDR_table.tsv")
  FDR_table12 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/ciliated/ciliated_WT_48hpi_adult_save/FDR_table.tsv")
    FDR_table13 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/ciliated/ciliated_WT_72hpi_adult_save/FDR_table.tsv")
      FDR_table14 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory/secretory_Alpha_adult_save/FDR_table.tsv")
        FDR_table15 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory/secretory_Alpha_child_save/FDR_table.tsv")
          FDR_table16 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory/secretory_WT_48hpi_adult_save/FDR_table.tsv")
            FDR_table17 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory/secretory_WT_48hpi_child_save/FDR_table.tsv")
              FDR_table18 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory/secretory_WT_72hpi_adult_save/FDR_table.tsv")
                FDR_table19 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory/secretory_WT_72hpi_child_save/FDR_table.tsv")
                  FDR_table20 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory-ciliated/secretory-ciliated_Alpha_adult_save/FDR_table.tsv")
                    FDR_table21 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory-ciliated/secretory-ciliated_Alpha_child_save/FDR_table.tsv")
                      FDR_table22 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory-ciliated/secretory-ciliated_WT_48hpi_adult_save/FDR_table.tsv")
                        FDR_table23 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory-ciliated/secretory-ciliated_WT_48hpi_child_save/FDR_table.tsv")
                          FDR_table24 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory-ciliated/secretory-ciliated_WT_72hpi_adult_save/FDR_table.tsv")
                            FDR_table25 <- read.delim("/data/gpfs/projects/punim1466/analysis/nanopore/isoform/DTU_new_no_sat_filt_0.5_ENST/Bystander_vs_control/secretory-ciliated/secretory-ciliated_WT_72hpi_child_save/FDR_table.tsv")
                            

merge <- append(FDR_table11$geneID, FDR_table12$geneID)
merge <- append(merge, FDR_table13$geneID)
merge <- append(merge, FDR_table14$geneID)
merge <- append(merge, FDR_table15$geneID)
merge <- append(merge, FDR_table16$geneID)
merge <- append(merge, FDR_table17$geneID)
merge <- append(merge, FDR_table18$geneID)
merge <- append(merge, FDR_table19$geneID)
merge <- append(merge, FDR_table20$geneID)
merge <- append(merge, FDR_table21$geneID)
merge <- append(merge, FDR_table22$geneID)
merge <- append(merge, FDR_table23$geneID)
merge <- append(merge, FDR_table24$geneID)
merge <- append(merge, FDR_table25$geneID)

                            
unimerge <- unique(merge)
                            
write.table(unimerge, file="bystander_vs_control_0.5_firststage_background_ENST.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)












