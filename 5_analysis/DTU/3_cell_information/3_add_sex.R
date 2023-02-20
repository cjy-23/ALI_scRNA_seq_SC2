#setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/")
setwd("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/")

scp_meta_2 <- read.delim("/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/cell_information/scp_meta_2.txt")

scp_meta_2$Sex <- ifelse((grepl("Sample1_Sample1|Sample3_Sample3|Sample6_Sample6",scp_meta_2$Patient)=="TRUE"), print("Male"), print("Female"))

write.table(scp_meta_2, file="scp_meta_3.txt")

scp_meta_2$Cell_Type3 <- ifelse((grepl("^Secretory",scp_meta_2$Cell_Type)=="TRUE"), "Secretory", scp_meta_2$Cell_Type)
scp_meta_2$Cell_Type3 <- ifelse((grepl("Secretory_Ciliated",scp_meta_2$Cell_Type)=="TRUE"), "Secretory_Ciliated", scp_meta_2$Cell_Type3)
scp_meta_2$Cell_Type3 <- ifelse((grepl("^Basal",scp_meta_2$Cell_Type)=="TRUE"), "Basal", scp_meta_2$Cell_Type3)
scp_meta_2$Cell_Type3 <- ifelse((grepl("^Ciliated",scp_meta_2$Cell_Type)=="TRUE"), "Ciliated", scp_meta_2$Cell_Type3)
scp_meta_2$Cell_Type3 <- ifelse((grepl("^Goblet",scp_meta_2$Cell_Type)=="TRUE"), "Goblet", scp_meta_2$Cell_Type3)
scp_meta_2$Cell_Type3 <- ifelse((grepl("^Goblet_Ionocyte",scp_meta_2$Cell_Type)=="TRUE"), "Goblet_Ionocyte", scp_meta_2$Cell_Type3)

scp_meta_2$group <- paste(scp_meta_2$Cell_Type3, scp_meta_2$Cohort, scp_meta_2$Cell_State, sep="_")


write.table(scp_meta_2, file="scp_meta_2_newcell_2.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) ## with celltype3
