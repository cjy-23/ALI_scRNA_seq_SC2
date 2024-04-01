setwd("analysis/velocity/")
library(purrr)

#Load all the data
barcodes_Long_read_UK_72hpi_adult <- read.csv("/velocity/barcodes_Long_read_UK_72hpi_adult.tsv", sep="")
barcodes_Long_read_uninfected_adult <- read.csv("/velocity/barcodes_Long_read_uninfected_adult.tsv", sep="")
barcodes_Long_read_VIC01_72hpi_adult <- read.csv("/velocity/barcodes_Long_read_VIC01_72hpi_adult.tsv", sep="")
barcodes_Long_read_VIC01_48hpi_adult <- read.csv("/velocity/barcodes_Long_read_VIC01_48hpi_adult.tsv", sep="")

barcodes_Short_read_UK_72hpi_adult <- read.csv("/velocity/barcodes_Short_read_UK_72hpi_adult.tsv", sep="")
barcodes_Short_read_uninfected_adult <- read.csv("/velocity/barcodes_Short_read_uninfected_adult.tsv", sep="")
barcodes_Short_read_VIC01_72hpi_adult <- read.csv("/velocity/barcodes_Short_read_VIC01_72hpi_adult.tsv", sep="")
barcodes_Short_read_VIC01_48hpi_adult <- read.csv("/velocity/barcodes_Short_read_VIC01_48hpi_adult.tsv", sep="")


barcodes_Long_read_UK_72hpi_child <- read.csv("/velocity/barcodes_Long_read_UK_72hpi_child.tsv", sep="")
barcodes_Long_read_uninfected_child <- read.csv("/velocity/barcodes_Long_read_uninfected_child.tsv", sep="")
barcodes_Long_read_VIC01_72hpi_child <- read.csv("/velocity/barcodes_Long_read_VIC01_72hpi_child.tsv", sep="")
barcodes_Long_read_VIC01_48hpi_child <- read.csv("/velocity/barcodes_Long_read_VIC01_48hpi_child.tsv", sep="")

barcodes_Short_read_UK_72hpi_child <- read.csv("/velocity/barcodes_Short_read_UK_72hpi_child.tsv", sep="")
barcodes_Short_read_uninfected_child <- read.csv("/velocity/barcodes_Short_read_uninfected_child.tsv", sep="")
barcodes_Short_read_VIC01_72hpi_child <- read.csv("/velocity/barcodes_Short_read_VIC01_72hpi_child.tsv", sep="")
barcodes_Short_read_VIC01_48hpi_child <- read.csv("/velocity/barcodes_Short_read_VIC01_48hpi_child.tsv", sep="")


##Check and change the suffices

nm <- c("barcodes_Long_read_UK_72hpi_adult",
        "barcodes_Long_read_uninfected_adult",
        "barcodes_Long_read_VIC01_72hpi_adult",
        "barcodes_Long_read_VIC01_48hpi_adult",
        
        "barcodes_Short_read_UK_72hpi_adult",
        "barcodes_Short_read_uninfected_adult",
        "barcodes_Short_read_VIC01_72hpi_adult",
        "barcodes_Short_read_VIC01_48hpi_adult",
        
        
        "barcodes_Long_read_UK_72hpi_child",
        "barcodes_Long_read_uninfected_child",
        "barcodes_Long_read_VIC01_72hpi_child",
        "barcodes_Long_read_VIC01_48hpi_child",
        
        "barcodes_Short_read_UK_72hpi_child",
        "barcodes_Short_read_uninfected_child",
        "barcodes_Short_read_VIC01_72hpi_child",
        "barcodes_Short_read_VIC01_48hpi_child" 
)



list <- list(barcodes_Long_read_UK_72hpi_adult,
        barcodes_Long_read_uninfected_adult,
        barcodes_Long_read_VIC01_72hpi_adult,
        barcodes_Long_read_VIC01_48hpi_adult,
        
        barcodes_Short_read_UK_72hpi_adult,
        barcodes_Short_read_uninfected_adult,
        barcodes_Short_read_VIC01_72hpi_adult,
        barcodes_Short_read_VIC01_48hpi_adult,
        
        
        barcodes_Long_read_UK_72hpi_child,
        barcodes_Long_read_uninfected_child,
        barcodes_Long_read_VIC01_72hpi_child,
        barcodes_Long_read_VIC01_48hpi_child,
        
        barcodes_Short_read_UK_72hpi_child,
        barcodes_Short_read_uninfected_child,
        barcodes_Short_read_VIC01_72hpi_child,
        barcodes_Short_read_VIC01_48hpi_child 
)

fun1 <- function(y) {
  
  y$x <- gsub("-.*", "", y$x)
  y$x <- paste0(y$x, "-1", sep = "")
  y
}


list2 <- lapply(list,fun1)
names(list2) <- nm


##Save the new barcodes


fun2 <- function(subset, name) {
  filename <- paste0("barcodes_fixed_", name, ".tsv")
  cat("Writing file:", filename, "\n")
  write.table(subset$x, file = filename, sep = "\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
}
walk2(list2, names(list2), fun2)

