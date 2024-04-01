#!/bin/bash
#SBATCH --job-name=run_velocyto  # Job name
#SBATCH --partition=cascade
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=4                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=150GB
#SBATCH --time=100:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE
#module load R/4.2.1
ml SAMtools
ml velocyto/0.17.17

velocyto run --bcfile "/barcodes_fixed/barcodes_fixed_barcodes_Long_read_UK_72hpi_adult.tsv.gz" "/velocyto_filtered/Long_read_UK_72hpi_adult/possorted_genome_bam.bam" "/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

echo "Complete!"


