#!/bin/bash
#SBATCH --job-name=bcl2fastq    # Job name
#SBATCH --partition=physical,snowy
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=8                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=120GB
#SBATCH --time=500:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE

cd "/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/"

source /usr/local/module/spartan_old.sh
module load bcl2fastq2/2-20-0

cellranger mkfastq --id=novaseq_cellranger --run="/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/" --csv="/sw/cellranger/cellranger-tiny-bcl-simple-1.2.0.csv"
echo "Complete!"


