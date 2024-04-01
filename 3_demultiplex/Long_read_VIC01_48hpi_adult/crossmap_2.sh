#!/bin/bash
#SBATCH --job-name=crossmap_long_adult_UK    # Job name
#SBATCH --partition=physical,snowy
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=2                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=120GB
#SBATCH --time=30:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE

module load python/3.8.2


CrossMap.py bam -a "/sw/hg38ToHg19.over.chain.gz" "/data/basecall_guppy_3.5.2/cellranger/Long_read_VIC01_48hpi_adult/outs/possorted_genome_bam.bam" possorted_genome_bam_hg19


echo "Complete!"


