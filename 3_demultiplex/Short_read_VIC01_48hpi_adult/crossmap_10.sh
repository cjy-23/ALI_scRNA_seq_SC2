#!/bin/bash
#SBATCH --job-name=crossmap_long_adult_UK    # Job name
#SBATCH --partition=physical,snowy
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=2                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=120GB
#SBATCH --time=30:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE

module load python/3.8.2


CrossMap.py bam -a "/data/gpfs/projects/punim1466/sw/hg38ToHg19.over.chain.gz" "/data/scratch/projects/punim1068/data/basecall_guppy_3.5.2/cellranger/Short_read_VIC01_48hpi_adult/outs/possorted_genome_bam.bam" possorted_genome_bam_hg19


echo "Complete!"


