#!/bin/bash
#SBATCH --job-name=demuxlet_long_adult_UK    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=4                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=100gb
#SBATCH --time=30:00:00               # Time limit hrs:min:sec
#SBATCH --partition=physical

echo "Running minimap on "`hostname`
 
#INPUT YOUR COMMAND HERE
module load gcc/8.3.0 
module load demuxlet/11022021

#shortread uninfected adult 



demuxlet --sam "/analysis/cellranger_new_run/run/Short_read_uninfected_adult/outs/possorted_genome_bam_hg19.sorted.bam" --tag-group CB --tag-UMI UB --vcf "/genotyping/adult/SNP_genotyping/plink_input/no_warnings/newfile_clean_3_chr_no0_lex.vcf.gz" --field GT --geno-error 0.01 --out Short_read_uninfected_adult_nowarning --group-list "/analysis/cellranger_new_run/run/Short_read_uninfected_adult/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"



echo "Complete!"


