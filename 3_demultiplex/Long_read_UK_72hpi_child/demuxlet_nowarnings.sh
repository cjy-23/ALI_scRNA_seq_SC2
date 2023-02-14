#!/bin/bash
#SBATCH --job-name=demuxlet_long_child_UK    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=4                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=100gb
#SBATCH --time=20:00:00               # Time limit hrs:min:sec
#SBATCH --partition=physical

echo "Running minimap on "`hostname`
 
#INPUT YOUR COMMAND HERE
module load gcc/8.3.0 
module load demuxlet/11022021

#shortread uninfected child 



demuxlet --sam "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/run/Long_read_UK_72hpi_child/outs/possorted_genome_bam_hg19.sorted.bam" --tag-group CB --tag-UMI UB --vcf "/data/gpfs/projects/punim1466/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr_no0_lex.vcf.gz" --field GT --geno-error 0.01 --out Long_read_UK_72hpi_child_nowarning --group-list "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/run/Long_read_UK_72hpi_child/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"



echo "Complete!"


