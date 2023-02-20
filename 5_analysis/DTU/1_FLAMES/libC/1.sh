#!/bin/bash
#SBATCH --job-name=flames_libC   # Job name
#SBATCH --partition=physical
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=100GB
#SBATCH --time=100:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE
#module load gcc/8.3.0
module load minimap2/2.17
#module load samtools/1.9
module load gcc/11.2.0

#usage: FLAMES [-h] -a GFF3 [-i INFQ] [-b INBAM] --outdir OUTDIR --genomefa
#             GENOMEFA [--minimap2_dir MINIMAP2_DIR]
#              [--config_file CONFIG_FILE]
#             [--downsample_ratio DOWNSAMPLE_RATIO]

#"/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames/flames/python/sc_long_pipeline.py" -a "/data/gpfs/projects/punim1466/db/refdata-gex-GRCh38-2020-A/genes/genes.gtf" -i "/data/scratch/projects/punim1068/ONT_promethion_scRNAseq/libA/fastq_5.0.11/merged_all_libA.fastq" --genomefa "/data/gpfs/projects/punim1466/db/refdata-gex-GRCh38-2020-A/fasta/genome.fa" --outdir "/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames/libA/"


#usage: <1.fastq folder> <2.output cell barcode statistics file> <3.fastq output reads that matched cell barcode> <4.barcode whitelist> <5.max edit distance> [6. UMI length (default: 10)]
"/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/FLAMES/src/bin/match_cell_barcode" "/data/scratch/projects/punim1068/ONT_promethion_scRNAseq/libC/" libC_match_cell_barcode C_match_cell_barcode_output "/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/libC/barcodes.tsv" 2 12
echo "Complete!"

