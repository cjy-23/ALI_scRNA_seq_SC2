#!/bin/bash
#SBATCH --job-name=flames_2_libF   # Job name
#SBATCH --partition=physical
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=500GB
#SBATCH --time=30-00:00:0               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE
#module load gcc/8.3.0
module load minimap2/2.17
#module load samtools/1.9
module load gcc/11.2.0

#usage: FLAMES [-h] -a GFF3 [-i INFQ] [-b INBAM] --outdir OUTDIR --genomefa
#             GENOMEFA [--minimap2_dir MINIMAP2_DIR]
#              [--config_file CONFIG_FILE]
#             [--downsample_ratio DOWNSAMPLE_RATIO]

"/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/FLAMES/python/sc_long_pipeline.py" -a "/data/gpfs/projects/punim1466/db/refdata-gex-GRCh38-2020-A/genes/genes.gtf" -i "/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/libF/F_match_cell_barcode_output" --genomefa "/data/gpfs/projects/punim1466/db/refdata-gex-GRCh38-2020-A/fasta/genome.fa" --outdir "/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/libF/" --config_file "/data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames_2/FLAMES/python/config_sclr_nanopore_default.json"


#minimap2 -ax splice -t 12 --junc-bed /data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames/libF_2/tmp.splice_anno.bed12 --junc-bonus 1 --splice-flank=no -k14 --secondary=no /data/gpfs/projects/punim1466/db/refdata-gex-GRCh38-2020-A/fasta/genome.fa /data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames/libF/H_match_cell_barcode_output | samtools view -bS -@ 4 -m 2G -o /data/gpfs/projects/punim1466/analysis/cellranger_hostonly_trim/flames/libF/tmp.align.bam
echo "Complete!"

