#!/bin/bash
#SBATCH --job-name=iseq_bcl2fastq    # Job name
#SBATCH --partition=physical,snowy
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=8                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=120GB
#SBATCH --time=500:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE


source /usr/local/module/spartan_old.sh
module load bcl2fastq2/2-20-0


##Short read adult
"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_UK_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_UK_72hpi_adult --r1-length 28 --r2-length 90 --expect-cells 10950

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_uninfected_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_uninfected_adult --r1-length 28 --r2-length 90 --expect-cells 10950

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_48hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_VIC01_48hpi_adult --r1-length 28 --r2-length 90 --expect-cells 10950 --no-bam

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_VIC01_72hpi_adult --r1-length 28 --r2-length 90 --expect-cells 10950

##Short read child

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_UK_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_UK_72hpi_child --r1-length 28 --r2-length 90 --expect-cells 10950

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_uninfected_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_uninfected_child --r1-length 28 --r2-length 90 --expect-cells 10950 --no-bam


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_48hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_VIC01_48hpi_child --r1-length 28 --r2-length 90 --expect-cells 10950 --no-bam

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Short_read_VIC01_72hpi_child --r1-length 28 --r2-length 90 --expect-cells 10950

##Long read adult

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_UK_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_UK_72hpi_adult --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam



"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_uninfected_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_uninfected_adult --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_48hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_VIC01_48hpi_adult --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_VIC01_72hpi_adult --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam




##Long read child


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_UK_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_UK_72hpi_child --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_uninfected_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_uninfected_child --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam



"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_48hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_VIC01_48hpi_child --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/analysis/illumina_fastq/fastq_path/" --sample Long_read_VIC01_72hpi_child --r1-length 28 --r2-length 90 --expect-cells 1930 --no-bam



echo "Complete!"


