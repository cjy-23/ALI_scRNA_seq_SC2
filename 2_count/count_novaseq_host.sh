#!/bin/bash
#SBATCH --job-name=iseq_bcl2fastq    # Job name
#SBATCH --partition=physical,snowy
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=8                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=120GB
#SBATCH --time=500:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE


source /usr/local/module/spartan_old.sh
module load bcl2fastq2/2-20-0


##Short read adult


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_uninfected_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_uninfected_adult --expect-cells 17850 --r1-length 28 --r2-length 90


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_48hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_VIC01_48hpi_adult --expect-cells 17850 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_VIC01_72hpi_adult --expect-cells 17850 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_UK_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_UK_72hpi_adult --expect-cells 17850 --r1-length 28 --r2-length 90



##Short read child

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_uninfected_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_uninfected_child --expect-cells 17850 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_48hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_VIC01_48hpi_child --expect-cells 17850 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_VIC01_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_VIC01_72hpi_child --expect-cells 17850 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Short_read_UK_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Short_read_UK_72hpi_child --expect-cells 17850 --r1-length 28 --r2-length 90


##Long read adult

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_uninfected_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_uninfected_adult --expect-cells 3150 --r1-length 28 --r2-length 90


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_48hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_VIC01_48hpi_adult --expect-cells 3150 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_VIC01_72hpi_adult --expect-cells 3150 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_UK_72hpi_adult --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_UK_72hpi_adult --expect-cells 3150 --r1-length 28 --r2-length 90



##Long read child

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_uninfected_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_uninfected_child --expect-cells 3150 --r1-length 28 --r2-length 90


"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_48hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_VIC01_48hpi_child --expect-cells 3150 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_VIC01_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_VIC01_72hpi_child --expect-cells 3150 --r1-length 28 --r2-length 90

"/data/gpfs/projects/punim1466/sw/cellranger/cellranger-6.1.1/cellranger" count --id Long_read_UK_72hpi_child --transcriptome "/data/gpfs/projects/punim1466/sw/cellranger/refdata-gex-GRCh38-2020-A/" --fastqs "/data/gpfs/projects/punim1466/sequencing_illumina/NovaSeq/CHA9252-NovaSeq/fastq/novaseq_cellranger/outs/fastq_path/" --sample Long_read_UK_72hpi_child --expect-cells 3150 --r1-length 28 --r2-length 90












echo "Complete!"

