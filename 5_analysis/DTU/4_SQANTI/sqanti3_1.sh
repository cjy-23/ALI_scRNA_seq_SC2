#!/bin/bash
#SBATCH --job-name=sqanti3    # Job name
#SBATCH --partition=physical
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=100GB
#SBATCH --time=1:00:00              # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE
module load gcc/8.3.0
module load minimap2/2.17
module load samtools/1.9


#python sqanti_qc2.py [-t cpus] [-n chunks]
 #    [--gtf] [--skipORF] 
 #    [-c shortread_STAR_junction_out] 
 #    [--cage_peak CAGE_PEAK_BED]
 #    [--polya_peak PolyA_PEAK_BED]
 #    [--polyA_motif_list POLYA_MOTIF_LIST]
 #    [--fl_count FL_COUNT]
 #    [--expression EXPRESSION]
 #    [--aligner_choice=minimap2,deSALT]
 #    [--is_fusion]
 #    <input_fasta> <annotation_gtf> <genome_fasta>

export PYTHONPATH=$PYTHONPATH:"/home/jessiejieyou/.conda/envs/SQANTI3_5.1.env/bin/cDNA_Cupcake/sequence/"
export PYTHONPATH=$PYTHONPATH:"/home/jessiejieyou/.conda/envs/SQANTI3_5.1.env/bin/cDNA_Cupcake/"


python "/data/gpfs/projects/punim1466/sw/SQANTI3-5.1/sqanti3_qc.py" -t 30 -n 20 "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/new_ref/refined_gtf.gtf" "/data/gpfs/projects/punim1466/db/refdata-gex-GRCh38-2020-A/genes/genes.gtf" "/data/gpfs/projects/punim1466/db/refdata-gex-GRCh38-2020-A/fasta/genome.fa" --aligner_choice=minimap2 --force_id_ignore --CAGE_peak "/data/gpfs/projects/punim1466/sw/SQANTI3-5.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed" --polyA_motif_list "/data/gpfs/projects/punim1466/sw/SQANTI3-5.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt" --polyA_peak "/data/gpfs/projects/punim1466/sw/SQANTI3-4.3/data/polyA_peak/atlas.clusters.2.0.GRCh38.96.bed" -o all_refined --saturation --report both




echo "Complete!"


