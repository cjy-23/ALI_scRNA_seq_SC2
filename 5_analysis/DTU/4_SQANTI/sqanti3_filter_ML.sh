#!/bin/bash
#SBATCH --job-name=sqanti3_filter    # Job name
#SBATCH --partition=physical
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=20GB
#SBATCH --time=3:00:00              # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE


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


python "/data/gpfs/projects/punim1466/sw/SQANTI3-5.1/sqanti3_filter.py" ML --isoforms "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/sqanti/all_refined_corrected.fasta" --gtf "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/sqanti/all_refined_corrected.gtf" --faa "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/sqanti/all_refined_corrected.faa" -j 0.5 -o filter_lib_ML all_refined_classification.txt

echo "Complete!"


