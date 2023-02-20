#!/bin/bash
#SBATCH --job-name=process   # Job name
#SBATCH --partition=physical
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=50GB
#SBATCH --time=1:00:00               # Time limit hrs:min:sec


#INPUT YOUR COMMAND HERE
module load r/4.2.0
Rscript "/data/gpfs/projects/punim1466/analysis/cellranger_new_run/DTU_15/inf_vs_uninf/secretory_Alpha_uninf_child/extra/dturtle_process_pseudobulk_final_1.R"
echo "Complete!"