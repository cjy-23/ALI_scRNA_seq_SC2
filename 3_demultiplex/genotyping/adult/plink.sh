#!/bin/bash
#SBATCH --job-name=plink_adult    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jessiejieyou@student.unimelb.edu.au     # Where to send mail
#SBATCH --ntasks=4                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=32000
#SBATCH --time=00:20:00               # Time limit hrs:min:sec
#SBATCH --partition=physical,snowy

echo "Running minimap on "`hostname`
 
#INPUT YOUR COMMAND HERE
module load plink/1.9b_6.21-x86_64


#make new file
plink --bfile HB00003100 --impute-sex --make-bed --out newfile 


#make clean new file_3
plink --bfile newfile --make-bed --set-hh-missing --out newfile_clean


#adult warning with codes
plink --bfile newfile_clean --recode vcf bgz --out newfile_clean_3 --snps-only just-acgt

echo "Complete!"


