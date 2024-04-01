#!/bin/bash
#SBATCH --job-name=plink_adult    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=8                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=32000
#SBATCH --time=100:00:00               # Time limit hrs:min:sec
#SBATCH --partition=physical,snowy

echo "Running minimap on "`hostname`
 
#INPUT YOUR COMMAND HERE


module load htslib/1.12
module load bcftools/1.12

sed -e 's/^/chr/' "/genotyping/adult/SNP_genotyping/plink_input/newfile_clean_3.vcf" > "/genotyping/adult/SNP_genotyping/plink_input/newfile_clean_3_chr.vcf"


#Then used nano to get rid of needless chr's


#Removed rows with "chr0"

sed '/chr0/d' "/genotyping/adult/SNP_genotyping/plink_input/newfile_clean_3_chr.vcf" > "/genotyping/adult/SNP_genotyping/plink_input/newfile_clean_3_chr_no0.vcf"


#Order chroms lexicographically

grep '^#' "/genotyping/adult/SNP_genotyping/plink_input/newfile_clean_3_chr_no0.vcf" > newfile_clean_3_chr_no0_lex.vcf && grep -v '^#' "/genotyping/adult/SNP_genotyping/plink_input/newfile_clean_3_chr_no0.vcf" | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> newfile_clean_3_chr_no0_lex.vcf

#bzgip lex 
bgzip -c "/genotyping/adult/SNP_genotyping/plink_input/no_warnings/newfile_clean_3_chr_no0_lex.vcf" > "/genotyping/adult/SNP_genotyping/plink_input/no_warnings/newfile_clean_3_chr_no0_lex.vcf.gz"

#tabix lex
tabix "/genotyping/adult/SNP_genotyping/plink_input/no_warnings/newfile_clean_3_chr_no0_lex.vcf.gz"

echo "Complete!"


