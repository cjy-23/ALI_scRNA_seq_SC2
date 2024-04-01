#!/bin/bash
#SBATCH --job-name=tabix    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                    # Run on a single CPU	
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --mem=32000
#SBATCH --time=1:00:00               # Time limit hrs:min:sec
#SBATCH --partition=physical,snowy

echo "Running minimap on "`hostname`
 
#INPUT YOUR COMMAND HERE


module load htslib/1.12
module load bcftools/1.12



#Add chr at each row

sed -e 's/^/chr/' "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3.vcf" > "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr.vcf"


#Then used nano to get rid of needless chr's


#Removed rows with "chr0"

sed '/chr0/d' "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr.vcf" > "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr_no0.vcf"


#Order chroms lexicographically

grep '^#' "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr_no0.vcf" > newfile_clean_3_chr_no0_lex.vcf && grep -v '^#' "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr_no0.vcf" | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> newfile_clean_3_chr_no0_lex.vcf



#bzgip lex 
bgzip -c "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr_no0_lex.vcf" > "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr_no0_lex.vcf.gz"

#tabix lex
tabix "/genotyping/child/SNP_genotyping/plink_input/newfile_clean_3_chr_no0_lex.vcf.gz"

echo "Complete!"


