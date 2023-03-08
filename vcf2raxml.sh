


## pooja.singh09@gmail.com
## code to make a phylogeny for whitefish from filtered SNPs in a vcf


##step1 keep vcf's zipped, I/O is faster

bsub gzip filter_all_basic_biallelic_minQ30_minDP2_maxDP50_minmeanDP2_maxmeanDP50_maf001_missing75.recode.vcf > filter_all_basic_biallelic_minQ30_minDP2_maxDP50_minmeanDP2_maxmeanDP50_maf001_missing75.recode.vcf.gz

##step2 ld prune your SNPs -> i think you have your own code for this antoine?

##step3 thin the full dataset to 1 site per kb; this helps the SNP calling run faster

bsub -W 120:00 -J "thin"$i -R "rusage[mem=3000]"  \
 "vcftools --gzvcf allGenomes.ldpruned.vcf.gz \
 --thin 1000 --stdout --recode | gzip > allGenomes.LDpruned.1kbThinned.vcf.gz"

#step4 filter SNPs: aeach site must have a homozygote of each allele
bsub < runFilter_raxmlsites.lsf

## runFilter_raxmlsites.lsf
#BSUB -J "raxmlfilter"
#BSUB -R "rusage[mem=8000]"
#BSUB -W 24:00
#BSUB -n 2
#BSUB -o log.%J.%I
#BSUB -e err.%J.%I

# load required modules
module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.8

bcftools view -i 'COUNT(GT="RR")>0 & COUNT(GT="AA")>0' \
 allGenomes.LDpruned.1kbThinned.vcf.gz | \
bgzip -c > allGenomes.LDpruned.1kbThinned.raxmlinput
tabix -fp vcf allGenomes.LDpruned.1kbThinned.raxmlinput

#count snps
bcftools view -H allGenomes.LDpruned.1kbThinned.raxmlinput.vcf.gz | wc -l

# Phylogeny

# Generate a phylogenetic tree: you can run this set of commands once with your thinned dataset and once without if you like but let's do it first with thinned SNPs
file=allGenomes.LDpruned.1kbThinned.raxmlinput

# Convert vcf to phylip

git clone https://github.com/edgardomortiz/vcf2phylip.git # this command you can run without bsub as it just downloads the vcf2phylip script from github

module load python/2.7
bsub "./vcf2phylip/vcf2phylip.py -i $file.vcf.gz -o $file.phylip"

# Run RAxML : multile outgroup individuals can be used with comma separations
# this command will help you understand the different parameters of raxml: raxmlHPC-PTHREADS-AVX -h 

module load gcc/4.8.2 gdc raxml/8.2.4
bsub -n 10 -W 120:00 -R 'rusage[mem=800]' "raxmlHPC-PTHREADS-AVX  \
 -m ASC_GTRGAMMA --asc-corr lewis -o outgroup_indivdual_ids -T 10 \
 -n ${file}.raxmlout.T10 -f a -p 12345 -x 12345 -s ${file}.phylip -N 100"

## download figtree and visualise the phylogeny on your computer

