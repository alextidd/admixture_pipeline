#!/bin/bash
cd ~/job_env_tidy/Data ;

#Order of pipe (data movement between folders):
##(1)ssh - raw VCF files downloaded from server
##(2)pop_filtered - VCFs > samples subsetted to populations of interest
##(3)biallelic_snps - VCFs > samples subsetted to populations of interest > biallelic SNPs only
##(4)AIMs - VCFs > samples subsetted to populations of interest > biallelic SNPs only > AIMs only
##(5)bcf - VCFs > samples subsetted to populations of interest > biallelic SNPs only > AIMs only > BCFs

function join_by { local IFS="$1"; shift; echo "$*"; } ;
populations=(ASW CEU YRI) ;
query_population=(ASW) ;
ref_populations=(CEU YRI) ;
populations_prefix=$(join_by _ "${populations[@]}") ;
query_population_prefix=$(join_by _ "${query_population[@]}") ;
ref_populations_prefix=$(join_by _ "${ref_populations[@]}") ;
lowerchr=(1) ;
upperchr=(22) ;

#Subset to populations of interest
echo "Subsetting to populations of interest..."
mkdir pop_filtered ;
../Code/subset_populations.R ${populations[@]} ALL ;
for chr in $(seq $lowerchr $upperchr); do
    vcftools --gzvcf ssh/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep "${populations_prefix}_ids.txt" --recode --stdout | gzip -c > pop_filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ;
done

#Filter to only biallelic SNPs. Extract rsIDs for AIM matching. Re-zip and index files.
echo "Filtering to only biallelic SNPs and extracting rsIDs for AIM matching. Then re-zipping and indexing the files..."
mkdir biallelic_snps rsIDs;
for chr in $(seq $lowerchr $upperchr); do
    bcftools view -Ov --max-alleles 2 --min-alleles 2 --types snps pop_filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf ;
    sed '/##/d' biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf | cut -f3 > rsIDs/rsIDs_chr"${chr}".txt ;
    bgzip -c biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf > biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ;
    tabix -p vcf biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ;
done

#Subsetting to only AIMs.
echo "Subsetting to only AIMs"
for chr in $(seq $lowerchr $upperchr); do
    ../Code/aims_matching.R ./rsIDs/rsIDs_chr"${chr}".txt ${chr} ;
    vcftools --vcf biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --snps AIMs/chr"${chr}"_all_matched_aims.txt --recode --recode-INFO-all --stdout > AIMs/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf ;
    bgzip -c AIMs/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf > AIMs/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ;
    tabix -p vcf AIMs/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ;
done
rm AIMs/chr*

#Create a karyogram showing genome coverage of the AIM sites
echo "Creating a karyogram showing genome coverage of the AIM sites"
mkdir matched_aims_genome_coverage ;
for chr in $(seq $lowerchr $upperchr); do
    sed '/##/d' AIMs/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf | cut -f2 > matched_aims_genome_coverage/chr"${chr}".txt ;
done
../Code/matched_aims_genome_coverage_karyogram.R ;

#Convert the 1000 Genomes files to BCF and index
echo "Converting the VCF files to BCF..."
mkdir bcf ;
for chr in $(seq $lowerchr $upperchr); do
   bcftools norm -m-any AIMs/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
   bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' > bcf/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
   bcftools index bcf/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done

#Convert the BCF files to PLINK format
echo "Converting the BCF files to PLINK format..."
mkdir plink ;
for chr in $(seq $lowerchr $upperchr); do
    plink --noweb --bcf bcf/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
    --keep-allele-order --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b37 no-fail --make-bed \
    --out plink/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done

##Identify and remove any duplicates
echo "Identifying and removing any duplicates..."
mkdir DupsRemoved ;
for chr in $(seq $lowerchr $upperchr); do
    plink --noweb --bfile plink/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --list-duplicate-vars ids-only ;
    plink --noweb --bfile plink/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --exclude plink.dupvar --make-bed \
    --out DupsRemoved/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
    rm plink.dupvar ;
done

#Prune variants from each chromosome (changed parameters for pruning)
echo "Pruning variants from each chromosome..."
mkdir Pruned ;
for chr in $(seq $lowerchr $upperchr); do
    plink --noweb --bfile DupsRemoved/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --maf 0.05 --indep 50 5 2 \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;

    plink --noweb --bfile DupsRemoved/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes \
    --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.prune.in --make-bed \
    --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes ;
done

#Get a list of all PLINK files (had to modify the sed command from the tutorial for mac)
echo "Getting a list of all PLINK files..."
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list.temporary ;
sed 's/.bim//g' ForMerge.list.temporary > ForMerge.list ;
rm ForMerge.list.temporary ;

#Merge all projects into a single PLINK file
echo "Merging all projects into a single PLINK file..."
plink --merge-list ForMerge.list --out Merge ;

#Perform PCA
echo "Performing PCA..."
plink --bfile Merge --pca ;

#Generate plot in R...
echo "Generating PCA plot in R..."
../Code/PCA_plot.R ; rm Rplots.pdf ;

echo "####RFMIX RUN####"
#Creating the subsetted vcf files for the QUERY and REFERENCE populations
echo "Creating the subsetted VCF files for the query and reference populations..."
../Code/subset_populations.R ${query_population[@]} ALL ;
../Code/subset_populations.R ${ref_populations[@]} ALL ;
mkdir query_pop reference_pop
for chr in $(seq $lowerchr $upperchr); do
    vcftools \
    --vcf biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \
    --keep "${query_population_prefix}_ids.txt" \
    --recode \
    --stdout > query_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf ;
    
    vcftools \
    --vcf biallelic_snps/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \
    --keep "${ref_populations_prefix}_ids.txt" \
    --recode \
    --stdout > reference_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf ;
done
rm *_ids.txt

#Writing the reference sample map. Reference sample map file for RFMIXv2.03 - The reference sample map used in the input for RFMIXv2.03 matches reference samples to their respective reference population. It consists of a tab-delimited text file with two columns: [Sample ID][Population]. The reference sample mape file is created during the run by wrangling the sample info file.
echo "Writing the reference sample map..."
../Code/reference_sample_map.R ${ref_populations[@]} ;

#Running RFMIX - make sure the RFMIXv2.03 executable is in the path. --crf-weight=2 bypasses slow "Generating internal simulation samples..." step of analysis, due to low number of reference samples. Remove for the final run / with more individuals.
echo "Running RFMIX"
for chr in $(seq $lowerchr $upperchr); do
    rfmix \
    --query-file=query_filtered_ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \
    --reference-file=reference_filtered_ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \
    --sample-map=filtered_${ref_populations_prefix}_reference_sample_map.txt \
    --genetic-map=./1000_Genomes_Phase3/reordered_rfmix2/genetic_map_chr"${chr}"_combined_b37.txt \
    -o RFMIX_output_chr"${chr}" \
    --chromosome=${chr} \
    --crf-weight=2 ;
done


#Plotting the RFMIX output on karyograms
echo "Plotting the RFMIX output on karyograms..."
../Code/chromomap_plot.R
