~/#!/bin/bash
cd ~/job_env/Data ;

#Final populations=(ASW GBR CEU TSI IBS GWD ESN MSL YRI LWK)
#Final query_population=(ASW)
#Final ref_populations=(GBR CEU TSI IBS GWD ESN MSL YRI LWK)

function join_by { local IFS="$1"; shift; echo "$*"; } ;
populations=(ASW GBR CEU TSI IBS GWD ESN MSL YRI LWK) ;
query_population=(ASW) ;
ref_populations=(GBR CEU TSI IBS GWD ESN MSL YRI LWK) ;
populations_prefix=$(join_by _ "${populations[@]}") ;
query_population_prefix=$(join_by _ "${query_population[@]}") ;
ref_populations_prefix=$(join_by _ "${ref_populations[@]}") ;
lowerchr=(1) ;
upperchr=(22) ;

#Creating subsetted populations ID files
../Code/subset_populations.R ${populations[@]} ALL ;
../Code/subset_populations.R ${query_population[@]} ALL ;
../Code/subset_populations.R ${ref_populations[@]} ALL ;

#Convert the 1000 Genomes files to BCF and index.
echo "Converting the VCF files to BCF..."
mkdir bcf ;
for chr in $(seq $lowerchr $upperchr); do
    bcftools view ssh/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Ob --threads 4 > bcf/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf   ;
    bcftools index bcf/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done

#Filtering. (1) Subset SAMPLES to populations of interest. (2) Filter SITES to biallelic SNPs. (3) Filter rsIDs to AIMs. (3) Index.
echo "Filtering (populations of interest, biallelic SNPs, AIMs)..."
mkdir filtered
for chr in $(seq $lowerchr $upperchr); do
    bcftools view --samples-file "${populations_prefix}_ids.txt" --force-samples \
    --max-alleles 2 --min-alleles 2 --types snps \
    --include 'ID=@AIMs/aims_rsIDs.txt' \
    bcf/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf | \
    bcftools norm --rm-dup both -Ob \
    > filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
    bcftools index filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done

#Create a karyogram showing genome coverage of the AIM sites
echo "Creating a karyogram showing genome coverage of the AIM sites"
mkdir coverage ;
for chr in $(seq $lowerchr $upperchr); do
    bcftools query -f '%POS\n' filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf > coverage/chr"${chr}".txt ;
done
../Code/coverage_karyogram.R ;

#Generate distinct variant IDs for PLINK wrangling. Sets the ID field to a unique value: CHROM:POS:REF:ALT.
mkdir distinct_ids
for chr in $(seq $lowerchr $upperchr); do
    bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \
    filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
    > distinct_ids/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
    bcftools index distinct_ids/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done


#Convert the BCF files to PLINK format
echo "Converting the BCF files to PLINK format..."
mkdir plink ;
for chr in $(seq $lowerchr $upperchr); do
    plink --noweb --bcf distinct_ids/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
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

#Get a list of all PLINK files
echo "Getting a list of all PLINK files..."
mkdir merged_plink pca
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list.temporary ;
sed 's/.bim//g' ForMerge.list.temporary > merged_plink/ForMerge.list ;
rm ForMerge.list.temporary ;

#Merge all projects into a single PLINK fileset
echo "Merging all projects into a single PLINK file..."
plink --merge-list merged_plink/ForMerge.list --out merged_plink/Merge ;
../Code/admixture_pop_file.R

#Perform smartPCA
awk '{print $1,$2,$3,$4,$5,1}' merged_plink/Merge.fam > pca/Merge.PCA.fam

echo genotypename: merged_plink/Merge.bed > pca/pca.par
echo snpname: merged_plink/Merge.bim >> pca/pca.par
echo indivname: pca/Merge.PCA.fam >> pca/pca.par
echo snpweightoutname: pca/pca.snpeigs >> pca/pca.par
echo evecoutname: pca/pca.eigs >> pca/pca.par
echo evaloutname: pca/pca.eval >> pca/pca.par
echo phylipoutname: pca/pca.fst >> pca/pca.par
echo numoutevec: 20 >> pca/pca.par
echo numoutlieriter: 0 >> pca/pca.par
echo outlieroutname: pca/pca >> pca/pca.par
echo altnormstyle: NO >> pca/pca.par
echo missingmode: NO >> pca/pca.par
echo nsnpldregress: 0 >> pca/pca.par
echo noxdata: YES >> pca/pca.par
echo nomalexhet: YES >> pca/pca.par

smartpca -p pca/pca.par

for chr in {1..22} ; do
    bcftools query -f '%POS\n' filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf | wc -l
done

echo "####RFMIX RUN####"
#Creating the subsetted bcf files for the QUERY and REFERENCE populations
echo "Creating the subsetted BCF files for the query and reference populations..."
mkdir query_pop reference_pop
for chr in $(seq $lowerchr $upperchr); do
    bcftools view --samples-file "${query_population_prefix}_ids.txt" --force-samples filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf -Ob > query_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
    bcftools index query_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
    
    bcftools view --samples-file "${ref_populations_prefix}_ids.txt" --force-samples filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf -Ob > reference_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
    bcftools index reference_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done

#Writing the reference sample map. Reference sample map file for RFMIXv2.03 - The reference sample map used in the input for RFMIXv2.03 matches reference samples to their respective reference population. It consists of a tab-delimited text file with two columns: [Sample ID][Population]. The reference sample mape file is created during the run by wrangling the sample info file.
echo "Writing the reference sample map for RFMIXv2..."
mkdir rfmixout
../Code/reference_sample_map.R ${ref_populations[@]} ;
cp filtered_"${ref_populations_prefix}"_reference_sample_map.txt ./rfmixout

#Running RFMIX.
echo "Running RFMIX"
for chr in $(seq $lowerchr $upperchr); do
    rfmix \
    --query-file=query_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
    --reference-file=reference_pop/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf \
    --sample-map=rfmixout/filtered_${ref_populations_prefix}_reference_sample_map.txt \
    --genetic-map=1000_Genomes_Phase3/reordered_rfmixv2/genetic_map_chr"${chr}"_combined_b37.txt \
    -o rfmixout/rfmix_output_chr"${chr}" \
    --chromosome=${chr} ;
done

#Plotting the RFMIX output on karyograms and pie charts can be done in the . Calculating global ancestry proportions. By summing up local ancestry proportions, global ancestry propportions can be obtained.


#Setting up admixture run directories & Generating population file for supervised ADMIXTURE run (.pop file specifying the ancestries of reference individuals. Each line of thhe .pop file corresponds to an individual listed on the same line number in the .fam file. If the individual is a population reference, the .pop file line should be a string designating the population. If the individual is of unknown ancestry, use '-' to indicate that the ancestry should be estimated.
echo "Prepping for admixture runs by setting up directories."
../Code/admixture_pop_file.R ;
mkdir admixture_supervised.K9 admixture_unsupervised.K9 admixture_unsupervised.K2 ;
cp merged_plink/Merge9.* admixture_supervised.K9 ; cp merged_plink/Merge.* admixture_supervised.K9 ;
cp merged_plink/Merge9.* admixture_unsupervised.K9 ; cp merged_plink/Merge.* admixture_unsupervised.K9 ;
cp merged_plink/Merge2.* admixture_unsupervised.K2 ; cp merged_plink/Merge.* admixture_unsupervised.K2 ;

#Colour assignment TSI IBS GBR CEU ASW YRI ESN MSL GWD LWK
#EF91CA,#E04B4B,#F0934E,#FCFB53,#B37855,#63BC6A,#6094C3,#A76BB2,#A4A4A4
GBR #E04B4B
IBS #F0934E
GWD #B37855
ESN #6094C3
MSL #63BC6A
CEU #EF91CA
YRI #A76BB2
LWK #A4A4A4
TSI #FCFB53



#Running ADMIXTURE and plotting - supervised, K=9
cd admixture_supervised.K9 ;
mv Merge9.pop Merge.pop ;
admixture --cv Merge.bed ${#ref_populations[@]} -j4 --supervised | tee log"${#ref_populations[@]}".log ;

sed 's/-/ASW/g' Merge.pop > ind2pop.txt
pop_order=(TSI IBS GBR CEU ASW YRI ESN MSL GWD LWK)

printf "%s\n" "${pop_order[@]}" > pop_order.txt
K=9
awk -v K=$K -v file=Merge 'BEGIN{ \
    printf("Supervised_Run_K9\t%d\t%s.%d.Q\n",K,file,K) }' > Merge.Qfilemap
pong -m Merge.Qfilemap -i ind2pop.txt --pop_names pop_order.txt --color_list ../pal_superpops.txt
cd .. ;

#Running ADMIXTURE - unsupervised, K=9

cd admixture_unsupervised.K9 ;
mv Merge9.pop Merge.pop ;
for i in $(seq 1 ${#ref_populations[@]}) ; do
    admixture --cv -j4 Merge.bed ${#ref_populations[@]} | tee log"${#ref_populations[@]}".out ;
done

sed 's/-/ASW/g' Merge.pop > ind2pop.txt
pop_order=(TSI IBS GBR CEU ASW YRI ESN MSL GWD LWK)
printf "%s\n" "${pop_order[@]}" > pop_order.txt
K=9
awk -v K=$K -v file=Merge 'BEGIN{ \
    printf("Unsupervised_Run_K9\t%d\t%s.%d.Q\n",K,file,K) }' > Merge.Qfilemap
pong -m Merge.Qfilemap -i ind2pop.txt --pop_names pop_order.txt --color_list ../pal_superpops.txt
cd ..

#Running ADMIXTURE - unsupervised, K=2
cd admixture_unsupervised.K2 ;
mv Merge2.pop Merge.pop

for i in $(seq 1 2) ; do
    admixture --cv -j4 Merge.bed 2 | tee log2.out ;
done

sed 's/-/ASW/g' Merge.pop > ind2pop.txt
pop_order=(EUR ASW AFR)
printf "%s\n" "${pop_order[@]}" > pop_order.txt
K=2
awk -v K=$K -v file=Merge 'BEGIN{ \
    printf("Unsupervised_Run_K2\t%d\t%s.%d.Q\n",K,file,K) }' > Merge.Qfilemap
pong -m Merge.Qfilemap -i ind2pop.txt --pop_names pop_order.txt --color_list ../pal_pops_reordered.txt
cd ..

#Running ADMIXTURE - unsupervised, considering different values of K (outputs CVs for each)
mkdir admixture_unsupervised.Kx ; cp merged_plink/Merge2.* admixture_unsupervised.Kx ; cp merged_plink/Merge.* admixture_unsupervised.Kx ; cd admixture_unsupervised.Kx
Klow=1 ;
Khigh=12 ;
prefix=K_search ;
for ((K=$Klow;K<=$Khigh;K++)); do
    admixture --cv Merge.bed $K | tee log.$prefix.${K}.out ;
done ;
echo '# CV results' > $prefix.CV.txt
for ((K=$Klow;K<=$Khigh;K++)); do
    awk -v K=$K '$1=="CV"{print K,$4}' log.$prefix.$K.out \
    >> $prefix.CV.txt
done

#5 runs for K of 2-12
prefix=K_search ; mkdir runs ;
for r in {1..5}; do for K in {2..12};
do
    admixture -s ${RANDOM} runs/${prefix}.bed $K
    mv ${prefix}.${K}.Q runs/${prefix}.K${K}r${r}.Q
done; done

cd runs
createQmap(){
local r=$1
local K=$2
awk -v K=$K -v r=$r -v file=${prefix}.K${K}r${r} 'BEGIN{ \
printf("K%dr%d\t%d\t%s.Q\n",K,r,K,file)
}' >> ${prefix}.multiplerun.Qfilemap
}
export -f createQmap
for K in {2..12}; do for r in {1..5}; do createQmap $r $K ; \
done; done

pong -m K_search.multiplerun.Qfilemap --greedy \
-s.95 -i ind2pop.txt


#Running ChromoPainterv2 and GLOBETROTTER (only for ASW versus CEU & YRI)
mkdir chromopainter_in ; mkdir chromopainter_in/individuals ;

bcftools concat -o chromopainter_in/merged_bcf filtered/*.bcf ; #creates concatenated, merged bcf

IFS=$'\r\n' GLOBIGNORE='*' command eval  'IND=($(cat ASW_CEU_YRI_ids.txt))' ;
for i in "${IND[@]}" ; do
bcftools view -s "${i}" --force-samples chromopainter_in/merged_bcf | bcftools query -f '%ID\t%CHROM\t%POS\t[%GT]\n' > "chromopainter_in/individuals/${i}.txt" ;
done

##Chromopainter input files
mkdir chromopainter_in ; cd chromopainter_in ;
#Run the R script to generate these files
cd chromopainter_in ;
ChromoPainterv2 -g haplotype_infile.txt -t label_infile.txt -r recom_rate_infile.txt -f population_list_infile.txt ;

./ChromoPainterv2 -g example/BrahuiYorubaSimulationChrom22.haplotypes -r example/BrahuiYorubaSimulationChrom22.recomrates -t example/BrahuiYorubaSimulation.idfile.txt -f example/BrahuiYorubaSimulation.poplist.txt 0 0 -o example/BrahuiYorubaSimulationChrom22

#TRACTS run
cd 2pops ;
mkdir tracts_in ;
#Create tracts input file in R
mv tracts_in 2pops

python2.7 ./ASW_tractlength_fix.py
python2.7 ./fancyplotting.py --name out out2 --population-tags European,African

cd .. ;


#Assortative mating run (geneset built in R, filter to those individuals and sites)
mkdir assortative_mating_analysis
for chr in $(seq $lowerchr $upperchr); do
    bcftools view -i 'ID=@assortative_mating_analysis/geneset_rsids.txt' bcf/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf | bcftools query -f '%ID\n' >> assortative_mating_analysis/matched.txt ;
done

echo "Filtering (populations of interest, pigmentation alleles)..."
mkdir pigment_filtered ;
for chr in $(seq $lowerchr $upperchr); do
    bcftools view --samples-file "${query_population_prefix}_ids.txt" --force-samples \
    --include 'ID=@assortative_mating_analysis/geneset_rsids.txt' \
    bcf/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf -Ou | \
    bcftools norm --rm-dup both -Ob \
    > pigment_filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
    bcftools index pigment_filtered/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf ;
done

bcftools concat -o pigment_filtered/merged_bcf pigment_filtered/*.bcf ;

bcftools query -f '%CHROM\t%POS\t%ID\t%EUR_AF\t%AFR_AF\t%REF\t%ALT\n' pigment_filtered/merged_bcf > assortative_mating_analysis/matched.txt;


