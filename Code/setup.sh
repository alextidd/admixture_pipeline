#!/bin/bash

#The job_env will contain several packages and input files that must be present before the run and several output files that will be produced by the run. This shell script explains how the input files were sourced, downloaded and possibly wrangled before they are usable in the run...
#All scripts should be run in the "Data" directory as they use relative paths.
cd ~/job_env_tidy/Data

#PACKAGES/SOFTWARE USED
##PLINK
##PLINK2 - used within ancestry_pipeline run to generate .haps and .sample files
##BCFTOOLS
##VCFTOOLS
##RFMIXv2.03 - used to generate local ancestry data for karyograms
##RFMIXv1.5.4 (RunRFMix.py) - used within the ancestry_pipeline run (requires Linux, run on VirtualBox Ubuntu)
##python2.7

#INPUT DATA

##VCFs and TBIs - 1KGP phased haplotype data VCF files and their index files per chromosome (GRCh37) were downloaded from the SSH server, using my login (ssh art4017@login.cx1.hpc.ic.ac.uk). This required TunnelBlick to switch to the Imperial College VPN during download. Commands (-o ServerAliveInterval=15 -o ServerAliveCountMax=3) were added in an attempt to keep the connection from becoming idle and terminating the download, although this did not work). Explicit filenames with .gz and .gz.tbi had to be specified rather than open-ended .gz* as some chr folders contained multiple copies (e.g. .gz,.gz.tbi,.gz1,.gz.tbi1).
cd ./ssh
lowerchr=(1) ;
upperchr=(22) ;
for chr in $(seq $lowerchr $upperchr); do
    sshpass -v -p 'Artartart99' scp -o ServerAliveInterval=15 -o ServerAliveCountMax=3 art4017@login.cx1.hpc.ic.ac.uk:/rds/general/project/human-popgen-datasets/live/1000genomes_phase3/GRCh37/vcf/"${chr}"/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz . ;
    sshpass -v -p 'Artartart99' scp -o ServerAliveInterval=15 -o ServerAliveCountMax=3 art4017@login.cx1.hpc.ic.ac.uk:/rds/general/project/human-popgen-datasets/live/1000genomes_phase3/GRCh37/vcf/"${chr}"/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi . ;
done

sshpass -v -p 'Artartart99' scp -o ServerAliveInterval=15 -o ServerAliveCountMax=3 art4017@login.cx1.hpc.ic.ac.uk:/rds/general/project/human-popgen-datasets/live/1000genomes_phase3/GRCh37/vcf/4/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz . ;

##AIMs - Ancestry informative markers are used in analysis to speed up the process. They were derived from two papers (Tandon, Patterson & Reich, 2015; Reich et al,. 2007). These panels (files:reich2007aims.csv,reich2015aims.xls) were downloaded from the Supplementary Materials of these studies and then wrangled (script:aims_panels_combining.R) to form a combined panel (aims_rsIDs.txt). Because the studies use different human genome builds, the rsIDs alone were extracted for compatibility to create the final panel to be used in AIMs matching during the run.
../Code/aims_panels_combining.R #creates aims_rsIDs.txt

##Sample info - Sample information for the phase 3 cohort (ID,population,family info) was obtained from the FTP server as an excel file. The first sheet of the file ('Sample info') was converted to a .csv file. This is used for creating ID lists for population-specific subsetting.
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx ;

##PED file - The PED file was downloaded from the 1KGP FTP server.
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;

##Genetic map files for ancestry_pipeline - The genetic map files were downloaded from a Google Drive folder: https://drive.google.com/drive/folders/1z1a961djGVH3zWaAkFithHLqHx8ZEnkt and used as is ('1000_Genomes_Phase3').

##Genetic map files for RFMIXv2.03 - This is a tab-delimited text file containing 3 columns: [Chromosome][Physical position (bp)][Genetic position (cM)]. This file is generated per chromosomes from the genetic map files in the folder '1000_Genomes_Phase3' by reordering columns and adding a chromosome column. The reordered files are saved in the subfolder '/1000_Genomes_Phase3/reordered_rfmixv2'
../Code/genetic_map_reorder.R

##Reference sample map file for RFMIXv2.03 - The reference sample map used in the input for RFMIXv2.03 matches reference samples to their respective reference population. It consists of a tab-delimited text file with two columns: [Sample ID][Population]. The reference sample map file is created during the run by wrangling the sample info file ('filtered_<reference_populations>_reference_sample_map.txt'). This is because it depends on the run-specific populations.

##Chromosome dimensions / coordinates for chromoMap - In order to create karyograms of Local Ancestry Inference from the RFMIX output, the R package chromoMap requires coordinates for the lengths and centromere positions of the chromosomes in the particular genome build. These were found online and compiled into a .txt file, with the columns: [Chromosome][Chromosome start (=1)][Chromosome stop][Centromere start].
chromomap_chromosome_coords.txt
