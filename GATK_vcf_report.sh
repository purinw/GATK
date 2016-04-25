#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script creats a tabular report of the variants within a VCF file that are related to genes specified by a given gene list file.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Run the command 'bash /path/to/GATK_vcf_report.sh file_name gene_list_file'
##-------------





##-------------
##Step0: Initialisation
##-------------
full_vcf_file=$1
gene_list=$2

prefix=$( echo ${full_vcf_file} | sed 's/.vcf$//' )



##-------------
##Step1: VCF File Reduction
##-------------
grep -v '^##' ${full_vcf_file} > ${prefix}.reduced.vcf



##-------------
##Step2: Allele Frequency Computation
##-------------
vcftools --vcf ${full_vcf_file} --counts --out ${prefix}.alleles



##-------------
##Step3: Tabulating Results (Via R)
##-------------
Rscript /mnt/data2/home2/purinw/Scripts/GATK_vcf_report.R ${prefix}.reduced.vcf ${prefix}.alleles.frq.count ${gene_list}



##-------------
##Step4: Removal of Temporary Files
##-------------
rm ${prefix}.reduced.vcf
rm ${prefix}.alleles.log
