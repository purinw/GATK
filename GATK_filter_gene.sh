#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script extracts variants within a VCF file that are related to genes specified by a given gene list file.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Run the command 'bash /path/to/GATK_filter_gene.sh file_name gene_list_file'
##-------------





##-------------
##Step0: Initialisation
##-------------
snpeff_dir=/mnt/data2/home2/purinw/snpEff
file_name=$1
gene_list=$2



##-------------
##Step1: Manipulating Arguments
##-------------
ARG=$(sed "s/\(.*\)/(ANN[*].GENE = '\1')/g" ${gene_list} | paste -sd '||')
prefix=$(echo ${file_name} | sed "s/.vcf$//")_$(echo ${gene_list} | sed "s/.txt$//" | sed "s/.*\///g")



##-------------
##Step2: Gene Filtration
##-------------
java -jar ${snpeff_dir}/SnpSift.jar filter "${ARG}" ${file_name} > ${prefix}.vcf
