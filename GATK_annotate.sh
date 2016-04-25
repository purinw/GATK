#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script annotates a given VCF file with the following databases:
##i>	dbSNP
##ii>	dbNSFP
##iii>	gwasCat
##iv>	PhastCons
##v>	ClinVar
##vi>	SnpEff
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Run the command 'bash /path/to/GATK_annotate.sh file_name'
##-------------





##-------------
##Step0: Initialisation
##-------------
snpeff_dir=/mnt/data2/home2/purinw/snpEff
DBSNP=/home/bhoom/data/hg19/gatk_bundle/dbsnp_138.hg19.vcf
java_mem=30g
file_name=$1
prefix=$(echo ${file_name} | sed 's/.vcf$//')



##-------------
##Step1: dbSNP
##-------------
echo '1/6 dbSNP Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar annotate ${DBSNP} ${file_name} > ${prefix}.dbsnp.vcf
echo '1/6 dbSNP Annotation Completed'



##-------------
##Step2: dbNSFP
##-------------
echo '2/6 dbNSFP Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar dbnsfp -v -db ${snpeff_dir}/data/dbNSFP.txt.gz ${prefix}.dbsnp.vcf > ${prefix}.dbsnp.dbnsfp.vcf
echo '2/6 dbNSFP Annotation Completed'



##-------------
##Step3: gwasCat
##-------------
echo '3/6 gwasCat Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar gwasCat -db ${snpeff_dir}/data/gwascatalog.txt ${prefix}.dbsnp.dbnsfp.vcf > ${prefix}.dbsnp.dbnsfp.gwascat.vcf
echo '3/6 gwasCat Annotation Completed'



##-------------
##Step4: PhastCons
##-------------
echo '4/6 PhastCons Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar phastCons ${snpeff_dir}/data/phastCons ${prefix}.dbsnp.dbnsfp.gwascat.vcf > ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.vcf
echo '4/6 PhastCons Annotation Completed'



##-------------
##Step5: ClinVar
##-------------
echo '5/6 ClinVar Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar annotate ${snpeff_dir}/data/clinvar.vcf ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.vcf > ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.clinvar.vcf
echo '5/6 ClinVar Annotation Completed'



##-------------
##Step6: SnpEff
##-------------
echo '6/6 SnpEff Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/snpEff.jar hg19 ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.clinvar.vcf > ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.clinvar.snpeff.vcf
echo '6/6 SnpEff Annotation Completed'



##-------------
##Step7: Export
##-------------
mv ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.clinvar.snpeff.vcf ${prefix}_annotated.vcf
rm ${prefix}.dbsnp*vcf
