#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script performs the entire joint variant-calling process of all the samples, following the Genome Analysis Toolkit (GATK)'s pipeline.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Check if the individual GATK process for every sample has been accomplished.
##iii>	Run the command 'bash /path/to/GATK_joint.sh'
##-------------





##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
fastq_dir=/mnt/data2/home2/purinw/WES/RawData
home_dir=/mnt/data2/home2/purinw/WES
ref_dir=/home/bhoom/data/hg19/gatk_bundle
out_dir=${home_dir}/Output
gatk_dir=/home/bhoom/bin/gatk3.3-0
samtools_dir=/home/bhoom/bin/samtools-0.1.19
resource_dir=/mnt/data2/home2/purinw/resource
sample_list_file=${home_dir}/sample_list.txt

##-------------
##Step0-2: References
##-------------
ref_genome=${ref_dir}/ucsc.hg19.fasta
indel_1=${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
indel_2=${ref_dir}/1000G_phase1.indels.hg19.sites.vcf
DBSNP=${ref_dir}/dbsnp_138.hg19.vcf

##-------------
##Step0-3: Other Parametres
##-------------
java_mem=30g

##-------------
##Step0-4: List of Samples
##-------------
sample_list=($(cat ${sample_list_file}))

##-------------
##Step0-5: Name of the joint file
##-------------
sample_name=COMBINED

##-------------
##Step0-6: Output Folders Creation
##-------------
mkdir -p ${out_dir}/${sample_name}
mkdir -p ${out_dir}/${sample_name}/{Script,LOG,GVCF,VCF,VQSR,VQSR+QC,VQSR+QC/FILTERED_ON_BAIT}



##-------------
##Step1: Combine GVCFs
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/01_${sample_name}_combine_vcfs.sh
#!/bin/bash
##-------------
##Step1-0: Create Variant Arguments
##-------------
ARG='$( sed 's/^\(.*\)/-V \1\/GVCF\/\1_GATK.vcf \\/g' ${home_dir}/sample_list.txt )'

##-------------
##Step1-1: Merge GVCFs
##-------------
cd ${out_dir}
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T CombineGVCFs \
-R ${ref_genome} \
\$ARG \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o ${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.vcf \
-log ${out_dir}/${sample_name}/LOG/01_combine_gvcfs.log
cd ${home_dir}

EOL



##-------------
##Step2: Joint Genotype
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/02_${sample_name}_joint_genotype.sh
#!/bin/bash
##-------------
##Step2: Joint Genotype
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/GVCF/COMBINED_GATK.vcf \
-nt 8 \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o  ${out_dir}/${sample_name}/VCF/${sample_name}_QC_RAW_OnBait.vcf \
-log ${out_dir}/${sample_name}/LOG/02_${sample_name}_genotype_gvcf.log

EOL


##-------------
##Step3: Variant Quality Score Recalibration
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/03_${sample_name}_recalibrate_variant.sh
#!/bin/bash
##-------------
##Step3-1-1: Select SNVs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
-selectType SNP \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-o ${out_dir}/${sample_name}/VQSR/${sample_name}_SELECT_SNV.vcf \
-log ${out_dir}/${sample_name}/LOG/03-1-1_${sample_name}_VQSR_select_SNV.log

##-------------
##Step3-1-2: Select Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VCF/${sample_name}_RAW.vcf \
-selectType INDEL \
-selectType MNP \
-selectType MIXED \
-selectType SYMBOLIC \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-o ${out_dir}/${sample_name}/VQSR/${sample_name}_SELECT_INDEL.vcf \
-log ${out_dir}/${sample_name}/LOG/03-1-2_${sample_name}_VQSR_select_INDEL.log

##-------------
##Step3-2-1: Recalibrate SNVs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref_genome} \
--input:VCF ${out_dir}/${sample_name}/VCF/${sample_name}_SELECT_SNV.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${resource_dir}/hapmap_3.3.hg19.sites.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 ${resource_dir}/1000G_omni2.5.hg19.sites.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${resource_dir}/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
-mode SNP \
-an QD \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.8 \
-tranche 99.7 \
-tranche 99.6 \
-tranche 99.5 \
-tranche 99.4 \
-tranche 99.3 \
-tranche 99.2 \
-tranche 99.1 \
-tranche 99.0 \
-tranche 98.0 \
-tranche 97.0 \
-tranche 96.0 \
-tranche 95.0 \
-tranche 90.0 \
-recalFile ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV.recal \
-tranchesFile ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV.tranches \
-rscriptFile ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV.R \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-log ${out_dir}/${sample_name}/LOG/03-1-1_${sample_name}_VQSR_snv_recalibration.log

##-------------
##Step3-2-2: Recalibrate Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref_genome} \
--input:VCF ${out_dir}/${sample_name}/VQSR/${sample_name}_SELECT_INDEL.vcf \
-resource:mills,known=true,training=true,truth=true,prior=12.0 ${indel_1} \
-mode INDEL \
--maxGaussians 4 \
-an QD \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.8 \
-tranche 99.7 \
-tranche 99.6 \
-tranche 99.5 \
-tranche 99.4 \
-tranche 99.3 \
-tranche 99.2 \
-tranche 99.1 \
-tranche 99.0 \
-tranche 98.0 \
-tranche 97.0 \
-tranche 96.0 \
-tranche 95.0 \
-tranche 90.0 \
-recalFile ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL.recal \
-tranchesFile ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL.tranches \
-rscriptFile ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL.R \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-log ${out_dir}/${sample_name}/LOG/03-2-2_${sample_name}_VQSR_indel_recalibration.log

##-------------
##Step3-3-1: Apply SNV Recalibration
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref_genome} \
--input:VCF ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV.vcf \
-ts_filter_level 99.9 \
-recalFile ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV.recal \
-tranchesFile ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV.tranches \
-mode SNP \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV_RECAL_APPLIED.vcf \
-log ${out_dir}/${sample_name}/LOG/03-3-1_${sample_name}_VQSR_apply_snv_recalibration.log

##-------------
##Step3-3-2: Apply Indels Recalibration
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref_genome} \
--input:VCF ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL.vcf \
-ts_filter_level 99.9 \
-recalFile ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL.recal \
-tranchesFile ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL.tranches \
-mode INDEL \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL_RECAL_APPLIED.vcf \
-log ${out_dir}/${sample_name}/LOG/03-3-2_${sample_name}_VQSR_apply_indel_recalibration.log

##-------------
##Step3-4: Combine SNVs + Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T CombineVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VQSR/${sample_name}_SNV_RECAL_APPLIED.vcf \
--variant ${out_dir}/${sample_name}/VQSR/${sample_name}_INDEL_RECAL_APPLIED.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--genotypemergeoption UNSORTED \
-o ${out_dir}/${sample_name}/VQSR/${sample_name}_FINAL_VQSR.vcf \
-log ${out_dir}/${sample_name}/LOG/03-4_${sample_name}_VQSR_combine_variants.log

EOL



##-------------
##Step4: SNP Quality Control
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/04_${sample_name}_SNV_quality_control.sh
#!/bin/bash
##-------------
##Step4-1-1: Extract SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VQSR/${sample_name}_FINAL_VQSR.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType SNP \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_RAW_SNV.vcf \
-log ${out_dir}/${sample_name}/LOG/04-1-1_${sample_name}_QC_select_snv.log

##-------------
##Step10-1-2: Extract Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VQSR/${sample_name}_final_VQSR.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType INDEL \
-selectType MNP \
-selectType MIXED \
-selectType SYMBOLIC \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_RAW_INDEL.vcf \
-log ${out_dir}/${sample_name}/LOG/04-1-2_${sample_name}_QC_select_INDEL.log

##-------------
##Step4-2-1: Annotate SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_QC_RAW_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--dbsnp ${DBSNP} \
-L ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_RAW_SNV.vcf \
-A GCContent \
-A VariantType \
-dt NONE \
-nt 1 \
-o  ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_RAW_SNV_ANNOTATED.vcf \
-log ${out_dir}/${sample_name}/LOG/04-2-1_${sample_name}_QC_snv_annotation.log

##-------------
##Step04-3-1: Filter SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_RAW_SNV_ANNOTATED.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--filterExpression 'QD < 2.0' \
--filterName 'QD' \
--filterExpression 'MQ < 30.0' \
--filterName 'MQ' \
--filterExpression 'FS > 40.0' \
--filterName 'FS' \
--filterExpression 'MQRankSum < -12.5' \
--filterName 'MQRankSum' \
--filterExpression 'ReadPosRankSum < -8.0' \
--filterName 'ReadPosRankSum' \
--filterExpression 'DP < 8.0' \
--filterName 'DP' \
--logging_level ERROR \
-o ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_FILTERED_SNV.vcf \
-log ${out_dir}/${sample_name}/LOG/04-3-1_${sample_name}_QC_filter_snv.log

##-------------
##Step4-4-1: Clean SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_FILTERED_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${sample_name}/VQSR+QC/FILTERED_ON_BAIT/${sample_name}_CLEAN_SNV.vcf \
-log ${out_dir}/${sample_name}/LOG/04-4-1_${sample_name}_QC_clean_snv.log

##-------------
##Step4-5: Combine SNVs + Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T CombineVariants \
-R ${ref_genome} \
--variant ${out_dir}/${sample_name}/VQSR+QC/FILTERED_ON_BAIT/${sample_name}_CLEAN_SNV.vcf \
--variant ${out_dir}/${sample_name}/VQSR+QC/${sample_name}_RAW_INDEL.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--genotypemergeoption UNSORTED \
-o ${out_dir}/${sample_name}/VQSR+QC/FILTERED_ON_BAIT/${sample_name}_CLEAN_SNV+INDEL.vcf \
-log ${out_dir}/${sample_name}/LOG/4-5_${sample_name}_QC_combine_variants.log

EOL





##-------------
##MASTER SCRIPT
##-------------
cat <<EOL > ${out_dir}/${sample_name}/Script/${sample_name}_GATK.sh
#!/bin/bash
##-------------
##Joint Vaiant Calling
##-------------
bash ${out_dir}/${sample_name}/Script/01_${sample_name}_combine_vcfs.sh
bash ${out_dir}/${sample_name}/Script/02_${sample_name}_joint_genotype.sh
bash ${out_dir}/${sample_name}/Script/03_${sample_name}_recalibrate_variant.sh
bash ${out_dir}/${sample_name}/Script/04_${sample_name}_SNV_quality_control.sh

EOL





##-------------
##EXECUTION
##-------------
bash ${out_dir}/${sample_name}/Script/${sample_name}_GATK.sh
