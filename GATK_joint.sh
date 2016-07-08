#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script performs the entire joint variant-calling process upon all samples, following the Genome Analysis Toolkit (GATK)'s pipeline.
##If no \"sample_name\" is specified in the command, all the folders in \"out_dir\" are automatically taken as the complete set of samples.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Check if the individual GATK process for every sample has been accomplished.
##iii>	Run the command 'bash /path/to/GATK_joint.sh [options] [sample_name]'
##-------------





##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
project_dir=/mnt/data2/home2/purinw/SUDS-12
ref_dir=/home/bhoom/data/hg19/gatk_bundle
resource_dir=/mnt/data2/home2/purinw/resource
bwa_dir=/home/bhoom/bin/bwa-0.7.5a
picard_dir=/home/bhoom/bin/picard-tools-1.119
gatk_dir=/home/bhoom/bin/gatk3.3-0
samtools_dir=/home/bhoom/bin/samtools-0.1.19
fastq_dir=${project_dir}/RawData
out_dir=${project_dir}/Output

##-------------
##Step0-2: References
##-------------
ref_genome=${ref_dir}/ucsc.hg19.fasta
indel_1=${ref_dir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
indel_2=${ref_dir}/1000G_phase1.indels.hg19.sites.vcf
DBSNP=${ref_dir}/dbsnp_138.hg19.vcf
exon_bed=/mnt/data2/home2/purinw/bed/exon_hg19.bed

##-------------
##Step0-3: Other Parametres
##-------------
java_mem=30g

##-------------
##Step0-4: Input Arguments
##-------------
while test $# -gt 0 ; do
        case "$1" in
                -h|--help)
				        echo ""
                        echo "Usage: bash $0 [options] [sample_name] [sample_name] [...]"
                        echo ""
                        echo "This script performs the entire joint variant-calling process upon all samples, following the Genome Analysis Toolkit (GATK)'s pipeline."
						echo "If no \"sample_name\" is specified in the command, all the folders in \"out_dir\" are automatically taken as the complete set of samples."
                        echo ""
                        echo "Options:"
                        echo "-h, --help				display this help and exit"
						echo "-v, --version				display version of this script and exit"
						echo "-XS, --no-summary			suppress the command summary before execution"
						echo "-XP, --no-prompt			suppress the user prompt before execution, only when the command summary is displayed"
						echo "-XX, --no-exec				suppress automatic execution, generating only script files"
						echo "-x, --exclude		SAMPLE_NAME	exclude SAMPLE_NAME from the (automatically generated) sample list"
                        echo "-e, --exome				call only exonic variants, drastically accelerating the Joint Genotype process"
						echo "-p, --prefix		PREFIX		specify the batch's name to be used, by default \"COMBINED\""
						echo "-s, --split				split the final joint file into single-sample files"
						echo ""
                        exit 0
                        ;;
				-v|--version)
						echo ""
						echo "GATK_joint.sh"
                        echo ""
						echo "Created MAR 2016"
						echo "Updated JUL 2016"
						echo "by"
						echo "PURIN WANGKIRATIKANT [purin.wan@mahidol.ac.th]"
						echo "Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)"
						echo "Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand"
						echo ""
						exit 0
						;;
				-XS|--no-summary)
						no_summary=1
						shift
						;;
				-XP|--no-prompt)
						no_prompt=1
						shift
						;;
				-XX|--no-exec)
						no_exec=1
						shift
						;;
				-x|--exclude)
						shift
						exclude_list+=( $1 )
						shift
						;;
                -e|--exome)
						seq_type='EXOME'
						bed_argument='-L '${exon_bed}
                        shift
                        ;;
				-p|--prefix)
                        shift
						joint_name=$1
                        shift
                        ;;
				-s|--split)
						sample_split='YES'
						shift
						;;
				*)
						sample_list+=( $1 )
                        shift
						;;
		esac
done

##-------------
##Step0-5: Default Value Setting
##-------------
if [[ ! -v seq_type ]] ; then
		seq_type='GENOME'
fi
if [[ ! -v joint_name ]] ; then
		joint_name='COMBINED'
fi
if [[ ! -v sample_split ]] ; then
		sample_split='NO'
fi
if [[ ! -v sample_list ]] ; then
		sample_list=($(ls ${out_dir}))
fi
for exclude_name in ${exclude_list[*]} ; do
		sample_list=(${sample_list[*]/"${exclude_name}"})
done

##-------------
##Step0-6: Input Verification
##-------------
for sample_name in ${sample_list[*]} ; do
		if [[ ! -e ${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.g.vcf ]] ; then
				echo
				echo 'Invalid SAMPLE NAME: '${sample_name}
				echo ${out_dir}/${sample_name}/GVCF/${sample_name}_GATK.g.vcf not found.
				echo 'Terminated.'
				echo
				exit 1
		fi
done

##-------------
##Step0-7: Summarisation & User's Confirmation Prompt
##-------------
if [[ ${no_summary} != 1 ]] ; then
		echo
		echo '---------------------------------------'
		echo 'JOINT VARIANT CALLING PROCESS'
		echo 'SAMPLE COUNT =			'${#sample_list[*]}
		echo 'SAMPLE NAMES =			'${sample_list[*]}
		echo 'SEQUENCED DATA =		'${seq_type}
		echo 'PREFIX =			'${joint_name}
		echo 'SAMPLE SPLIT =			'${sample_split}
		echo '---------------------------------------'
		echo

		if [[ ${no_prompt} != 1 ]] ; then
				while true ; do
						read -p "Are all the input arguments correct? (Y/N): " confirm
						case ${confirm} in
								Y|y)
										echo "Confirmed. Initiating..."
										echo
										break
										;;
								N|n)
										echo "Terminated."
										echo
										exit 1
										;;
								* )
										echo "Please enter Y or N."
										echo
										;;
						esac
				done
		fi
fi

##-------------
##Step0-8: Output Folders Creation
##-------------
mkdir -p ${out_dir}/${joint_name}
mkdir -p ${out_dir}/${joint_name}/{Script,LOG,GVCF,VCF,VQSR,QC,QC/FILTERED}



##-------------
##Step1: Combine GVCFs
##-------------
cat <<EOL > ${out_dir}/${joint_name}/Script/1_${joint_name}_combine_vcfs.sh
#!/bin/bash
##-------------
##Step1-0: Create Variant Arguments
##-------------
samples_argument="$( echo ${sample_list[*]} | sed 's/ /\n/g' | sed 's/^\(.*\)/-V \1\/GVCF\/\1_GATK.g.vcf \\/g' )"

##-------------
##Step1-1: Merge GVCFs
##-------------
present_dir=\$(pwd)
cd ${out_dir}
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T CombineGVCFs \
-R ${ref_genome} \
\${samples_argument} \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o ${out_dir}/${joint_name}/GVCF/${joint_name}_GATK.g.vcf \
-log ${out_dir}/${joint_name}/LOG/1_combine_gvcfs.log
cd \${present_dir}

EOL



##-------------
##Step2: Joint Genotype
##-------------
cat <<EOL > ${out_dir}/${joint_name}/Script/2_${joint_name}_joint_genotype.sh
#!/bin/bash
##-------------
##Step2: Joint Genotype
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/GVCF/${joint_name}_GATK.g.vcf \
-nt 8 \
${bed_argument} \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o  ${out_dir}/${joint_name}/VCF/${joint_name}_RAW.vcf \
-log ${out_dir}/${joint_name}/LOG/2_${joint_name}_genotype_gvcf.log

EOL


##-------------
##Step3: Variant Quality Score Recalibration
##-------------
cat <<EOL > ${out_dir}/${joint_name}/Script/3_${joint_name}_recalibrate_variant.sh
#!/bin/bash
##-------------
##Step3-1-1: Select SNVs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/VCF/${joint_name}_RAW.vcf \
-selectType SNP \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-o ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.vcf \
-log ${out_dir}/${joint_name}/LOG/3-1-1_${joint_name}_VQSR_select_SNV.log

##-------------
##Step3-1-2: Select Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/VCF/${joint_name}_RAW.vcf \
-selectType INDEL \
-selectType MNP \
-selectType MIXED \
-selectType SYMBOLIC \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-o ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.vcf \
-log ${out_dir}/${joint_name}/LOG/3-1-2_${joint_name}_VQSR_select_INDEL.log

##-------------
##Step3-2-1: Recalibrate SNVs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref_genome} \
--input:VCF ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.vcf \
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
-recalFile ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.recal \
-tranchesFile ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.tranches \
-rscriptFile ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.R \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-log ${out_dir}/${joint_name}/LOG/3-1-1_${joint_name}_VQSR_snv_recalibration.log

##-------------
##Step3-2-2: Recalibrate Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref_genome} \
--input:VCF ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.vcf \
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
-recalFile ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.recal \
-tranchesFile ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.tranches \
-rscriptFile ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.R \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-log ${out_dir}/${joint_name}/LOG/3-2-2_${joint_name}_VQSR_indel_recalibration.log

##-------------
##Step3-3-1: Apply SNV Recalibration
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref_genome} \
--input:VCF ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.vcf \
-ts_filter_level 99.9 \
-recalFile ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.recal \
-tranchesFile ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV.tranches \
-mode SNP \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV_RECAL_APPLIED.vcf \
-log ${out_dir}/${joint_name}/LOG/3-3-1_${joint_name}_VQSR_apply_snv_recalibration.log

##-------------
##Step3-3-2: Apply Indels Recalibration
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref_genome} \
--input:VCF ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.vcf \
-ts_filter_level 99.9 \
-recalFile ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.recal \
-tranchesFile ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL.tranches \
-mode INDEL \
-dt NONE \
--disable_auto_index_creation_and_locking_when_reading_rods \
-o ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL_RECAL_APPLIED.vcf \
-log ${out_dir}/${joint_name}/LOG/3-3-2_${joint_name}_VQSR_apply_indel_recalibration.log

##-------------
##Step3-4: Combine SNVs + Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T CombineVariants \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/VQSR/${joint_name}_SNV_RECAL_APPLIED.vcf \
--variant ${out_dir}/${joint_name}/VQSR/${joint_name}_INDEL_RECAL_APPLIED.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--genotypemergeoption UNSORTED \
-o ${out_dir}/${joint_name}/VQSR/${joint_name}_FINAL_VQSR.vcf \
-log ${out_dir}/${joint_name}/LOG/3-4_${joint_name}_VQSR_combine_variants.log

EOL



##-------------
##Step4: SNP Quality Control
##-------------
cat <<EOL > ${out_dir}/${joint_name}/Script/4_${joint_name}_SNV_quality_control.sh
#!/bin/bash
##-------------
##Step4-1-1: Extract SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/VQSR/${joint_name}_FINAL_VQSR.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType SNP \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${joint_name}/QC/${joint_name}_RAW_SNV.vcf \
-log ${out_dir}/${joint_name}/LOG/4-1-1_${joint_name}_QC_select_snv.log

##-------------
##Step4-1-2: Extract Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/VQSR/${joint_name}_FINAL_VQSR.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
-selectType INDEL \
-selectType MNP \
-selectType MIXED \
-selectType SYMBOLIC \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${joint_name}/QC/${joint_name}_RAW_INDEL.vcf \
-log ${out_dir}/${joint_name}/LOG/4-1-2_${joint_name}_QC_select_INDEL.log

##-------------
##Step4-2-1: Annotate SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/QC/${joint_name}_RAW_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--dbsnp ${DBSNP} \
-L ${out_dir}/${joint_name}/QC/${joint_name}_RAW_SNV.vcf \
-A GCContent \
-A VariantType \
-dt NONE \
-nt 1 \
-o  ${out_dir}/${joint_name}/QC/${joint_name}_RAW_SNV_ANNOTATED.vcf \
-log ${out_dir}/${joint_name}/LOG/4-2-1_${joint_name}_QC_snv_annotation.log

##-------------
##Step4-3-1: Filter SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/QC/${joint_name}_RAW_SNV_ANNOTATED.vcf \
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
-o ${out_dir}/${joint_name}/QC/${joint_name}_FILTERED_SNV.vcf \
-log ${out_dir}/${joint_name}/LOG/4-3-1_${joint_name}_QC_filter_snv.log

##-------------
##Step4-4-1: Clean SNPs
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/QC/${joint_name}_FILTERED_SNV.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--excludeFiltered \
-nt 1 \
-o  ${out_dir}/${joint_name}/QC/FILTERED/${joint_name}_CLEAN_SNV.vcf \
-log ${out_dir}/${joint_name}/LOG/4-4-1_${joint_name}_QC_clean_snv.log

##-------------
##Step4-5: Combine SNVs + Indels
##-------------
java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
-T CombineVariants \
-R ${ref_genome} \
--variant ${out_dir}/${joint_name}/QC/FILTERED/${joint_name}_CLEAN_SNV.vcf \
--variant ${out_dir}/${joint_name}/QC/${joint_name}_RAW_INDEL.vcf \
--disable_auto_index_creation_and_locking_when_reading_rods \
--genotypemergeoption UNSORTED \
-o ${out_dir}/${joint_name}/QC/FILTERED/${joint_name}_CLEAN_SNV+INDEL.vcf \
-log ${out_dir}/${joint_name}/LOG/4-5_${joint_name}_QC_combine_variants.log

EOL



##-------------
##Step5: Sample Spliting
##-------------
cat <<EOL > ${out_dir}/${joint_name}/Script/5_${joint_name}_split_samples.sh
#!/bin/bash
##-------------
##Step5-0: Create Sample List
##-------------
columns=(\$(grep '^#[^#]' ${out_dir}/${joint_name}/QC/FILTERED/${joint_name}_CLEAN_SNV+INDEL.vcf))
sample_list=(\${columns[*]:9})

##-------------
##Step5-1: Split File
##-------------
for sample_name in \${sample_list[*]} ; do
		mkdir -p ${out_dir}/\${sample_name}/QC/FILTERED
		java -Xmx${java_mem} -jar ${gatk_dir}/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R ${ref_genome} \
		-sn \${sample_name} \
		--variant ${out_dir}/${joint_name}/QC/FILTERED/${joint_name}_CLEAN_SNV+INDEL.vcf \
		-o ${out_dir}/\${sample_name}/QC/FILTERED/\${sample_name}_CLEAN_SNV+INDEL.vcf \
		-log ${out_dir}/${joint_name}/LOG/5_${joint_name}_split_sample_\${sample_name}.log
done

EOL





##-------------
##MASTER SCRIPT
##-------------
cat <<EOL > ${out_dir}/${joint_name}/Script/${joint_name}_GATK.sh
#!/bin/bash
##-------------
##Joint Vaiant Calling
##-------------
bash ${out_dir}/${joint_name}/Script/1_${joint_name}_combine_vcfs.sh
bash ${out_dir}/${joint_name}/Script/2_${joint_name}_joint_genotype.sh
bash ${out_dir}/${joint_name}/Script/3_${joint_name}_recalibrate_variant.sh
bash ${out_dir}/${joint_name}/Script/4_${joint_name}_SNV_quality_control.sh

if [[ ${sample_split} == 'YES' ]] ; then
		bash ${out_dir}/${joint_name}/Script/5_${joint_name}_split_samples.sh
fi

EOL





##-------------
##EXECUTION
##-------------
if [[ ${no_exec} != 1 ]] ; then
		bash ${out_dir}/${joint_name}/Script/${joint_name}_GATK.sh
fi
