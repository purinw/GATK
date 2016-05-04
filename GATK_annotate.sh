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
##ii>	Run the command 'bash /path/to/GATK_annotate.sh [options] input_file'
##-------------





##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
snpeff_dir=/mnt/data2/home2/purinw/snpEff

##-------------
##Step0-2: References
##-------------
DBSNP=/home/bhoom/data/hg19/gatk_bundle/dbsnp_138.hg19.vcf
DBNSFP=${snpeff_dir}/data/dbNSFP.txt.gz
GWASCATALOG=${snpeff_dir}/data/gwascatalog.txt
PHASTCONS=${snpeff_dir}/data/phastCons
CLINVAR=${snpeff_dir}/data/clinvar.vcf 

##-------------
##Step0-3: Other Parametres
##-------------
java_mem=30g

##-------------
##Step0-4: Input Arguments
##-------------
while test $# -gt 0; do
        case "$1" in
                -h|--help)
                        echo "Usage: bash /path/to/GATK_individual.sh [options] input_file"
                        echo "This script performs the entire variant-calling process upon one sample, following the Genome Analysis Toolkit (GATK)'s pipeline."
                        echo ""
                        echo "Options:"
                        echo "-h, --help				display this help and exit"
						echo "-O, --out-dir"
                        # echo "-o, --out-file				"
                        exit 0
                        ;;
                # -o|--out-file)
                        # shift
						# if [[ $( echo $1 | grep -c '\/' ) -gt 0 ]]
								# then out_dir=$( echo $1 | sed 's/\/[^/]*$//' )
								# else out_dir=$( pwd )
						# fi
						# out_file=$( echo $1 | sed 's/.*\///' )
                        # shift
                        # ;;
                -O|--out-dir)
                        shift
						out_dir=$( echo $1 | sed 's/\/$//' )
                        shift
                        ;;
				*)
						in_file=$1
						shift
						;;
		esac
done

if [[ ${out_dir} == '' ]]
		then out_dir=$( echo ${in_file} | sed 's/\/[^/]*$//' )
fi
if [[ ${out_file} == '' ]]
		then out_file=$( echo ${in_file} | sed 's/.*\/\([^/]*\).vcf$/\1_annotated.vcf/' )
fi
out_full=${out_dir}/${out_file}
prefix=${out_full}.temp

##-------------
##Step0-5: Sample Verification
##-------------
if [ ! -e ${in_file} ] ; then
		echo 'Invalid INPUT FILE: '${in_file}'. Terminated.'
		echo
		exit 1
fi

##-------------
##Step0-5: Summarisation
##-------------
echo
echo '---------------------------------------'
echo 'VARIANT ANNOTATION'
echo 'INPUT FILE =			'${in_file}
echo 'OUTPUT FILE =			'${out_full}
echo '---------------------------------------'
echo

##-------------
##Step0-6: User's Confirmation Prompt
##-------------
while true; do
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

##-------------
##Step0-7: Output Folders Creation
##-------------
mkdir -p ${out_dir}



##-------------
##Step1: dbSNP
##-------------
echo '1/6 dbSNP Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar annotate ${DBSNP} ${in_file} > ${prefix}.dbsnp.vcf
echo '1/6 dbSNP Annotation Completed'



##-------------
##Step2: dbNSFP
##-------------
echo '2/6 dbNSFP Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar dbnsfp -v -db ${DBNSFP} ${prefix}.dbsnp.vcf > ${prefix}.dbsnp.dbnsfp.vcf
echo '2/6 dbNSFP Annotation Completed'



##-------------
##Step3: gwasCat
##-------------
echo '3/6 gwasCat Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar gwasCat -db ${GWASCATALOG} ${prefix}.dbsnp.dbnsfp.vcf > ${prefix}.dbsnp.dbnsfp.gwascat.vcf
echo '3/6 gwasCat Annotation Completed'



##-------------
##Step4: PhastCons
##-------------
echo '4/6 PhastCons Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar phastCons ${PHASTCONS} ${prefix}.dbsnp.dbnsfp.gwascat.vcf > ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.vcf
echo '4/6 PhastCons Annotation Completed'



##-------------
##Step5: ClinVar
##-------------
echo '5/6 ClinVar Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar annotate ${CLINVAR} ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.vcf > ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.clinvar.vcf
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
mv ${prefix}.dbsnp.dbnsfp.gwascat.phastcons.clinvar.snpeff.vcf ${out_full}
rm ${prefix}.dbsnp*vcf
