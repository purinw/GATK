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
while test $# -gt 0 ; do
        case "$1" in
                -h|--help)
                        echo ""
						echo "Usage: bash /path/to/GATK_individual.sh [options] input_file"
                        echo ""
                        echo "This script annotates a given VCF file with the following databases:"
                        echo "i>	dbSNP"
                        echo "ii>	dbNSFP"
                        echo "iii>	gwasCat"
                        echo "iv>	PhastCons"
                        echo "v>	ClinVar"
                        echo "vi>	SnpEff"
                        echo ""
                        echo "Options:"
                        echo "-h, --help				display this help and exit"
						echo "-v, --version				display version of this script"
						echo "-XS, --no-summary			suppress the command summary before execution"
						echo "-XP, --no-prompt			suppress the user prompt before execution, only when the command summary is displayed"
						# echo "-O, --out-dir	OUT_DIR			specify output directory, the same as input's by default"
                        echo "-o, --out-file		OUT_FILE	specify output file, by default the input file's name with the suffix \"annotated\" appended"
						echo "-r, --replace				instruct the output annotated file to replace the input file, treating OUT_FILE as the temporary file if both activated"
						echo ""
                        exit 0
                        ;;
				-v|--version)
						echo ""
						echo "GATK_annotate.sh"
                        echo ""
						echo "Updated JUN 2016"
						echo "by"
						echo "PURIN WANGKIRATIKANT"
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
                -o|--out-file)
                        shift
						out_file=$1
						out_dir=$( echo $1 | sed 's/\/[^\/]*$//')/.
                        shift
                        ;;
				-r|--replace)
						replacement=YES
						shift
						;;
                # -O|--out-dir)
                        # shift
						# out_dir=$( echo $1 | sed 's/\/$//' )
                        # shift
                        # ;;
				*)
						in_file=$1
						shift
						;;
		esac
done

##-------------
##Step0-5: Default Value Setting
##-------------
if [[ ! -v out_dir ]] ; then
		out_dir=$( pwd )
fi
if [[ ! -v out_file ]] ; then
		out_file=$( echo ${in_file} | sed 's/.vcf$/_annotated.vcf/' )
fi
if [[ ! -v replacement ]] ; then
		replacement=NO
fi

##-------------
##Step0-6: Input Verification
##-------------
if [[ ! -e ${in_file} ]] ; then
		echo
		echo 'Invalid INPUT FILE: '${in_file}
		echo ${in_file} not found.
		echo 'Terminated.'
		echo
		exit 1
fi

##-------------
##Step0-7: Summarisation & User's Confirmation Prompt
##-------------
if [[ ${no_summary} != 1 ]] ; then
		echo
		echo '---------------------------------------'
		echo 'VARIANT ANNOTATION'
		echo 'INPUT FILE =			'${in_file}
		echo 'OUTPUT FILE =			'${out_file}
		echo 'REPLACEMENT =			'${replacement}
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
##Step0-8: Output Folder Creation
##-------------
mkdir -p ${out_dir}
cp ${in_file} ${out_file}



##-------------
##Step1: dbSNP
##-------------
echo '1/6 dbSNP Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar annotate ${DBSNP} ${out_file} > ${out_file}.temp.vcf
mv ${out_file}.temp.vcf ${out_file}
echo '1/6 dbSNP Annotation Completed'



##-------------
##Step2: dbNSFP
##-------------
echo '2/6 dbNSFP Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar dbnsfp -v -db ${DBNSFP} ${out_file} > ${out_file}.temp.vcf
mv ${out_file}.temp.vcf ${out_file}
echo '2/6 dbNSFP Annotation Completed'



##-------------
##Step3: gwasCat
##-------------
echo '3/6 gwasCat Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar gwasCat -db ${GWASCATALOG} ${out_file} > ${out_file}.temp.vcf
mv ${out_file}.temp.vcf ${out_file}
echo '3/6 gwasCat Annotation Completed'



##-------------
##Step4: PhastCons
##-------------
echo '4/6 PhastCons Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar phastCons ${PHASTCONS} ${out_file} > ${out_file}.temp.vcf
mv ${out_file}.temp.vcf ${out_file}
echo '4/6 PhastCons Annotation Completed'



##-------------
##Step5: ClinVar
##-------------
echo '5/6 ClinVar Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/SnpSift.jar annotate ${CLINVAR} ${out_file} > ${out_file}.temp.vcf
mv ${out_file}.temp.vcf ${out_file}
echo '5/6 ClinVar Annotation Completed'



##-------------
##Step6: SnpEff
##-------------
echo '6/6 SnpEff Annotation Started'
java -Xmx${java_mem} -jar ${snpeff_dir}/snpEff.jar hg19 ${out_file} > ${out_file}.temp.vcf
mv ${out_file}.temp.vcf ${out_file}
echo '6/6 SnpEff Annotation Completed'



##-------------
##Step7: Replacement
##-------------
if [[ ${replacement} == 'YES' ]] ; then
		mv ${out_file} ${in_file}
fi
