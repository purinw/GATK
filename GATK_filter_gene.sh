#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script extracts variants within a VCF file that are related to genes specified by given gene files or given gene names.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Check if the annotation process for the input file has been accomplished.
##iii>	Run the command 'bash /path/to/GATK_filter_gene.sh [options] input_file'
##-------------





##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
snpeff_dir=/mnt/data2/home2/purinw/snpEff

##-------------
##Step0-2: Other Parametres
##-------------
java_mem=30g

##-------------
##Step0-3: Input Arguments
##-------------
while test $# -gt 0 ; do
        case "$1" in
                -h|--help)
                        echo ""
						echo "Usage: bash /path/to/GATK_filter_gene.sh [options] input_file"
                        echo ""
                        echo "This script extracts variants within a VCF file that are related to genes specified by given gene files or given gene names."
                        echo ""
                        echo "Options:"
                        echo "-h, --help				display this help and exit"
						echo "-v, --version				display version of this script"
						echo "-XS, --no-summary			suppress the command summary before execution"
						echo "-XP, --no-prompt			suppress the user prompt before execution, only when the command summary is displayed"
						# echo "-O, --out-dir	OUT_DIR			specify output directory, the same as input's by default"
                        echo "-o, --out-file		OUT_FILE	specify output file, by default the input file's name with the suffix \"filtered\" appended"
						echo "-r, --replace				instruct the output filtered file to replace the input file, treating OUT_FILE as the temporary file if both activated"
						echo "-G, --gene-file		GENE_FILE	specify gene file containing gene names to be used in the filtration"
						echo "-g, --gene-name		GENE_NAME	specify gene names to be used in filtration"
						echo ""
                        exit 0
                        ;;
				-v|--version)
						echo ""
						echo "GATK_filter_gene.sh"
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
				-G|--gene-file)
						shift
						gene_file_list+=( $1 )
						shift
						;;
				-g|--gene-name)
						shift
						gene_list+=( $1 )
						shift
						;;
				*)
						in_file=$1
						shift
						;;
		esac
done

##-------------
##Step0-4: Default Value Setting
##-------------
if [[ ! -v out_dir ]] ; then
		out_dir=$( pwd )
fi
if [[ ! -v out_file ]] ; then
		out_file=$( echo ${in_file} | sed 's/.vcf$/_filtered.vcf/' )
fi
if [[ ! -v replacement ]] ; then
		replacement=NO
fi

##-------------
##Step0-5: Input Verification
##-------------
if [[ ! -e ${in_file} ]] ; then
		echo
		echo 'Invalid INPUT FILE: '${in_file}
		echo 'Terminated.'
		echo
		exit 1
fi
for gene_file in ${gene_file_list[*]} ; do
		if [[ ! -e ${gene_file} ]] ; then
				echo
				echo 'Invalid GENE FILE: '${gene_file}
				echo 'Terminated.'
				echo
				exit 1
		fi
done

##-------------
##Step0-6: Summarisation & User's Confirmation Prompt
##-------------
if [[ ${no_summary} != 1 ]] ; then
		echo
		echo '---------------------------------------'
		echo 'GENE FILTRATION'
		echo 'INPUT FILE =			'${in_file}
		echo 'OUTPUT FILE =			'${out_file}
		echo 'REPLACEMENT =			'${replacement}
		echo 'GENE FILE COUNT =		'${#gene_file_list[*]}
		echo 'GENE FILES =			'${gene_file_list[*]}
		echo 'GENE NAME COUNT =		'${#gene_list[*]}
		echo 'GENE NAMES =			'${gene_list[*]}
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
##Step0-7: Output Folder Creation
##-------------
mkdir -p ${out_dir}



# ##-------------
# ##Step1: SnpEff Annotation Verificiation
# ##-------------
# if [[ $( grep "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' \">" ${in_file} | wc -l ) -eq 0 ]] ; then
		# java -Xmx${java_mem} -jar ${snpeff_dir}/snpEff.jar hg19 ${in_file} > ${out_file}.temp.vcf
# fi


##-------------
##Step1: Manipulating Arguments
##-------------
for gene_file in ${gene_file_list[*]} ; do
		gene_list+=($( cat ${gene_file} ))
done
gene_argument=$(echo ${gene_list[*]} | sed 's/ /\n/g' | sed "s/\(.*\)/(ANN[*].GENE = '\1')/g" | paste -sd '||')



##-------------
##Step2: Gene Filtration
##-------------
java -jar ${snpeff_dir}/SnpSift.jar filter "${gene_argument}" ${in_file} > ${out_file}



##-------------
##Step4: Replacement
##-------------
if [[ ${replacement} == 'YES' ]] ; then
		mv ${out_file} ${in_file}
fi
