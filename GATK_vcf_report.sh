#!/bin/bash
##-------------
##Purin Wangkiratikant, MAR 2016
##Clinical Database Centre, Institute of Personalised Genomics and Gene Therapy (IPGG)
##Faculty of Medicine Siriraj Hospital, Mahidol University, Bangkok, Thailand
##-------------
##This script creats a tabular report of the variants within a VCF file.
##
##How to run:
##i>	Change all the directories and files within Step0 of this script accordingly.
##ii>	Check if the annotation process for the input file has been accomplished.
##iii>	Run the command 'bash /path/to/GATK_vcf_report.sh [options] input_file'
##-------------





##-------------
##Step0: Initialisation
##-------------

##-------------
##Step0-1: Directories
##-------------
snpeff_dir=/mnt/data2/home2/purinw/snpEff
script_dir=/mnt/data2/home2/purinw/Scripts
gene_filtration_script=${script_dir}/GATK_filter_gene.sh
r_report_script=${script_dir}/GATK_vcf_report.R

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
						echo "Usage: bash /path/to/GATK_vcf_report.sh [options] input_file"
                        echo ""
                        echo "This script creats a tabular report of the variants within a VCF file."
                        echo ""
                        echo "Options:"
                        echo "-h, --help				display this help and exit"
						echo "-v, --version				display version of this script"
						echo "-XS, --no-summary			suppress the command summary before execution"
						echo "-XP, --no-prompt			suppress the user prompt before execution, only when the command summary is displayed"
						# echo "-O, --out-dir	OUT_DIR			specify output directory, the same as input's by default"
                        echo "-p, --prefix		PREFIX		specify output prefix to be used, by default the input file's name with the suffix \"report\" appended"
						echo "-G, --gene-file		GENE_FILE	specify gene file containing gene names to be used in the filtration"
						echo "-f, --filter				include only variants within the specified gene files"
						echo "-w, --wide				output only the wide report, with all the samples' genotype in each row"
						echo "-l, --long				output only the long report, with only one sample's genotype in each row"
						echo "-ph, --phenotype			include \"PHENOTYPE\" field in the report(s), retrieved from each gene file's name"
						echo "-t, --tier				include \"TIER\" field in the long report, calculated by using clinical evidences, genotype, and background frequencies"
						echo ""
                        exit 0
                        ;;
				-v|--version)
						echo ""
						echo "GATK_vcf_report.sh"
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
                -p|--prefix)
                        shift
						prefix=$1
						out_dir=$( echo $1 | sed 's/\/[^\/]*$//')/.
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
				-f|--filter)
						filtration=YES
						shift
						;;
				-w|--wide)
						wide_flag=1
						shift
						;;
				-l|--long)
						long_flag=1
						shift
						;;
				-ph|--phenotype)
						phenotype=YES
						shift
						;;
				-t|--tier)
						tier=YES
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
if [[ ! -v prefix ]] ; then
		prefix=$( echo ${in_file} | sed 's/.vcf$//' )
fi
if [[ ! -v filtration ]] ; then
		filtration=NO
fi
if [[ -v wide_flag || ! -v long_flag ]] ; then
		wide=YES
		else
		wide=NO
fi
if [[ -v long_flag || ! -v wide_flag ]] ; then
		long=YES
		else
		long=NO
fi
if [[ ! -v phenotype ]] ; then
		phenotype=NO
fi
if [[ ! -v tier ]] ; then
		tier=NO
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
if [[ ${filtration} == 'YES' && ${gene_file_list} == '' ]] ; then
		echo
		echo 'Invalid INPUT: Filter without gene files'
		echo 'Terminated.'
		echo
		exit 1
fi
if [[ ${phenotype} == 'YES' && ${gene_file_list} == '' ]] ; then
		echo
		echo 'Invalid INPUT: Phenotype without gene files'
		echo 'Terminated.'
		echo
		exit 1
fi
if [[ ${tier} == 'YES' && ${long} == 'NO' ]] ; then
		echo
		echo 'Invalid INPUT: Tier without long output'
		echo 'Terminated.'
		echo
		exit 1
fi


##-------------
##Step0-6: Summarisation
##-------------
if [[ ${no_summary} != 1 ]] ; then
		echo
		echo '---------------------------------------'
		echo 'VARIANT TABULAR REPORT'
		echo 'INPUT FILE =			'${in_file}
		echo 'OUTPUT PREFIX =			'${prefix}
		echo 'GENE FILE COUNT =		'${#gene_file_list[*]}
		echo 'GENE FILES =			'${gene_file_list[*]}
		echo 'GENE FILTRATION =		'${filtration}
		echo 'WIDE OUTPUT REPORT =		'${wide}
		echo 'LONG OUTPUT REPORT =		'${long}
		echo 'PHENOTYPE FIELD =		'${phenotype}
		echo 'TIER FIELD =			'${tier}
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



##-------------
##Step1: Gene Filtration
##-------------
if [[ ${filtration} == 'YES' ]] ; then
		gene_file_argument="$( echo ${gene_file_list[*]} | sed 's/ /\n/g' | sed 's/^\(.*\)/-G \1 /g' )"
		bash ${gene_filtration_script} ${in_file} -o ${prefix}_filtered.vcf -XS ${gene_file_argument}
		in_file=${prefix}_filtered.vcf
fi



##-------------
##Step2: VCF File Reduction
##-------------
( grep -v '^##' ${in_file} | sed 's/^#//' ) > ${prefix}.reduced.vcf



##-------------
##Step3: Allele Frequency Computation
##-------------
vcftools --vcf ${in_file} --counts --out ${prefix}.alleles



##-------------
##Step4: Creating Gene-Phenotype File
##-------------
(for gene_file in ${gene_file_list[*]} ; do awk '{print $1"\t"FILENAME}' ${gene_file} ; done ) | sed 's/\t\(.*\)\..*/\t\1/g' | sed 's/\t.*\/\(.*\)/\t\1/g' > ${prefix}_gene_pheno.txt



##-------------
##Step4: Tabulating Results (Via R)
##-------------
Rscript ${r_report_script} ${prefix}.reduced.vcf ${prefix}.alleles.frq.count ${wide} ${long} ${phenotype} ${tier} ${filtration} ${prefix}_gene_pheno.txt



##-------------
##Step5: Removal of Temporary Files
##-------------
rm ${prefix}.reduced.vcf
rm ${prefix}.alleles.log
rm ${prefix}_gene_pheno.txt