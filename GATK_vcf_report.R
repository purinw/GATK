input_files <- commandArgs(trailingOnly = TRUE)

vcf_file <- input_files[1]
count_file <- input_files[2]
gene_file <- input_files[3]
output_prefix <- gsub(pattern=".reduced.vcf$",replacement = "",vcf_file)

IGNORED_FUNC <- c("intron_variant","synonymous_variant","intergenic_region","downstream_gene_variant","upstream_gene_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant")
IGNORED_GENO <- c("0/0","./.","0/.")

ANNOTATION_FIELD <- list(
  DBNSFP_FIELD = c("dbNSFP_1000Gp1_AF","dbNSFP_1000Gp1_AFR_AF","dbNSFP_1000Gp1_AMR_AF","dbNSFP_1000Gp1_ASN_AF","dbNSFP_1000Gp1_EUR_AF","dbNSFP_ExAC_AF","dbNSFP_ExAC_Adj_AF","dbNSFP_ExAC_AFR_AF","dbNSFP_ExAC_AMR_AF","dbNSFP_ExAC_EAS_AF","dbNSFP_ExAC_FIN_AF","dbNSFP_ExAC_NFE_AF","dbNSFP_ExAC_SAS_AF","dbNSFP_ESP6500_AA_AF","dbNSFP_ESP6500_EA_AF","dbNSFP_SIFT_pred","dbNSFP_Polyphen2_HDIV_pred","dbNSFP_Polyphen2_HVAR_pred","dbNSFP_MutationTaster_pred","dbNSFP_MutationAssessor_pred","dbNSFP_FATHMM_pred","dbNSFP_LRT_pred","dbNSFP_CADD_phred","dbNSFP_GERP++_RS","dbNSFP_GERP++_NR","dbNSFP_phastCons100way_vertebrate","dbNSFP_PROVEAN_pred","dbNSFP_Interpro_domain","dbNSFP_Uniprot_acc","dbNSFP_MetaSVM_pred"),
  CLINVAR_FIELD = c("CLNALLE","CLNSIG","CLNDBN","CLNREVSTAT","CLNACC","CLNSRC","CLNSRCID","CLNDSDB","CLNDSDBID"),
  GWASCAT_FIELD = c("GWASCAT_TRAIT","GWASCAT_OR_BETA","GWASCAT_PUBMED_ID","GWASCAT_P_VALUE")
)

VCF <- read.delim(vcf_file, stringsAsFactors=FALSE)

n_field <- count.fields(count_file, sep = "\t")
n_row <- length(n_field)-1
COUNTS <- read.table(count_file, stringsAsFactors=FALSE,fill=TRUE,col.names=1:max(n_field),skip=1)
GENES <- readLines(gene_file)

for(i in (min(n_field)+1):max(n_field)){
  COUNTS[(COUNTS[,i]=="")|is.na(COUNTS[,i]),i] <- "X:0"
}

ALT <- as.matrix(COUNTS[,(min(n_field)+1):max(n_field)])
N_CHROM_ALT <- sapply(1:n_row,function(N){temp <- ALT[N,]; alt_chrom <- as.numeric(unlist(strsplit(temp,split=":"))) ; return(sum(alt_chrom,na.rm=TRUE))})

n_allele <- dim(ALT)[2]
for(i in 1:n_allele){
  ALT[(ALT[,i]=="X:0"),i] <- ""
}

INFO <- VCF$INFO
SNPEFF_ANN <- strsplit(gsub(pattern = ".*ANN=",replacement="",INFO),split=",")

ALLELES <- character(n_row)
CANDIDATE_GENES <- character(n_row)
FUNCTIONS <- character(n_row)
ANNOTATIONS <- matrix(".",ncol=length(unlist(ANNOTATION_FIELD)),nrow=n_row)

ELIGIBLE <- logical(n_row)

for(i in 1:n_row){
  scatter <- lapply(SNPEFF_ANN[[i]],FUN=function(one_SNPEFF){unlist(strsplit(one_SNPEFF,split="|",fixed=T))})
  N_scatter <- length(scatter)
  allele_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][1]})
  gene_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][4]})
  func_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][2]})
  
  relevant_ind <- (gene_scatter %in% GENES)&!(func_scatter %in% IGNORED_FUNC)
  allele_scatter <- allele_scatter[relevant_ind]
  gene_scatter <- gene_scatter[relevant_ind]
  func_scatter <- func_scatter[relevant_ind]
  scatter <- strsplit(unique(paste(allele_scatter,gene_scatter,func_scatter,sep=":")),split=":")
  N_scatter <- length(scatter)
  
  if(N_scatter>0){
    ELIGIBLE[i] <- TRUE
    allele_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][1]})
    gene_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][2]})
    func_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][3]})
    ALLELES[i] <- paste(allele_scatter,collapse=";")
    CANDIDATE_GENES[i] <- paste(gene_scatter,collapse=";")
    FUNCTIONS[i] <- paste(func_scatter,collapse=";")
    info_scatter <- unlist(strsplit(INFO[i],split=";"))
    ANNOTATIONS[i,] <- sapply(unlist(ANNOTATION_FIELD),FUN=function(x){ind <- grepl(info_scatter,pattern=paste(x,"=",sep=""),fixed=T) ; if(any(ind)){return(gsub(pattern="[^=]+=","",info_scatter[ind]))}else{return(".")}})
  }
}

ANNOTATION_FIELD$DBNSFP_FIELD <- gsub(pattern="^dbNSFP_",replacement="",ANNOTATION_FIELD$DBNSFP_FIELD)
ANNOTATION_FIELD$CLINVAR_FIELD <- gsub(pattern="^CLN",replacement = "ClinVar_",ANNOTATION_FIELD$CLINVAR_FIELD)
dimnames(ANNOTATIONS)[[2]] <- unlist(ANNOTATION_FIELD)
dimnames(ALT)[[2]] <- paste("ALLELE.FREQ",1:n_allele,sep=".")

OUTPUT <- cbind(VCF[,1:5],N_CHROM_ALT,ALT,ALLELES,CANDIDATE_GENES,FUNCTIONS,ANNOTATIONS,VCF[,6:7],stringsAsFactors=FALSE)
OUTPUT <- OUTPUT[ELIGIBLE,]
INFO <- INFO[ELIGIBLE]

for(i in rev(grep("^ALLELE.FREQ.",dimnames(OUTPUT)[[2]]))){
  if(all(OUTPUT[,i]=="")) OUTPUT[,i] <- NULL}

GENO <- gsub(":.*$","",as.matrix(VCF[ELIGIBLE,-(1:9)]))

sample_names <- dimnames(GENO)[[2]]
dimnames(GENO)[[2]] <- paste("SampleID",sample_names,sep="=")
write.table(cbind(OUTPUT,GENO,INFO),paste(output_prefix,"report_wide.tsv",sep="_"),sep="\t",row.names=F,quote=F)  

GENO <- data.frame(SampleID=rep(sample_names,each=sum(ELIGIBLE)),GENOTYPE=as.character(GENO),stringsAsFactors = F)
OUTPUT <- cbind(OUTPUT,GENO,INFO)

OUTPUT <- OUTPUT[!(OUTPUT$GENOTYPE %in% IGNORED_GENO),]

write.table(OUTPUT,paste(output_prefix,"report_long.tsv",sep="_"),sep="\t",row.names=F,quote=F)
