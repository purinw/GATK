input_arguments <- commandArgs(trailingOnly = TRUE)

vcf_file <- input_arguments[1]
count_file <- input_arguments[2]
wide_flag <- input_arguments[3]
long_flag <- input_arguments[4]
pheno_flag <- input_arguments[5]
tier_flag <- input_arguments[6]
filter_flag <- input_arguments[7]
if(pheno_flag=='YES' | filter_flag=='YES')
  gene_file <- input_arguments[8]

output_prefix <- gsub(pattern=".reduced.vcf$",replacement = "",vcf_file)

IGNORED_FUNC <- c("intron_variant","synonymous_variant","intergenic_region","downstream_gene_variant","upstream_gene_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_exon_variant","splice_region_variant&intron_variant","splice_region_variant&synonymous_variant")
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
##GENES <- readLines(gene_file)
if(pheno_flag=='YES' | filter_flag=='YES'){
  candidate_genes_pheno <- read.delim(gene_file, header=FALSE, stringsAsFactors=FALSE)
  
  if(pheno_flag=='YES')
    pheno <- sapply(unique(candidate_genes_pheno[,1]),function(gene){return(paste(candidate_genes_pheno[candidate_genes_pheno[,1]==gene,2],collapse=";"))})
  
  if(filter_flag=='YES')
    GENES <- unique(candidate_genes_pheno[,1])
}

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
PHENOTYPE <- character(n_row)
ANNOTATIONS <- matrix(".",ncol=length(unlist(ANNOTATION_FIELD)),nrow=n_row)

ELIGIBLE <- logical(n_row)

for(i in 1:n_row){
  scatter <- lapply(SNPEFF_ANN[[i]],FUN=function(one_SNPEFF){unlist(strsplit(one_SNPEFF,split="|",fixed=T))})
  N_scatter <- length(scatter)
  allele_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][1]})
  gene_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][4]})
  func_scatter <- sapply(1:N_scatter,FUN=function(x){scatter[[x]][2]})
  
  relevant_ind <- !(func_scatter %in% IGNORED_FUNC)
  if(filter_flag=='YES')
    relevant_ind <- relevant_ind & (gene_scatter %in% GENES)
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
    if(pheno_flag=='YES')
      PHENOTYPE[i] <- paste(sapply(gene_scatter,FUN=function(x){paste(pheno[x],collapse=",")}),collapse=";")
    FUNCTIONS[i] <- paste(func_scatter,collapse=";")
    info_scatter <- unlist(strsplit(INFO[i],split=";"))
    ANNOTATIONS[i,] <- sapply(unlist(ANNOTATION_FIELD),FUN=function(x){ind <- grepl(info_scatter,pattern=paste(x,"=",sep=""),fixed=T) ; if(any(ind)){return(gsub(pattern="[^=]+=","",info_scatter[ind]))}else{return(".")}})
  }
}

ANNOTATION_FIELD$DBNSFP_FIELD <- gsub(pattern="^dbNSFP_",replacement="",ANNOTATION_FIELD$DBNSFP_FIELD)
ANNOTATION_FIELD$CLINVAR_FIELD <- gsub(pattern="^CLN",replacement = "ClinVar_",ANNOTATION_FIELD$CLINVAR_FIELD)
dimnames(ANNOTATIONS)[[2]] <- unlist(ANNOTATION_FIELD)
dimnames(ALT)[[2]] <- paste("ALLELE.FREQ",1:n_allele,sep=".")
if(pheno_flag=='YES'){
  OUTPUT <- cbind(VCF[,1:5],N_CHROM_ALT,ALT,ALLELES,CANDIDATE_GENES,PHENOTYPE,FUNCTIONS,ANNOTATIONS,VCF[,6:7],stringsAsFactors=FALSE)
}else{
  OUTPUT <- cbind(VCF[,1:5],N_CHROM_ALT,ALT,ALLELES,CANDIDATE_GENES,FUNCTIONS,ANNOTATIONS,VCF[,6:7],stringsAsFactors=FALSE)
}
OUTPUT <- OUTPUT[ELIGIBLE,]
INFO <- INFO[ELIGIBLE]

for(i in rev(grep("^ALLELE.FREQ.",dimnames(OUTPUT)[[2]]))){
  if(all(OUTPUT[,i]=="")) OUTPUT[,i] <- NULL
}

GENO <- gsub(":.*$","",as.matrix(VCF[ELIGIBLE,-(1:9)]))

sample_names <- colnames(VCF)[-(1:9)]
colnames(GENO) <- paste("SampleID",sample_names,sep="=")
OUTPUT_WIDE <- cbind(OUTPUT,GENO,INFO)
if(wide_flag=='YES'){
  write.table(OUTPUT_WIDE,paste(output_prefix,"report_wide.tsv",sep="_"),sep="\t",row.names=F,quote=F)  
}
if(long_flag=='YES'){
  GENO_LONG <- data.frame(SampleID=rep(sample_names,each=sum(ELIGIBLE)),GENOTYPE=as.character(GENO),stringsAsFactors = F)
  if(tier_flag=='YES'){
    OUTPUT_LONG <- cbind(OUTPUT,GENO_LONG)
    
    # Tier Calculation
    clinvar.db <- grepl(pattern="cardiomyopathy",OUTPUT_LONG$ClinVar_DBN,ignore.case = T)|grepl(pattern="arrhythmia",OUTPUT_LONG$ClinVar_DBN,ignore.case = T)|grepl(pattern=".",OUTPUT_LONG$ClinVar_DBN,ignore.case = T,fixed=T)
    clinvar.sig <- (grepl(pattern="[45]",OUTPUT_LONG$ClinVar_SIG)&!grepl(pattern="255",OUTPUT_LONG$ClinVar_SIG))|grepl(pattern=".",OUTPUT_LONG$ClinVar_SIG,ignore.case = T,fixed=T)
    clinvar.nosig <- !grepl(pattern="[45]",OUTPUT_LONG$ClinVar_SIG)|grepl(pattern="255",OUTPUT_LONG$ClinVar_SIG)
    population.10pc <- apply(OUTPUT_LONG[,gsub(pattern="^dbNSFP_",replacement="",ANNOTATION_FIELD$DBNSFP_FIELD[1:15])],1,function(pop_vec){pop_vec <- as.numeric(pop_vec) ; return(any(is.na(pop_vec)|(pop_vec<=.1)))})
    
    heteroz <- OUTPUT_LONG$GENOTYPE %in% c("1/0","0/1")
    homoz <- OUTPUT_LONG$GENOTYPE %in% c("1/1")
    
    TIER <- character(dim(OUTPUT_LONG)[1])
    
    TIER[clinvar.db&clinvar.sig&homoz] <- "Tier3"
    TIER[clinvar.db&clinvar.nosig&population.10pc&homoz] <- "Tier2"
    TIER[clinvar.db&population.10pc&heteroz] <- "Tier1"
    
    OUTPUT_LONG <- cbind(OUTPUT_LONG,TIER,INFO)}
  else{
    OUTPUT_LONG <- cbind(OUTPUT,GENO_LONG,INFO)
  }
  
  OUTPUT_LONG <- OUTPUT_LONG[!(OUTPUT_LONG$GENOTYPE %in% IGNORED_GENO),]
  
  write.table(OUTPUT_LONG,paste(output_prefix,"report_long.tsv",sep="_"),sep="\t",row.names=F,quote=F)
}
