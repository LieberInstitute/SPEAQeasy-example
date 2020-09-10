library(VariantAnnotation)

# load genotype data
genotyped = readVcf("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/Sample_file_check_vcffile/LIBD_Brain_merged_topmed_RNASeq_annotated_variants_043020.vcf",
	genome="hg38")


#load brains pd sheet with rna from
load('/dcl01/lieber/ajaffe/lab/SPEAQeasy-example/sample_selection/pd_example.Rdata')

#load brain sentrix to get ID's
brain_sentrix<- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/brain_sentrix.csv")

#ad sentrix_id to pd table
pd_example<-merge(pd_example,brain_sentrix)


#subset for samples that match those brains
snpsGeno_example<-genotyped[,colnames(genotyped) %in% pd_example$ID]
dt<-rowRanges(snpsGeno_example)
rownames(snpsGeno_example)<-paste0(as.character(seqnames(dt)),":",as.character(start(dt)),"_",as.character(dt$REF),"/",as.character(unlist(dt$ALT)))
save(snpsGeno_example,file="snpsGeno_example.RData")



