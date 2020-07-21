###
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(limma)
library(jaffelab)
library(RColorBrewer)

dir.create("pdfs")
dir.create("tables")

### load data
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda")
load('/dcl01/lieber/ajaffe/lab/SPEAQeasy-example/sample_selection/pd_example.Rdata')
rm(getRPKM)
## keep subset of sample
Index = pd_example$SAMPLE_ID
rse_gene = rse_gene[,colData(rse_gene)$SAMPLE_ID %in% Index]

## cell PCs
cellPca = prcomp(as.data.frame(colData(rse_gene)[,49:58]))
rse_gene$cellPC = cellPca$x[,1]
getPcaVars(cellPca)[1] # 91.8
round(cellPca$rot[,1],3)
# fetal quiescent and adult neuron increase




## filter for expressed
rse_gene = rse_gene[rowMeans(getRPKM(rse_gene,"Length")) > 0.2,]

##############
## metrics ###

## check if ratios of cell changed by batch
pdf(file = "pdfs/Region_Race_cellcheck.pdf")
boxplot(rse_gene$rRNA_rate ~ rse_gene$BrainRegion,xlab="")
boxplot(rse_gene$mitoRate ~ rse_gene$BrainRegion,xlab="")
boxplot(rse_gene$gene_Assigned ~ rse_gene$BrainRegion,xlab="")
boxplot(rse_gene$mitoRate ~ rse_gene$Race,las=3,xlab="")
boxplot(rse_gene$gene_Assigned ~ rse_gene$Race,las=3,xlab="")
dev.off()

#### explore human
geneExprs = log2(getRPKM(rse_gene,"Length")+1)
pca = prcomp(t(geneExprs))
pca_vars = getPcaVars(pca)
pca_vars_lab = paste0("PC", seq(along=pca_vars), ": ",
	pca_vars, "% Var Expl")


##########
pdf("pdfs/PCA_plotsExprs.pdf",w=9)
par(mar=c(8,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4,"Dark2"))

## pc1 vs pc2
plot(pca$x, pch=21, bg=factor(rse_gene$PrimaryDx),cex=1.2,
	xlab = pca_vars_lab[1], ylab = pca_vars_lab[2])
legend("bottomleft", levels(rse_gene$PrimaryDx), col=1:2, pch=15,cex=2)


## by line
for(i in 1:10) {
	boxplot(pca$x[,i] ~ rse_gene$Sex,
		ylab=pca_vars_lab[i], las = 3,xlab="Sex",outline=FALSE)
	points( pca$x[,i] ~ jitter(as.numeric(factor(rse_gene$Sex))),
		pch = 21, bg = rse_gene$PrimaryDx,cex=1.2)
}

## by experiment
for(i in 1:10) {
	boxplot(pca$x[,i] ~ rse_gene$Race,
		ylab=pca_vars_lab[i], las = 3,xlab="Race",outline=FALSE)
	points( pca$x[,i] ~ jitter(as.numeric(factor(rse_gene$Race))),
		pch = 21, bg = rse_gene$PrimaryDx,cex=1.2)
}
dev.off()


####################
## modeling ########
####################

dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)

## mean-variance
mod = model.matrix(~PrimaryDx + cellPC + AgeDeath + BrainRegion,
	data=colData(rse_gene))
pdf(file = "pdfs/vGene.pdf")
vGene = voom(dge,mod,plot=TRUE)
dev.off()

gene_dupCorr = duplicateCorrelation(vGene$E, mod,
	block=colData(rse_gene)$SAMPLE_ID)
save(gene_dupCorr, file = "rdas/gene_dupCorr_neurons.rda")

fitGeneDupl = lmFit(vGene,
	correlation=gene_dupCorr$consensus.correlation,
	block=colData(rse_gene)$SAMPLEID)

ebGeneDupl = eBayes(fitGeneDupl)
outGeneDupl = topTable(ebGeneDupl,coef=2,
	p.value = 1,number=nrow(rse_gene),sort="none")

pdf(file = "pdfs/hist_pval.pdf")
hist(outGeneDupl$P.Value)
dev.off()
table(outGeneDupl$adj.P.Val < 0.05)
table(outGeneDupl$adj.P.Val < 0.1)

sigGeneDupl =  topTable(ebGeneDupl,coef=2,
	p.value = 0.1,number=nrow(rse_gene))

sigGeneDupl[,c("Symbol","logFC", "P.Value","AveExpr")]
sigGeneDupl[sigGeneDupl$logFC > 0,c("Symbol","logFC", "P.Value")]
sigGeneDupl[sigGeneDupl$logFC <  0,c("Symbol","logFC", "P.Value")]

write.csv(outGeneDupl, file = "tables/de_stats_allExprs.csv")
write.csv(sigGeneDupl, file = "tables/de_stats_fdr10_sorted.csv")

###################
## check plots ####
###################

exprs = vGene$E[rownames(sigGeneDupl),]
#exprsClean = cleaningY(exprs, mod, 2)


### make boxplots
# cleanGeneExprs = cleaningY(geneExprs, mod[,!is.na(eBGene$p.value[1,])], P=3)


pdf("pdfs/DE_boxplots_byDiagnosis.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = exprs[i,]
	boxplot(yy ~ rse_gene$PrimaryDx, outline=FALSE,
		ylim=range(yy), ylab="Normalized log2 Exprs", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene$PrimaryDx)),
		pch = 21, bg= rse_gene$PrimaryDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()

#### RPKM #####
e = geneExprs[rownames(sigGeneDupl),]

pdf("pdfs/DE_boxplots_byGenome_log2RPKM.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = e[i,]
	boxplot(yy ~ rse_gene$PrimaryDx, las=3,outline=FALSE,
		ylim=range(yy), ylab="log2(RPKM+1)", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene$PrimaryDx)),
		pch = 21, bg= rse_gene$PrimaryDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()


#################### no rat astrocyte differences

### gene ontology
library(clusterProfiler)
library(org.Hs.eg.db)

## get significant genes by sign
sigGene = outGeneDupl[outGeneDupl$P.Value < 0.005,]
sigGeneList = split(as.character(sigGene$EntrezID), sign(sigGene$logFC))
sigGeneList = lapply(sigGeneList, function(x) x[!is.na(x)])
geneUniverse = as.character(outGeneDupl$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

## do GO and KEGG
goBP_Adj <- compareCluster(sigGeneList, fun = "enrichGO",
	universe = geneUniverse, OrgDb = org.Hs.eg.db,
	ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 1,
	qvalueCutoff  = 1,	readable= TRUE)

goMF_Adj <- compareCluster(sigGeneList, fun = "enrichGO",
	universe = geneUniverse, OrgDb = org.Hs.eg.db,
	ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = 1,
	qvalueCutoff  = 1,	readable= TRUE)

goCC_Adj <- compareCluster(sigGeneList, fun = "enrichGO",
	universe = geneUniverse, OrgDb = org.Hs.eg.db,
	ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = 1,
	qvalueCutoff  = 1,	readable= TRUE)

kegg_Adj <- compareCluster(sigGeneList, fun = "enrichKEGG",
	universe = geneUniverse,  pAdjustMethod = "BH",
	pvalueCutoff  = 1, qvalueCutoff  = 1)

dir.create("rdas")
save(goBP_Adj, goCC_Adj, goMF_Adj, kegg_Adj,
	file = "rdas/gene_set_objects_p005.rda")

goList = list(BP = goBP_Adj, MF = goMF_Adj, CC = goCC_Adj, KEGG = kegg_Adj)
goDf = dplyr::bind_rows(lapply(goList, as.data.frame), .id = "Ontology")
goDf = goDf[order(goDf$pvalue),]

write.csv(goDf, file = "tables/geneSet_output.csv", row.names=FALSE)

options(width=130)
goDf[goDf$p.adjust < 0.05, c(1:5,7)]

################################################
#make heatmap of differentially expressed genes#
################################################
library("pheatmap")

exprs_heatmap = vGene$E[rownames(sigGene),]

df <- as.data.frame(colData(rse_gene)[,c("PrimaryDx")])
rownames(df) <- colnames(exprs_heatmap)
colnames(df)<-"diagnosis"

pdf(file="pdfs/de_heatmap.pdf")
pheatmap(exprs_heatmap, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()
