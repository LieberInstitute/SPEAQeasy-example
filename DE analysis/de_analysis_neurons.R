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
load("count_data/stemcell_pipeline_rse_gene_annotated_n106.Rdata")

## keep subset of samples
neuronIndex = which(rse_gene$DIVgroup %in% c("56","70") &
	rse_gene$SensitivitySubset & !grepl("SCP$", rse_gene$SAMPLE_ID))
rse_gene = rse_gene[,neuronIndex]

## cell PCs
cellPca = prcomp(as.data.frame(colData(rse_gene)[,49:58]))
rse_gene$cellPC = cellPca$x[,1]
getPcaVars(cellPca)[1] # 84%
round(cellPca$rot[,1],3)
# fetal quiescent and adult neuron increase
# NPC and endothelial decrease


### split by human and rat
rse_gene_hs = rse_gene[rowData(rse_gene)$gene_species == "human",]
rse_gene_rn = rse_gene[rowData(rse_gene)$gene_species == "rat",]

## filter for expressed
rse_gene_hs = rse_gene_hs[rowMeans(getRPKM(rse_gene_hs,"Length")) > 0.2,]
rse_gene_rn = rse_gene_rn[rowMeans(getRPKM(rse_gene_rn,"Length")) > 0.2,]

##############
## metrics ###

## check if ratios of cell changed by batch
boxplot(rse_gene$rn6_fraction ~ rse_gene$Batch,xlab="")
boxplot(rse_gene$r_geneRate ~ rse_gene$Batch,xlab="")
boxplot(rse_gene$h_geneRate ~ rse_gene$Batch,xlab="")
boxplot(rse_gene$rn6_fraction ~ rse_gene$Exp,las=3,xlab="")
boxplot(rse_gene$r_geneRate ~ rse_gene$Exp,las=3,xlab="")
boxplot(rse_gene$hg38_fraction ~ rse_gene$Exp,las=3,xlab="")
boxplot(rse_gene$h_geneRate ~ rse_gene$Exp,las=3,xlab="")


#### explore human
geneExprs_hs = log2(getRPKM(rse_gene_hs,"Length")+1)
pca_hs = prcomp(t(geneExprs_hs))
pca_hs_vars = getPcaVars(pca_hs)
pca_hs_vars_lab = paste0("PC", seq(along=pca_hs_vars), ": ",
	pca_hs_vars, "% Var Expl")


##########
pdf("pdfs/PCA_plots_hsExprs.pdf",w=9)
par(mar=c(8,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4,"Dark2"))

## pc1 vs pc2
plot(pca_hs$x, pch=21, bg=factor(rse_gene_hs$RealDx),cex=1.2,
	xlab = pca_hs_vars_lab[1], ylab = pca_hs_vars_lab[2])
legend("bottomleft", levels(rse_gene_hs$RealDx), col=1:2, pch=15,cex=2)

# pc1 is neuronal fraction
plot(pca_hs$x[,1] ~ rse_gene_hs$Neurons,
	ylab = pca_hs_vars_lab[1],  xlab="Adult Neuron RNA Fraction",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))
summary(lm(pca_hs$x[,1] ~ rse_gene_hs$Neurons))

plot(pca_hs$x[,1] ~ rse_gene_hs$NPC,
	ylab = pca_hs_vars_lab[1],  xlab="NPC RNA Fraction",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))

plot(pca_hs$x[,1] ~ rse_gene_hs$Fetal_quiescent,
	ylab = pca_hs_vars_lab[1],  xlab="Fetal Quiescent RNA Fraction",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))

plot(pca_hs$x[,1] ~ rse_gene_hs$Endothelial,
	ylab = pca_hs_vars_lab[1],  xlab="Endothelial RNA Fraction",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))

plot(pca_hs$x[,1] ~ rse_gene_hs$cellPC,
	ylab = pca_hs_vars_lab[1],  xlab="PC1 of RNA Fractions",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))

## PC2 - human fraction?
plot(pca_hs$x[,2] ~ rse_gene_hs$hg38_fraction,
	ylab = pca_hs_vars_lab[2],  xlab="Human Alignment Fraction",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))

## PC3?
boxplot(pca_hs$x[,3] ~ rse_gene_hs$Batch,
	ylab = pca_hs_vars_lab[3],  xlab="Color Batch",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))
points( pca_hs$x[,3] ~ jitter(as.numeric(factor(rse_gene_hs$Batch))),
		pch = 21, bg = rse_gene$RealDx,cex=1.2)

## pc4 is gene assignment
plot(pca_hs$x[,4] ~ rse_gene_hs$hg38_fraction,
	ylab = pca_hs_vars_lab[4],  xlab="Human Alignment Rate",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))
plot(pca_hs$x[,4] ~ rse_gene_hs$h_geneRate,
	ylab = pca_hs_vars_lab[4],  xlab="Human Gene Assignment Rate",
	pch = 21,cex=1.2,bg=factor(rse_gene_hs$RealDx))

## by line
for(i in 1:10) {
	boxplot(pca_hs$x[,i] ~ rse_gene_hs$RealCode,
		ylab=pca_hs_vars_lab[i], las = 3,xlab="",outline=FALSE)
	points( pca_hs$x[,i] ~ jitter(as.numeric(factor(rse_gene_hs$RealCode))),
		pch = 21, bg = rse_gene$RealDx,cex=1.2)
}

## by experiment
for(i in 1:10) {
	boxplot(pca_hs$x[,i] ~ rse_gene_hs$Exp,
		ylab=pca_hs_vars_lab[i], las = 3,xlab="",outline=FALSE)
	points( pca_hs$x[,i] ~ jitter(as.numeric(factor(rse_gene_hs$Exp))),
		pch = 21, bg = rse_gene$RealDx,cex=1.2)
}
dev.off()

#########################3
#### explore rat ########
geneExprs_rn = log2(getRPKM(rse_gene_rn,"Length")+1)
pca_rn = prcomp(t(geneExprs_rn))
pca_rn_vars = getPcaVars(pca_rn)
pca_rn_vars_lab = paste0("PC", seq(along=pca_rn_vars), ": ",
	pca_rn_vars, "% Var Expl")

## 10 most expressed sequences
exInd = order(rowMeans(geneExprs_rn),decreasing=TRUE)[1:10]
colSums(assays(rse_gene_rn)$counts[exInd,])/1e6
colSums(assays(rse_gene_rn)$counts[-exInd,])/1e6

#######
pdf("pdfs/PCA_plots_rnExprs.pdf",w=9)
par(mar=c(8,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4,"Dark2"))

## pc1 vs pc2
plot(pca_rn$x, pch=21, bg=factor(rse_gene_rn$RealDx),cex=1.2,
	xlab = pca_rn_vars_lab[1], ylab = pca_rn_vars_lab[2])
legend("bottomleft", levels(rse_gene_rn$RealDx), col=1:3, pch=15,cex=1.2)

# pc1 is fraction rat
plot(pca_rn$x[,1] ~ rse_gene_rn$rn6_fraction,
	ylab = pca_rn_vars_lab[1],  xlab="Rat Alignment Rate",
	pch = 21,cex=1.2,bg=factor(rse_gene_rn$RealDx))

# pc2 is batch - do we know astrocyte batches?
boxplot(pca_rn$x[,2] ~ rse_gene_rn$Batch, outline=FALSE,
	ylab = pca_rn_vars_lab[2],  xlab="Rat Alignment Rate")
points(pca_rn$x[,2] ~ jitter(as.numeric(factor(rse_gene_rn$Batch)),amount=0.15),
	pch = 21,cex=1.2,bg=factor(rse_gene_rn$RealDx))

boxplot(pca_rn$x[,2] ~ rse_gene_rn$Exp, outline=FALSE,las=3,
	ylab = pca_rn_vars_lab[2],  xlab="")
points(pca_rn$x[,2] ~ jitter(as.numeric(rse_gene_rn$Exp),amount=0.15),
	pch = 21,cex=1.2,bg=factor(rse_gene_rn$RealDx))

## PC3 - cell type?
plot(pca_rn$x[,3] ~ rse_gene_rn$NPC,
	ylab = pca_rn_vars_lab[3],  xlab="NPC RNA Fraction",
	pch = 21,cex=1.2,bg=factor(rse_gene_rn$RealDx))

plot(pca_rn$x[,3] ~ rse_gene_rn$Neurons,
	ylab = pca_rn_vars_lab[3],  xlab="Adult Neuron RNA Fraction",
	pch = 21,cex=1.2,bg=factor(rse_gene_rn$RealDx))

## by line
for(i in 1:10) {
	boxplot(pca_rn$x[,i] ~ rse_gene_rn$RealCode,
		ylab=pca_rn_vars_lab[i], las = 3,xlab="",outline=FALSE)
	points( pca_rn$x[,i] ~ jitter(as.numeric(factor(rse_gene_rn$RealCode))),
		pch = 21, bg = rse_gene$RealDx,cex=1.2)
}

## by experiment
for(i in 1:10) {
	boxplot(pca_rn$x[,i] ~ rse_gene_rn$Exp,
		ylab=pca_rn_vars_lab[i], las = 3,xlab="",outline=FALSE)
	points( pca_rn$x[,i] ~ jitter(as.numeric(factor(rse_gene_rn$Exp))),
		pch = 21, bg = rse_gene$RealDx,cex=1.2)
}
dev.off()

####################
## modeling ########
####################

dge = DGEList(counts = assays(rse_gene_hs)$counts,
	genes = rowData(rse_gene_hs))
dge = calcNormFactors(dge)

## mean-variance
mod = model.matrix(~RealDx + cellPC + hg38_fraction + h_geneRate + DIVgroup,
	data=colData(rse_gene_hs))
vGene = voom(dge,mod,plot=TRUE)

gene_dupCorr = duplicateCorrelation(vGene$E, mod,
	block=colData(rse_gene)$RealCode)
save(gene_dupCorr, file = "rdas/gene_dupCorr_neurons.rda")

fitGeneDupl = lmFit(vGene,
	correlation=gene_dupCorr$consensus.correlation,
	block=colData(rse_gene)$RealCode)

ebGeneDupl = eBayes(fitGeneDupl)
outGeneDupl = topTable(ebGeneDupl,coef=2,
	p.value = 1,number=nrow(rse_gene_hs),sort="none")

hist(outGeneDupl$P.Value)
table(outGeneDupl$adj.P.Val < 0.05)
table(outGeneDupl$adj.P.Val < 0.1)

sigGeneDupl =  topTable(ebGeneDupl,coef=2,
	p.value = 0.1,number=nrow(rse_gene_hs))

sigGeneDupl[,c("Symbol","logFC", "P.Value","AveExpr")]
sigGeneDupl[sigGeneDupl$logFC > 0,c("Symbol","logFC", "P.Value")]
sigGeneDupl[sigGeneDupl$logFC <  0,c("Symbol","logFC", "P.Value")]

write.csv(outGeneDupl, file = "tables/de_stats_allExprs.csv")
write.csv(sigGeneDupl, file = "tables/de_stats_fdr10_sorted.csv")

###################
## check plots ####
###################

exprs = vGene$E[rownames(sigGeneDupl),]
exprsClean = cleaningY(exprs, mod, 2)


### make boxplots
# cleanGeneExprs = cleaningY(geneExprs_hs, mod[,!is.na(eBGene$p.value[1,])], P=3)
pdf("pdfs/DE_boxplots_byGenome.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = exprs[i,]
	boxplot(yy ~ rse_gene_hs$RealCode, las=3,outline=FALSE,
		ylim=range(yy), ylab="Normalized log2 Exprs", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene_hs$RealCode)),
		pch = 21, bg= rse_gene_hs$RealDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()

pdf("pdfs/DE_boxplots_byGenome_clean.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = exprsClean[i,]
	boxplot(yy ~ rse_gene_hs$RealCode, las=3,outline=FALSE,
		ylim=range(yy), ylab="Cleaned log2 Exprs", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene_hs$RealCode)),
		pch = 21, bg= rse_gene_hs$RealDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()

pdf("pdfs/DE_boxplots_byDiagnosis.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = exprs[i,]
	boxplot(yy ~ rse_gene_hs$RealDx, outline=FALSE,
		ylim=range(yy), ylab="Normalized log2 Exprs", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene_hs$RealDx)),
		pch = 21, bg= rse_gene_hs$RealDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()

pdf("pdfs/DE_boxplots_byDiagnosis_clean.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = exprsClean[i,]
	boxplot(yy ~ rse_gene_hs$RealDx, outline=FALSE,
		ylim=range(yy), ylab="Cleaned log2 Exprs", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene_hs$RealDx)),
		pch = 21, bg= rse_gene_hs$RealDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()

#### RPKM #####
e = geneExprs_hs[rownames(sigGeneDupl),]
ec = cleaningY(e, mod, P=2)

pdf("pdfs/DE_boxplots_byGenome_log2RPKM.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = e[i,]
	boxplot(yy ~ rse_gene_hs$RealCode, las=3,outline=FALSE,
		ylim=range(yy), ylab="log2(RPKM+1)", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene_hs$RealCode)),
		pch = 21, bg= rse_gene_hs$RealDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()

pdf("pdfs/DE_boxplots_byGenome_clean_log2RPKM.pdf",w=10)
par(mar=c(8,6,4,2),cex.axis=1.8,cex.lab=1.8, cex.main=1.8)
palette(brewer.pal(4,"Dark2"))
for(i in 1:nrow(sigGeneDupl)) {
	yy = ec[i,]
	boxplot(yy ~ rse_gene_hs$RealCode, las=3,outline=FALSE,
		ylim=range(yy), ylab="Cleaned log2(RPKM+1)", xlab="",
		main = paste(sigGeneDupl$Symbol[i], "-", sigGeneDupl$gencodeID[i]))
	points(yy ~ jitter(as.numeric(rse_gene_hs$RealCode)),
		pch = 21, bg= rse_gene_hs$RealDx,cex=1.3)
	ll = ifelse(sigGeneDupl$logFC[i] > 0, "topleft", "topright")
	legend(ll, paste0("p=", signif(sigGeneDupl$P.Value[i],3)), cex=1.3)
}
dev.off()


##################
### RAT
dge_rn = DGEList(counts = assays(rse_gene_rn)$counts,
	genes = rowData(rse_gene_rn))
dge_rn = calcNormFactors(dge_rn)

## mean-variance
vGene_rn = voom(dge_rn,mod,plot=TRUE)

gene_dupCorr_rat = duplicateCorrelation(vGene_rn$E, mod,
	block=colData(rse_gene)$RealCode)
save(gene_dupCorr,gene_dupCorr_rat, file = "rdas/gene_dupCorr_neurons.rda")


fitGeneDupl_rat = lmFit(vGene_rn,
	correlation=gene_dupCorr_rat$consensus.correlation,
	block=colData(rse_gene)$RealCode)

ebGeneDupl_rat = eBayes(fitGeneDupl_rat)
outGeneDupl_rat = topTable(ebGeneDupl_rat,coef=2,
	p.value = 1,number=nrow(rse_gene_hs),sort="none")

min(outGeneDupl_rat$adj.P.Val)

sum(outGeneDupl_rat$adj.P.Val < 0.05)
sum(outGeneDupl_rat$adj.P.Val < 0.2)
colSums(apply(ebGeneDupl_rat$p.value,2,p.adjust,"fdr") < 0.05)
## no rat astrocyte differences

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

#############################
## compare to postmortem ####
#############################

## read back in stats
outGeneDupl = read.csv("tables/de_stats_allExprs.csv", row.names=1, as.is=TRUE)

## TWAS
load("/dcl01/ajaffe/data/lab/dg_hippo_paper/rdas/tt_objects_gene.Rdata",verbose=TRUE)
tt = as.data.frame(tt)

## split by region
ttList = split(tt, tt$region)
names(ttList)[1] = "DG-GCL"

## put in same order
g = unique(tt$geneid)
ttList = lapply(ttList, function(x) {
	x = as.data.frame(x[match(g, x$geneid),])
	rownames(x) = g
	x})

## get out matrices
tMat = sapply(ttList, "[[", "TWAS.Z")
pMat = sapply(ttList, "[[", "TWAS.P")
pMat[is.na(pMat)] = 1
fdrMat = sapply(ttList, "[[", "TWAS.FDR")
fdrMat[is.na(fdrMat)] = 1
bonfMat = sapply(ttList, "[[", "TWAS.Bonf")
bonfMat[is.na(bonfMat)] = 1
rownames(tMat) = rownames(pMat) = rownames(fdrMat) = rownames(bonfMat) = g

## match up
outGeneDupl_tt = outGeneDupl[match(rownames(tMat), rownames(outGeneDupl)),]
outGeneDupl_tt$TWAS_Z_DLPFC = tMat[,"DLPFC"]
outGeneDupl_tt$TWAS_P_DLPFC = pMat[,"DLPFC"]
outGeneDupl_tt$TWAS_FDR_DLPFC = fdrMat[,"DLPFC"]

## plot DLPFC
pdf("pdfs/TWAS_vs_StemDE.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(t ~ TWAS_Z_DLPFC, data = outGeneDupl_tt,pch=21,bg="grey",
	xlab = "TWAS SCZD Z-score (DLPFC)",
	ylab = "SCZD IPSC Neuron Z-score")

sigIndex = which(outGeneDupl_tt$adj.P.Val < 0.1)
plot(t ~ TWAS_Z_DLPFC, data = outGeneDupl_tt[sigIndex,],
	pch=21,bg="grey",ylim = range(outGeneDupl_tt$t,na.rm=TRUE),
	xlim = range(outGeneDupl_tt$TWAS_Z_DLPFC,na.rm=TRUE),
	xlab = "TWAS SCZD Z-score (DLPFC)",
	ylab = "SCZD IPSC Neuron Z-score")
abline(h=0,v=0,lty=2)
dev.off()

cor(tMat[,"DLPFC"], outGeneDupl_tt$t ,use="comp")

cor(tMat[sigIndex,"DLPFC"], outGeneDupl_tt$t[sigIndex] ,use="comp")
which(tMat[,"DLPFC"] > 7)

# PCDHA5
outGeneDupl_tt[which(outGeneDupl_tt$t > 7),]

## CNTN4
outGeneDupl_tt[which(outGeneDupl_tt$t > 4 & outGeneDupl_tt$TWAS_Z_DLPFC > 4),]
outGeneDupl_tt[which(outGeneDupl_tt$t > 4 & outGeneDupl_tt$TWAS_Z_DLPFC > 3),]
outGeneDupl_tt[which(outGeneDupl_tt$t < -3 & outGeneDupl_tt$TWAS_Z_DLPFC < -2),]

outGeneDupl_tt[which(outGeneDupl_tt$t < -3 & tMat[,"DLPFC"] < -2),]

#### #meta??
outGeneDupl_tt$zMeta_it = (outGeneDupl_tt$t + outGeneDupl_tt$TWAS_Z_DLPFC)/sqrt(2)
outGeneDupl_tt$pMeta_it = 2*pnorm(-abs(outGeneDupl_tt$zMeta_it))
bonfMeta_it = p.adjust(outGeneDupl_tt$pMeta_it,"bonf")
table(bonfMeta_it < 0.05)
sigGene_it = outGeneDupl_tt[order(pMeta_it)[1:sum(bonfMeta_it < 0.05,na.rm=TRUE)],]
write.csv(sigGene_it, file ="tables/metaAnalysis_ipscAndTwas_bonf05.csv")

## hippo
plot(tMat[,"HIPPO"], outGeneDupl_tt$t )
cor(tMat[,"HIPPO"], outGeneDupl_tt$t ,use="comp")
sigIndex = which(outGeneDupl_tt$adj.P.Val < 0.1)
plot(tMat[sigIndex,"HIPPO"], outGeneDupl_tt$t[sigIndex])
abline(h=0,v=0,lty=2)
cor(tMat[sigIndex,"HIPPO"], outGeneDupl_tt$t[sigIndex] ,use="comp")

###################
##### load DE #####
###################
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda",verbose=TRUE)
outGeneD = outGene
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda",verbose=TRUE)
outGeneH = outGene

outGeneDupl_de = outGeneDupl[match(rownames(outGeneD), rownames(outGeneDupl)),]
outGeneDupl_de$DE_t_DLPFC = outGeneD$t
outGeneDupl_de$DE_p_DLPFC = outGeneD$P.Value


pdf("pdfs/DlpfcDE_vs_StemDE.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(t ~ DE_t_DLPFC, data = outGeneDupl_de,
	pch=21,bg="grey",
	xlab = "SCZD DE Z-score (DLPFC)",
	ylab = "SCZD DE Z-score (Neuron)")

cor(outGeneDupl_de$DE_t_DLPFC, outGeneDupl_de$t ,use="comp")
sigIndex = which(outGeneDupl_de$adj.P.Val < 0.1)
plot(t ~ DE_t_DLPFC, data = outGeneDupl_de[sigIndex,],
	pch=21,bg="grey",ylim = range(outGeneDupl_de$t,na.rm=TRUE),
	xlim = range(outGeneDupl_de$DE_t_DLPFC,na.rm=TRUE),
	xlab = "SCZD DE Z-score (DLPFC)",
	ylab = "SCZD DE Z-score (Neuron)")
abline(h=0,v=0,lty=2)
dev.off()

outGeneDupl_de[which(outGeneDupl_de$t > 3 & outGeneDupl_de$DE_t_DLPFC > 3),]

#### #meta??
outGeneDupl_de$zMeta_id = (outGeneDupl_de$t + outGeneDupl_de$DE_t_DLPFC)/sqrt(2)
outGeneDupl_de$pMeta_id = 2*pnorm(-abs(outGeneDupl_de$zMeta_id))
bonfMeta_id = p.adjust(outGeneDupl_de$pMeta_id,"bonf")
table(bonfMeta_id < 0.05)
sigGene_id = outGeneDupl_de[order(pMeta_id)[1:sum(bonfMeta_id < 0.05,na.rm=TRUE)],]
write.csv(sigGene_id, file ="tables/metaAnalysis_ipscAndDE_bonf05.csv")

####################
## all three meta
outGeneDupl_de$TWAS_Z_DLPFC = tMat[match(rownames(outGeneDupl_de), rownames(tMat)),"DLPFC"]
outGeneDupl_de$TWAS_P_DLPFC = pMat[match(rownames(outGeneDupl_de), rownames(pMat)),"DLPFC"]
outGeneDupl_de$zMeta_idt = (outGeneDupl_de$t + outGeneDupl_de$DE_t_DLPFC +outGeneDupl_de$TWAS_Z_DLPFC )/sqrt(3)
outGeneDupl_de$pMeta_idt = 2*pnorm(-abs(outGeneDupl_de$zMeta_idt))
bonfMeta_idt = p.adjust(pMeta_idt,"bonf")
table(bonfMeta_idt < 0.05)
sigGene_idt = outGeneDupl_de[order(pMeta_idt)[1:sum(bonfMeta_idt < 0.05,na.rm=TRUE)],]

write.csv(sigGene_idt, file ="tables/metaAnalysis_ipscAndDEandTwas_bonf05.csv")

######################
# vs npcs
load("/dcl01/lieber/ajaffe/Brady/ipsc/hiler_NPC_SCZD_timecourse/analysis_MTN/rdas/DEstats_n83-Sep2019-iteration_MTN23Sep2019.Rdata")

mm = match(rownames(outGeneDupl), rownames(outGene.anja.add))
outGeneDupl$t_NPC = outGene.anja.add$t[mm]
outGeneDupl$logFC_NPC = outGene.anja.add$logFC[mm]
outGeneDupl$P.Value_NPC = outGene.anja.add$P.Value[mm]
outGeneDupl$adj.P.Val_NPC = outGene.anja.add$adj.P.Val[mm]

pdf("pdfs/NPCvsNeuron_DE.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(outGeneDupl$t_NPC, outGeneDupl$t,
	pch = 21, bg="grey", xlab = "SCZD iPSC NPC Z-score",
	ylab = "SCZD iPSC Neuron Z-score")
cc = cor(outGeneDupl$t_NPC, outGeneDupl$t,use="comp")
legend("topleft", paste0("r=",signif(cc,3)),cex=1.5)

sigIndexNeuron = which(outGeneDupl$adj.P.Val < 0.1)
plot(outGeneDupl$t_NPC[sigIndexNeuron], outGeneDupl$t[sigIndexNeuron],
	xlab = "SCZD iPSC NPC Z-score",pch = 21, bg="grey",
	ylab = "SCZD iPSC Neuron Z-score",
	ylim = range(outGeneDupl$t,na.rm=TRUE),
	xlim = range(outGeneDupl$t_NPC,na.rm=TRUE))

sigIndexNpc = which(outGeneDupl$adj.P.Val_NPC < 0.1)
plot(outGeneDupl$t_NPC[sigIndexNpc], outGeneDupl$t[sigIndexNpc],
	xlab = "SCZD iPSC NPC Z-score",pch = 21, bg="grey",
	ylab = "SCZD iPSC Neuron Z-score",
	ylim = range(outGeneDupl$t,na.rm=TRUE),
	xlim = range(outGeneDupl$t_NPC,na.rm=TRUE))
dev.off()

########################
## vs brennand
prev =  read.csv("genelist_JAMMA.csv",as.is=TRUE)
outGeneDupl_prev = outGeneDupl[match(prev$Ensembl, outGeneDupl$ensemblID),]
outGeneDupl_prev$logFC_prev = prev$logFC
outGeneDupl_prev$t_prev = prev$t

pdf("pdfs/usVsBrennand_Neuron_DE.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(outGeneDupl_prev$t_prev, outGeneDupl_prev$t,
	pch = 21, bg="grey", xlab = "Brennand SCZD Neuron Z-score",
	ylab = "LIBD SCZD Neuron Z-score")
cc = cor(outGeneDupl_prev$t_prev, outGeneDupl_prev$t,use="comp")
legend("topleft", paste0("r=",signif(cc,3)),cex=1.5)

sigIndexNeuron = which(outGeneDupl_prev$adj.P.Val < 0.1)
plot(outGeneDupl_prev$t_prev[sigIndexNeuron], outGeneDupl_prev$t[sigIndexNeuron],
	xlab = "Brennand SCZD Neuron Z-score",pch = 21, bg="grey",
	ylab = "LIBD SCZD Neuron Z-score",
	ylim = range(outGeneDupl_prev$t,na.rm=TRUE),
	xlim = range(outGeneDupl_prev$t_prev,na.rm=TRUE))
	abline(h=0,v=0,lty=2)

dev.off()

outGeneDupl$prevSig = outGeneDupl$ensemblID %in% prev$Ensembl
tt = table(outGeneDupl$prevSig, outGeneDupl$adj.P.Val < 0.1)
getOR(tt)
chisq.test(tt)

######################
## ron
load("/dcl01/ajaffe/data/lab/libd_stem_timecourse/tc_analysis/scz_vs_cnt/de_stats_szControl.rda")
outGeneDupl_old = outGeneDupl[match(rownames(voomStats), rownames(outGeneDupl)),]

outGeneDupl_old$t_old = voomStats$t_NEURON

pdf("pdfs/usVsOldUs_Neuron_DE.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(outGeneDupl_old$t_old, outGeneDupl_old$t,
	pch = 21, bg="grey", xlab = "Old LIBD SCZD Neuron Z-score",
	ylab = "LIBD SCZD Neuron Z-score")
cc = cor(outGeneDupl_old$t_old, outGeneDupl_old$t,use="comp")
legend("topleft", paste0("r=",signif(cc,3)),cex=1.5)

sigIndexNeuron = which(outGeneDupl_old$adj.P.Val < 0.1)
plot(outGeneDupl_old$t_old[sigIndexNeuron], outGeneDupl_old$t[sigIndexNeuron],
	xlab = "Old LIBD SCZD Neuron Z-score",pch = 21, bg="grey",
	ylab = "LIBD SCZD Neuron Z-score",
	ylim = range(outGeneDupl_old$t,na.rm=TRUE),
	xlim = range(outGeneDupl_old$t_old,na.rm=TRUE))
abline(h=0,v=0,lty=2)

dev.off()

######################
## other brennand ####

library(readxl)

prev2 = read_excel("41467_2017_2330_MOESM11_ESM.xlsx", sheet = "Neuron")
prev2 = as.data.frame(prev2)

outGeneDupl_prev2 = outGeneDupl[match(prev2$geneName, outGeneDupl$ensemblID),]
outGeneDupl_prev2$t_prev = prev2$t

pdf("pdfs/usVsBrennandNatComm_Neuron_DE.pdf")
par(mar=c(5,6,2,2),cex.axis=2,cex.lab=2)
plot(outGeneDupl_prev2$t_prev, outGeneDupl_prev2$t,
	pch = 21, bg="grey", xlab = "Brennand NatComm SCZD Neuron Z-score",
	ylab = "LIBD SCZD Neuron Z-score")
cc = cor(outGeneDupl_prev2$t_prev, outGeneDupl_prev2$t,use="comp")
legend("topleft", paste0("r=",signif(cc,3)),cex=1.5)

sigIndexNeuron = which(outGeneDupl_prev2$adj.P.Val < 0.1)
plot(outGeneDupl_prev2$t_prev[sigIndexNeuron], outGeneDupl_prev2$t[sigIndexNeuron],
	xlab = "Brennand NatComm SCZD Neuron Z-score",pch = 21, bg="grey",
	ylab = "LIBD SCZD Neuron Z-score",
	ylim = range(outGeneDupl_prev2$t,na.rm=TRUE),
	xlim = range(outGeneDupl_prev2$t_prev,na.rm=TRUE))
abline(h=0,v=0,lty=2)
dev.off()

cor(outGeneDupl_prev2$t_prev[sigIndexNeuron],
	outGeneDupl_prev2$t[sigIndexNeuron],use="comp")

####################
## ion checks ###

library(readxl)
ions =read_excel("Ion Channel Gene List.xlsx",col_names=FALSE)
colnames(ions) = "Symbol"
ions$GeneSymbol = toupper(ions$Symbol)
ions$GeneSymbol[ions$Symbol == "Scn2a1"] = "SCN2A"
table(ions$GeneSymbol[1:89] %in%  rowData(rse_gene)$Symbol)
table(outGeneDupl$Symbol %in% ions$GeneSymbol[1:89])
outGeneDupl$isIon = outGeneDupl$Symbol %in% ions$GeneSymbol[1:89]

table(outGeneDupl$isIon, outGeneDupl$adj.P.Val < 0.1)
table(outGeneDupl$isIon, outGeneDupl$P.Value < 0.05)
cat(outGeneDupl[which(outGeneDupl$isIon & outGeneDupl$P.Value < 0.05),"Symbol"], sep= "\n")
outGeneDupl[which(outGeneDupl$isIon & outGeneDupl$P.Value < 0.1),]
