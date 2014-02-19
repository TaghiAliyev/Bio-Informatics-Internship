getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R";
setwd(workingDir);
library(WGCNA)
options(stringsAsFactors = FALSE);

# Enables multi threading
enableWGCNAThreads()
# Loading of file
lnames = load(file = "FemaleLiver-01-dataInput.RData");

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],xlab = "Soft Threshold Power", 
	ylab = "Scale Free Topology Model Fit, signed R^2", type="n",main=paste("Scale Independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],labels = powers, cex = cex1, col="red")

abline(h=0.90, col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold(power)"
	,ylab="Mean Connectivity",type="n",main = paste("Mean Connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

bwnet = blockwiseModules(datExpr, maxBlockSize = 2000, power = 6, TOMType = "unsigned",
		minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
		saveTOMs = TRUE, saveTOMFileBase = "femaleMouseTOM-blockwise", verbose = 3)
getwd()
load(file = "FemaleLiver-02-networkConstruction-auto.RData")

bwLabels = matchLabels(bwnet$colors, moduleLabels)

bwModuleColors = labels2colors(bwLabels)

table(bwLabels)

sizeGrWindow(6,6)

plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
		"Module Colors", main = "Gene dendrogram and module colors in block 1",
		dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
			"Module Colors", main = "Gene dendrogram and module colors in block 2",
			dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

