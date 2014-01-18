getwd();
workingDir = "C:/Users/Taghi/Desktop/Internship";
setwd(workingDir);
library(WGCNA)
options(stringsAsFactors = FALSE);

# Enables multi threading
enableWGCNAThreads()
# Loading of file
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames

# Generating Weighted network using 1-step network construction
powers = c(c(1:10),seq(from = 12,to = 20, by = 2))
sft = pickSoftThreshold(datExpr,powerVector = powers , verbose = 5)

sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab = "Soft Threshold (power)", 
	ylab = "Scale Free Topology Model Fit, signed R^2",
	type = "n", main = paste("Scale Independence"));

Taghi = 10
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1],sft$fitIndices[,5], xlab = "Soft Threshold(power)",ylab="Mean Connectivity",
	type = "n", main=paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],powers, cex1, col="red")

#Constructing actual network

network = blockwiseModules(datExpr,power = 6, TOMType = "unsigned", minModuleSize = 30,
		reassignThreshold = 0, mergeCutHeight = 0.25, 
		numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, 
		saveTOMFileBase="femaleMouseTOM",verbose = 3)


sizeGrWindow(12, 9)
mergedColors = labels2colors(network$colors)
plotDendroAndColors(network$dendrograms[[1]],mergedColors[network$blockGenes[[1]]],
		"Module colors",
		dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

moduleLabels = network$colors
moduleColors = labels2colors(network$colors)
MEs = network$MEs;
geneTree = network$dendrograms[[1]];
save(MEs,moduleLabels,moduleColors,geneTree,file="FemaleLive-02-networkConstruction-auto.RData")
