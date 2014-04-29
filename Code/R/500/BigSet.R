getwd();
wd = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/500";
setwd(wd);
library(WGCNA);
library(cluster)
enableWGCNAThreads()
options(stringsAsFactors = FALSE);

bigSet = read.csv("mRNA_norm_combined.csv");
write.table(bigSet,"bigSetReady.txt", sep = " ");

library(survival)
args(coxph)

# Check the coxph method. That is the one needed for cox analysis. Just think of the formula

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(bigSet, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

oneStepNetwork = blockwiseModules(bigSet, power = 8, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
		numericLabels = TRUE, saveTOMs = TRUE, saveTOMFileBase = "oneStepTOM", verbose = 3)
network = blockwiseModules(bigSet, maxBlockSize = 2500, power = 8, TOMType = "unsigned",
		minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
		numericLabels = TRUE, saveTOMs = TRUE, saveTOMFileBase = "bigSetTom", verbose = 3)

colors = network$colors
moduleColors = labels2colors(colors)
table(colors)

blockToDraw = 10;

plotDendroAndColors(network$dendrograms[[blockToDraw]], moduleColors[network$blockGenes[[blockToDraw]]], "Module Colors", main = "Gene dendrogram and module colors in block 1",
		dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

smallSample = bigSet[,1:12000]
smaller = bigSet[,1:10]

j = softConnectivity(datE = smallSample, power = 8)

ADJ1=abs(cor(smaller,use="p"))
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(j)
scaleFreePlot(j, main="Check scale free topology\n")