getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/secondTutorial";
setwd(workingDir);
library(WGCNA);
options(stringsAsFactors = FALSE);
load("Simulated-NetworkConstruction.RData");
detach(ModuleEigengeneNetwork1)
attach(ModuleEigengeneNetwork1)

datME = moduleEigengenes(datExpr,colorh1)$eigengenes
signif(cor(datME,use = "p"),2)

dissimME = (1-cor(datME,method = "p"))/2
hclustdatME = hclust(as.dist(dissimME), method = "average")

par(mfrow = c(1,1))
plot(hclustdatME, main = "Clustering tree based of the module eigengenes")

sizeGrWindow(8,9)
plotMEpairs(datME,  y = y)

signif(cor(datME, ModuleEigengeneNetwork1[,-1]),2)

sizeGrWindow(8,9)
par(mfrow = c(3,1), mar = c(1,2,4,1))
which.module = "turquoise";
plotMat(t(scale(datExpr[,colorh1 == which.module])),nrgcols = 30, rlabels = T,
			clabels = T, rcols = which.module, title = which.module)

which.module = "blue";
plotMat(t(scale(datExpr[,colorh1 == which.module])), nrgcols = 30, rlabels = T,clabels = T,
			rcols = which.module, title = which.module)

which.module = "brown";
plotMat(t(scale(datExpr[,colorh1 == which.module])),nrgcols = 30,rlabels = T, clabels = T,
		rcols = which.module, title = which.module)

sizeGrWindow(8,7);
which.module = "green";
ME = datME[,paste("ME",which.module,sep = "")]
par(mfrow = c(2,1), mar = c(0.3,5.5,3,2))
plotMat(t(scale(datExpr[,colorh1 == which.module])),nrgcols = 30, rlabels = F, rcols = which.module,
			main = which.module, cex.main = 2)
par(mar = c(5,4.2,0,0.7))
barplot(ME,col = which.module, main = "", cex.main = 2, ylab = "eigengene expression", xlab = "array sample")

signif(cor(y,datME, use = "p"),2)

cor.test(y,datME$MEbrown)

p.values = corPvalueStudent(cor(y,datME, use = "p"), nSamples = length(y))

GS1 = as.numeric(cor(y,datExpr,use = "p"))
GeneSignificance = abs(GS1)

ModuleSignificance = tapply(GeneSignificance, colorh1, mean, na.rm = T)
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance, colorh1)
collectGarbage()
save.image("Simulated-RelatingToExt.RData")