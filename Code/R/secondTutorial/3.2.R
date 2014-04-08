getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/secondTutorial";
setwd(workingDir);
library(WGCNA);
options(stringsAsFactors = FALSE);

load("Simulated-dataSimulation.RData");
attach(ModuleEigengeneNetwork1)

meanExpressionByArray = apply(datExpr,1,mean,na.rm = T)
NumberMissingByArray = apply(is.na(data.frame(datExpr)),1,sum)

sizeGrWindow(9,5)
barplot(meanExpressionByArray, xlab = "Sample", ylab = "Mean Expression",
			main = "Mean expression across samples",
			names.arg = c(1:50), cex.names = 0.7)

NumberMissingByArray

KeepArray = NumberMissingByArray < 500
table(KeepArray)
datExpr = datExpr[KeepArray,]
y = y[KeepArray]
ArrayName[KeepArray]

NumberMissingByGene = apply(is.na(data.frame(datExpr)),2,sum)

summary(NumberMissingByGene)

variancedatExpr = as.vector(apply(as.matrix(datExpr),2,var,na.rm=T))
no.presentdatExpr = as.vector(apply(!is.na(as.matrix(datExpr)),2,sum))

table(no.presentdatExpr)

KeepGenes = variancedatExpr > 0 & no.presentdatExpr >= 4
table(KeepGenes)

datExpr = datExpr[,KeepGenes]
GeneName = GeneName[KeepGenes]

sizeGrWindow(9,5)
plotClusterTreeSamples(datExpr = datExpr, y = y)