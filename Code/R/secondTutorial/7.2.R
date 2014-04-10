getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/secondTutorial";
setwd(workingDir);
library(WGCNA);
library(cluster);
options(stringsAsFactors = FALSE);

load("Simulated-RelatingToExt.RData");

ADJ1 = abs(cor(datExpr,use = "p"))^6
Alldegrees1 = intramodularConnectivity(ADJ1,colorh1)

head(Alldegrees1)

colorLevels = unique(colorh1)
sizeGrWindow(9,6)
par(mfrow = c(2,as.integer(0.5 + length(colorLevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorLevels)))
{
	whichmodule = colorLevels[[i]];
	restrict1 = (colorh1 == whichmodule);
	verboseScatterplot(Alldegrees1$kWithin[restrict1],
			GeneSignificance[restrict1], col = colorh1[restrict1],
			main = whichmodule,
			xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

# Overall connectivity of all genes to modules not just intramodular connectiviy
datKME = signedKME(datExpr, datME, outputColumnName = "MM.")

head(datKME)

FilterGenes = abs(GS1) > .2 & abs(datKME$MM.brown) > .8
table(FilterGenes)

dimnames(data.frame(datExpr))[[2]][FilterGenes]

sizeGrWindow(8,6)
par(mfrow = c(2,2))

which.color = "turquoise";
restrictGenes = colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[restrictGenes],
			(datKME[restrictGenes, paste("MM.",which.color,sep = "")])^6,
			col = which.color, xlab = "Intramodular Connectivity",
			ylab = "(Module Membership)^6")

which.color = "blue";
restrictGenes = colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[restrictGenes],
			(datKME[restrictGenes, paste("MM.",which.color,sep = "")])^6,
			col = which.color, xlab = "Intramodular Connectivity",
			ylab = "(Module Membership)^6")

which.color = "brown";
restrictGenes = colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[restrictGenes],
			(datKME[restrictGenes, paste("MM.",which.color,sep = "")])^6,
			col = which.color, xlab = "Intramodular Connectivity",
			ylab = "(Module Membership)^6")

which.color = "green";
restrictGenes = colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[restrictGenes],
			(datKME[restrictGenes, paste("MM.",which.color,sep = "")])^6,
			col = which.color, xlab = "Intramodular Connectivity",
			ylab = "(Module Membership)^6")

NS1=networkScreening(y=y, datME=datME, datExpr=datExpr,
oddPower=3, blockSize=1000, minimumSampleSize=4,
addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)

# network screening analysis
mean(NoiseGeneIndicator[rank(NS1$p.Weighted,ties.method="first")<=100])
# standard analysis based on the correlation p-values (or Student T test)
mean(NoiseGeneIndicator[rank(NS1$p.Standard,ties.method="first")<=100])


topNumbers=c(10,20,50,100)
for (i in c(1:length(topNumbers)) )
{
print(paste("Proportion of noise genes in the top", topNumbers[i], "list"))
WGCNApropNoise=mean(NoiseGeneIndicator[rank(NS1$p.Weighted,ties.method="first")<=topNumbers[i]])
StandardpropNoise=mean(NoiseGeneIndicator[rank(NS1$p.Standard,ties.method="first")<=topNumbers[i]])
print(paste("WGCNA, proportion of noise=", WGCNApropNoise,
", Standard, prop. noise=", StandardpropNoise))
if (WGCNApropNoise< StandardpropNoise) print("WGCNA wins")
if (WGCNApropNoise==StandardpropNoise) print("both methods tie")
if (WGCNApropNoise>StandardpropNoise) print("standard screening wins")
}

rm(dissTOM); collectGarbage();
