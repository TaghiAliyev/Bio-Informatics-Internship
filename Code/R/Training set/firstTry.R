getwd();
wd = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/Training set";
setwd(wd);
library(WGCNA);
library(cluster)
enableWGCNAThreads()
options(stringsAsFactors = FALSE);

# Reading csv file that contains data from training set.
trainingSet = read.csv("extractedData.csv");

dim(trainingSet);
names(trainingSet);
head (trainingSet);

datExpr = as.data.frame(t(trainingSet));
names(datExpr) = trainingSet$ID_REF;
rownames(datExpr) =  names(trainingSet);
datExpr = datExpr[-c(1),];
dim(datExpr)
datExpr0 = mat.or.vec(11,20000)
datExprToFile = mat.or.vec(11,54675)

for (i in 1:11)
{
	for (j in 1:54675)
	{
		datExprToFile[i,j] = as.double(datExpr[i,j]);
}
}

write.table(datExprToFile, file = "datExpr.txt", sep = " ")
datExprToFile = read.table(file = "datExpr.txt", sep = " ")
rownames(datExprToFile) = rownames(datExpr);
colnames(datExprToFile) = trainingSet$ID_REF;

meanExpressionByArray = apply(datExpr0,1,mean, na.rm = T)
NumberMissingByArray = apply(is.na(data.frame(datExpr0)),1,sum)

meanExpressionByArray
NumberMissingByArray

sizeGrWindow(9,5)
barplot(meanExpressionByArray, xlab = "Sample", ylab = "Mean expression",
		main = "Mean expression across samples", names.arg = c(1:11),
		cex.names = 0.7)

sizeGrWindow(9,5)
plotClusterTreeSamples(datExprToFile)

# Getting beta value/power for network construction
powers = c(c(1:20),seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)


sizeGrWindow(9,5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
		xlab = "Soft threshold(power)", ylab = "Scale free topology model fit, signed R^2",
		type = "n", main = paste("Scale Independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
	labels = powers, cex = cex1, col = "red");
abline(h=0.90, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold(power)",
	ylab = "Mean connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

collectGarbage()
Adjacency = abs(cor(datExpr0, use = "p"))^22
k = softConnectivity(datE = datExpr0, power = 22)
sizeGrWindow(10,5)
par(mfrow = c(1,2))
hist(k)
scaleFreePlot(k, main = "Check scale free topology\n")

# Only 4000 most connected genes are kept
restrictedDatExpr = datExprToFile[,rank(-k,ties.method="first") <= 4000]
dim(restrictedDatExpr)


# Using only 4000 of most expressed genes out of 20000
network = blockwiseModules(restrictedDatExpr, power = 22, TOMType = "unsigned",
		minModuleSize = 25, reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE,
	saveTOMFileBase = "trainingSet4000TOM", verbose = 3)

# Block-wise creation that is suited for big datasets. Using whole dataset
blockWiseNetwork = blockwiseModules(datExprToFile, maxBlockSize = 2000, power = 16,
	TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, saveTOMs = TRUE, saveTOMFileBase = "trainingSetBlockWiseTOM",
	verbose = 3)


dissimME = (1-t(cor(datExprToFile,method = "p")))/2
hclustdat = hclust(as.dist(dissimME), method = "average")

par(mfrow = c(1,1))
plot(hclustdat, main = "Clustering tree based of the module eigengenes")

# Now that I have expression profiles, might be a bit more interesting to look into differentially expressed guys
# I have stable plaques and ruptured plaques.
# So, checking if something is differentially expressed between those is probably a good idea

# TODO : check pairwise t-test and see how it actually works. Then for each gene compute it and see if you find some very small p values
# Take the ones where p < 0.05 and then see how many there are and also, take the smallest as well
# See if I can get same results as other guys.

pValues = mat.or.vec(1,54675);
foldChanges = mat.or.vec(1,54675);
for (j in 1:54675)
{
rupturedValues = mat.or.vec(1,5)
stableValues = mat.or.vec(1,6)
for (i in 1:5)
{
	rupturedValues[1,i] = datExprToFile[i,j]
}
for (i in 1:6)
{
stableValues[1,i] = datExprToFile[5+i,j]
}

pValues[1,j] = t.test(rupturedValues,stableValues, paired = TRUE)$p.value
foldChanges[1,j] = mean(rupturedValues)/mean(stableValues);
}
minValue = 1000
minIndex = -1
less = 0

for (i in 1:54675)
{
	if (pValues[1,i] < minValue)
	{
		minValue = pValues[1,i]	
		minIndex = i
	}
	if (pValues[1,i] < 0.05)
	{
		less = less + 1
	}
}
print(minValue)
print(minIndex)
pValues = read.table("pValues.txt", sep = " ")
write.table(pValues,"pValues.txt", sep = " ")
write.table(pValues,"pValuesPaired.txt", sep = " ")

top10 = mat.or.vec(1,10)
top10[1] = 21659 # SULF1
top10[2] = 14829 # leucine rich repeat containing 17. LRRC17
top10[3] = 4402 # phosphoinositide-3-kinase, regulatory subunit 6. PIK3R6
top10[4] = 1256 # 5'-nucleotidase, ecto (CD73)
top10[5] = 20569 # calcium/calmodulin-dependent serine protein kinase (MAGUK family)
top10[6] = 19362 # fibroblast activation protein, alpha
top10[7] = 18240 # cyclin-dependent kinase 16
top10[8] = 39588 # Transcribed locus
top10[9] = 14568 # sarcoglycan, beta (43kDa dystrophin-associated glycoprotein)
top10[10] = 54355 # serine/threonine kinase 10
for (i in 1:10)
{
	print(colnames(datExprToFile)[top10[i]])
}


# TODO : Look into softwares and do something man.
# Hierarchical clustering is important.
less = 0
for (i in 1:54675)
{
	if (adjustedP[i] < 0.05)
	{
		less = less + 1
	}
}
trainingSet = trainingSet[,-c(1)]
limit = 100
copy = trainingSet[1:limit,1:11]
r = c(t(as.matrix(copy)))

f1 = c("CEAr17Norma", "CEAr33Norma", "CEAr35Norma", "CEAr57Norma","CEAr100Norm", "CEAs14Norma", "CEAs25Norma", "CEAs36Norma"
 , "CEAs41Norma", "CEAs45Norma", "CEAs81Norma")
f2 = names(datExprToFile)
k1 = 11
k2 = limit
tm1 = gl(k1,1,k1*k2,factor(f2))
tm2 = gl(k2,k1,k1*k2,factor(f1))

av = aov(r ~ tm1*tm2)
