getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/secondTutorial";
setwd(workingDir);
library(WGCNA);
options(stringsAsFactors = FALSE);

load("Simulated-StandardScreening.RData");
attach(ModuleEigengeneNetwork1)

ADJ1 = abs(cor(datExpr,use="p"))^6

k = as.vector(apply(ADJ1,2,sum,na.rm=T))

sizeGrWindow(10,5)
par(mfrow = c(1,2))
hist(k)
scaleFreePlot(k,main = "Check scale free topology\n")

datExpr = datExpr[, rank(-k,ties.method = "first")<=3600]

dissADJ = 1 - ADJ1

dissTOM = TOMdist(ADJ1)
collectGarbage();

pam4 = pam(as.dist(dissADJ), 4)
pam5 = pam(as.dist(dissADJ), 5)
pam6 = pam(as.dist(dissADJ), 6)

table(pam4$clustering,truemodule)
table(pam5$clustering,truemodule)
table(pam6$clustering,truemodule)


pamTOM4 = pam(as.dist(dissTOM),4)
pamTOM5 = pam(as.dist(dissTOM),5)
pamTOM6 = pam(as.dist(dissTOM),6)

table(pamTOM4$clustering,truemodule)
table(pamTOM5$clustering,truemodule)
table(pamTOM6$clustering,truemodule)

hierADJ = hclust(as.dist(dissADJ), method = "average")
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ,colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03,
				main = "Gene hierarchical clustering dendrogram and simulated module colors")

colorStaticADJ = as.character(cutreeStaticColor(hierADJ, cutHeight = 0.99, minSize = 20))

sizeGrWindow(10,5)
plotDendroAndColors(hierADJ, colors = data.frame(truemodule,colorStaticADJ), dendroLabels = FALSE, abHeight = 0.99,
			main = "Gene dendrogram and module colors")

branch.number = cutreeDynamic(hierADJ, method = "tree")

colorDynamicADJ = labels2colors(branch.number)

colorDynamicHybridADJ = labels2colors(cutreeDynamic(hierADJ,distM = dissADJ, cutHeight = 0.998,
			deepSplit = 2, pamRespectsDendro = FALSE))

sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ, colors = data.frame(truemodule, colorStaticADJ, colorDynamicADJ, colorDynamicHybridADJ),
			dendroLabels = FALSE, marAll = c(0.2,8, 2.7, 0.2),
		main = "Gene dendrogram and module colors")


hierTOM = hclust(as.dist(dissTOM),method = "average");

colorStaticTOM = as.character(cutreeStaticColor(hierTOM,cutHeight = .99, minSize = 20))
colorDynamicTOM = labels2colors(cutreeDynamic(hierTOM,method = "tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM,distM = dissTOM, cutHeight = 0.998,
				deepSplit = 2, pamRespectsDendro = FALSE))

sizeGrWindow(10,5)
plotDendroAndColors(hierTOM, colors = data.frame(truemodule, colorStaticTOM, 
				colorDynamicTOM,colorDynamicHybridTOM),
			dendroLabels = FALSE, marAll = c(1,8,3,1),
			main = "Gene dendrogram and module colors, TOM adjacency")

tabStaticADJ = table (colorStaticADJ, truemodule)
tabStaticTOM = table(colorStaticTOM, truemodule)
tabDynamicADJ = table(colorDynamicADJ,truemodule)
tabDynamicTOM = table(colorDynamicTOM, truemodule)
tabDynamicHybridADJ = table(colorDynamicHybridADJ,truemodule)
tabDynamicHybridTOM = table(colorDynamicHybridTOM,truemodule)

randIndex(tabStaticADJ, adjust = F)
randIndex(tabStaticTOM, adjust = F)
randIndex(tabDynamicADJ, adjust = F)
randIndex(tabDynamicTOM, adjust = F)
randIndex(tabDynamicHybridADJ, adjust = F)
randIndex(tabDynamicHybridTOM, adjust = F)

tabDynamicHybridTOM

colorh1 = colorDynamicHybridTOM

rm(ADJ1); rm(dissADJ);
collectGarbage();
save.image("Simulated-NetworkConstruction.RData");