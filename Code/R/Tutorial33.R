getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R";
setwd(workingDir);
library(WGCNA)
options(stringsAsFactors = FALSE);

# Loading of file
lnames = load(file = "FemaleLiver-01-dataInput.RData");
lnames = load(file = "FemaleLive-02-networkConstruction-auto.RData");

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6);

textMatrix = paste(signif(moduleTraitCor,2), "\n(", signif(moduleTraitPvalue,1),")",sep = "");

dim(textMatrix) = dim(moduleTraitCor);
par(mar = c(6,8.5,3,3));

labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs),
		colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5,
		zlim = c(-1,1), main = paste("Module-trait relationships"))

weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"

modNames = substring(names(MEs),3)

geneModuleMembership = as.data.frame(cor(datExpr,MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples));

names(geneModuleMembership) = paste("MM",modNames,sep = "");
names(MMPvalue) = paste("p.MM",modNames,sep="");

geneTraitSignificance = as.data.frame(cor(datExpr,weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples));

names(geneTraitSignificance) = paste("GS.",names(weight),sep = "");
names(GSPvalue) = paste("p.GS.",names(weight),sep = "");


module = "salmon"
column = match(module, modNames);
moduleGenes = moduleColors == module;

sizeGrWindow(7,7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes,1]),
		xlab = paste("Module Membership in", module, "module"),
		ylab = "Gene significance for body weight",
		main = paste("Module membership vs. gene significance \n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes,annot$substanceBXH);
sum(is.na(probes2annot))


geneInfo0 = data.frame(substanceBXH = probes, geneSymbol = annot$gene_symbol[probes2annot],
			LocusLinkID = annot$LocusLinkID[probes2annot], moduleColor = moduleColors,
			geneTraitSignificance, GSPvalue);

modOrder = order(-abs(cor(MEs,weight,use = "p")));

for (mod in 1:ncol(geneModuleMembership))
{
	oldNames = names(geneInfo0)
	geneInfo0 = data.frame(geneInfo0,geneModuleMembership[,modOrder[mod]], MMPvalue[,modOrder[mod]]);
	names(geneInfo0) = c(oldNames,paste("MM.", modNames[modOrder[mod]],sep = ""),paste("p.MM.",modNames[modOrder[mod]],sep=""));
}

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder,]

write.csv(geneInfo, file = "geneInfo.csv")