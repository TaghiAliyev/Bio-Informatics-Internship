getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/secondTutorial";
setwd(workingDir);
library(WGCNA);
options(stringsAsFactors = FALSE);

datGeneSummary = read.csv("GeneSummaryTutorial.csv")
datTraits = read.csv("TraitsTutorial.csv")
datMicroarrays = read.csv("MicroarrayDataTutorial.csv")

datMicroarrays[1:5,1:4]

ArrayName = names(data.frame(datMicroarrays[,-1]))

GeneName = datMicroarrays$GeneName

datExpr = data.frame(t(datMicroarrays[,-1]))
names(datExpr) = datMicroarrays[,1]
dimnames(datExpr)[[1]] = names(data.frame(datMicroarrays[,-1]))

truemodule = datGeneSummary$truemodule
rm(datMicroarrays)
collectGarbage()

datExpr[1:5,1:5]

table(dimnames(datExpr)[[1]] == datTraits$ArrayName)

y = datTraits$y