getwd();
workingDir = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/secondTutorial";
setwd(workingDir);
library(WGCNA);
options(stringsAsFactors = FALSE);

load("Simulated-dataSimulation.RData");
attach(ModuleEigengeneNetwork1)

GS1 = as.numeric(cor(y,datExpr,use = "p"))

p.Standard = corPvalueFisher(GS1,nSamples = length(y) )
p.Standard2 = p.Standard

p.Standard2[is.na(p.Standard)] = 1
q.Standard = qvalue(p.Standard2)$qvalues

StandardGeneScreeningResults = data.frame(GeneName,PearsonCorrelation = GS1, p.Standard, q.Standard)

head(StandardGeneScreeningResults)

NoiseGeneIndicator = is.element(truemodule,c("turquoise","blue","yellow","grey"))+.0
SignalGeneIndicator = 1 - NoiseGeneIndicator

mean(NoiseGeneIndicator[rank(p.Standard)<=20])
mean(NoiseGeneIndicator[rank(p.Standard)<=200])
mean(NoiseGeneIndicator[rank(p.Standard)<=100])

table(q.Standard <.20)

mean(NoiseGeneIndicator[q.Standard<=0.20])

save.image(file = "Simulated-StandardScreening.RData")