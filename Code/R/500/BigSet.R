getwd();
wd = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/500";
setwd(wd);
library(WGCNA);
library(cluster)
enableWGCNAThreads()
options(stringsAsFactors = FALSE);

bigSet = read.csv("mRNA_norm_combined.csv");

library(survival)
args(coxph)

# Check the coxph method. That is the one needed for cox analysis. Just think of the formula

