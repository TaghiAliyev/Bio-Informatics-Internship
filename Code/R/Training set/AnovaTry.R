getwd();
wd = "C:/Users/Taghi/eclipseWorkspace/School/Project/Bio-Informatics-Internship/Code/R/Training set";
setwd(wd);
library(WGCNA);
library(cluster)
enableWGCNAThreads()
options(stringsAsFactors = FALSE);

r = t(datExprToFile)
k1 = 11

f1 = c("CEAr17Norma", "CEAr33Norma", "CEAr35Norma", "CEAr57Norma","CEAr100Norm", 
"CEAs14Norma", "CEAs25Norma", "CEAs36Norma"
 , "CEAs41Norma", "CEAs45Norma", "CEAs81Norma")

sdata<-read.table("sandberg-sampledata.txt", header=T, row.names=1)
strain <- gl(2,12,24, label=c("129","bl6"))
region <- gl(6,2,24, label=c("ag", "cb", "cx", "ec", "hp", "mb"))
# 2 groups : stable and ruptured. Summation of these is what is interesting
groups <- gl(2,5,11, label = c("r","s"))
numbers <- gl(5,1,11, label = c("1","2","3","4","5","6"))
aof <- function(x) {
	m <- data.frame(groups,numbers,x);
	anova(aov(x ~ groups * numbers,m))
}

anovaresults2 <- apply(r,1,aof)

temp.1 <- unlist(lapply(anovaresults,function(x) { x["Pr(>F)"][1:2,]}))
temp.2 <- matrix(temp.1,length(anovaresults),2,byrow = T)
dimnames(temp.2) <- list(names(anovaresults), dimnames(anovaresults[[1]])[[1]][1:2])
pvalues <- data.frame(t(temp.2))

reg.hi.p <- t(data.frame(pvalues[2,pvalues[2,] < 0.0001 & pvalues[3,] > 0.1]))
reg.hi.pdata <- sdata[row.names(reg.hi.p),]

# FDR
bh <- function(x, fdr) {
  thresh <- F;
  crit<-0;
  len<-length(x)
  answer <- array(len);
  first <- T;
  for(i in c(len:0)) {
    crit<-fdr*i/len;
    if (x[i] < crit || thresh == T) {
      answer[i]<-T
      thresh <- T
      if (first) {
        cat(i ,"genes selected at FDR =", fdr ,"\n")
        first = F;
      }
    } else {
      answer[i]<-F
    }
  }
  answer
}
regionp <- sort(t(pvalues[2,])) # the p values must be sorted in increasing order!
fdr.result <- bh(regionp, 0.1) # reports that 192 genes are selected
bhthresh <- cbind(regionp, fdr.result)