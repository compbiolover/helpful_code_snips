#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")

library(limma)
library(sva)
setwd("D:\\biowolf\\82metabolism\\10.intersect")

rt = read.table("tcgaMetabExp.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
metab=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
metab=avereps(metab)

rt = read.table("geoMatrix.txt",header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)
geo=log2(geo+1)

sameGene=intersect(row.names(metab),row.names(geo))
metabOut=metab[sameGene,]
geoOut=geo[sameGene,]

all=cbind(metabOut,geoOut)
batchType=c(rep(1,ncol(metabOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType,par.prior=TRUE)
metabOut=outTab[,colnames(metabOut)]
metabOut=rbind(ID=colnames(metabOut),metabOut)
write.table(metabOut,file="tcgaMetabExp.share.txt",sep="\t",quote=F,col.names=F)
geoOut=outTab[,colnames(geoOut)]
geoOut=rbind(ID=colnames(geoOut),geoOut)
write.table(geoOut,file="geoMetabExp.share.txt",sep="\t",quote=F,col.names=F)


#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

library("limma")
library(pheatmap)

setwd("D:\\biowolf\\metabolism\\11.diff")               
inputFile="tcgaMetabExp.share.txt"                                  
fdrFilter=0.05                                                    
logFCfilter=0.5                                                    
conNum=32                                                         
treatNum=375                                                         

outTab=data.frame()
grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data[data<0]=0

for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

write.table(outTab,file="all.xls",sep="\t",row.names=F,quote=F)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.xls",sep="\t",row.names=F,quote=F)

heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="tcgaDiffMetabExp.txt",sep="\t",col.names=F,quote=F)

pdf(file="vol.pdf",height=5,width=5)
xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

hmExp=data[as.vector(outDiff[,1]),]
hmExp=log2(hmExp+0.1)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=7,width=10)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=3,
         fontsize_col=10)
dev.off()
#install.packages('survival')

library(survival)
pFilter=0.05                                                        
setwd("D:\\biowolf\\metabolism\\14.uniCox")         
rt=read.table("tcgaExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)  

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outTab,file="tcgaUniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="tcgaUniSigExp.txt",sep="\t",row.names=F,quote=F)


rt <- read.table("tcgaUniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

pdf(file="forest.pdf", width = 6,height = 4.5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#install.packages("glmnet")
#install.packages("survival")


library("glmnet")
library("survival")

setwd("D:\\biowolf\\metabolism\\15.lasso")                
rt=read.table("tcgaUniSigExp.txt",header=T,sep="\t",row.names=1)           
rt$futime=rt$futime/365

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)

coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="tcgaRisk.txt",sep="\t",quote=F,row.names=F)

rt=read.table("geoExpTime.txt",header=T,sep="\t",row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
rt$futime=rt$futime/365
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="geoRisk.txt",sep="\t",quote=F,row.names=F)

#install.packages("survival")
#install.packages("survminer")

library(survival)
library("survminer")
setwd("D:\\biowolf\\metabolism\\16.survival")             

bioSurvival=function(inputFile=null,outFile=null){
  rt=read.table(inputFile,header=T,sep="\t")                   
  diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=TRUE,
                     pval=paste0("p=",pValue),
                     pval.size=4,
                     risk.table=TRUE,
                     legend.labs=c("High risk", "Low risk"),
                     legend.title="Risk",
                     xlab="Time(years)",
                     break.time.by = 1,
                     risk.table.title="",
                     palette=c("red", "blue"),
                     risk.table.height=.25)
  pdf(file=outFile,onefile = FALSE,width = 6.5,height =5.5)
  print(surPlot)
  dev.off()
}
bioSurvival(inputFile="tcgaRisk.txt",outFile="tcgaRisk.pdf")
bioSurvival(inputFile="geoRisk.txt",outFile="geoRisk.pdf")


#install.packages("pheatmap")

library(pheatmap)
setwd("D:\\biowolf\\metabolism\\17.riskPlot")            

bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null,heatmapFile=null){
  rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)   
  rt=rt[order(rt$riskScore),]                                          
  
  
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  line=rt[,"riskScore"]
  line[line>10]=10
  pdf(file=riskScoreFile,width = 10,height = 3.5)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("green",lowLength),rep("red",highLength)) )
  abline(h=median(rt$riskScore),v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  dev.off()
  
  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="green"
  pdf(file=survStatFile,width = 10,height = 3.5)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=log2(rt1+1)
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=heatmapFile,width = 10,height = 3.5)
  pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           fontsize_row=11,
           show_colnames = F,
           fontsize_col=3,
           color = colorRampPalette(c("green", "black", "red"))(50) )
  dev.off()
}
bioRiskPlot(inputFile="geoRisk.txt",riskScoreFile="geo.riskScore.pdf",survStatFile="geo.survStat.pdf",heatmapFile="geo.heatmap.pdf")
bioRiskPlot(inputFile="tcgaRisk.txt",riskScoreFile="tcga.riskScore.pdf",survStatFile="tcga.survStat.pdf",heatmapFile="tcga.heatmap.pdf")

library(survival)
setwd("D:\\biowolf\\metabolism\\18.tcgaIndep")                       
risk=read.table("tcgaRisk.txt",header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table("tcgaClinical.txt",sep="\t",check.names=F,header=T,row.names=1)    
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])

uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="tcga.uniCox.txt",sep="\t",row.names=F,quote=F)

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="tcga.multiCox.txt",sep="\t",row.names=F,quote=F)


bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  pdf(file=forestFile, width = 6.3,height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}

bioForest(coxFile="tcga.uniCox.txt",forestFile="tcga.uniForest.pdf",forestCol="green")
bioForest(coxFile="tcga.multiCox.txt",forestFile="tcga.multiForest.pdf",forestCol="red")

#install.packages('survival')

library(survival)
setwd("D:\\biowolf\\metabolism\\19.geoIndep")                        
risk=read.table("geoRisk.txt",header=T,sep="\t",check.names=F,row.names=1)          
cli=read.table("geoClinical.txt",sep="\t",check.names=F,header=T,row.names=1)      
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])

uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="geo.uniCox.txt",sep="\t",row.names=F,quote=F)

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="geo.multiCox.txt",sep="\t",row.names=F,quote=F)

bioForest=function(coxFile=null,forestFile=null,forestCol=null){
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  
  pdf(file=forestFile, width = 6.3,height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
  axis(1)
  dev.off()
}

bioForest(coxFile="geo.uniCox.txt",forestFile="geo.uniForest.pdf",forestCol="green")
bioForest(coxFile="geo.multiCox.txt",forestFile="geo.multiForest.pdf",forestCol="red")

#install.packages("survivalROC")

library(survivalROC)
setwd("D:\\biowolf\\metabolism\\20.ROC")                         

bioROC=function(riskFile=null,cliFile=null,outFile=null){
  risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
  rocCol=rainbow(ncol(rt)-2)
  aucText=c()
  
  pdf(file=outFile,width=6,height=6)
  par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
  roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, predict.time =1, method="KM")
  plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[1], 
       xlab="False positive rate", ylab="True positive rate",
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  aucText=c(aucText,paste0("risk score"," (AUC=",sprintf("%.3f",roc$AUC),")"))
  abline(0,1)
  
  j=1
  for(i in colnames(rt[,3:(ncol(rt)-1)])){
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i], predict.time =1, method="KM")
    j=j+1
    aucText=c(aucText,paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
    lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j],lwd = 2)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
  dev.off()
}
bioROC(riskFile="tcgaRisk.txt",cliFile="tcgaClinical.txt",outFile="tcgaROC.pdf")
bioROC(riskFile="geoRisk.txt",cliFile="geoClinical.txt",outFile="geoROC.pdf")




library(rms)
setwd("C:\\Users\\lexb4\\Desktop\\21.Nomogram")                         

riskFile="tcgaRisk.txt"
cliFile="tcgaClinical.txt"
outFile="tcga.Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)         
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
dd <- datadist(rt)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~ ., x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)

nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
                lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
                maxscale=100, 
                fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
pdf(file=outFile,height=7.5,width=11)
plot(nom)
dev.off()

riskFile="geoRisk.txt"
cliFile="geoClinical.txt"
outFile="geo.Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
dd <- datadist(rt)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~ ., x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)

nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
                lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
                maxscale=100, 
                fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
pdf(file=outFile,height=6,width=9)
plot(nom)
dev.off()



#install.packages("ggplot2")

library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("D:\\biowolf\\metabolism\\24.multipleGSEA")           
files=grep(".xls",dir(),value=T)                                        
data = lapply(files,read.delim)                                          
names(data) = files

dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".xls","",dataSet$.id)                           

gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,fill=pathway,group=pathway))+
  geom_point(shape=21) + scale_fill_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits =c(min(dataSet$RUNNING.ES-0.02), max(dataSet$RUNNING.ES+0.02))) +   
  theme_bw() + theme(panel.grid =element_blank()) + theme(panel.border = element_blank()) + 
  theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) + guides(fill=guide_legend(title = NULL)) + 
  theme(legend.background = element_blank()) + theme(legend.key = element_blank())
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "High risk<----------->Low risk", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

pdf('multipleGSEA.pdf',      
    width=9,                
    height=5)               
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()
