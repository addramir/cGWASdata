
library(corrplot)
library(ppcor)
par(mfrow=c(1,2))

xb=read.table("data_table.txt",header=T,stringsAsFactors = F)
setwd("...")
load("data_for_plot.RData")

par(mfrow=c(1,2))
#FADS1
i=which(xb[,"Gene"]=="FADS1")
snp=xb[i,"SNP"]
trait=xb[i,"trait"]
cvrts=unlist(strsplit(xb[i,"bd_cvrts"],";"))
dta[is.na(dta[,snp]),snp]=mean(dta[,snp],na.rm=T)
sdta=cbind(snp=dta[,snp],phedata[,c(cvrts,trait)])

Mpc=pcor(sdta)$estimate
Mpc=Mpc[colnames(M),colnames(M)]

M <- cor(sdta[,c(trait,"snp",cvrts)])
M[lower.tri(M)]=0
M[2:(length(cvrts)+2),1]=lm(paste(trait,"~snp+",paste(cvrts,collapse="+")),sdta)$coef[2:(length(cvrts)+2)]

colnames(M)[2]=rownames(M)[2]="SNP"
colnames(M)[1]=rownames(M)[1]="lysoPC a C20:4"
colnames(M)[3]=rownames(M)[3]="lysoPC a C20:3"
diag(M)=0
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M,cl.pos="b",na.label="X",method = "square",
         addCoef.col = "black",col = col(200),number.cex = 1, tl.cex=2,tl.col=c("blue","darkblue","red"))

#ETFDH
i=which(xb[,"Gene"]=="ETFDH")
snp=xb[i,"SNP"]
trait=xb[i,"trait"]
cvrts=unlist(strsplit(xb[i,"bd_cvrts"],";"))
dta[is.na(dta[,snp]),snp]=mean(dta[,snp],na.rm=T)
sdta=cbind(snp=dta[,snp],phedata[,c(cvrts,trait)])

M <- cor(sdta[,c(trait,"snp",cvrts)])
M[lower.tri(M)]=0
M[2:(length(cvrts)+2),1]=lm(paste(trait,"~snp+",paste(cvrts,collapse="+")),sdta)$coef[2:(length(cvrts)+2)]

colnames(M)[2]=rownames(M)[2]="SNP"
diag(M)=0
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M,cl.pos="b",na.label="X",method = "square",
         addCoef.col = "black",col = col(200),number.cex = 1, tl.cex=2,tl.col=c("blue","darkblue","red","red"))

