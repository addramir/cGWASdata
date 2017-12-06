#####
##### Butiful example - figure 2 
#####
library(corrplot)
library(ppcor)
par(mfrow=c(1,2))

xb=read.table("20161110_big.txt",header=T,stringsAsFactors = F)
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
corrplot(M, method="circle",cl.pos="b")

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
corrplot(M, method="circle",cl.pos="b")

