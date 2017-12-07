########### Functions for exact estimations of cGWAS from uGWAS
########### Yakov Tsepilov
########### last maintained: 10.11.2016

# This file consist several fucntions used to calculate cGWAS results. 


#core functions
	#betas
	.b=function(response,pred,S){
		out=solve(S[pred,pred])%*%S[response,pred]
		return(out)
	}
	#standart errors of betas
	.seb=function(response,pred,S,N){
		sigma_joint=(S[response,response]-t(.b(response,pred,S=S))%*%S[response,pred])
		#sigma_joint=S[response,response]
		out=as.vector(sigma_joint)*solve(S[pred,pred])/(N-length(pred)-1)
		
		out=diag(out)
		out=sqrt(out)
		return(out)
	}
	.manova=function(response,pred,S,N){
		sigma_joint=(S[response,response]-t(.b(response,pred,S=S))%*%S[response,pred])
		Fst=((S[response,response]-as.vector(sigma_joint))/as.vector(sigma_joint))*((N-length(pred)-1)/length(pred))
		return(Fst)
	}
	#T statistics (not T^2)
	.t=function(response,pred,S,N){
		out=.b(response=response,pred=pred,S=S)/.seb(response=response,pred=pred,S=S,N=N)
		return(out)
	}

ReadCovariates=function(covariates,path2data,suffix=".txt"){
	Nprd=length(covariates)
	prd=covariates[1]
	for (prd in covariates){
		if (!exists(prd)){
			eval(text=paste(prd,"<-read.table(paste(path2data,prd,suffix,sep=\"\"),header=T,stringsAsFactors=F)",sep=""))
		}
	}
}

# The main function. It uses uGWAS results stored in .RData format. 
# CovM - covariance matrix of response vairable and covariates
# all_varg - the vector of variances of studied SNPs
# response - the name of response variable
# N - sample size
# covariates - the list of covariates
# cn_* - column names
# output_threshold - the threshold for p-value filter for storing of outputs
# gut_snps - the list of SNPs that should be used in anlisys (in case you want to filter them)
# all_CR - call rate for each SNP
# correction - should be GC correction applied?
# path_uGWAS - path to uGWAS in .RData format

exact_cGWAS=function(CovM,all_varg,response,N,covariates=NULL,
	cn_b="beta_SNP",cn_snp="SNP",cn_se="se_SNP",
	output_threshold=1e-6,gut_snps=NULL,all_CR=NULL,
	correction=TRUE,path_uGWAS){
	
	#checking the data
	cat("Load data...","\n")
	Nprd=length(covariates)
	#eval(text=paste("snps=",response",[,cn_snp]",sep=""))
	snps=names(all_varg)
	resp_pred<-c(response,covariates)
	
	prd=resp_pred[1]
	for (prd in resp_pred){
		load(paste(path_uGWAS,prd,".RData",sep=""))
	}
	
	cat("Filtering for NAs and searching overlapping SNPs...","\n")
	prd=resp_pred[1]
	for (prd in resp_pred){
		eval(parse(text=paste(prd,"=",prd,"[!is.na(",prd,"[,cn_b]),]",sep="")))
	}
	
	prd=resp_pred[1]
	for (prd in resp_pred){
		eval(parse(text=paste("snps_prd=",prd,"[,cn_snp]",sep="")))
		snps=intersect(snps,snps_prd)
	}
	
	if (!is.null(gut_snps))	snps=intersect(gut_snps,snps)
	if (!is.null(all_CR)) snps=intersect(names(all_CR),snps)
	
	if (is.null(snps)|length(snps)==0){
		warning("No overlapping SNPs between inputs")
		opt <- options(show.error.messages=FALSE) 
		on.exit(options(opt)) 
		stop() 
	#stop("No overlapping SNPs between inputs")
	}	
	
	snps=unique(snps)
	Nsnps=length(snps)
		
	all_varg=all_varg[snps]
	all_CR=all_CR[snps]
	prd=resp_pred[1]
	for (prd in resp_pred){
		eval(parse(text=paste("index=match(snps,",prd,"[,cn_snp])",sep="")))
		eval(parse(text=paste(prd,"=",prd,"[index,]",sep="")))
	}
	
	cat(Nsnps,"overlapping SNPs in total...","\n")
	
	if (is.null(covariates)|length(covariates)==0){
	
		cat("No covariates were specified. No cGWAS was performed...","\n")
		
		out=array(NA,c(Nsnps,4))
		colnames(out)=c("b","se","chi2","Pval")
		rownames(out)=snps
		
		prd=resp_pred[1]
		eval(parse(text=paste("out[,\"b\"]<-",prd,"[,cn_b]",sep="")))
		eval(parse(text=paste("out[,\"se\"]<-",prd,"[,cn_se]",sep="")))
		out[,"chi2"]<-(out[,"b"]/out[,"se"])^2
		
		snps=snps[!is.na(out[,"chi2"])]
		out=out[!is.na(out[,"chi2"]),]	
		
	} else{
	
		#creating of S matrix
		S=array(NA,c(Nprd+2,Nprd+2))
		colnames(S)=rownames(S)=c(response,covariates,"g")
		S[resp_pred,resp_pred] <- CovM[resp_pred,resp_pred]
		
		#calculation
		out=array(NA,c(Nsnps,4))
		colnames(out)=c("b","se","chi2","Pval")
		rownames(out)=snps
		
		#beta matrix
		betas=array(NA,c(Nsnps,length(resp_pred)))
		rownames(betas)=snps
		colnames(betas)=resp_pred
		
		for (prd in resp_pred){
			eval(parse(text=paste("betas[,prd]<-",prd,"[,cn_b]",sep="")))
		}
		cat("Starting cGWAS...","\n")
		
		i=1
		pb <- txtProgressBar(style=3)
		for (i in (1:Nsnps)){
			setTxtProgressBar(pb, value=i/Nsnps)
			S["g","g"]<-varg<-all_varg[i]
			prd=resp_pred[1]
			for (prd in resp_pred){
				b1=betas[i,prd]
				S["g",prd]<-S[prd,"g"]<-b1*varg
			}
			b<-.b(response=1,pred=2:(Nprd+2),S=S)["g",]
			se<-.seb(response=1,pred=2:(Nprd+2),S=S,N=N*all_CR[i])["g"]
			out[i,"b"]<-b
			out[i,"se"]<-se
			#out[i,"chi2"]<-(b/se)^2
		}
		out[,"chi2"]<-(out[,"b"]/out[,"se"])^2
		close(pb)
		
		#out=as.data.frame(out)
		#out=cbind(SNP=snps,out)
		snps=snps[!is.na(out[,"chi2"])]
		out=out[!is.na(out[,"chi2"]),]	
		cat("cGWAS was calculated for",dim(out)[1],"SNPs","\n")
	}	
	
	#output
	output=list()
	
	output$gc_lambda=median(out[,"chi2"],na.rm=T)/qchisq(0.5,1,lower.tail=F)
	cat("GC Lambda of all",dim(out)[1],"SNPs is",output$gc_lambda,"\n")
	
	if (correction){
		cat("Correcting...","\n")
		out[,"chi2"]=out[,"chi2"]/output$gc_lambda
		output$status="corrected for GC lambda"
	} else{
		output$status="Not corrected for GC lambda"
	}
	
	out[,"Pval"]=pchisq(out[,"chi2"],1,lower.tail=F)
		
	output$response=response
	output$covariates=paste(covariates,collapse=",")
	
	cat("Filtering results using threshold =",output_threshold,"\n")
	
	snps=snps[out[,"Pval"]<=output_threshold]
	out=out[out[,"Pval"]<=output_threshold,]
	
	cat(dim(out)[1],"SNPs with p-value <=",output_threshold,"\n")
	
	#colnames(out)=c("b","se","chi2","Pval")
	output$results=out
	output$snps=snps
	
	cat("Finished. Enjoy=)","\n")
	return(output)
}


###### Technical functions bellow

####clumping
function_for_making_full_table_without_gcv = function(pathRData,thr=5e-8,delta=2.5e5,list_of_phe=NULL){
	if (is.null(list_of_phe)){
		list_of_phe=system(paste("ls",pathRData,sep=" "),intern=T)
		list_of_phe=unlist(strsplit(list_of_phe,".RData"))
	}
	pheindex=1
	locus_table=as.data.frame(NULL)
	j=0
	for (pheindex in (1:length(list_of_phe))){	
		phe_name=list_of_phe[pheindex]
		file_name=paste(pathRData,phe_name,".RData",sep="")
		load(file_name)
		snps=Zx$snps
		Zx=Zx$results
		Zx=rbind(Zx,rep(1,4))
		Zx=as.data.frame(Zx)
		Zx=Zx[-dim(Zx)[1],]
		rownames(Zx)=snps
		colnames(Zx)=c("beta_SNP","se_SNP","Chi2","P-value")
		Zx=cbind(SNP=snps,Zx)
		
		Zx=Zx[Zx[,"P-value"]<=thr,]
		
		snps=Zx[,"SNP"]
		snps=unique(snps)
		snps=intersect(snps,gut_snps)
		index=match(snps,Zx[,"SNP"])
		Zx=Zx[index,]
		
		if (length(snps)>0){
			index=match(snps,rownames(snp_info))
			Chromosome=snp_info[index,"chr"]
			Position=snp_info[index,"pos"]
			Zx=cbind(Zx,Chromosome,Position)
			Zx=cbind(Zx,trait=phe_name)
			Zx=Zx[order(Zx$Chromosome,Zx$Position),]
			
			if (length(snps)>1){
				i=2
				while (i<=(dim(Zx)[1])){
					if ((abs(Zx$Position[i]-Zx$Position[i-1])<=delta)&(Zx$Chromosome[i]==Zx$Chromosome[i-1])){
						if (Zx[i,"P-value"]<=Zx[i-1,"P-value"]){
							Zx=Zx[-(i-1),]
						} else{
							Zx=Zx[-i,]
						}
					} else{
						i=i+1
					}
				}
				locus_table=rbind(locus_table,Zx)
			} else {
				locus_table=rbind(locus_table,Zx)
			}
		}
	}
	locus_table[,"P-value"]=as.numeric(locus_table[,"P-value"])
	locus_table[,"Position"]=as.numeric(locus_table[,"Position"])
	locus_table=locus_table[order(locus_table$Chromosome,locus_table$Position),]
	return(locus_table)
}

function_for_shlop_24_10_2013=function(locus_table,p_value="imp_P_value",pos="imp_pos",snp="imp_SNP",delta=5e5,chr="chr"){
	locus_table[,p_value]=as.numeric(locus_table[,p_value])
	locus_table[,pos]=as.numeric(locus_table[,pos])
	Zx <-locus_table
	Zx=Zx[order(Zx[,chr],Zx[,pos]),]
	n_traits=1
	Zx=cbind(Zx,n_traits)
	i=2
	while (i<=(dim(Zx)[1])){
		if ((abs(Zx[i,pos]-Zx[i-1,pos])<=delta)&(Zx[i,chr]==Zx[i-1,chr])){
			if (Zx[i,p_value]<=Zx[i-1,p_value]){
				Zx=Zx[-(i-1),]
				Zx[i-1,"n_traits"]=Zx[i-1,"n_traits"]+1
			} else{
				Zx=Zx[-i,]
				Zx[i-1,"n_traits"]=Zx[i-1,"n_traits"]+1
			}
		} else{
			i=i+1
		}
	}
	locus_table<-Zx
	rownames(locus_table)=as.character(locus_table[,snp])
	return(locus_table)
}

checking_function=function(pathRData,CovM,thr){
	list_of_phe=system(paste("ls",pathRData,sep=" "),intern=T)
	list_of_phe=unlist(strsplit(list_of_phe,".RData"))
	pheindex=1
	Nsnps_1e_6=0
	Nsnps_5e_8=0
	Nsnps_thr=0
	AGRH=FALSE
	for (pheindex in (1:length(list_of_phe))){	
		phe_name=list_of_phe[pheindex]
		file_name=paste(pathRData,phe_name,".RData",sep="")
		load(file_name)
		phe=Zx$response
		cvrts=Zx$covariates
		cvrts=unlist(strsplit(cvrts,","))
		
		cvrts_CovM=colnames(CovM)[which(CovM[phe,]>0)]
		if (length(cvrts)==length(cvrts_CovM)){
			if (length(cvrts)>0){
				if (sum(cvrts%in%cvrts_CovM)!=length(cvrts)){
					AGRH=TRUE
					cat(phe,"\n")
				}
			}
		} else {
			AGRH=TRUE
			cat(phe,"\n")
		}	
		snps=Zx$snps
		Zx=Zx$results
		Zx=rbind(Zx,rep(1,4))
		Zx=as.data.frame(Zx)
		Zx=Zx[-dim(Zx)[1],]
		rownames(Zx)=snps
		colnames(Zx)=c("beta_SNP","se_SNP","Chi2","P-value")
		Nsnps_1e_6=Nsnps_1e_6+length(which(Zx[,"P-value"]<=1e-6))
		Nsnps_5e_8=Nsnps_5e_8+length(which(Zx[,"P-value"]<=5e-8))
		Nsnps_thr=Nsnps_thr+length(which(Zx[,"P-value"]<=thr))
	}
	
	if (AGRH) cat("Something Wrong!","\n") else cat("Seems all right","\n")
	cat("SNPs with p-value <= 1e-6 :",Nsnps_1e_6,"\n")
	cat("SNPs with p-value <= 5e-8 :",Nsnps_5e_8,"\n")
	cat("SNPs with p-value <=",thr,":",Nsnps_thr,"\n")
}
