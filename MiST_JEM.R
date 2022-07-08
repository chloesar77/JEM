#!/usr/bin/Rscript

# Script to perform MiST & JEM as well as famSKAT and single linear regression
# Written by Chloe Sarnowski (chloesar@bu.edu)
# 04/06/2020
# Note: use R > 3.2.0

args = commandArgs(trailingOnly=TRUE)

###
#args[1]=N simulations
#args[2]=N SNPs
#args[3]="LD" or "noLD"
#args[4]=r2 value if LD
#args[5]=annotation: "binary", "quant_random", "quant_perfect","quant_2distrib", "noannot"
#args[6]=N causal SNPs
#args[7]=N causal SNPs with misspecified annotation 1
#args[8]=N causal SNPs with misspecified annotation 2
#args[9]=fixed effect for the annotation (pi), "fixed" for equal effect size for all causal SNPs, or "weight" for unequal effect size for causal SNPs.
#args[10]=maximum value for the weight (weight will range from 1 to maximum value, need to be proportional to N causal SNPs, must be different from 0)
#args[11]=R2, proportion of variance explained
#args[12]="H0" or "H1"
#args[13]="test.ped" haplotype file from HapGen with simulated family structure (5 first columns are family id, individual id,father id,mother id,sex)
#args[14]="test.map" map file (4 columns for each SNP: chromosome, rsID, pos in cM, pos in hg19)
#args[15]= N individuals to extract randomly at each simulation
#args[16]=N other SNPS with misspecfied annotation 1
#args[17]=N other SNPS with misspecfied annotation 2
#args[18]=influence of frequency of SNPs: "mixed","rare" ,"lowfreq","common"
###

# Example of command line to run the script
#Rscript JEM.R 1 50 "noLD" 0 "binary" 25 5 5 "fixed" 0 0.02 "H1" extract_test.ped extract_test.map 2000 0 0 "mixed"
#Rscript JEM.R 1 50 "noLD" 0 "binary" 25 5 5 "weight" 5 0.02 "H1" extract_test.ped extract_test.map 2000 0 0 "mixed"
#Rscript JEM.R 1 50 "LD" 0.9 "binary" 25 5 5 "fixed" 0 0.02 "H1" extract_test.ped extract_test.map 2000 0 0 "mixed"
#Rscript JEM.R 1 50 "noLD" 0 "binary" 25 5 5 "fixed" 0 0.05 "H1" extract_test.ped extract_test.map 2000 0 0 "mixed"

require(gaston)
library('data.table')
library('MASS')
library('kinship2')
library('seqMeta')
library('coxme')

#open the ped file
print("reading the ped file")
ped=fread(args[13], head=F, sep=" ", data.table=FALSE)
ped <- as.matrix(ped)

#open the map file
print("reading the map file")
map=fread(args[14], head=F, sep="\t", data.table=FALSE)
ped_info=ped[,1:5]
ped_info=as.matrix(ped_info)
colnames(ped_info)=c("famid","id","fa","mo","sex")

#determine the number of individuals
Nind=dim(ped_info)[1]

#recode the genotypes from haplotypes
geno=matrix( NA, nrow = Nind, ncol = dim(map)[1])
i<- 1:dim(map)[1]
temp1=5+2*(i-1)+1
temp2=5+2*(i-1)+2
geno=ped[,temp1]+ped[,temp2]
colnames(geno) = map[,2]
checkfreq <- (.colSums(geno==2, m=nrow(geno), n=ncol(geno)) + .colSums(geno==1, m=nrow(geno), n=ncol(geno))/2 )/Nind
# code 2 for minor allele, thus freq of allele should be less than 0.50
geno[,checkfreq >=0.5] <- 2 - geno[,checkfreq >=0.5]
checkfreq <- (.colSums(geno==2, m=nrow(geno), n=ncol(geno)) + .colSums(geno==1, m=nrow(geno), n=ncol(geno))/2 )/Nind

#determine the number of SNPs
q=as.numeric(args[2])

if(args[18]=="rare"){
	#SNPs=sample(which(checkfreq>=0.005 & checkfreq<0.01),q)
	#the number of SNPS is fixed in the simulations and we also use the same SNPs
	SNPs=c(20897,55683,21966,6793,77266,101990,50199,101018,61285,49180,24282,104164,27171,20327,57091,25335,29445,54562,34405,12808,98312,6378,21635,25515,72660,37268,51115,38131,95666,50607,10462,8742,79413,55019,92071,49008,45893,16186,105476,48646,64269,90383,66243,13977,54356,58135,3294,19741,48273,64025)
}

if(args[18]=="common"){
	#SNPs=sample(which(checkfreq>=0.05),q)
	#the number of SNPS is fixed in the simulations and we also use the same SNPs
	SNPs=c(86139,26209,15295,38204,10991,64927,6847,69406,25111,5400,87656,23287,73023,76552,27566,26872,100297,27632,92888,114266,72446,114629,114007, 69276, 		43134,109601,103924, 65855, 85901, 81539,26299,15261,88898,27970,11751,105919,13880,20432,66650,43046,67760,80595,96744,59578,82393,20529,54691,82681,79459,71583)

}

if(args[18]=="lowfreq"){
	#SNPs=sample(which(checkfreq>=0.01 & checkfreq<0.05),q)
	#the number of SNPS is fixed in the simulations and we also use the same SNPs
			SNPs=c(6099,20226,50974,15146,20355,4684,114539,42521,65777,12935,78066,81386,65442,114698,99419,110036,109346,48045,100070,78702,43186,15148,111632,100854,69445,8727,13655,71350,61373,42106,7854,1380,88341,53022,112265,85462,6509,61190,89608,113052,112139,111711,11433,54282,38229,41131,102746,19598,97878,98575)
}

if(args[18]=="mixed"){
	#SNPs=sample(which(checkfreq>=0.005),q)
	#the number of SNPS is fixed in the simulations and we also use the same SNPs
	SNPs=c(9825,52465,93957,82238,31811,7947,73294,100751,75758,54179,72800,81272,10541,111694,43423,75684,15960,31134,106329,73471,3653,43085,62296,61817,19118,62251,105926,114409,3154,76747,44916,74823,43003,60009,13216,28801,9591,50837,18405,7568,114655,49322,59454,61228,10461,26383,90627,100266,19983,28973)
}

#fix the number of simulations
Nsim <- as.numeric(args[1])

res <- matrix(NA,Nsim,4*q+5)
NSNPs_effect=as.numeric(args[6])
NSNPs_noneffect=q-NSNPs_effect
NSNPs_effect_wrongannot1=as.numeric(args[7])
NSNPs_noneffect_wrongannot1=as.numeric(args[8])
NSNPs_effect_wrongannot2=as.numeric(args[16])
NSNPs_noneffect_wrongannot2=as.numeric(args[17])

#determine the number of individuals
Nindu=as.numeric(args[15])
vect2=rep(1,NSNPs_effect)

#generate the annotations
if(args[5]=="binary"){
	vect=c(rep(0,NSNPs_effect_wrongannot1),rep(1,NSNPs_effect-NSNPs_effect_wrongannot1),rep(1,NSNPs_noneffect_wrongannot1),rep(0,NSNPs_noneffect-NSNPs_noneffect_wrongannot1))
	vect1=c(rep(1,NSNPs_effect_wrongannot2),rep(0,NSNPs_effect-NSNPs_effect_wrongannot2),rep(0,NSNPs_noneffect_wrongannot2),rep(1,NSNPs_noneffect-NSNPs_noneffect_wrongannot2))
} else if(args[5]=="noannot"){
	vect=c(rep(1,NSNPs_effect_wrongannot1),rep(1,NSNPs_effect-NSNPs_effect_wrongannot1),rep(1,NSNPs_noneffect_wrongannot1),rep(1,NSNPs_noneffect-NSNPs_noneffect_wrongannot1))		
}

if(args[9]=="fixed"){
	Z2=matrix(rep(1,length(vect2)),ncol=1,nrow=length(vect2))
} else if(args[9]=="weight"){
	weights=as.numeric(args[10])
	Z2=matrix(rep(c(1:weights),NSNPs_effect/weights),ncol=1,nrow=NSNPs_effect)		
}

R2=as.numeric(args[11])

#for each simulation
for (k in 1:Nsim){
	print(paste("Simulation",k,sep=""))
	ind=sample(dim(geno)[1],Nindu)
	ped_info2=ped_info[which(ped_info[,2]%in%ind),]
  
	#extract q variants from the genotype matrix
	x2=geno[ind,SNPs]
	checkfreq2=checkfreq[SNPs]
  
  	# simulate LD separately for causal and non causal SNPs
	if(args[3]=="LD"){
		print("generating LD")
    		VarX_simu1=matrix(sqrt(as.numeric(args[4])),nrow=NSNPs_effect,NSNPs_effect)
		diag(VarX_simu1)=1
		mat1=VarX_simu1
		VarX_simu2=matrix(sqrt(as.numeric(args[4])),nrow=NSNPs_noneffect,NSNPs_noneffect)
		diag(VarX_simu2)=1
		mat2=VarX_simu2
		final_mat=matrix(0,q,q)
		final_mat[1:(NSNPs_effect),1:(NSNPs_effect)]=mat1
		final_mat[(NSNPs_effect+1):q,(NSNPs_effect+1):q]=mat2
		VarX_simu=final_mat
		SNPs_simu=mvrnorm(as.numeric(Nindu),Sigma=VarX_simu,mu=rep(0,q))
		MAFf=sample(seq(0.01,0.05,by=0.001),q,replace=T)
		q1=qnorm((1-MAFf)*(1-MAFf))
		q2=qnorm((1-MAFf)*(1-MAFf) + 2*MAFf*(1-MAFf))
		newSNPs=SNPs_simu
		colnames(newSNPs)=c(1:q)
    
    		for(i in 1:q){
			newSNPs[which(newSNPs[,i] <= q1[i]),i] = 0
			newSNPs[which(newSNPs[,i] > q1[i] & newSNPs[,i] <= q2[i]),i] = 1
			newSNPs[which(newSNPs[,i] > q2[i]),i] = 2
		}
    
    		geno_std=scale(newSNPs)
		geno_std3=newSNPs
  	}
  
  
	if(args[3]=="all_in_LD"){
   		print("generating LD")
		VarX_simu=matrix(sqrt(as.numeric(args[4])),nrow=q,ncol=q)
		diag(VarX_simu)=1
		SNPs_simu=mvrnorm(as.numeric(Nindu),Sigma=VarX_simu,mu=rep(0,q))
		MAFf=sample(seq(0.01,0.05,by=0.001),q,replace=T)
		q1=qnorm((1-MAFf)*(1-MAFf))
		q2=qnorm((1-MAFf)*(1-MAFf) + 2*MAFf*(1-MAFf))
		newSNPs=SNPs_simu
		colnames(newSNPs)=c(1:q)
    
		for(i in 1:q){
			newSNPs[which(newSNPs[,i] <= q1[i]),i] = 0
			newSNPs[which(newSNPs[,i] > q1[i] & newSNPs[,i] <= q2[i]),i] = 1
			newSNPs[which(newSNPs[,i] > q2[i]),i] = 2
    		}
    
    		geno_std=scale(newSNPs)
		geno_std3=newSNPs
	  }

  
	if(args[3]=="noLD"){
		geno_std=scale(x2)
		geno_std3=x2
	}
    
	if(args[5]=="quant_random"){
		vect <- c(sample(seq(from=0,to=30,by=0.1),q,replace=T))
	} else if(args[5]=="quant_perfect"){
		vect <- c(sample(seq(from=15,to=20,by=0.1),NSNPs_effect,replace=T),sample(seq(from=0,to=5,by=0.1),NSNPs_noneffect,replace=T))
	} else if(args[5]=="quant_2distrib"){
		vect=c(sample(rnorm(Nind,22,5.7) ,NSNPs_effect,replace=T),sample(rnorm(Nind,3,2.8) ,NSNPs_noneffect,replace=T))
	} 

	geno_std2=geno_std[,1:NSNPs_effect]
  	  
	#generate the matrix of annotations Z	
	if(args[16]!=0){
		Z=matrix(c(vect,vect1),ncol=2,nrow=q)
	} else if(args[16]==0){
		Z=matrix(vect,ncol=1,nrow=q)
	}
 	# the first column of the annotation matrix is a column of 1 to include an intercept
	Z=cbind(rep(1,q),Z)	
	# by default the fixed effect of the annotation is fixed to 0.1
	pi=matrix(c(0,0.1),ncol=2,nrow=1)
	Z2=cbind(rep(1,NSNPs_effect),Z2)
  
	# define the matrix XZ, genotype matrix multiplied by annotations
	XZ= geno_std %*% Z
	XZ2= geno_std2 %*% Z2
  
	# define K1 matrix
	K1 = geno_std %*% t(geno_std)
  
	#define variance of residuals
	W=var(XZ2%*%t(pi))   
	vare = (W*(1-R2))/R2
  	e=rnorm(dim(geno_std2)[1], sd =sqrt(vare))
  
	# generate the phenotype y
	if(args[12]=="H1"){
		# Generation of y under the alternate model (Y depends on X)
		y <- XZ2%*%t(pi) + e
		y=y[,1]
	} else if(args[12]=="H0"){
		# Generation of y under the null model (Y is independent of X)
		y <- rnorm(Nindu)
  	}  

	# MiST: joint test #
  	print("performing MiST")
	mu=mean(y)
	Yu=y-mu	
	X1=rep(1,Nindu)
	
	deltamoins1=diag(var(y),nrow=Nindu,ncol=Nindu)
	delta=solve(deltamoins1)
	deltasqrt=sqrt(deltamoins1)
	H1=deltasqrt %*% X1 %*% solve(X1%*%deltamoins1%*%X1) %*% t(X1) %*% deltasqrt
	IN=diag(1,nrow=Nindu,ncol=Nindu)

	if(args[16]!=0){
		Z=matrix(c(vect,vect1),ncol=2,nrow=q)
	} else if(args[16]==0){
		Z=matrix(vect,ncol=1,nrow=q)
	}
  
	B= geno_std %*% Z
	V=t(B) %*% deltamoins1 %*% (IN-H1) %*% delta %*% (IN-H1) %*% deltamoins1 %*% B
	U=t(Yu) %*% B %*% solve(V) %*% t(B) %*% Yu
	score_p1=1-pchisq(U,dim(Z)[2])

	Z=cbind(rep(1,q),Z)
	B= geno_std %*% Z

	if(args[5]!="noannot"){
		V=t(B) %*% deltamoins1 %*% (IN-H1) %*% delta %*% (IN-H1) %*% deltamoins1 %*% B
		U=t(Yu) %*% B %*% solve(V) %*% t(B) %*% Yu
		score_p2=1-pchisq(U,dim(Z)[2]) 
	}
  
  	if(args[5]=="noannot"){
		score_p2=1 
	}

	#famSKAT (to compare power with MiST)
	print("performing famSKAT")
	SNPInfo=matrix(NA,q,3)
	SNPInfo[,1]=colnames(geno_std3)
	SNPInfo[,2]=rep(1,q)
	SNPInfo[,3]=Z[,2]
	colnames(SNPInfo)=c("Name","gene","Annot")
	o=prepScores(geno_std3, y~1,SNPInfo,family = stats::gaussian(), kins = makekinship(ped_info2[,1], ped_info2[,2], ped_info2[,3], ped_info2[,4]), sparse = TRUE, verbose = FALSE)
	p=skatMeta(o, SNPInfo=SNPInfo, wts = "Annot", method = "saddlepoint", mafRange = c(0, 0.5), verbose = FALSE)
	p1=skatMeta(o, SNPInfo=SNPInfo, wts = function(maf){1/(maf*(1-maf))}, method = "saddlepoint", mafRange = c(0, 0.5), verbose = FALSE) 
	p2=skatMeta(o, SNPInfo=SNPInfo, wts = function(maf){stats::dbeta(maf,1,25)}, method = "saddlepoint", mafRange = c(0, 0.5), verbose = FALSE) 
  
  	skat_pvalue=p$p
	skat_pvalue1=p1$p  
	skat_pvalue2=p2$p  
  
	#JEM
	print("performing JEM")
	estimates <- lmm.aireml(y, X = XZ, K = K1, verbose = FALSE, get.P=TRUE)
	tau=estimates$tau
	py = estimates$Py
	varpi=estimates$varbeta
	pi_est=estimates$BLUP_beta
	delta2=diag(rep(tau,q)) %*% t(geno_std) %*% py
	vardelta2=(diag(rep(tau,q)) %*% t(geno_std) %*% estimates$P %*% geno_std %*% diag(rep(tau,q)))
	varbeta = Z%*%varpi%*%t(Z) + vardelta2
	beta_est=as.vector(Z%*%pi_est) + delta2
	sdbeta=as.vector(sqrt(diag(varbeta)))
	Zbeta=beta_est/sdbeta
	Chi2=Zbeta*Zbeta
	pbeta = 2* (1-pnorm(abs(Zbeta)))
	results=cbind(Z[,2],beta_est,delta2,sdbeta,Zbeta,Chi2,pbeta)
	colnames(results)=c("Annot","Beta_est","Delta_est","SD","Z","Chi2","P")
  
	#simple linear regression (to compare effect sizes with JEM)
	kmat= makekinship(famid=ped_info2[,1], id=ped_info2[,2], mother.id=ped_info2[,3], father.id=ped_info2[,4])
	beta_est_lin=c()
	pbeta_lin=c()
	for(i in 1:dim(geno_std)[2]){
		s=lmekin(y~geno_std[,i]+(1|ped_info2[,2]), varlist=list(kmat), na.action=na.omit)
		beta <- s$coefficients$fixed
		nvar <- length(beta)
		nfrail <- nrow(s$var) - nvar
		se <- sqrt(diag(s$var)[nfrail + 1:nvar])
		z<- round(beta/se, 2)
		p<- signif(1 - pchisq((beta/se)^2, 1), 2)
		beta_est_lin=c(beta_est_lin,beta[2])
		pbeta_lin=c(pbeta_lin,p[2])
	}

  	res[k,]=c(as.vector(t(beta_est)),as.vector(t(pbeta)),as.vector(t(beta_est_lin)),as.vector(t(pbeta_lin)),score_p1,score_p2,skat_pvalue,skat_pvalue1,skat_pvalue2)
	
}

# write the results
print("writing the results")
write.table(res,paste("Simu_",args[12],"_",q,"SNPs_",args[3],"_",args[4],"_",args[5],"_",args[6],"eff_",args[7],"effwr_",args[8],"noeffwr_",args[9],"_R",args[11],"_",args[16],"effwr_",args[17],"noeffwr_",args[18],".csv",sep=""),sep=",",row.names=F,col.names=F,quote=F)

