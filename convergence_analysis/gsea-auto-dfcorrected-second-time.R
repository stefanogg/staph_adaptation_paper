# Get relative position and ordering of core SNPs in a collection of velvet assemblies
help = paste(
"gsea-auto.R Gene set enrichment analysis testing for enrichment of substitutions in some ontological classes",
"Daniel Wilson (2015)",
"",
"Usage: Rscript gsea-auto.R gene_info.txt nsubs_by_gene.txt gene_ontologies.txt results-prefix 1> log_file.txt 2> error_log_file.txt",
sep="\n")

# Defaults
COUTPUT = FALSE
# For a range of dispersion parameters
DISP = 10^(-2:0); # This is now (one over) the shape parameter of the gamma prior IN THE ALTERNATIVE MODEL ONLY
nsim = 10000
# Prior probs?

# Read options from command line
args = commandArgs(trailingOnly = TRUE)
if(length(args)!=4) {
	cat(help,sep="\n")
	stop("\nIncorrect usage\n")
}
geneinfo_filename = args[1]
nsubs_filename = args[2]
ontologies_filename = args[3]
results_prefix = args[4]

# Functions
# Now using Negative Binomial distributions rather than Multinomial and Multinomial Dirichlet distributions
dnb = function(x,L,a,b,log=FALSE) {
	ret = x*log(L)+a*log(b)+lgamma(x+a)-lfactorial(x)-lgamma(a)-(x+a)*log(L+b)
	if(log) return(ret)
	return(exp(ret))
}
dmultinb = function(x,L,a,b,log=FALSE) {
	ret = sum(x*log(L))+a*log(b)+lgamma(sum(x)+a)-sum(lfactorial(x))-lgamma(a)-(sum(x)+a)*log(sum(L)+b)
	if(log) return(ret)
	return(exp(ret))
}
# Define the probability only up to a constant of proportionality to allow improper priors (i.e. a=0 or b=0)
dmultinb0 = function(x,L,a,b,log=FALSE) {
	ret = sum(x*log(L))+lgamma(sum(x)+a)-sum(lfactorial(x))-(sum(x)+a)*log(sum(L)+b)
	if(log) return(ret)
	return(exp(ret))
}
# Log mean of logged values
logmean = function(x,na.rm=TRUE) { mx=max(x,na.rm=na.rm); mx+log(mean(exp(x-mx),na.rm=na.rm)) }

geneinfo = read.delim(geneinfo_filename,as.is=T,h=T,sep="\t",comment.char='')
bipclass = read.delim(nsubs_filename,as.is=T,h=T,sep="\t",comment.char='')
ontologies = as.matrix(read.delim(ontologies_filename,as.is=T,h=T,sep="\t",comment.char='',check.names=FALSE))

# Check
if(!all(rownames(bipclass)==geneinfo$Name)) stop("Row names of nsubs_by_gene.txt do not match Names column in gene_info.txt")
if(!all(rownames(ontologies)==geneinfo$Name)) stop("Row names of gene_ontologies.txt do not match Names column in gene_info.txt")

y = as.numeric(bipclass[,1]); x1 = as.numeric(geneinfo$Length)
if(any(is.na(y))) stop("Non-numeric values found in first column of nsubs_by_gene.txt")
if(any(is.na(x1))) stop("Non-numeric values found in gene_info.txt Length column")

FAC = levels(factor(ontologies)); NFAC = length(FAC)
if(NFAC>100) stop("Very large number of levels found in ontologies (>10) may produce slow computation and/or low power")
if(NFAC>10) warning("Large number of levels found in ontologies (>10) may produce slow computation and/or low power")

# For general (up/down/un-regulated) analysis
ALPHA = 1/DISP
# Find the MAP estimate of lambda0 under the null
MAP = sum(y)/sum(x1)
BF = matrix(NA,ncol(ontologies),length(DISP))
# Now an uninformative improper log-uniform distribution is used for the rate parameter under the null
# and the loglikelihood under the null (which is the same regardless of ALPHA) is defined up to a constant only
NL = dmultinb0(y,x1,0,0,log=TRUE)
fit0 = glm(y~0+x1,family=poisson(link="identity"))
# Number of substitutions and total length of genes in each gene class
NSUB = array(NA,c(ncol(ontologies),NFAC))
TLEN = array(NA,c(ncol(ontologies),NFAC))
ELAMBDA0 = matrix(NA,ncol(ontologies),length(DISP))
VLAMBDA0 = matrix(NA,ncol(ontologies),length(DISP))
PCLASS = rep(NA,ncol(ontologies))
INTEGRAND = list()
for(i in 1:ncol(ontologies)) {
#	Choose an experiment. For each gene define up/down/unknown
	updown = ontologies[,i]; x2 = factor(updown)
	mtch = match(levels(x2),FAC)
	nsub = xtabs(y~x2); tlen = xtabs(x1~x2)
	NSUB[i,mtch] = nsub
	TLEN[i,mtch] = tlen
	for(j in 1:length(DISP)) {
		lognonintegrand = function() {
			nfac = length(nsub)
			const = sum(y*log(x1))-sum(lfactorial(y))+(nfac*ALPHA[j])*log(ALPHA[j])-nfac*lgamma(ALPHA[j])+sum(lgamma(nsub+ALPHA[j]))
			const - NL
		}
		get_integrand = function(nsub,tlen,j,offset=0)
			Vectorize(function(lambda0) {
				nfac = length(nsub)
				exp(offset + (-nfac*ALPHA[j]-1)*log(lambda0)+sum(-(nsub+ALPHA[j])*log(tlen+ALPHA[j]/lambda0)))
			})
		integrand0 = get_integrand(nsub,tlen,j,0)
		OFFSET = lognonintegrand()
#		integrand = Vectorize(function(lambda0) integrand0(lambda0,OFFSET))
		integrand = get_integrand(nsub,tlen,j,OFFSET)
		BF[i,j] = log(integrate(integrand,0,Inf)$value)
#		Approximate using gamma distribution with specified mean and variance
#		integrand = Vectorize(function(lambda0) integrand0(lambda0,OFFSET-BF[i,j]))
		integrand = get_integrand(nsub,tlen,j,OFFSET-BF[i,j])
		integrate(function(x) integrand(x),0,Inf,stop.on.error=F)$value
		ELAMBDA0[i,j] = integrate(function(x) x*integrand(x),0,Inf,stop.on.error=F,subdivisions=1001)$value
		VLAMBDA0[i,j] = integrate(function(x) x*x*integrand(x),0,Inf,stop.on.error=F,subdivisions=1001)$value - ELAMBDA0[i,j]^2
		INTEGRAND[[paste0(i,":",j)]] = integrand
	}
	fit1 = glm(y~0+x1:x2,family=poisson(link="identity"),start=rep(coef(fit0),sum(table(x2)>0)))
	PCLASS[i] = pchisq(-fit1$deviance+fit0$deviance,nlevels(x2)-1,low=F)
	if(COUTPUT & (i%%10)==0) cat("Done",i,"of",ncol(ontologies),"\n")
}


# Simulate from the posterior with the following probabilities
gd = 1:length(DISP)
mnBF = apply(BF[,gd],1,logmean,na.rm=FALSE)
pexpt = exp(mnBF); pexpt = pexpt/sum(pexpt)
#barplot(pexpt[order(pexpt)],names=substr(colnames(ontologies)[order(pexpt)],1,5),las=3,cex.names=.5)
#data.frame("Expt"=colnames(ontologies),"Post.prob"=pexpt,"BF"=exp(mnBF))[order(pexpt,decreasing=TRUE)[1:30],]

rexpt.mat = sample.int(length(BF[,gd]),nsim,replace=TRUE,prob=exp(BF[,gd]))
rexpt = matrix(1:nrow(BF),nrow(BF),ncol(BF))[rexpt.mat]
wdsim = matrix(gd,nrow(BF),ncol(BF),byrow=T)[rexpt.mat]

# Partial Bayes factors should behave correctly, but computationally intensive
# Use cross-validation, which is anti-conservative because the data are used twice, once in fitting and once in assessing fit (although this is not the purpose here)
#LAMBDA0 = rgamma(nsim,ELAMBDA0[cbind(rexpt,wdsim)]^2/VLAMBDA0[cbind(rexpt,wdsim)],ELAMBDA0[cbind(rexpt,wdsim)]/VLAMBDA0[cbind(rexpt,wdsim)]); #LAMBDA0[LAMBDA0<1e-20] = 1e-20
# The gamma distribution approximation is poor in some cases. Instead, just plug in posterior mean estimate.
LAMBDA0 = ELAMBDA0[cbind(rexpt,wdsim)]
# Although the effect on the hyperparameter cannot be removed, ameliorate the double use of the data (there is some modest change of ordering)
ELAMBDA = rep(0.0,nrow(ontologies)); VLAMBDA = rep(0.0,nrow(ontologies)); xvBF = matrix(0.0,nsim,nrow(ontologies))
xvBF0 = (y*log(x1)+(sum(y)-y)*log(sum(x1)-x1)+lgamma(sum(y))-lfactorial(y)-lgamma(sum(y)-y)-(sum(y))*log(sum(x1)))
for(i in 1:nsim) {
	ix = rexpt[i]
	updown = ontologies[,ix]; x2 = factor(updown)
	mtch = match(levels(x2),FAC)
	APOST.loc = NSUB[ix,mtch[x2]]+ALPHA[wdsim[i]]-y
	BPOST.loc = TLEN[ix,mtch[x2]]+ALPHA[wdsim[i]]/LAMBDA0[i]-x1
	ELAMBDA = ELAMBDA + APOST.loc/BPOST.loc
	VLAMBDA = VLAMBDA + APOST.loc/BPOST.loc^2
	xvBF[i,] = y*log(x1)+APOST.loc*log(BPOST.loc)+lgamma(y+APOST.loc)-lfactorial(y)-lgamma(APOST.loc)-(y+APOST.loc)*log(x1+BPOST.loc) - xvBF0
}
# Post-process
ELAMBDA = ELAMBDA/nsim
VLAMBDA = VLAMBDA/nsim
xvBF = apply(xvBF,2,logmean,na.rm=FALSE)

# Approximate the parameters per category, conditional on ontology
# No: Approximate on the log-scale by assuming multiplication of independent gamma distributions, and manipulating the mean and variance on the log scale, before back-transforming
# No: E = digamma(a) - log(b?); V = trigamma(a)
# Delta approximation:
# Second order gives negative results so simply use first order
if(FALSE) {
ELAMBDAFAC = t(sapply(1:ncol(ontologies), function(i) {
	pr = exp(BF[i,]-mnBF[i]-log(length(DISP)))
	sapply(1:NFAC,function(j) {
		f = function(lambda0) (ALPHA+NSUB[i,j])/(ALPHA/lambda0+TLEN[i,j])
#		ddf = function(lambda0) -2*(ALPHA+NSUB[i,j])*ALPHA*TLEN[i,j]/(ALPHA+TLEN[i,j]*lambda0)^3
		   ECOND = f(ELAMBDA0[i,]); #+0.5*ddf(ELAMBDA0[i,])*VLAMBDA0[i,]
		sum(ECOND*pr)
	})
}))
VLAMBDAFAC = t(sapply(1:ncol(ontologies), function(i) {
	pr = exp(BF[i,]-mnBF[i]-log(length(DISP)))
	sapply(1:NFAC,function(j) {
		   f = function(lambda0) (ALPHA+NSUB[i,j])/(ALPHA/lambda0+TLEN[i,j])^2
#		   ddf = function(lambda0) -2*ALPHA*(ALPHA+NSUB[i,j])*(2*TLEN[i,j]*lambda0-ALPHA)/(ALPHA+TLEN[i,j]*lambda0)^4
		   ECOND = f(ELAMBDA0[i,]); #+0.5*ddf(ELAMBDA0[i,])*VLAMBDA0[i,]
		   sum(ECOND*pr)
		   })
}))
plot(ELAMBDAFAC,VLAMBDAFAC)
}

ELAMBDAFAC = matrix(NA,ncol(ontologies),NFAC); VLAMBDAFAC = matrix(NA,ncol(ontologies),NFAC)
for(i in 1:ncol(ontologies)) {
	pr = exp(BF[i,]-mnBF[i]-log(length(DISP)))
	for(k in 1:NFAC) if(!is.na(NSUB[i,k])) {
		get_fun = function(i) {
			wh = paste0(i,":",1:length(DISP))
			Vectorize(function(lambda0) {
				sapply(wh,function(WH) INTEGRAND[[WH]](lambda0))
			})
		}
		fun = get_fun(i)
		get_integrand = function(FUN,i,pr)
			Vectorize(function(lambda0) {
#				sum(pr*(ALPHA+NSUB[i,k])/(ALPHA/lambda0+TLEN[i,k])*dgamma(lambda0,ELAMBDA0[i,]^2/VLAMBDA0[i,],ELAMBDA0[i,]/VLAMBDA0[i,]))
#					  sum(pr*dgamma(lambda0,ELAMBDA0[i,]^2/VLAMBDA0[i,],ELAMBDA0[i,]/VLAMBDA0[i,]))
					sum(pr*(ALPHA+NSUB[i,k])/(ALPHA/lambda0+TLEN[i,k])*FUN(lambda0))
#					sum(pr*FUN(lambda0))
			})
		integrand = get_integrand(fun,i,pr)
		ELAMBDAFAC[i,k] = integrate(integrand,0,Inf,stop.on.error=F)$val
		get_integrand = function(FUN)
		Vectorize(function(lambda0) {
#				sum(pr*(ALPHA+NSUB[i,k])/(ALPHA/lambda0+TLEN[i,k])^2*dgamma(lambda0,ELAMBDA0[i,]^2/VLAMBDA0[i,],ELAMBDA0[i,]/VLAMBDA0[i,]))
				  sum(pr*(ALPHA+NSUB[i,k])/(ALPHA/lambda0+TLEN[i,k])^2*FUN(lambda0))
				  })
		integrand = get_integrand(fun)
		VLAMBDAFAC[i,k] = integrate(integrand,0,Inf,stop.on.error=F)$val
	}
	if(COUTPUT & (i%%10)==0) cat("Done",i,"of",ncol(ontologies),"\n")	
}
if(FALSE) {
plot(NSUB/TLEN,ELAMBDAFAC,ylim=range(ELAMBDAFAC-2*sqrt(VLAMBDAFAC),ELAMBDAFAC+2*sqrt(VLAMBDAFAC),na.rm=T))
arrows(NSUB/TLEN,ELAMBDAFAC-2*sqrt(VLAMBDAFAC),NSUB/TLEN,ELAMBDAFAC+2*sqrt(VLAMBDAFAC),len=0,col=rep(1:3,each=nrow(ELAMBDAFAC)))
abline(0,1,col="grey",lty=3)

1->k; plot(NSUB[,k]/TLEN[,k],ELAMBDAFAC[,k],ylim=range(ELAMBDAFAC-2*sqrt(VLAMBDAFAC),ELAMBDAFAC+2*sqrt(VLAMBDAFAC),na.rm=T))
arrows(NSUB[,k]/TLEN[,k],ELAMBDAFAC[,k]-2*sqrt(VLAMBDAFAC[,k]),NSUB[,k]/TLEN[,k],ELAMBDAFAC[,k]+2*sqrt(VLAMBDAFAC[,k]),len=0,col=k)
abline(0,1,col="grey",lty=3)
}

# Order loci by cross-validation partial Bayes factor
od = order(xvBF,decreasing=TRUE)
out = data.frame("Locus"=geneinfo$Name,"nsub"=y,"kb"=x1,"enrich"=y*sum(x1)/sum(y)/x1,xvBF,ELAMBDA,"SDLAMBDA"=sqrt(VLAMBDA))
#head(out[od,],n=30)
write.table(out,paste0(results_prefix,".results.gene.txt"),row=F,col=T,quote=F,sep="\t",na="-")

# Order ontologies by Bayes factor
od = order(pexpt,decreasing=TRUE)
colnames(NSUB) = FAC; colnames(TLEN) = FAC; colnames(ELAMBDAFAC) = FAC; colnames(VLAMBDAFAC) = FAC
out = data.frame("Ontology"=colnames(ontologies),"Post.prob"=pexpt,"mnBF"=exp(mnBF),"nsub"=NSUB,"tlen"=TLEN,"ELAMBA"=ELAMBDAFAC,"SDLAMBDA"=sqrt(VLAMBDAFAC),"glm.pval"=PCLASS,check.names=F)
#head(out[od,],n=30)
write.table(out,paste0(results_prefix,".results.ontology.txt"),row=F,col=T,quote=F,sep="\t",na="-")
