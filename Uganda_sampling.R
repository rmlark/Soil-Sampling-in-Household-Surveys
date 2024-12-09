#########################################################################
#  
# This script uses the variance component set from the Uganda data 
# in order to compute the standard error of a District mean computed as
# the combination of the BLUE of the overall mean and the BLUP of the 
# random effect of the sampled District.
#
# The output is shown in Figure 4 of the paper. But note that for 
# SOC the variance components are for log of SOC%, and the published
# graphs shows the width of +1 one standard error based on a mean SOC
# of 1.5% and 1.2% for topsoil and subsoil respectively.
#
# To run this script, first run the section which defines functions
# at the end.  Note that these functions are set up for the specific
# case where each district contains EA and HH within EA, where one
# sample is collected per HH, and where the supplied variance components
# are between-district, between-EA within-district and residual (between
# HH within EA).
#
#
#########################################################################
#
#
# Read in the variance components

VC<-read.table("VC_Uganda.dat",header=T)

# Which variable? Edit vname below to select the variable of interest
# and Vname to adjust figure labels

colnames(VC)[2:5]
vname<-"log_SOC_t"     
Vname<-"Topsoil SOC" # for use in Figures
var.index<-which(colnames(VC)==vname)

#######################################################################

n.District<-2 # We assume that two Districts are sampled

#  We set up a vector of fixed total sample sizes to consider
#  and colours to use in subsequent plots

TSS.vec<-c(100,200,400,600) 
TSS.col<-c("black","red","blue","green")

# We now loop over the TSS values

for (j in 1:4){
TSS<-TSS.vec[j]      # extract the total sample size
SWD<-TSS/n.District  # sample size within district 

noEA<-(FAC(SWD))$pos # find all factors of the within-district sample size,
			   # each to be treated as number of EA, with a corresponding
             	   # number of HH per EA given fixed sample size.	 

n.options<-length(noEA)

op<-matrix(0,nrow=n.options,ncol=4)
colnames(op)<-c("TSS","EA.per.District","HH.per.EA","SE")

for(i in 1:n.options){ # Loop over all the options for a given TSS
nEA.w.D<-noEA[i]
nHH.w.EA<-SWD/nEA.w.D
op[i,4]<-SE.BLUP.district(n.District,nEA.w.D,nHH.w.EA,VC,var.index)
op[i,1:3]<-c(TSS,nEA.w.D,nHH.w.EA)
}

if(j==1){
OP<-op
plot.col<-rep(TSS.col[j],n.options)
}else{
plot.col<-c(plot.col,rep(TSS.col[j],n.options))
OP<-rbind(OP,op)}
}


#
# The script below will plot the standard errors of the predictions
#
plot(OP[,2],OP[,4],pch=16,xlab="Number of EA per District",
ylim=c(0,(max(OP[,4]))),ylab="Standard error",col=plot.col,main=Vname,
log="x")
legend("topright",legend=TSS.vec,pch=16,col=TSS.col,title="Total sample size")
#lines(c(0.05,10000),c(0.05,0.05),lty=5)

#
# Below for making SOC plot with back-transformation.  Mean SOC is soc.mean set
# to 1.5 for topsoil and 1.2 for subsoil

soc.mean<-1.5

Width.SE<-exp(log(soc.mean)+OP[,4])-exp(log(soc.mean)-OP[,4])

plot(OP[,2],Width.SE,pch=16,xlab="Number of EA per District",
ylim=c(0,(max(Width.SE))),ylab="Width of ± one standard error",col=plot.col,main=Vname,
log="x")
legend("topright",legend=TSS.vec,pch=16,col=TSS.col,title="Total sample size")
#lines(c(0.001,10000),c(0.2,0.2),lty=5)




#########################################################################
#########################################################################
#
#  FUNCTIONS
#
#########################################################################
#########################################################################

#########################################################################
#
#  This function computes the SE for a district BLUP, given the number
#  of districts (this is not varied, but affects SE via the BLUE for the
#  overall mean), the number of EA within each district, the number of HH
#  within each EA, the variance components and the variable of interest. 

SE.BLUP.district<-function(n.District,nEA.w.D,nHH.w.EA,VC,var.index){

n.EA<-n.District*nEA.w.D
n.HH<-n.EA*nHH.w.EA
n<-n.HH
r<-n.District+n.EA

#########################################################################
#
#  Make design matrix for random effects, Z
#

Z<-matrix(0,nrow=n,ncol=r)

District.code<-seq(1,n.District)

EA.code<-seq(1,n.EA)

HH.code<-seq(1,n.HH)

# unified code is from 1 - r

codes<-(cbind(EA.code,District.code))[rep(1:n.EA,times=nHH.w.EA),]
add.term<-t(as.vector(c(0,n.EA)))

unified.codes<-sweep(codes,2,add.term, "+")

for (i in 1:n){
Z[i,unified.codes[i,]]<-1
}
Zt<-t(Z)


#########################################################################
#
#  Make design matrix for fixed effects, X
#

X<-matrix(1,nrow=n,ncol=1)
Xt<-t(X)

# Make R structure

R<-diag(n)*VC[nrow(VC),var.index]
Rinv<-solve(R)

# Make G structure

G<-diag(c(rep(VC[2,var.index],n.EA),rep(VC[1,var.index],n.District)))
Ginv<-solve(G)

#

##############################################################################################
#
#  Now make error covariance matrix for BLUEs and BLUPs, Ci
#

W<-rbind(Xt,Zt)
p<-nrow(W)-nrow(G)
Opp<-matrix(0,nrow=p,ncol=p)
Orp<-matrix(0,nrow=r,ncol=p)
Opr<-matrix(0,nrow=p,ncol=r)

Gstar<-rbind(cbind(Opp,Opr),cbind(Orp,solve(G)))

C<-W%*%Rinv%*%t(W)+Gstar

Ci<-solve(C)

#########################################################################
#
#  Make a vector, lam, to combine the BLUEs and BLUPs into the BLUP for
#  a district mean, then find its standard error.  This is returned.

m<-r+1
lam<-matrix(0,nrow=m,ncol=1)
lam[1]<-1
lam[m]<-1

SE<-as.numeric(sqrt(t(lam)%*%Ci%*%lam))

return(SE)
}

#########################################################################
FAC<- function(x) {
    x <- as.integer(x)
    div <- seq_len(abs(x))
    factors <- div[x %% div == 0L]
    factors <- list(neg = -factors, pos = factors)
    return(factors)
}



