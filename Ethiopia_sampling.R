#######################################################################
#  
# This script uses the variance component set from the Ethiopia data 
# in order to compute the standard error of a Zone mean computed as
# the combination of the BLUE of the overall mean and the BLUP of the 
# random effect of the sampled Zone.
#
# The output is shown in Figure 5 of the paper.
#
# Read in the variance components from VC_Ethiopia.dat
#
# In the proposed set of sample designs we consider just one
# field per parcel, so the between-parcel variance and residual
# variance (the last two rows in VC_Ethiopa.data) are merged 
# into a common residual.
#
# To run this script, first run the section which defines functions
# at the end.  Note that these functions are set up for the specific
# case where eachzone contains EA and HH within EA, where one
# two samples are collected per HH (different parcels).
#
#
#########################################################################
#
#
#  Read in the variance components, and merge the last two rows

VC.prov<-read.table("VC_Ethiopia.dat",header=T)
VC<-VC.prov[-4,2:5]
VC[4,1:4]<-VC.prov[4,2:5]+VC.prov[5,2:5]

# Which variable? Edit vname below to select the variable of interest
# and Vname to adjust figure labels

colnames(VC)[1:4]
vname<-"SOC_t"
Vname<-"Topsoil SOC"
var.index<-which(colnames(VC)==vname)

#######################################################################

n.Zone<-2  	# We assume that two Zones are sampled
nP.w.HH<-2	# number of parcels within each HH 

#  We set up a vector of fixed total sample sizes to consider
#  and colours to use in subsequent plots

TSS.vec<-c(100,200,400,600)
TSS.col<-c("black","red","blue","green")

# We now loop over the TSS values

for (j in 1:4){
TSS<-TSS.vec[j]      		# extract the total sample size
EA.HH<-TSS/(n.Zone*nP.w.HH) 	# number of EA per HH

noEA<-FAC(EA.HH)$pos # find all factors of the within-district sample size,
			   # each to be treated as number of EA, with a corresponding
             	   # number of HH per EA given fixed sample size.

n.options<-length(noEA)

op<-matrix(0,nrow=n.options,ncol=5)
colnames(op)<-c("TSS","EA.per.Zone","HH.per.EA","Parcel.per.HH","SE")

for(i in 1:n.options){

nEA.w.Z<-noEA[i] # number of EA within each district
nHH.w.EA<-EA.HH/nEA.w.Z # number of HH within each EA
op[i,1:4]<-c(TSS,nEA.w.Z,nHH.w.EA,nP.w.HH)
op[i,5]<-SE.BLUP.zone(n.Zone,nEA.w.Z,nHH.w.EA,nP.w.HH,VC,var.index)
}

if(j==1){
OP<-op
plot.col<-rep(TSS.col[j],n.options)
}else{
plot.col<-c(plot.col,rep(TSS.col[j],n.options))
OP<-rbind(OP,op)}
}

plot(OP[,2],OP[,5],pch=16,xlab="Number of EA per Zone",
ylim=c(0,(max(OP[,5]))),ylab="Standard error",col=plot.col,main=Vname,
log="x")
legend("topright",legend=TSS.vec,pch=16,col=TSS.col,title="Total sample size")
#lines(c(0.001,10000),c(0.075,0.075),lty=5)



#########################################################################
#########################################################################
#
#  FUNCTION
#
#########################################################################
#########################################################################

#########################################################################
#
#  This function computes the SE for a zone BLUP, given the number
#  of zones (this is not varied, but affects SE via the BLUE for the
#  overall mean), the number of EA within each zone, the number of HH
#  within each EA, the number of parcels within each HH the variance 
#  components and the variable of interest. 

SE.BLUP.zone<-function(n.Zone,nEA.w.Z,nHH.w.EA,nP.w.HH,VC,var.index){

n.EA<-n.Zone*nEA.w.Z
n.HH<-n.EA*nHH.w.EA
n.P<-n.HH*nP.w.HH
n<-n.P
r<-n.Zone+n.EA+n.HH


#########################################################################
#
#  Make design matrix for random effects, Z
#

Z<-matrix(0,nrow=n,ncol=r)

Zone.code<-seq(1,n.Zone)

EA.code<-seq(1,n.EA)

HH.code<-seq(1,n.HH)


# unified code is from 1 - r

codes<-(cbind(HH.code,EA.code,Zone.code))[rep(1:n.HH,times=nP.w.HH),]
add.term<-t(as.vector(c(0,n.HH,n.HH+n.EA)))

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

G<-diag(c(rep(VC[3,var.index],n.HH),rep(VC[2,var.index],n.EA),
rep(VC[1,var.index],n.Zone)))
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




