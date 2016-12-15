rank(c(5,3,2,3,2,1))
? rank

max(rank(c(5,3,2,3,2,1)))+1-rank(c(5,3,2,3,2,1))

library(igraph)
library(BioNet)
set.seed(45)
g <- erdos.renyi.game(10, 0.5, directed = TRUE)
get.adjacency(g)
#10 x 10 sparse Matrix of class "dgCMatrix"

[1,] . 1 . 1 . . 1 1 1 1
[2,] . . . 1 . 1 1 1 1 .
[3,] 1 1 . 1 . 1 . . 1 .
[4,] 1 1 . . . . 1 1 1 1
[5,] . 1 . . . . 1 . . 1
[6,] . 1 1 1 1 . 1 . . .
[7,] 1 . 1 . . 1 . . 1 .
[8,] . 1 1 1 . 1 . . 1 .
[9,] 1 . 1 1 . 1 . . . 1
[10,] 1 . 1 . 1 1 . . 1 .

edges <- as.data.frame(get.edgelist(g))
edges <- edges[order(edges$V1, edges$V2), ]
names.edges<-paste(as.character(edges$V1),as.character(edges$V2),sep="-")
names(edges)<-paste(as.character(edges$V1),as.character(edges$V2),sep="-")
#input: the data of the correlations,the names of two nodes for this correlation must be specified
#Re-order the names of the nodes in the network for the correlations
#undirected graph, re-order the names of the names of first node and second node by combining the first node and the second node by "-"
#so that the names of the edges become unique
#specify the names of two nodes for the correlations
#order the names of first node and the second node, and if the names of two nodes for two different correlations are the same
#there must some error

getwd()

D:/Research project 2015/Protein network analysis/Real data/CATHGEN/

#set up the dataset from Cathgen project

#pair group id
library(survival)
setwd("D:/Research project/Protein network analysis/Real data/CATHGEN/")
#response dataset
dat.resp<-read.table("ResponseVector.txt")
dim(dat.resp)
#covariates dataset
dat.var<-read.table("Data.txt")
n<-dim(dat.var)[1]
p<-dim(dat.var)[2]
corr.true<-cor(dat.var, use="complete.obs")
rownames(corr.true)
thre<-0.1
mat.edge<-NA
names.mat.edge<-NA
if(length(which(is.na(rownames(corr.true))))==0 && length(which(is.na(colnames(corr.true))))==0 ){
for (i in 1:(p-1)){
  for (j in (i+1):p ){
    if(abs(corr.true[i,j])>thre){
      mat.edge<-c(mat.edge,corr.true[i,j])


      names.mat.edge<-c(names.mat.edge,paste(rownames(corr.true)[i],colnames(corr.true)[j],sep="-"))
      
    }}}}
if(is.null(rownames(corr.true)) && is.null(colnames(corr.true)) ){
for (i in 1:(p-1)){
  for (j in (i+1):p ){
    if(abs(corr.true[i,j])>thre){
      mat.edge<-c(mat.edge,corr.true[i,j])
      
      
      names.mat.edge<-c(names.mat.edge,paste(paste("V",as.character(i),sep=""),paste("V",as.character(j),sep=""),sep="-"))
      
    }}}}
edge.corr.thre<-mat.edge
names(edge.corr.thre)<-names.mat.edge




#calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
#inputs:dat-the dataset containing the variables for correlations,nsim: the number of permutation samples
pval.perm.corr.submat<-function(dat,nsim){
  n<-dim(dat)[1]
  #the number of variales for correlatins
  p<-dim(dat)[2]
  #the number of correlations from p variables
  pp<-p*(p-1)/2
  mat.corr.perm<-matrix(NA,ncol=pp,nrow=nsim)
  
  for ( ii in 1:nsim){

    #permuation random samples for further correlation
    var.perm<-apply(dat,2,function(x){return(sample(x,size=n))})
    #correlaton from nsim permutation samples
    corr.perm<-cor(var.perm, use="complete.obs") 
    mat.corr.perm[ii,]<-trans.mat.vec(corr.perm)}
  
  #Pearson's correlation from real data
  corr.true<-trans.mat.vec(cor(dat, use="complete.obs"))
  #calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
  p.val.perm<-c()
  for (j in 1:pp){
    p.val.perm<-c(p.val.perm,(sum(abs(mat.corr.perm[,j])>abs(corr.true[j])))/nsim)}
  from<-c()
  to<-c()
  for (i in 1:(dim(dat)[2]-1)){
    from<-c(from,rep(names(dat)[i],length((i+1):dim(dat)[2])))
    to<-c(to,names(dat)[-c(1:i)])}
  names(p.val.perm)<-paste(from,to,sep="_")
  return(p.val.perm)}



pval.perm.corr<-function(dat,nsim){
#permutation test 
#sample size(the number of observations-n
  n<-dim(dat)[1]
#the number of variales for correlatins
p<-dim(dat)[2]
n.submat<-ceiling(p/20)
if(n.submat==1){
  return(pval.perm.corr.submat(dat,nsim))}
if(n.submat>1){
  pval.perm.corr.sub<-c()
  pval.perm.corr.cross<-c()
  for(i in 1:n.submat){
    pval.perm.corr.sub<-c(pval.perm.corr.sub,pval.perm.corr.submat(dat=dat[,c(((i-1)*20+1):min(p,20*i))],nsim=nsim))}
  for(i in 1:(n.submat-1)){
     for (j in (i+1):n.submat){
    pval.perm.corr.cross<-c(pval.perm.corr.cross,pval.perm.corr.cross(dat1=dat[,c(((i-1)*20+1):min(p,20*i))],dat2=dat[,c(((j-1)*20+1):min(p,20*j))],nsim=nsim))}}
 return(c(pval.perm.corr.sub,pval.perm.corr.cross))
}}
 
  dat.var[,c(4:6)]
ppval.test<-pval.perm.corr(dat=dat.var[,1:32],nsim=15)

length(ppval.test)
ppval.test[1:100]
min(4,5)

  pval.perm.corr<-function(dat,nsim){
    n<-dim(dat)[1]
    #the number of variales for correlatins
    p<-dim(dat)[2]
#the number of correlations from p variables
pp<-p*(p-1)/2

mat.corr.perm<-matrix(NA,ncol=pp,nrow=nsim)

for ( i in 1:nsim){
  dat<-dat.var
#permuation random samples for further correlation
var.perm<-apply(dat,2,function(x){return(sample(x,size=n))})
#correlaton from nsim permutation samples
corr.perm<-cor(var.perm, use="complete.obs") 
mat.corr.perm[i,]<-trans.mat.vec(corr.perm)[,3]}
dim(trans.mat.vec(corr.perm))
#Pearson's correlation from real data
corr.true<-trans.mat.vec(cor(dat, use="complete.obs"))
#calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
p.val.perm<-c()
for (j in 1:pp){
p.val.perm<-c(p.val.perm,sum(abs(mat.corr.perm[,j])>abs(corr.true[j]))/length(corr.true))}
from<-c()
to<-c()
for (i in 1:(dim(dat)[2]-1)){
  from<-c(from,rep(names(dat)[i],length((i+1):dim(dat)[2])))
    to<-c(to,names(dat)[-c(1:i)])}
names(p.val.perm)<-paste(from,to,sep="_")
return(p.val.perm)}

#partion the size of 472 variables into 30 groups with relatively smaller size
p.val.perm.var<-pval.perm.corr(dat.var,nsim=100)







pval.perm.corr.cross<-function(dat1,dat2,nsim){
  n<-dim(dat1)[1]
  #the number of variales for correlatins
  p1<-dim(dat1)[2]
  p2<-dim(dat2)[2]
  #the number of correlations from p variables
  pp<-p1*p2
  
  mat.corr.perm<-matrix(NA,ncol=pp,nrow=nsim)
  dat<-cbind(dat1,dat2)
  for ( i in 1:nsim){
    #permuation random samples for further correlation
    var.perm<-apply(dat,2,function(x){return(sample(x,size=n))})
    #correlaton from nsim permutation samples
    corr.perm<-cor(var.perm[,1:ncol(dat1)],var.perm[,(1+ncol(dat1)):(ncol(dat1)+ncol(dat2))], use="complete.obs") 
    mat.corr.perm[i,]<-trans.mat.vec.cross(corr.perm)}
  #Pearson's correlation from real data
  corr.true<-trans.mat.vec.cross(cor(dat1,dat2, use="complete.obs"))
  #calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
  p.val.perm<-c()
  for (j in 1:pp){
    p.val.perm<-c(p.val.perm,sum(abs(mat.corr.perm[,j])>abs(corr.true[j]))/nsim)}
    from<-rep(names(dat1),ncol(dat2))
    to<-rep(names(dat2),each=ncol(dat1))
  names(p.val.perm)<-paste(from,to,sep="_")
  return(p.val.perm)}



? rep(c("a","b"),2)
rep(c("a","b"),each=2)



corr.x<-pval.perm.corr.cross(dat.var[,1:20],dat[,21:40], nsim=15)
corr.x














sam.order<-sam[order(sam)]
if(length(which((sam.order-3.5)==0))>=1)
  p.val<-mean(which((sam.order-3.5)==0)))/length(sam.order)
if()
if(length(which(abs(sam.order-3.5)==min(abs(sam.order-3.5)))>=1 && (sam.order-3.5)>=0)
  p.val<-which(abs(sam.order[(sam.order-3.5)>=0,]-3.5)==min(abs(sam.order[(sam.order-3.5)>=0]-3.5)))]/length(sim.order)
and (sam-3.5)>=0)
rank(sam)/length(sam)

sam<-sample(1:10,size=10)

sum(sam<=3.5)/length(sam)





#the function-trans.mat.vec,transform the pXp correlation matrix to the vector with the length-p(p-1)/2
#input:x-correlation matrix
#output:the vector of all correlations among p variables with length-p(p-1)/2
trans.mat.vec<-function(x){
  corr.pearson.vec<-c()
for (i in 1:(dim(x)[1]-1)){
  corr.pearson.vec<-c(corr.pearson.vec,x[-c(1:i),i])
}
return(corr.pearson.vec)
}


trans.mat.vec.cross<-function(x){
  corr.pearson.vec<-c()
  for (i in 1:dim(x)[1]){
    corr.pearson.vec<-c(corr.pearson.vec,x[i,])
  }
  return(corr.pearson.vec)
}


























data.var.nozero<-dat.var[,which(apply(dat.var,2,function(x){length(which(x==0))})!=0)]
any(dat.var[,1]==0)

length(which(zerocount.by.var


dat.mat<-cbind(dat.resp,dat.pairid,dat.var)
dat.tot<-as.data.frame(dat.mat)
varnames<-paste(rep("Var",472),as.character(1:472),sep="")
varnames[1:20]
colnames(dat.tot)<-c("Response","Pair_id",varnames)
head(dat.tot)[,1:20]

cor.mat<-cor(dat.var, use="complete.obs") 
#Correlation matrix of 472 log transformed covariates
cor.mat.log<-cor(dat.var.log,use="complete.obs") 
#transform the correlation matrix to the vector with the length-p(p-1)/2
corr.pearson.log.vec<-c()
from<-c()
to<-c()
for (i in 1:(dim(cor.mat.log)[1]-1)){
  from<-c(from,rep(paste("V",as.character(i),sep=""),length((i+1):dim(cor.mat.log)[1])))
  if(length(cor.mat.log[-c(1:i),i])!= 1){
    to<-c(to,names(cor.mat.log[-c(1:i),i]))}
  if(length(cor.mat.log[-c(1:i),i])== 1){
    to<-c(to,names(cor.mat.log[i:(i+1),i])[2])}
  corr.pearson.log.vec<-c(corr.pearson.log.vec,cor.mat.log[-c(1:i),i])
  
}
#fisher's Z transform
fisher.z.corr.log<-0.5*log((1+corr.pearson.log.vec)/(1-corr.pearson.log.vec))
hist(fisher.z.corr.log,density=200, breaks=200,prob=TRUE, ,col="skyblue")
m.f.l<-mean(fisher.z.corr.log)
std.f.l<-sqrt(var(fisher.z.corr.log))
curve(dnorm(x, mean=m.f.l, sd=std.f.l), 
      col="red", lwd=2, add=TRUE, yaxt="n")


fisher.z.corr.log[which(fisher.z.corr.log>3)]

dim(cor.mat)[1]
#the function-trans.mat.vec,transform the pXp correlation matrix to the vector with the length-p(p-1)/2
#input:x-correlation matrix
#output:the vecotor of correlation among p variables with length-p(p-1)/2
trans.mat.vec<-function(x)
corr.pearson.vec<-c()
from<-c()
to<-c()
for (i in 1:(dim(x)[1]-1)){
  from<-c(from,rep(paste("V",as.character(i),sep=""),length((i+1):dim(x)[1])))
  if(length(x[-c(1:i),i])!= 1){
  to<-c(to,names(x[-c(1:i),i]))}
  if(length(x[-c(1:i),i])== 1){
    to<-c(to,names(x[i:(i+1),i])[2])}
  corr.pearson.vec<-c(corr.pearson.vec,x[-c(1:i),i])
  return(corr.pearson.vec)
 
}
#fisher's Z transform
fisher.z.corr<-0.5*log((1+corr.pearson.vec)/(1-corr.pearson.vec))
hist(fisher.z.corr,density=200, breaks=200,prob=TRUE, ,col="skyblue")
m.f<-mean(fisher.z.corr)
std.f<-sqrt(var(fisher.z.corr))
curve(dnorm(x, mean=m.f, sd=std.f), 
      col="red", lwd=2, add=TRUE, yaxt="n")


#test statistic based on fisher's Z transform
tstat.fisher.z.corr<-sqrt(n-3)*abs(fisher.z.corr.log)

tstat.fisher.z.corr1<-sqrt(n-3)*fisher.z.corr
hist(tstat.fisher.z.corr1,density=200, breaks=200,prob=TRUE, col="skyblue")
m<-mean(tstat.fisher.z.corr1)
std<-sqrt(var(tstat.fisher.z.corr1))
curve(dnorm(x, mean=m, sd=std), 
      col="red", lwd=2, add=TRUE, yaxt="n")

#P value based on test statistic of fisher's Z transformation of correlation
pvalue.fisher.z.corr<-1-pnorm(tstat.fisher.z.corr)
length(which(pvalue.fisher.z.corr==0))/length(corr.pearson.vec)

#Robust multiple testing method on Tony TT and Liu 2016 JASA 
kurt.est.vec<-apply(dat.var.log,2,function(x){return(sum((x-mean(x,na.rm=T))^4,na.rm=T)/(sum((x-mean(x,na.rm=T))^2,na.rm=T)^2))})
kurt.est.corr<-n*sum(kurt.est.vec)/(3*p)
corr.pearson.threshold<-(abs(corr.pearson.log.vec)>2*sqrt(log(p)/n))
tstat.kurt.corr<-corr.pearson.threshold*sqrt(n)*abs(corr.pearson.log.vec)/(sqrt(kurt.est.corr*(1-corr.pearson.log.vec^2)))
names(tstat.kurt.corr)<-paste(from,to,sep="-")

length(which(tstat.kurt.corr1<0))
tstat.kurt.corr1<-corr.pearson.threshold*sqrt(n)*corr.pearson.log.vec/(sqrt(kurt.est.corr*(1-corr.pearson.log.vec^2)))
tstat.kurt.corr.threshold1<-tstat.kurt.corr1[tstat.kurt.corr1!=0]
hist(tstat.kurt.corr.threshold1,density=200, breaks=200,prob=TRUE, col="skyblue")
m1<-mean(tstat.kurt.corr.threshold1)
std1<-sqrt(var(tstat.kurt.corr.threshold1))
curve(dnorm(x, mean=m1, sd=std1), 
      col="red", lwd=2, add=TRUE, yaxt="n")

? mean

tstat.kurt.corr.threshold<-tstat.kurt.corr[tstat.kurt.corr!=0]
length(tstat.kurt.corr.threshold)/length(tstat.kurt.corr)
p.value.tstat.kurt.corr<-1-pnorm(tstat.kurt.corr.threshold)
length(which(corr.pearson.vec>0.9))
length(which(corr.pearson.vec>0.8))
0.5*log((1+0.9)/(1-0.9))
0.5*log((1+0.8)/(1-0.8))
sqrt(n-3)
corr.pearson.vec[which(corr.pearson.vec>0.9)]

length(which(tstat.kurt.corr!=0))
pnorm(2)

length(corr.pearson.vec)


length(corr.pearson.vec)==p*(p-1)/2
length(from)==p*(p-1)/2
length(to)==p*(p-1)/2

to[111152]
i<-dim(cor.mat)[1]-1
rep(paste("V",as.character(i),sep=""),length(c(i+1,dim(cor.mat)[1]))
rep("v1",1)

rep("v1",5)
paste("V",as.character(i),sep="")
names(corr.pearson.vec)[1:10]
dat.cathgen<-data.frame(resp=dat.resp,matchpairid=dat.pairid,var1
~472=dat.var)


cor.mat[-c(1:5),1]

dim(dat.cathgen)
clogit(case ~ spontaneous + induced + strata(stratum), data=infert)## End(Not run)
# Correlations with significance levels
library(Hmisc)
rcorr(x, type="pearson") # type can be pearson or spearman

#mtcars is a data frame 
rcorr(as.matrix(mtcars))

x<-rnorm(10000)
y<-(exp(2*x)-1)/(exp(2*x)+1)
mean(y)
ynew<-sample(y)
mean(ynew)

.libPaths();

index = sample(y, size=length(y)/2, replace=F)
median(index)

qbeta(0.025,11,91)
qbeta(0.975,11,91)
qnorm(0.975)

#test paired data permutation 

#simulate the dependent observations based on bivariate normal distribution with correlation rho

dat.dep<-function(nsim,rho){
  library(MASS)
  n.simm<-nsim/2
  Sigma <- matrix(c(1,rho,rho,1),2,2)
dat<-mvrnorm(n = n.simm, rep(0, 2), Sigma)
dat.obs<-c(dat[,1],dat[,2])
return(dat.obs)}
dat.test<-dat.dep(500,0.5)

c(c(1,2),c(3,4))

pairid

for( i in 1:99){
  dat.test <-cbind(dat.test,dat.dep(500,0.5))  
}
pairid<-rep(1:250,2)
dat.test.all<-data.frame(dat.test,pairid)
colnames(dat.test.all)<-c(paste(rep("V",100),as.character(1:100),sep=""),"Pairid")
dat.test<-data.frame(dat.test)
names(dat.test)<-paste(rep("V",100),as.character(1:100),sep="")
summary(dat.test)
dim(dat.test)
setwd("D:/Research project 2015/Protein network analysis/Real data/CATHGEN/")
source("pval_perm_corr.R")
#p-values on permutation based on block pairid
pval.sim.pair<-pval.perm.corr(dat=dat.test,nsim=100000,bvar=pairid)
write.table(pval.sim.pair,"p_value_test_pair.txt",row.names=T,col.names=F)
#p-values on permutation based on all the obsevations
pval.sim<-pval.perm.corr(dat=dat.test,nsim=100000,bvar=NA)
write.table(pval.sim,"p_value_test_full.txt",row.names=T,col.names=F)
length(pval.sim.pair)
length(pval.sim)
sum(pval.sim<0.05)/length(pval.sim)

sum(pval.sim.pair<0.05)/length(pval.sim)
pval.pair.nopair.perm<-c(0.01010101,0.07959596)




dat.test.0.9<-dat.dep(500,0.9)
for( i in 1:99){
  dat.test.0.9 <-cbind(dat.test.0.9,dat.dep(500,0.9))  
}

dat.test.all.0.9<-data.frame(dat.test.0.9,pairid)

setwd("D:/Research project 2015/Protein network analysis/Real data/CATHGEN/")
source("pval_perm_corr.R")
#p-values on permutation based on block pairid
pval.sim.pair.9<-pval.perm.corr(dat=dat.test.0.9,nsim=100000,bvar=pairid)
write.table(pval.sim.pair.9,"p_value_test_pair_9.txt",row.names=T,col.names=F)
#p-values on permutation based on all the obsevations
pval.sim.9<-pval.perm.corr(dat=dat.test.0.9,nsim=100000,bvar=NA)
write.table(pval.sim.9,"p_value_test_full_9.txt",row.names=T,col.names=F)
length(pval.sim.pair.9)
length(pval.sim.9)
pval.pair.nopair.perm.0.9<-c(sum(pval.sim.pair.9<0.05)/length(pval.sim.pair.9),sum(pval.sim.9<0.05)/length(pval.sim.9))


dat.test.0.1<-dat.dep(500,0.1)
for( i in 1:99){
  dat.test.0.1 <-cbind(dat.test.0.1,dat.dep(500,0.1))  
}

dat.test.all.0.1<-data.frame(dat.test.0.1,pairid)

setwd("D:/Research project 2015/Protein network analysis/Real data/CATHGEN/")
source("pval_perm_corr.R")
#p-values on permutation based on block pairid
pval.sim.pair.1<-pval.perm.corr(dat=dat.test.0.1,nsim=100000,bvar=pairid)
write.table(pval.sim.pair.1,"p_value_test_pair_1.txt",row.names=T,col.names=F)
#p-values on permutation based on all the obsevations
pval.sim.1<-pval.perm.corr(dat=dat.test.0.1,nsim=100000,bvar=NA)
write.table(pval.sim.1,"p_value_test_full_1.txt",row.names=T,col.names=F)
length(pval.sim.pair.1)
length(pval.sim.1)
pval.pair.nopair.perm.0.1<-c(sum(pval.sim.pair.1<0.05)/length(pval.sim.pair.1),sum(pval.sim.1<0.05)/length(pval.sim.1))

pval.perm.test.all<-cbind(pval.pair.nopair.perm.0,pval.pair.nopair.perm.0.1,pval.pair.nopair.perm.0.5,pval.pair.nopair.perm.0.9)
colnames(pval.perm.test.all)<-c("Corr_0","Corr_0.1","Corr_0.5","Corr_0.9")
rownames(pval.perm.test.all)<-c("Perm_paired","Perm_nopair")
write.table(pval.perm.test.all,"p_value_sim_perm_pair.txt",row.names=T,col.names=T)
sum(pval.sim.pair<0.05)/length(pval.sim)
99*50



dat.test.0<-dat.dep(500,0)
for( i in 1:99){
  dat.test.0 <-cbind(dat.test.0,dat.dep(500,0))  
}


#p-values on permutation based on block pairid
pval.sim.pair.0<-pval.perm.corr(dat=dat.test.0,nsim=100000,bvar=pairid)
write.table(pval.sim.pair.0,"p_value_test_pair_0.txt",row.names=T,col.names=F)
#p-values on permutation based on all the obsevations
pval.sim.0<-pval.perm.corr(dat=dat.test.0,nsim=100000,bvar=NA)
write.table(pval.sim.0,"p_value_test_full_0.txt",row.names=T,col.names=F)

pval.pair.nopair.perm.0<-c(sum(pval.sim.pair.0<0.05)/length(pval.sim.pair.0),sum(pval.sim.0<0.05)/length(pval.sim.0))

sum(pval.sim.pair.0<0.05)/length(pval.sim.0)

length(pval.sim.0[pval.sim.0<0.05])/length(pval.sim.0)

all(is.na(c(NA,NA)))==T


? write.table



dim(dat.test)
head(dat.test)
summary(dat.test)
#response dataset
dat.resp<-read.table("ResponseVector.txt")
dim(dat.resp)
#covariates dataset
dat.var<-read.table("Data.txt")
dim(dat.var)
#pair group id
dat.pairid<-read.table("MatchedPairIndicator.txt")
dim(dat.pairid)
length(unique(dat.pairid))
table(dat.pairid)
dat.var[dat.pairid==5,]
length(unique(pairid))
any(is.na(pairid))==F
?rep

if(is.na(dat.pairid)==F && any(is.na(dat.pairid)))
{warning("NA in paired group id")
  break}
if(is.na(dat.pairid)==F && any(is.na(dat.pairid))==F){
  size.b<-length(unique(dat.pairid))
}
dat.pairid<-unlist(dat.pairid)
dat.pairid[[1]]

-dat.var[,1:2]
dat<-dat.var
var.perm<-apply(dat,2,function(x){return(x[rep(sample(unique(dat.pairid),size=size.b),each=n/size.b)])})
dim(var.perm)
apply(dat.pairid,2,function(x){return(length(unique(unlist(x))))})
dat.var[dat.pairid==3,1:2]

sample(unique(bvar),size=size.b)

library(igraph)
## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))


g <- graph_from_data_frame(relations, directed=F, vertices=actors)
V(g)$age
get.vertex.attribute(g,"age")
print(g, e=TRUE, v=TRUE)
is(g,"graphNEL")




g1 <- graph.ring(10)
is(g1,"igraph")
g1 <- set.graph.attribute(g, "name", "RING")
get.graph.attribute(g1,"name")
# It is the same as
g1$name <- "dfasf"
g$name











node<-data.frame(name=colnames(dat.var),pval=pval.node)
edge<-data.frame(from=str)

## The opposite operation
as_data_frame(g, what="vertices")
as_data_frame(g, what="edges")
datt<-c(50.1,70.1,137,166.9,170.5,152.8,80.5,123.5,112.6,148.5,160,125.4)
min(datt)*length(datt)
mean(datt)


setwd("D:/Research project 2015/Protein network analysis/Real data/CATHGEN/")
source("edge_score_func.R")

p.val.corr.perm<-read.table("p_value_test_full_9.txt",row.names=1,header=T)

p.val.corr.perm.order<-p.val.corr.perm[order(p.val.corr.perm),1]
names(p.val.corr.perm.order)<-rownames(p.val.corr.perm)

edge.score<-uniform.beta.edge.score(pval=p.val.corr.perm.order,fdr=0.01)

edge.score.order.name<-edge.score[order(names(edge.score))]

edge.score.order.name[1:20]


max(edge.score)
min(edge.score)
length(which(edge.score>0))
length(which(edge.score<0))



typeof(unique(from.name))
unique((to.name))
unique(names(dat.test))

length(from.name)

from.name[order(from.name)[1]]==name[order(name)[1]]
dim(dat.test)
dim(dat.var)
names(dat.test)
pval.sim.test<-pval.perm.corr(dat=dat.test,nsim=10,bvar=NA)
pval.sim.test[1:10]
names(pval.sim.test)
g<-induced.graph.data.frame(dat=dat.test,edge.score=pval.sim.test,node.weight=NA,edge.weight=NA)
print(g)

from.name<-c()
for(i in 1:length(colnames(dat))){
  from.name<-c(from.name,rep(order(names(dat))[i],(length(names(dat))-i)))
}
to.name<-c()
for(j in 1:(length(colnames(dat))-1)){
  to.name<-c(to.name,order(names(dat))[(j+1):length(names(dat))])
}
? read.table


#Construct the test data for the new adapted RunFastHeitz_adapted function
setwd("D:/Research project 2015/Protein network analysis/Real data/CATHGEN/")
source("pval_perm_corr.R")

#response dataset
dat.resp<-read.table("ResponseVector.txt")
dim(dat.resp)
#covariates dataset
dat.var<-read.table("Data.txt")
dim(dat.var)
#pair group id
dat.pairid<-read.table("MatchedPairIndicator.txt")
dim(dat.pairid)
length(unique(dat.pairid))
table(dat.pairid)










pval.edge<-read.table("pvalue_Permutation.txt",head=T)
pval.edge.1<-pval.edge$pval
names(pval.edge.1)<-pval.edge$name
source("edge_score_func.R")
edge.score<-uniform.beta.edge.score(pval=pval.edge.1,fdr=0.01)
#P-value from conditional likelihood for paired case-control study
#Contruct the node scores from the functions in BioNet libraryd
dat.obs<-data.frame(dat.resp,dat.pairid,dat.var)
colnames(dat.obs)<-c("Response","Pair_id",colnames(dat.var))
library(survival)
pval.node<-c()
for(i in 3:dim(dat.obs)[2]){
  logit.v<-clogit(Response ~dat.obs[,i]  + strata(Pair_id), dat.obs)  
  pval.node<-c(pval.node,summary(logit.v)$logtest["pvalue"])
}
names(pval.node)<-colnames(dat.obs)[3:dim(dat.obs)[2]]
#contruct the test network with the format-igraph based on the dataset of the covariates and edge scores

da.igraph.test<-induced.graph.data.frame(dat=dat.var,edge.score=edge.scores,node.weight=NA,edge.weight=NA)
min(E(da.igraph.test)$score)
min(get.edge.attribute(da.igraph.test,"score"))
list.edge.attributes(da.igraph.test)
#get the node score based on the network and p.values of the nodes in the network
node.scores<-node.score(da.igraph=da.igraph.test,pval=pval.node,fdr=0.20)
list.test<-list(da.igraph=da.igraph.test,node.scores=node.scores,edge.scores=edge.score)
length(edge.scores[edge.scores>0])/length(edge.scores)
length(node.scores[node.scores>0])

max(edge_attr(da.igraph.test,"score"))

#if Y=aX+(1-a)N(0,1),X~N(0,1),corr(X,Y)=rho,a=a.rho(rho),0<=a<1,-1<rho<1
a.rho<-function(rho){
  if(rho==(sqrt(2)/2)){
    a=2
  }
  if(rho==0){
    a=0
  }
  if(rho!=(sqrt(2)/2) && rho>0){
 a=(rho^2-rho*sqrt(1-rho^2))/(2*rho^2-1)
  }
  if(rho!=(sqrt(2)/2) && rho<0){
    a=(rho^2+rho*sqrt(1-rho^2))/(2*rho^2-1)
  }
 if(!(a<=1 && a>=0)){
   warning("a should be betweeen 0 and 1")
 } 
return(a)}
x<-seq(-1,1,by=0.001)
a.r<-c()
for (a in seq(-1,1,by=0.001))
{a.r<-c(a.r,a.rho(a))}
plot(x,a.r,col=2)
abline(h=1,col=4)
a.rho(-0.1)