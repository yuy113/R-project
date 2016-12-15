#Construct the test data for the new adapted RunFastHeitz_adapted function
setwd("D:/Research project/Protein network analysis/Real data/CATHGEN/")
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
###################################################################################################
pval.edge<-read.table("pvalue_Permutation.txt",head=T)


pval.edge.1<-pval.edge$pval


names(pval.edge.1)<-pval.edge$name
source("edge_score_func.R")
edge.score<-uniform.beta.edge.score(pval=pval.edge.1,fdr=0.1)
######################################################################################################
##generate the node scores for the network test
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

######################################################################################################
##correlation threshold######
#####################################################################################################
#pair group id
library(survival)
setwd("D:/Research project/Protein network analysis/Real data/CATHGEN/")
#response dataset
dat.resp<-read.table("ResponseVector.txt")
dim(dat.resp)
#covariates dataset
dat.var<-read.table("Data.txt")
############################################################################################
##Generate the edges in the network after removing the edges with correlation below threshold-thre

edgecorrthre<-function(dat.var,thre){
n<-dim(dat.var)[1]
p<-dim(dat.var)[2]
corr.true<-cor(dat.var, use="complete.obs")
rownames(corr.true)

mat.edge<-NULL
names.mat.edge<-NULL
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
return(edge.corr.thre)}

edge.corr.thre<-edgecorrthre(dat.var=dat.var,thre=0.1)
edge.corr.thre2<-edgecorrthre(dat.var=dat.var,thre=0.2)
edge.corr.thre3<-edgecorrthre(dat.var=dat.var,thre=0.3)
####################################################################################################################
#generate the network with the edges after correlation thresholding

dat=dat.var
node.score=NA
edge.score=edge.corr.thre3
node.weight=NA
edge.weight=NA
library(igraph)
## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
node <- data.frame(name=colnames(dat)[order(colnames(dat))],weight=node.weight,score=node.score)
from.name<-unlist(strsplit(names(edge.score),"-"))[seq(1,2*length(names(edge.score)),by=2)]
to.name<-unlist(strsplit(names(edge.score),"-"))[seq(2,2*length(names(edge.score)),by=2)]
relations <- data.frame(from=from.name,
                          to=to.name,
                          score=edge.score,
                          weight=edge.weight)
network.thre3 <- graph_from_data_frame(relations, directed=F, vertices=node)
source("edge_score_func.R")
#get the node score based on the network and p.values of the nodes in the network
node.scores<-node.score(da.igraph=network.thre,pval=pval.node,fdr=0.2) 
node.scores3<-node.score(da.igraph=network.thre3,pval=pval.node,fdr=0.2) 
module.thre<-runFastHeinz(network=network.thre, scores=node.scores)
module.thre2<-runFastHeinz(network=network.thre2, scores=node.scores2)
module.thre3<-runFastHeinz(network=network.thre3, scores=node.scores3)
return(module.thre)}

plotModule(module.thre2)
V(module.thre2)
V(module.thre)
V(module.thre3)
node.scores[node.scores>0]











which(is.na(names(edge.corr.thre)))

names(edge.corr.thre)[1:100]

mat.edge<-NA
mat.edge<-c(mat.edge,c(2,3))