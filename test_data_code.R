#Construct the test data for the new adapted RunFastHeitz_adapted function
setwd("~/Documents/Network analysis/R code/")
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
edge.score<-uniform.beta.edge.score(pval=pval.edge.1,fdr=0.1)
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
? clogit

sum(edge.score>-1)/length(edge.score)

sum(edge.score>0)/length(edge.score)

edge.scores<-edge.score[edge.score>-1]
#contruct the test network with the format-igraph based on the dataset of the covariates and edge scores
da.igraph.test<-induced.graph.data.frame(dat=dat.var,node.score=NA,edge.score=edge.scores,node.weight=NA,edge.weight=NA)
#get the node score based on the network and p.values of the nodes in the network
node.scores<-node.score(da.igraph=da.igraph.test,pval=pval.node,fdr=0.1)

network<-induced.graph.data.frame(dat=dat.var,node.score=NA,edge.score=edge.scores,node.weight=NA,edge.weight=NA)
V(network)$score<-node.scores













list.test<-list(da.igraph=da.igraph.test,node.scores=node.scores,edge.scores=edge.scores)

length(E(network)$score)
length(edge.scores)
length(names(edge.scores))
vcount(network)
length(E(network)$name)