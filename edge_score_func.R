uniform.beta.edge.score<-function(pval,fdr){
  #fit mixture of uniform-beta model to p-values of the edges
  #input,parameters,pval-the vector of the p-values of the edges in the graph
  #input,fdr-pre-specified FDR for the multiple testing framework
  #a-the proportion of the noise
  #b-the signal distribution-Beta(b,1)
  #pval-the vector of p-values
  log.like.edge<-function(par){
    a<- par[1]
    b<-par[2]
    a1<-exp(a)/(1+exp(a))
    b1<-exp(b)/(1+exp(b))
    like.edge<-a1+(1-a1)*b1*pval^(b1-1)
    log.like.edge1<-sum(-log(like.edge))
    return(log.like.edge1)}
  logit.inv<-function(x){
    return(exp(x)/(1+exp(x)))
  }
  init<-runif(2,0.1,0.9)
  a.est<-logit.inv(optim(par=init,log.like.edge)$par[2])
  lambda.est<-logit.inv(optim(par=init,log.like.edge)$par[1])
#p(p-1)/2=length(pval)
  p<-sqrt(1+2*length(pval))
  fdr2<-fdr/2
  pi.est<-lambda.est+(1-lambda.est)*a.est
  tau.fdr<-((pi.est-fdr2*lambda.est)/(fdr2*(1-lambda.est)))^(1/(a.est-1))
  edge.score<-(a.est-1)*(log(pval)-log(tau.fdr))
  names(edge.score)<-names(pval)
  return(edge.score)
}

#Contruct the igraph object from the names of the covariates in the studies
#incoporate the node scores as the attributes of the nodes in the network 
#incoporate the edge scores as the attributes of the edges in the network 
#use the function-graph_from_data_frame from the library-igraph
#input:dat-the dataset of the covariates in the study, without any nonmissing names of the covariates
#input:edge.score-the edge scores from the function-edge_score_func
#input:node.weight, the weights or scores for the covariates as the nodes in the network
#input:edge.weight, the weights for the relationships among the covariates in the dataset-dat
#the inputs, dat and edge.score are requried for this function
#other input variables are optional
induced.graph.data.frame<-function(dat,node.score,edge.score,node.weight=NA,edge.weight=NA){
library(igraph)
## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
node <- data.frame(name=colnames(dat)[order(colnames(dat))],weight=node.weight,score=node.score)
  from.name<-unlist(strsplit(names(edge.score),"_"))[seq(1,2*length(names(edge.score)),by=2)]
  to.name<-unlist(strsplit(names(edge.score),"_"))[seq(2,2*length(names(edge.score)),by=2)]
relations <- data.frame(from=from.name,
                        to=to.name,
                        score=edge.score,
                       weight=edge.weight)


g <- graph_from_data_frame(relations, directed=F, vertices=node)

return(g)}

#get node scores from the functions in BioNet library
#input, da.igraph-the graph object with the format-igraph from igraph library
#input:the vector of p-values of the nodes , 
#and the name of the vector must equal to names of nodes in the network
#adjusted node score based on p-values 
#split the overall FDR-fdr based on the number of tests for the nodes and the edges
#p-the number of the nodes in the network
#fdr1 for the edge-fdr1=(p-1)fdr/(p+1)
#fdr2 for the node-fdr2=2fdr/(p+1)
node.score<-function(da.igraph,pval,fdr){
p<-length(pval)
fdr1<-fdr/(2)
library(BioNet)
fb.bm.node<-fitBumModel(pval,plot=F)
score.bm<-scoreNodes(da.igraph,fb=fb.bm.node,fdr=fdr1)
return(score.bm)
}

