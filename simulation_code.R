#Simulation setting 1
#detect consistency of the optimization method to detect functional module
#Assuming p-values uniformly distributed between 0 and 1
nsim<-1000
cluster.size.sim<-c()
for (i in 1:nsim){
  pval.sim<-runif(vcount(network))
  names(pval.sim)<-V(network)$name
  node.scores<-node.score(da.igraph=da.igraph.test,pval=pval.sim,fdr=0.2)
  
  cluster.size.sim<-c(cluster.size.sim,vcount(runFastHeinz.e(network=network, node.scores=node.scores,edge.scores=edge.scores)))}
percentage.pos.module<-sum(cluster.size.sim>0)/length(cluster.size.sim)

#p-percentage of signal component distribution 
#a-parameter of signal component distribution-B(a,1)
p<-0.5
x<-rbinom(vcount(network),1,p)
pval.sim.pos<-runif(vcount(network))[x==1]+rbeta(vcount(network),a,1)[x==0]
names(pval.sim.pos)<-V(network)$name

#test statistic-
#nodes of signal normal distribution with N(mu,1) mu=1,2,3
#nodes of noise- standard normal distribution-N(0,1)
#p-the proportion of signal distribution
p<-0.1
mu<-2
x<-rbinom(vcount(network),1,p)
pval.sim.pos<-rnorm(vcount(network),mean=mu,1)[x==1]+rnorm(vcount(network),a,1)[x==0]
#expected number of significant nodes
length(x==1)
#
cluster.size.sim<-c(cluster.size.sim,vcount(runFastHeinz.e(network=network, node.scores=node.scores,edge.scores=edge.scores)))}



#create random graph by ER graph or degree preference random graph
vcount(network)
ecount(network)
range(E(network)$score)
range(V(network)$score)



n.size=20
p=0.4


#generate the random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of positive scoring cluster3,node scores in positive scoring clusters-uniform(0.1,3)
#total 3 positive scoring clusters-within first one edge scores all positive-uniform(0.1,3)
#within another positive scoring cluster-edge scores all negative-uniform(-3,-0.1)
#within last positive scoring cluster-edge scores-uniform(-3,3)
#node scores in other nodes in the output random graph-uniform(-3,3)
#edge scores in other edges in the output random graph-uniform(-3,3)
rand.pos.network<-function(network,p,n.size){
library(igraph)
library(BioNet)
g <- erdos.renyi.game(vcount(network), p)
V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
V(g)$score<-runif(vcount(g),min=-2,max=2)
E(g)$score<-runif(ecount(g),min=-4,max=4)
E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
names(E(g)$score)<-E(g)$name
from<-get.edgelist(g,names=T)[,1]
to<-get.edgelist(g,names=T)[,2]
#construct the positive scoring cluster using multivariate normal distribution
#n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
#p, the percentage of connecting edges 
library(MASS)
ind.cluster1<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(1,n.size,1))))==2)))
new.edgelist.cluster1<-get.edgelist(g,names=T)[ind.cluster1==1,]
edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
n.edge1=length(edge.name.cluster1)
score.edge1<-runif(n.edge1,min=0.1,max=4)
names(score.edge1)<-edge.name.cluster1
score.node.cluster1<-runif(n.size,min=0.1,max=2)
names(score.node.cluster1)<-V(g)$name[1:n.size]
V(g)$score[1:n.size]<-score.node.cluster1
E(g)$score[ind.cluster1==1]<-score.edge1
ind.cluster3<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
new.edgelist.cluster3<-get.edgelist(g,names=T)[ind.cluster3==1,]
edge.name.cluster3<-paste(new.edgelist.cluster3[,1],new.edgelist.cluster3[,2],sep="_")
n.edge3=length(edge.name.cluster3)
score.edge3<-runif(n.edge3,min=-4,max=-0.1)
names(score.edge3)<-edge.name.cluster3
score.node.cluster3<-runif(n.size,min=0.1,max=2)
names(score.node.cluster3)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]

V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster3
E(g)$score[ind.cluster3==1]<-score.edge3

n.start.cluster2<-floor(vcount(g)/2)
ind.cluster2<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(n.start.cluster2,n.start.cluster2+n.size-1,1))))==2)))
new.edgelist.cluster2<-get.edgelist(g,names=T)[ind.cluster2==1,]
edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
n.edge2=length(edge.name.cluster2)
score.edge2<-runif(n.edge2,min=-4,max=4)
names(score.edge2)<-edge.name.cluster2
score.node.cluster2<-runif(n.size,min=0.1,max=2)
names(score.node.cluster2)<-V(g)$name[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]
V(g)$score[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]<-score.node.cluster2
E(g)$score[ind.cluster2==1]<-score.edge2
return(g)}















#generate the random erdos renyi graph model using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of positive scoring cluster3,node scores in positive scoring clusters-uniform(0.1,2)
#total 3 isolated positive scoring clusters-within first one edge scores all positive-uniform(0.1,4)
#within another positive scoring cluster-edge scores all negative-uniform(-4,-0.1)
#within last positive scoring cluster-edge scores-uniform(-4,4)
#connecting with other nodes in the random graph with negative edge scores
#node scores in other nodes in the output random graph-uniform(-2,2)
#edge scores in other edges in the output random graph-uniform(-4,4)
rand.pos.network2<-function(network,p,n.size){
  library(igraph)
  library(BioNet)
  g <- erdos.renyi.game(vcount(network), p)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-2,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  #the first isolated positive scoring cluster
  g.cluster1<-erdos.renyi.game(n.size, p)
  V(g.cluster1)$name<-paste("Var",as.character(1:vcount(g.cluster1)),sep="")

  new.edgelist.cluster1<-get.edgelist(g.cluster1,names=T)
  edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
  n.edge1=  ecount(g.cluster1)
  score.edge1<-runif(n.edge1,min=0.1,max=4)
  names(score.edge1)<-edge.name.cluster1
  score.node.cluster1<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster1)<-V(g)$name[1:n.size]
  V(g.cluster1)$score<-score.node.cluster1
  E(g.cluster1)$score<-score.edge1
  E(g.cluster1)$name<- edge.name.cluster1
  #the third isolated positive scoring cluster
  g.cluster3<-erdos.renyi.game(n.size, p)
  V(g.cluster3)$name<-paste("Var",as.character(seq((vcount(g)-n.size+1),vcount(g),1)),sep="")
  
  new.edgelist.cluster3<-get.edgelist(g.cluster3,names=T)
  edge.name.cluster3<-paste(new.edgelist.cluster3[,1],new.edgelist.cluster3[,2],sep="_")
  n.edge3=  ecount(g.cluster3)
  score.edge3<-runif(n.edge3,min=-4,max=-0.1)
  names(score.edge3)<-edge.name.cluster3
  score.node.cluster3<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster3)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g.cluster3)$score<-score.node.cluster3
  E(g.cluster3)$score<-score.edge3
  E(g.cluster3)$name<- edge.name.cluster3
  #the second isolated positive scoring cluster
  g.cluster2<-erdos.renyi.game(n.size, p)
  n.start.cluster2<-floor(vcount(g)/2)
  V(g.cluster2)$name<-paste("Var",as.character(seq(n.start.cluster2,n.start.cluster2+n.size-1,1)),sep="")
  
  new.edgelist.cluster2<-get.edgelist(g.cluster2,names=T)
  edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
  n.edge2=  ecount(g.cluster2)
  score.edge2<-runif(n.edge2,min=-4,max=4)
  names(score.edge2)<-edge.name.cluster2
  score.node.cluster2<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster2)<-V(g)$name[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]
  V(g.cluster2)$score<-score.node.cluster2
  E(g.cluster2)$score<-score.edge2
  E(g.cluster2)$name<- edge.name.cluster2 
  

  
nod.num.other<-c(seq((n.size+1),(n.start.cluster2-1),1),seq((n.start.cluster2+n.size),(vcount(g)-n.size),1))
n.p=floor(p*length(nod.num.other)*n.size)
from.to.tot.2<-cbind(rep(seq(n.start.cluster2,n.start.cluster2+n.size-1,1),each=length(nod.num.other)),rep(nod.num.other,n.size))
ind.p.2= sample(1:(length(nod.num.other)*n.size),size=n.p)
from.to.p.2<-from.to.tot.2[ind.p.2,]
from.to.tot.1<-cbind(rep(1:n.size,each=length(nod.num.other)),rep(nod.num.other,n.size))
ind.p.1= sample(1:(length(nod.num.other)*n.size),size=n.p)
from.to.p.1<-from.to.tot.1[ind.p.1,]
from.to.tot.3<-cbind(rep(seq((vcount(g)-n.size+1),vcount(g),1),each=length(nod.num.other)),rep(nod.num.other,n.size))
ind.p.3= sample(1:(length(nod.num.other)*n.size),size=n.p)
from.to.p.3<-from.to.tot.3[ind.p.1,]
from.to.other.cluster<-rbind(from.to.p.1,from.to.p.2,from.to.p.3)
from.other.cluster<-from.to.other.cluster[,1]
to.other.cluster<-from.to.other.cluster[,2]


edge.score.cluster.other<-runif(length(to.other.cluster),min=-4,max=-2.1)
names(edge.score.cluster.other)<-paste(V(g)$name[from.other.cluster],V(g)$name[to.other.cluster],sep="_")
#conbine three isolated positive scoring clusters and other parts together to form a new graph
 node.score<-c(score.node.cluster1,V(g)$score[seq((n.size+1),(n.start.cluster2-1),1)],
               score.node.cluster2,V(g)$score[seq((n.start.cluster2+n.size),(vcount(g)-n.size),1)],score.node.cluster3)
names(node.score)<-V(g)$name
ind.other<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,nod.num.other))==2))))
new.edgelist.other<-get.edgelist(g,names=T)[ind.other==1,]
edge.name.other<-paste(new.edgelist.other[,1],new.edgelist.other[,2],sep="_")
edge.score.tot<-E(g)$score
names(edge.score.tot)<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
edge.score.other<-edge.score.tot[edge.name.other]
edge.score.tot.new<-c(score.edge1,score.edge2,score.edge3,edge.score.other,edge.score.cluster.other) 
from.new<-unlist(strsplit(names(edge.score.tot.new),"_"))[seq(1,(2*length(edge.score.tot.new)),2)]
to.new<-unlist(strsplit(names(edge.score.tot.new),"_"))[seq(2,(2*length(edge.score.tot.new)),2)]
node <- data.frame(name=names(node.score),weight=NA,score=node.score)

relations <- data.frame(from=from.new,
                        to=to.new,
                        score=edge.score.tot.new,
                        weight=NA)
g.new <- graph_from_data_frame(relations, directed=F, vertices=node)
E(g.new)$name<-paste(from.new,to.new,sep="_")
return(g.new)}














n.size=20
p=0.2
get.edgelist(network.test,name=F)[1:1000,]
E(network.test)$score[1:1000]
sim.test.pos.cluster2<-function(network,n.size,p,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  fnr.e2<-c()
  tpr.n<-c()
  fnr.n1<-c()
  fnr.n2<-c()
  for (i in 1:n.sim){
network.test<-rand.pos.network2(network=network,n.size=n.size,p=p)
node.score<-V(network.test)$score
names(node.score)<-V(network.test)$name

edge.score<-E(network.test)$score
names(edge.score)<-E(network.test)$name


module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
module.n<-runFastHeinz(network=network.test, scores=node.score)

node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
tpr.e=c(tpr.e,sum(node.num.e<=n.size)/n.size)

n.start.cluster2<-floor(vcount(network.test)/2)
fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
)/(n.size))
fnr.e2<-c(fnr.e2,length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/(n.size))


node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
tpr.n=c(tpr.n,sum(node.num.n<=n.size)/n.size)
n.start.cluster2<-floor(vcount(network.test)/2)
fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
            )/(n.size*2))

fnr.n2<-c(fnr.n2,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1))))/(n.size))

}
return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,fnr_edge2=fnr.e2,tpr_node=tpr.n,fnr_node1=fnr.n1,fnr_node2=fnr.n2))}
#Explore the TPR and FNR for 3 isolated clusters in random graph with different cluster size and connection probability
result.sim21<-sim.test.pos.cluster2(network,n.size=20,p=0.3,n.sim=100)

result.sim211<-sim.test.pos.cluster2(network,n.size=20,p=0.3,n.sim=100)


result.sim11<-sim.test.pos.cluster2(network,n.size=10,p=0.3,n.sim=100)
result.sim22<-sim.test.pos.cluster2(network,n.size=20,p=0.1,n.sim=100)

result.sim12<-sim.test.pos.cluster2(network,n.size=10,p=0.1,n.sim=100)

result.sim31<-sim.test.pos.cluster2(network,n.size=50,p=0.3,n.sim=100)
result.sim32<-sim.test.pos.cluster2(network,n.size=50,p=0.1,n.sim=100)

fnr=sum()
to.node.num<-as.numeric(unlist(strsplit(to,"r"))[seq(2,2*length(to),2)])
sim.result.summary.32<-apply(result.sim32,2,mean)
sim.result.summary.31<-apply(result.sim31,2,mean)
sim.result.summary.21<-apply(result.sim21,2,mean)
sim.result.summary.22<-apply(result.sim22,2,mean)
sim.result.summary.12<-apply(result.sim12,2,mean)
sim.result.summary.111<-apply(result.sim11,2,mean)


sim.result.tot<-rbind(sim.result.summary.11,sim.result.summary.12,sim.result.summary.22,sim.result.summary.31,
                      sim.result.summary.32)

rownames(sim.result.tot)<-c("cluster_size_10_p_0.3","cluster_size_10_p_0.1",
                            "cluster_size_20_p_0.1","cluster_size_50_p_0.3",
                            "cluster_size_50_p_0.1")
write.table(sim.result.tot,"sim_result_1000.txt",row.names=T,col.names=T)



V(module.n)$name

V(module.e)$name

plotModule(module.e)
vcount(module.e)
plotModule(module.n)
vcount(module.n)
E(network.test)$score[1:1000]














library(igraph)
vcount(network)
library(BioNet)
writeHeinzEdges(network,file="D:/Research project 2015/Protein network analysis/Real data/CATHGEN/edge.txt" ,edge.scores=0, use.score=T)
writeHeinzNodes(network, file="D:/Research project 2015/Protein network analysis/Real data/CATHGEN/node.txt", node.scores=0, use.score=T)

runHeinz(heinz.folder="D:\\Research project 2015\\Protein network analysis\\Real data\\CATHGEN\\heinz-master\\script\\", 
         heinz.e.file="D:\\Research project 2015\\Protein network analysis\\Real data\\CATHGEN\\edge.txt", 
         heinz.n.file="D:\\Research project 2015\\Protein network analysis\\Real data\\CATHGEN\\node.txt", N=TRUE, E=T, diff=-1, n=1)

runHeinz <- function (heinz.folder = "", heinz.e.file, heinz.n.file, N = TRUE, 
                      E = FALSE, diff = -1, n = 1) 
{
  if (heinz.folder != "") {
    heinz.folder <- file.path(heinz.folder, "", fsep = .Platform$file.sep)
  }
  N = if (N) 
    "True"
  else "False"
  E = if (E) 
    "True"
  else "False"
  command <- paste("python ", heinz.folder, ";", "heinz.py -e ", heinz.e.file, " -n ", heinz.n.file, " -N ", N, " -E ", E, " -d ", diff, " -s ", n, sep = "")
  system2(command)
}

system('python -c "a = \'hello world\' ; print a; import pandas"')
system("python C:\\Users\\Name\\Desktop\\my.py")
system("python -c 'print \\"hello world\\"'")

system('python -c "a = 2 + 2; print a"')



network.test<-rand.pos.network(network,n.size=20,p=0.6)
E(network.test)$name[1:1000]








#assign random node score and edge scores to other remaining nodes and edges 
V(g)$score[c((n.size+1):(n.start.cluster2-1),(n.start.cluster2+n.size):(vcount(g)-n.size))]<-runif((vcount(g)-3*n.size),min=-3,max=3)
name.edge.tot<-paste(get.edgelist(g,names=T)[,1],get.edgelist(g,names=T)[,2],sep="_")

node.num.other<-c((n.size+1):(n.start.cluster2-1),(n.start.cluster2+n.size):(vcount(g)-n.size))

ind.node.other<-intersect(intersect(which(ind.cluster1==0),which(ind.cluster2==0)),which(ind.cluster3==0))
new.edgelist.other<-get.edgelist(g,names=T)[ind.node.other,]
edge.name.other<-paste(new.edgelist.other[,1],new.edgelist.other[,2],sep="_")
E(g)$score[]
sum(ind.cluster3==0)
length(edge.name.cluster2)+length(edge.name.cluster1)+length(edge.name.cluster3)+length(edge.name.other)
length(E(g))
dim(get.edgelist(g,names=T))

6>5 || 0>=0

length(edge.name.cluster2)+length(edge.name.cluster1)+length(edge.name.cluster3)+length(intersect(intersect(which(ind.cluster1==0),which(ind.cluster2==0)),which(ind.cluster3==0)))

edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
intersect(union(get.edgelist(g,names=F)[1,]),seq(1,20,1))==2

? union

from.node.num.cluster.1<-from.node.num[from.node.num<=n.size]
sig.mat<-diag(1,n.size,n.size)
for (i in 2:n.size){
  
  sig.mat[c(1:(i-1)),i]<-sig.mat[i,c(1:(i-1))]<-cor.edge[((i-2)*(i-1)/2+1):(i*(i-1)/2)]
}


vcount(network.test)*(vcount(network.test)-1)/2













sapply(1:n.size, function(i) 
sig.mat[-c(1:19),20]

for (i in 1:n.size){
  for (j in 1:n.size){
    if(j!=i)
  }
}

cor.mat<-function(y,n.size){
  return(sapply(1:n.size,cor.y,y))
}
sigma.mat<-sapply(1:n.size,cor.mat,n.size=n.size)
mu.vec<-runif(n.size,min=0.5,max=4)
cluster.pos1<-mvrnorm(1,mu.vec,sig.mat)
sig.mat[1,20]
library(clusterGeneration)
rcorrmatrix(n.size,alphad=1)

exp(-1)
exp(1)
# VINE METHOD to generate random correlation matrices
#with all partial correlations distributed ~ beta(betaparam,betaparam)
#rescaled to [-1, 1]
betaparam<-10
function<-vineBeta(d, betaparam){
P = matrix(0,d,d)           #storing partial correlations
S = diag(1,d,d)

for (k in 1:(d-1)){
  for (i in (k+1):d){
    P[k,i] = rbeta(1,betaparam,betaparam); # sampling from beta
    P[k,i] = (P[k,i] -0.5)*2   # linearly shifting to [-1, 1]
    p = P[k,i]
   for (l in seq(k-1,1,-1) ) {#converting partial correlation to raw correlation
      p = p * sqrt((1-P[l,i]^2)*(1-P[l,k]^2)) + P[l,i]*P[l,k]}

 S[k,i] = p
 S[i,k] = p
}}
seq(5,1,-1)
#permuting the variables to make the distribution permutation-invariant
permutation = sample(1:d,d)
S = S[permutation, permutation]
return(S)}


1000*999/2

g <- barabasi.game(vcount(network.test),m=80,power=1,directed=F)
degree.distribution(g)
vcount(g)
ecount(g)



barabasi.game(n, power = 1, m = NULL, out.dist = NULL, out.seq = NULL, 
              out.pref = FALSE, zero.appeal = 1, directed = TRUE,
              algorithm = c("psumtree", "psumtree-multiple", "bag"),
              start.graph = NULL)

vcount(network.test)*(vcount(network.test)-1)/2*0.1


g2 <- barabasi.game(vcount(network.test),m=24,power=1,directed=F)
vcount(g2)
ecount(g2)
vcount(network.test)*(vcount(network.test)-1)/2*0.1


g1 <- barabasi.game(vcount(network.test),m=77,power=0.5,directed=F)
vcount(g1)
ecount(g1)
vcount(network.test)*(vcount(network.test)-1)/2*0.3


#generate the random scale-free graphs according to the Barabasi-Albert model 
#using the number of nodes in the network
#m:Numeric constant, the number of edges to add in each time step 
#m=24,approximately genearte the same number of edge connections as expected of connecting probability 0.1 in erdos renyi graph model
#m=77,,approximately genearte the same number of edge connections as expected of connecting probability 0.3 in erdos renyi graph model
#power: The power of the preferential attachment, distribution of degrees of the nodes in Barabasi-Albert model 
#n.size-the size of positive scoring cluster3,node scores in positive scoring clusters-uniform(0.1,2)
#total 3 isolated positive scoring clusters-within first one edge scores all positive-uniform(0.1,4)
#within another positive scoring cluster-edge scores all negative-uniform(-4,-0.1)
#within last positive scoring cluster-edge scores-uniform(-4,4)
#connecting with other nodes in the random graph with negative edge scores
#node scores in other nodes in the output random graph-uniform(-2,2)
#edge scores in other edges in the output random graph-uniform(-4,4)
m=24
power=1
n.size=50
g.test<-randba.pos.network2(network=network,m=24,power=1,n.size=50)
randba.pos.network2<-function(network,m,power,n.size){

  library(igraph)
  library(BioNet)
  g <- barabasi.game(vcount(network),m=m,power=power,directed=F)
  p<-2*ecount(g)/(vcount(g)*(vcount(g)-1))
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-2,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  #the first isolated positive scoring cluster
  g.cluster1<-erdos.renyi.game(n.size, p)
  V(g.cluster1)$name<-paste("Var",as.character(1:vcount(g.cluster1)),sep="")
  
  new.edgelist.cluster1<-get.edgelist(g.cluster1,names=T)
  edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
  n.edge1=  ecount(g.cluster1)
  score.edge1<-runif(n.edge1,min=0.1,max=4)
  names(score.edge1)<-edge.name.cluster1
  score.node.cluster1<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster1)<-V(g)$name[1:n.size]
  V(g.cluster1)$score<-score.node.cluster1
  E(g.cluster1)$score<-score.edge1
  E(g.cluster1)$name<- edge.name.cluster1
  #the third isolated positive scoring cluster
  g.cluster3<-erdos.renyi.game(n.size, p)
  V(g.cluster3)$name<-paste("Var",as.character(seq((vcount(g)-n.size+1),vcount(g),1)),sep="")
  
  new.edgelist.cluster3<-get.edgelist(g.cluster3,names=T)
  edge.name.cluster3<-paste(new.edgelist.cluster3[,1],new.edgelist.cluster3[,2],sep="_")
  n.edge3=  ecount(g.cluster3)
  score.edge3<-runif(n.edge3,min=-4,max=-0.1)
  names(score.edge3)<-edge.name.cluster3
  score.node.cluster3<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster3)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g.cluster3)$score<-score.node.cluster3
  E(g.cluster3)$score<-score.edge3
  E(g.cluster3)$name<- edge.name.cluster3
  #the second isolated positive scoring cluster
  g.cluster2<-erdos.renyi.game(n.size, p)
  n.start.cluster2<-floor(vcount(g)/2)
  V(g.cluster2)$name<-paste("Var",as.character(seq(n.start.cluster2,n.start.cluster2+n.size-1,1)),sep="")
  
  new.edgelist.cluster2<-get.edgelist(g.cluster2,names=T)
  edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
  n.edge2=  ecount(g.cluster2)
  score.edge2<-runif(n.edge2,min=-4,max=4)
  names(score.edge2)<-edge.name.cluster2
  score.node.cluster2<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster2)<-V(g)$name[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]
  V(g.cluster2)$score<-score.node.cluster2
  E(g.cluster2)$score<-score.edge2
  E(g.cluster2)$name<- edge.name.cluster2 
  
  
  
  nod.num.other<-c(seq((n.size+1),(n.start.cluster2-1),1),seq((n.start.cluster2+n.size),(vcount(g)-n.size),1))
  n.p=floor(p*length(nod.num.other)*n.size)
  from.to.tot.2<-cbind(rep(seq(n.start.cluster2,n.start.cluster2+n.size-1,1),each=length(nod.num.other)),rep(nod.num.other,n.size))
  ind.p.2= sample(1:(length(nod.num.other)*n.size),size=n.p)
  from.to.p.2<-from.to.tot.2[ind.p.2,]
  from.to.tot.1<-cbind(rep(1:n.size,each=length(nod.num.other)),rep(nod.num.other,n.size))
  ind.p.1= sample(1:(length(nod.num.other)*n.size),size=n.p)
  from.to.p.1<-from.to.tot.1[ind.p.1,]
  from.to.tot.3<-cbind(rep(seq((vcount(g)-n.size+1),vcount(g),1),each=length(nod.num.other)),rep(nod.num.other,n.size))
  ind.p.3= sample(1:(length(nod.num.other)*n.size),size=n.p)
  from.to.p.3<-from.to.tot.3[ind.p.1,]
  from.to.other.cluster<-rbind(from.to.p.1,from.to.p.2,from.to.p.3)
  from.other.cluster<-from.to.other.cluster[,1]
  to.other.cluster<-from.to.other.cluster[,2]
  
  
  edge.score.cluster.other<-runif(length(to.other.cluster),min=-4,max=-2.1)
  names(edge.score.cluster.other)<-paste(V(g)$name[from.other.cluster],V(g)$name[to.other.cluster],sep="_")
  #conbine three isolated positive scoring clusters and other parts together to form a new graph
  node.score<-c(score.node.cluster1,V(g)$score[seq((n.size+1),(n.start.cluster2-1),1)],
                score.node.cluster2,V(g)$score[seq((n.start.cluster2+n.size),(vcount(g)-n.size),1)],score.node.cluster3)
  names(node.score)<-V(g)$name
  ind.other<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,nod.num.other))==2))))
  new.edgelist.other<-get.edgelist(g,names=T)[ind.other==1,]
  edge.name.other<-paste(new.edgelist.other[,1],new.edgelist.other[,2],sep="_")
  edge.score.tot<-E(g)$score
  names(edge.score.tot)<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  edge.score.other<-edge.score.tot[edge.name.other]
  edge.score.tot.new<-c(score.edge1,score.edge2,score.edge3,edge.score.other,edge.score.cluster.other) 
  from.new<-unlist(strsplit(names(edge.score.tot.new),"_"))[seq(1,(2*length(edge.score.tot.new)),2)]
  to.new<-unlist(strsplit(names(edge.score.tot.new),"_"))[seq(2,(2*length(edge.score.tot.new)),2)]
  node <- data.frame(name=names(node.score),weight=NA,score=node.score)
  
  relations <- data.frame(from=from.new,
                          to=to.new,
                          score=edge.score.tot.new,
                          weight=NA)
  g.new <- graph_from_data_frame(relations, directed=F, vertices=node)
  E(g.new)$name<-paste(from.new,to.new,sep="_")
  return(g.new)}


sim.test.pos.cluster3<-function(network,m,power,n.size,n.sim){

  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  fnr.e2<-c()
  tpr.n<-c()
  fnr.n1<-c()
  fnr.n2<-c()
  for (i in 1:n.sim){
    network.test<-randba.pos.network2(network=network,m=m,power=power,n.size=n.size)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
    tpr.e=c(tpr.e,sum(node.num.e<=n.size)/n.size)
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
    fnr.e2<-c(fnr.e2,length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/(n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
    tpr.n=c(tpr.n,sum(node.num.n<=n.size)/n.size)
    n.start.cluster2<-floor(vcount(network.test)/2)
    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size*2))
    
    fnr.n2<-c(fnr.n2,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1))))/(n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,fnr_edge2=fnr.e2,tpr_node=tpr.n,fnr_node1=fnr.n1,fnr_node2=fnr.n2))}


#Explore the TPR and FNR for 3 isolated clusters in random Barabasi-Albert model graph 
#with different cluster size and powers, and number of edges
#number of edges approximately 0.1 connecting probability in Erdos-Renyi model graph,m=24
#number of edges approximately 0.3 connecting probability in Erdos-Renyi model graph,m=77
resultab.sim21<-sim.test.pos.cluster3(network,m=24,power=1,n.size=20,n.sim=100)
resultab.sim22<-sim.test.pos.cluster3(network,m=24,power=3,n.size=20,n.sim=100)
resultab.sim23<-sim.test.pos.cluster3(network,m=24,power=1,n.size=50,n.sim=100)
resultab.sim24<-sim.test.pos.cluster3(network,m=24,power=3,n.size=50,n.sim=100)
resultab.sim25<-sim.test.pos.cluster3(network,m=24,power=1,n.size=10,n.sim=100)
resultab.sim26<-sim.test.pos.cluster3(network,m=24,power=3,n.size=10,n.sim=100)

resultab.sim31<-sim.test.pos.cluster3(network,m=77,power=1,n.size=20,n.sim=100)
resultab.sim32<-sim.test.pos.cluster3(network,m=77,power=3,n.size=20,n.sim=100)
resultab.sim33<-sim.test.pos.cluster3(network,m=77,power=1,n.size=50,n.sim=100)
resultab.sim34<-sim.test.pos.cluster3(network,m=77,power=3,n.size=50,n.sim=100)
resultab.sim35<-sim.test.pos.cluster3(network,m=77,power=1,n.size=10,n.sim=100)
resultab.sim36<-sim.test.pos.cluster3(network,m=77,power=3,n.size=10,n.sim=100)


sim.result.summary.32<-apply(resultab.sim32,2,mean)
sim.result.summary.31<-apply(resultab.sim31,2,mean)
sim.resultab.summary.21<-apply(resultab.sim21,2,mean)
sim.resultab.summary.22<-apply(resultab.sim22,2,mean)
sim.resultab.summary.23<-apply(resultab.sim23,2,mean)
sim.resultab.summary.24<-apply(resultab.sim24,2,mean)
sim.resultab.summary.25<-apply(resultab.sim25,2,mean)
sim.resultab.summary.26<-apply(resultab.sim26,2,mean)

sim.resultab.summary.31<-apply(resultab.sim31,2,mean)
sim.resultab.summary.32<-apply(resultab.sim32,2,mean)
sim.resultab.summary.33<-apply(resultab.sim33,2,mean)
sim.resultab.summary.34<-apply(resultab.sim34,2,mean)
sim.resultab.summary.35<-apply(resultab.sim35,2,mean)
sim.resultab.summary.36<-apply(resultab.sim36,2,mean)






resultab.sim211<-sim.test.pos.cluster2(network,n.size=20,p=0.3,n.sim=100)


resultab.sim11<-sim.test.pos.cluster2(network,n.size=10,p=0.3,n.sim=100)
resultab.sim22<-sim.test.pos.cluster2(network,n.size=20,p=0.1,n.sim=100)

resultab.sim12<-sim.test.pos.cluster2(network,n.size=10,p=0.1,n.sim=100)

resultab.sim31<-sim.test.pos.cluster2(network,n.size=50,p=0.3,n.sim=100)
resultab.sim32<-sim.test.pos.cluster2(network,n.size=50,p=0.1,n.sim=100)

fnr=sum()
to.node.num<-as.numeric(unlist(strsplit(to,"r"))[seq(2,2*length(to),2)])
sim.result.summary.32<-apply(result.sim32,2,mean)
sim.result.summary.31<-apply(result.sim31,2,mean)
sim.result.summary.21<-apply(result.sim21,2,mean)
sim.result.summary.22<-apply(result.sim22,2,mean)
sim.result.summary.12<-apply(result.sim12,2,mean)
sim.result.summary.111<-apply(result.sim11,2,mean)
