#generate the random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of 3 positive scoring cluster,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 3 positive scoring clusters-within first high edge scores all positive-uniform(2,4) 
#with low node scores-uniform(0.1,1)
#within another positive scoring cluster-low edge scores-uniform(0.1,1) and high node scores-uniform(2,4) 
#within last positive scoring cluster-edge scores-low edge scores-uniform(0.1,1) and low node scores-uniform(0.1,1) 
#node scores in other nodes in the output random graph-uniform(-4,-0.1)
#edge scores in other edges in the output random graph-uniform(-4,-0.1)
rand.pos.network3<-function(network,p,n.size){
  library(igraph)
  library(BioNet)
  g <- erdos.renyi.game(vcount(network), p)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  ##################################################################################
  #moderate node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
  new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
  edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
  n.edge4=length(edge.name.cluster4)
  score.edge4<-runif(n.edge4,min=0.1,max=1)
  names(score.edge4)<-edge.name.cluster4
  score.node.cluster4<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster4)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster4
  E(g)$score[ind.cluster4==1]<-score.edge4
  
  #high node score,moderate edge score cluster
  ind.cluster3<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-3*n.size+1),(vcount(g)-2*n.size),1))))==2)))
  new.edgelist.cluster3<-get.edgelist(g,names=T)[ind.cluster3==1,]
  edge.name.cluster3<-paste(new.edgelist.cluster3[,1],new.edgelist.cluster3[,2],sep="_")
  n.edge3=length(edge.name.cluster3)
  score.edge3<-runif(n.edge3,min=0.1,max=1)
  names(score.edge3)<-edge.name.cluster3
  score.node.cluster3<-runif(n.size,min=2,max=4)
  names(score.node.cluster3)<-V(g)$name[seq((vcount(g)-3*n.size+1),(vcount(g)-2*n.size),1)]
  V(g)$score[seq((vcount(g)-3*n.size+1),(vcount(g)-2*n.size),1)]<-score.node.cluster3
  E(g)$score[ind.cluster3==1]<-score.edge3
  ##moderate node score,high edge score cluster
  n.start.cluster2<-floor(vcount(g)/2)
  ind.cluster2<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(n.start.cluster2,n.start.cluster2+n.size-1,1))))==2)))
  new.edgelist.cluster2<-get.edgelist(g,names=T)[ind.cluster2==1,]
  edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
  n.edge2=length(edge.name.cluster2)
  score.edge2<-runif(n.edge2,min=2,max=4)
  names(score.edge2)<-edge.name.cluster2
  score.node.cluster2<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster2)<-V(g)$name[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]
  V(g)$score[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]<-score.node.cluster2
  E(g)$score[ind.cluster2==1]<-score.edge2
  
  
  
  
  
  return(g)}


p=0.3
n.size=50
sim.test.pos.cluster3<-function(network,n.size,p,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  fnr.e2<-c()
  tpr.n<-c()
  fnr.n1<-c()
  fnr.n2<-c()
  for (i in 1:n.sim){
    network.test<-rand.pos.network3(network=network,n.size=n.size,p=p)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    
    
    
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
   
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.e<-c(tpr.e,(length(intersect(node.num.e,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq((vcount(network.test)-3*n.size+1),(vcount(network.test)-2*n.size),1)))/n.size))
    
    fnr.e2<-c(fnr.e2,(length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
  
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.n<-c(tpr.n,(length(intersect(node.num.n,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
    
    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq((vcount(network.test)-3*n.size+1),(vcount(network.test)-2*n.size),1)))/n.size))
    
    
    
    
    fnr.n2<-c(fnr.n2,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,fnr_edge2=fnr.e2,tpr_node=tpr.n,fnr_node1=fnr.n1,fnr_node2=fnr.n2))}
#Explore the TPR and FNR for 2 isolated clusters in random graph with different cluster size and connection probability


result31.sim21<-sim.test.pos.cluster3(network,n.size=50,p=0.3,n.sim=100)



result31.sim22<-sim.test.pos.cluster3(network,n.size=50,p=0.1,n.sim=100)

result31.sim121<-sim.test.pos.cluster3(network,n.size=10,p=0.3,n.sim=100)
result31.sim122<-sim.test.pos.cluster3(network,n.size=10,p=0.1,n.sim=100)
result31.sim221<-sim.test.pos.cluster3(network,n.size=20,p=0.3,n.sim=100)
result31.sim222<-sim.test.pos.cluster3(network,n.size=20,p=0.1,n.sim=100)




sim.result31.summary.21<-apply(result31.sim21,2,mean)
sim.result31.summary.22<-apply(result31.sim22,2,mean)
sim.result31.summary.121<-apply(result31.sim121,2,mean)
sim.result31.summary.122<-apply(result31.sim122,2,mean)
sim.result31.summary.221<-apply(result31.sim221,2,mean)
sim.result31.summary.222<-apply(result31.sim222,2,mean)


sim3.result.tot<-rbind(sim.result31.summary.21,sim.result31.summary.22,sim.result31.summary.121,sim.result31.summary.122,
                        sim.result31.summary.221,sim.result31.summary.222)

rownames(sim3.result.tot)<-c("3cluster_size_50_p_0.3","3cluster_size_50_p_0.1",
                            "3cluster_size_10_p_0.3","3cluster_size_10_p_0.1",
                            "3cluster_size_20_p_0.1")
write.table(sim3.result.tot,"sim_3cluster_result1_100.txt",row.names=T,col.names=T)



########################simulation scenario 3,2 positive scoring clusters#########################################################

#generate the random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of 2 positive scoring cluster,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 2 positive scoring clusters
#within another positive scoring cluster-low edge scores-uniform(0.1,1) and high node scores-uniform(2,4) 
#within last positive scoring cluster-edge scores-low edge scores-uniform(0.1,1) and low node scores-uniform(0.1,1) 
#node scores in other nodes in the output random graph-uniform(-4,-0.1)
#edge scores in other edges in the output random graph-uniform(-4,-0.1)
rand.pos.network21<-function(network,p,n.size){
  library(igraph)
  library(BioNet)
  g <- erdos.renyi.game(vcount(network), p)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  ##################################################################################
  #moderate node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
  new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
  edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
  n.edge4=length(edge.name.cluster4)
  score.edge4<-runif(n.edge4,min=0.1,max=1)
  names(score.edge4)<-edge.name.cluster4
  score.node.cluster4<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster4)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster4
  E(g)$score[ind.cluster4==1]<-score.edge4
  
  #high node score,moderate edge score cluster
  n.start.cluster2<-floor(vcount(g)/2)
  ind.cluster2<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(n.start.cluster2,n.start.cluster2+n.size-1,1))))==2)))
  new.edgelist.cluster2<-get.edgelist(g,names=T)[ind.cluster2==1,]
  edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
  n.edge2=length(edge.name.cluster2)
  score.edge2<-runif(n.edge2,min=0.1,max=1)
  names(score.edge2)<-edge.name.cluster2
  score.node.cluster2<-runif(n.size,min=2,max=4)
  names(score.node.cluster2)<-V(g)$name[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]
  V(g)$score[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]<-score.node.cluster2
  E(g)$score[ind.cluster2==1]<-score.edge2
  return(g)}


p=0.3
n.size=50
sim.test.pos.cluster21<-function(network,n.size,p,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
 
  tpr.n<-c()
  fnr.n1<-c()
 
  for (i in 1:n.sim){
    network.test<-rand.pos.network21(network=network,n.size=n.size,p=p)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
    
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.e<-c(tpr.e,(length(intersect(node.num.e,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
    
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.n<-c(tpr.n,(length(intersect(node.num.n,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,tpr_node=tpr.n,fnr_node1=fnr.n1))}
#Explore the TPR and FNR for 4 isolated clusters in random graph with different cluster size and connection probability


result21.sim21<-sim.test.pos.cluster21(network,n.size=50,p=0.3,n.sim=100)



result21.sim22<-sim.test.pos.cluster21(network,n.size=50,p=0.1,n.sim=100)

result21.sim121<-sim.test.pos.cluster21(network,n.size=10,p=0.3,n.sim=100)
result21.sim122<-sim.test.pos.cluster21(network,n.size=10,p=0.1,n.sim=100)
result21.sim221<-sim.test.pos.cluster21(network,n.size=20,p=0.3,n.sim=100)
result21.sim222<-sim.test.pos.cluster21(network,n.size=20,p=0.1,n.sim=100)




sim.result21.summary.21<-apply(result21.sim21,2,mean)
sim.result21.summary.22<-apply(result21.sim22,2,mean)
sim.result21.summary.121<-apply(result21.sim121,2,mean)
sim.result21.summary.122<-apply(result21.sim122,2,mean)
sim.result21.summary.221<-apply(result21.sim221,2,mean)
sim.result21.summary.222<-apply(result21.sim222,2,mean)


sim21.result.tot<-rbind(sim.result21.summary.21,sim.result21.summary.22,sim.result21.summary.121,sim.result21.summary.122,
                       sim.result21.summary.221,sim.result21.summary.222)

rownames(sim21.result.tot)<-c("3cluster_size_50_p_0.3","3cluster_size_50_p_0.1",
                             "3cluster_size_10_p_0.3","3cluster_size_10_p_0.1",
                             "3cluster_size_20_p_0.1")
write.table(sim21.result.tot,"sim_3cluster_result1_100.txt",row.names=T,col.names=T)



########################simulation scenario 4,2 positive scoring clusters#########################################################

#generate the random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of 2 positive scoring cluster,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 2 positive scoring clusters
#within one positive scoring cluster-low node scores-uniform(0.1,1) and high edge scores-uniform(2,4) 
#within last positive scoring cluster-edge scores-low edge scores-uniform(0.1,1) and low node scores-uniform(0.1,1) 
#node scores in other nodes in the output random graph-uniform(-4,-0.1)
#edge scores in other edges in the output random graph-uniform(-4,-0.1)
rand.pos.network22<-function(network,p,n.size){
  library(igraph)
  library(BioNet)
  g <- erdos.renyi.game(vcount(network), p)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  ##################################################################################
  #moderate node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
  new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
  edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
  n.edge4=length(edge.name.cluster4)
  score.edge4<-runif(n.edge4,min=0.1,max=1)
  names(score.edge4)<-edge.name.cluster4
  score.node.cluster4<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster4)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster4
  E(g)$score[ind.cluster4==1]<-score.edge4
  
  #low node score,high edge score cluster
  n.start.cluster2<-floor(vcount(g)/2)
  ind.cluster2<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(n.start.cluster2,n.start.cluster2+n.size-1,1))))==2)))
  new.edgelist.cluster2<-get.edgelist(g,names=T)[ind.cluster2==1,]
  edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
  n.edge2=length(edge.name.cluster2)
  score.edge2<-runif(n.edge2,min=2,max=4)
  names(score.edge2)<-edge.name.cluster2
  score.node.cluster2<-runif(n.size,min=0.1,max=2)
  names(score.node.cluster2)<-V(g)$name[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]
  V(g)$score[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]<-score.node.cluster2
  E(g)$score[ind.cluster2==1]<-score.edge2
  return(g)}

sim.test.pos.cluster22<-function(network,n.size,p,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  
  tpr.n<-c()
  fnr.n1<-c()
  
  for (i in 1:n.sim){
    network.test<-rand.pos.network22(network=network,n.size=n.size,p=p)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
    
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.e<-c(tpr.e,(length(intersect(node.num.e,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
   
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.n<-c(tpr.n,(length(intersect(node.num.n,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,tpr_node=tpr.n,fnr_node1=fnr.n1))}
#Explore the TPR and FNR for 4 isolated clusters in random graph with different cluster size and connection probability


result22.sim21<-sim.test.pos.cluster22(network,n.size=50,p=0.3,n.sim=100)
result22.sim22<-sim.test.pos.cluster22(network,n.size=50,p=0.1,n.sim=100)
result22.sim121<-sim.test.pos.cluster22(network,n.size=10,p=0.3,n.sim=100)
result22.sim122<-sim.test.pos.cluster22(network,n.size=10,p=0.1,n.sim=100)
result22.sim221<-sim.test.pos.cluster22(network,n.size=20,p=0.3,n.sim=100)
result22.sim222<-sim.test.pos.cluster22(network,n.size=20,p=0.1,n.sim=100)

sim.result22.summary.21<-apply(result22.sim21,2,mean)
sim.result22.summary.22<-apply(result22.sim22,2,mean)
sim.result22.summary.121<-apply(result22.sim121,2,mean)
sim.result22.summary.122<-apply(result22.sim122,2,mean)
sim.result22.summary.221<-apply(result22.sim221,2,mean)
sim.result22.summary.222<-apply(result22.sim222,2,mean)


sim22.result.tot<-rbind(sim.result22.summary.21,sim.result22.summary.22,sim.result22.summary.121,sim.result22.summary.122,
                        sim.result22.summary.221,sim.result22.summary.222)

rownames(sim22.result.tot)<-c("3cluster_size_50_p_0.3","3cluster_size_50_p_0.1",
                              "3cluster_size_10_p_0.3","3cluster_size_10_p_0.1",
                              "3cluster_size_20_p_0.1")
write.table(sim22.result.tot,"sim_3cluster_result1_100.txt",row.names=T,col.names=T)

















#generate the Barabasi-Albert scale free random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of positive scoring cluster3,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 3 positive scoring clusters-within first high edge scores all positive-uniform(2,4) 
#with low node scores-uniform(0.1,1)
#within another positive scoring cluster-low edge scores-uniform(0.1,1) and high node scores-uniform(2,4) 
#within last positive scoring cluster-edge scores-low edge scores-uniform(0.1,1) and low node scores-uniform(0.1,1) 
#node scores in other nodes in the output random graph-uniform(-4,-0.1)
#edge scores in other edges in the output random graph-uniform(-4,-0.1)
network=network
m=24
power=3
n.size=10
randba.pos3.network4<-function(network,m,power,n.size){
  library(igraph)
  library(BioNet)
  g <- barabasi.game(vcount(network),m=m,power=power,directed=F)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  ind.cluster1<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(1,n.size,1))))==2)))
  if(sum(ind.cluster1)>=1){
    new.edgelist.cluster1<-get.edgelist(g,names=T)[ind.cluster1==1,]
    edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
    n.edge1=length(edge.name.cluster1)
    score.edge1<-runif(n.edge1,min=2,max=4)
    names(score.edge1)<-edge.name.cluster1
    E(g)$score[ind.cluster1==1]<-score.edge1}
  
  score.node.cluster1<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster1)<-V(g)$name[1:n.size]
  V(g)$score[1:n.size]<-score.node.cluster1
  
  ##################################################################################
  
  
  #moderate node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
  if(sum(ind.cluster4)>=1){
    
    new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
    edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
    n.edge4=length(edge.name.cluster4)
    score.edge4<-runif(n.edge4,min=0.1,max=1)
    names(score.edge4)<-edge.name.cluster4
    E(g)$score[ind.cluster4==1]<-score.edge4
  }
  score.node.cluster4<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster4)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster4
  
  
  
  
  ##high node score,moderate edge score cluster
  n.start.cluster2<-floor(vcount(g)/2)
  ind.cluster2<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(n.start.cluster2,n.start.cluster2+n.size-1,1))))==2)))
  if(sum(ind.cluster2)>=1){
    new.edgelist.cluster2<-get.edgelist(g,names=T)[ind.cluster2==1,]
    edge.name.cluster2<-paste(new.edgelist.cluster2[,1],new.edgelist.cluster2[,2],sep="_")
    n.edge2=length(edge.name.cluster2)
    score.edge2<-runif(n.edge2,min=0.1,max=1)
    names(score.edge2)<-edge.name.cluster2
    E(g)$score[ind.cluster2==1]<-score.edge2}
  
  score.node.cluster2<-runif(n.size,min=2,max=4)
  names(score.node.cluster2)<-V(g)$name[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]
  V(g)$score[seq(n.start.cluster2,n.start.cluster2+n.size-1,1)]<-score.node.cluster2
  
  return(g)}

length(intersect(c(1,3),1:3))


sim.test.pos3.cluster4<-function(network,n.size,m,power,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  fnr.e2<-c()
  tpr.n<-c()
  fnr.n1<-c()
  fnr.n2<-c()

  for (i in 1:n.sim){
    network.test<-randba.pos3.network4(network=network,m=m,power=power,n.size=n.size)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
  
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.e<-c(tpr.e,(length(intersect(node.num.e,seq(1,n.size,1)))
    )/(n.size))
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq(n.start.cluster2,n.start.cluster2+n.size-1,1)))
    )/(n.size))
   
    fnr.e2<-c(fnr.e2,(length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
 
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.n<-c(tpr.n,(length(intersect(node.num.n,seq(1,n.size,1))))/(n.size))
    
    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq(n.start.cluster2,n.start.cluster2+n.size-1,1))))/(n.size))
    
    fnr.n2<-c(fnr.n2,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,fnr_edge2=fnr.e2,tpr_node=tpr.n,fnr_node1=fnr.n1,fnr_node2=fnr.n2))}
#Explore the TPR and FNR for 4 isolated clusters in random graph with different cluster size and connection probability
resultab3.sim11<-sim.test.pos3.cluster4(network=network,n.size=50,m=24,power=3,n.sim=100)
resultab3.sim12<-sim.test.pos3.cluster4(network=network,n.size=50,m=24,power=1,n.sim=100)
resultab3.sim21<-sim.test.pos3.cluster4(network=network,n.size=50,m=77,power=3,n.sim=100)
resultab3.sim22<-sim.test.pos3.cluster4(network=network,n.size=50,m=77,power=1,n.sim=100)



resultab3.sim211<-sim.test.pos3.cluster4(network=network,n.size=20,m=24,power=3,n.sim=100)
resultab3.sim212<-sim.test.pos3.cluster4(network=network,n.size=20,m=24,power=1,n.sim=100)
resultab3.sim221<-sim.test.pos3.cluster4(network=network,n.size=20,m=77,power=3,n.sim=100)
resultab3.sim222<-sim.test.pos3.cluster4(network=network,n.size=20,m=77,power=1,n.sim=100)




sim.resultab3.summary.21<-apply(resultab3.sim21,2,mean)
sim.resultab3.summary.22<-apply(resultab3.sim22,2,mean)
sim.resultab3.summary.12<-apply(resultab3.sim12,2,mean)
sim.resultab3.summary.11<-apply(resultab3.sim11,2,mean)

sim.resultab3.summary.211<-apply(resultab3.sim211,2,mean)
sim.resultab3.summary.212<-apply(resultab3.sim212,2,mean)
sim.resultab3.summary.221<-apply(resultab3.sim221,2,mean)
sim.resultab3.summary.222<-apply(resultab3.sim222,2,mean)


#####################simulation scenario 5,2 positive scoring clusters############################
#generate the Barabasi-Albert scale free random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of positive scoring cluster3,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 2 positive scoring clusters-within first high node scores all positive-uniform(2,4) 
#with low edge scores-uniform(0.1,1)
#within last positive scoring cluster-edge scores-low edge scores-uniform(0.1,1) and low node scores-uniform(0.1,1) 
#node scores in other nodes in the output random graph-uniform(-4,-0.1)
#edge scores in other edges in the output random graph-uniform(-4,-0.1)
network=network
m=24
power=3
n.size=10
randba.pos21.network4<-function(network,m,power,n.size){
  library(igraph)
  library(BioNet)
  g <- barabasi.game(vcount(network),m=m,power=power,directed=F)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  ind.cluster1<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(1,n.size,1))))==2)))
  if(sum(ind.cluster1)>=1){
    new.edgelist.cluster1<-get.edgelist(g,names=T)[ind.cluster1==1,]
    edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
    n.edge1=length(edge.name.cluster1)
    score.edge1<-runif(n.edge1,min=0.1,max=1)
    names(score.edge1)<-edge.name.cluster1
    E(g)$score[ind.cluster1==1]<-score.edge1}
  
  score.node.cluster1<-runif(n.size,min=2,max=4)
  names(score.node.cluster1)<-V(g)$name[1:n.size]
  V(g)$score[1:n.size]<-score.node.cluster1
  
  ##################################################################################
  
  
  #moderate node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
  if(sum(ind.cluster4)>=1){
    
    new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
    edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
    n.edge4=length(edge.name.cluster4)
    score.edge4<-runif(n.edge4,min=0.1,max=1)
    names(score.edge4)<-edge.name.cluster4
    E(g)$score[ind.cluster4==1]<-score.edge4
  }
  score.node.cluster4<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster4)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster4
  
  


  
  return(g)}

length(intersect(c(1,3),1:3))


a=c(1,5)

b=1:4

any(a %in% b)

sim.test.pos21.cluster4<-function(network,n.size,m,power,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  tpr.n<-c()
  fnr.n1<-c()
  
  for (i in 1:n.sim){
    network.test<-randba.pos21.network4(network=network,m=m,power=power,n.size=n.size)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
    
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.e<-c(tpr.e,(length(intersect(node.num.e,seq(1,n.size,1)))
    )/(n.size))
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.n<-c(tpr.n,(length(intersect(node.num.n,seq(1,n.size,1)))
    )/(n.size))

    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,tpr_node=tpr.n,fnr_node1=fnr.n1))}
#Explore the TPR and FNR for 4 isolated clusters in random graph with different cluster size and connection probability
resultab21.sim11<-sim.test.pos21.cluster4(network=network,n.size=50,m=24,power=3,n.sim=100)
resultab21.sim12<-sim.test.pos21.cluster4(network=network,n.size=50,m=24,power=1,n.sim=100)
resultab21.sim21<-sim.test.pos21.cluster4(network=network,n.size=50,m=77,power=3,n.sim=100)
resultab21.sim22<-sim.test.pos21.cluster4(network=network,n.size=50,m=77,power=1,n.sim=100)


resultab21.sim13<-sim.test.pos21.cluster4(network=network,n.size=50,m=77,power=2,n.sim=100)
resultab21.sim23<-sim.test.pos21.cluster4(network=network,n.size=50,m=24,power=2,n.sim=100)
resultab21.sim223<-sim.test.pos21.cluster4(network=network,n.size=20,m=77,power=2,n.sim=100)
resultab21.sim213<-sim.test.pos21.cluster4(network=network,n.size=20,m=24,power=2,n.sim=100)
resultab21.sim214<-sim.test.pos21.cluster4(network=network,n.size=20,m=24,power=1.5,n.sim=100)

resultab21.sim214<-sim.test.pos21.cluster4(network=network,n.size=20,m=24,power=1.5,n.sim=2)
resultab21.sim215<-sim.test.pos21.cluster4(network=network,n.size=20,m=77,power=1.5,n.sim=100)

resultab21.sim211<-sim.test.pos21.cluster4(network=network,n.size=20,m=24,power=3,n.sim=100)
resultab21.sim212<-sim.test.pos21.cluster4(network=network,n.size=20,m=24,power=1,n.sim=100)
resultab21.sim221<-sim.test.pos21.cluster4(network=network,n.size=20,m=77,power=3,n.sim=100)
resultab21.sim222<-sim.test.pos21.cluster4(network=network,n.size=20,m=77,power=1,n.sim=100)


sim.resultab21.summary.23<-apply(resultab21.sim23,2,mean)
sim.resultab21.summary.13<-apply(resultab21.sim13,2,mean)

sim.resultab21.summary.21<-apply(resultab21.sim21,2,mean)
sim.resultab21.summary.22<-apply(resultab21.sim22,2,mean)
sim.resultab21.summary.12<-apply(resultab21.sim12,2,mean)
sim.resultab21.summary.11<-apply(resultab21.sim11,2,mean)

sim.resultab21.summary.211<-apply(resultab21.sim211,2,mean)
sim.resultab21.summary.212<-apply(resultab21.sim212,2,mean)
sim.resultab21.summary.221<-apply(resultab21.sim221,2,mean)
sim.resultab21.summary.222<-apply(resultab21.sim222,2,mean)


#####################simulation scenario 6,2 positive scoring clusters############################
#generate the Barabasi-Albert scale free random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of positive scoring cluster3,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 2 positive scoring clusters-within first high edge scores all positive-uniform(2,4) 
#with low node scores-uniform(0.1,1)
#within last positive scoring cluster-edge scores-low edge scores-uniform(0.1,1) and low node scores-uniform(0.1,1) 
#node scores in other nodes in the output random graph-uniform(-4,-0.1)
#edge scores in other edges in the output random graph-uniform(-4,-0.1)
network=network
m=24
power=3
n.size=10
randba.pos22.network4<-function(network,m,power,n.size){
  library(igraph)
  library(BioNet)
  g <- barabasi.game(vcount(network),m=m,power=power,directed=F)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  ind.cluster1<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(1,n.size,1))))==2)))
  if(sum(ind.cluster1)>=1){
    new.edgelist.cluster1<-get.edgelist(g,names=T)[ind.cluster1==1,]
    edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
    n.edge1=length(edge.name.cluster1)
    score.edge1<-runif(n.edge1,min=2,max=4)
    names(score.edge1)<-edge.name.cluster1
    E(g)$score[ind.cluster1==1]<-score.edge1}
  
  score.node.cluster1<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster1)<-V(g)$name[1:n.size]
  V(g)$score[1:n.size]<-score.node.cluster1
  
  ##################################################################################
  
  
  #large node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
  if(sum(ind.cluster4)>=1){
    
    new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
    edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
    n.edge4=length(edge.name.cluster4)
    score.edge4<-runif(n.edge4,min=0.1,max=1)
    names(score.edge4)<-edge.name.cluster4
    E(g)$score[ind.cluster4==1]<-score.edge4
  }
  score.node.cluster4<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster4)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster4
  
  return(g)}

length(intersect(c(1,3),1:3))


sim.test.pos22.cluster4<-function(network,n.size,m,power,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  tpr.n<-c()
  fnr.n1<-c()
  
  for (i in 1:n.sim){
    network.test<-randba.pos22.network4(network=network,m=m,power=power,n.size=n.size)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
    
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.e<-c(tpr.e,(length(intersect(node.num.e,seq(1,n.size,1)))
    )/(n.size))
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.n<-c(tpr.n,(length(intersect(node.num.n,seq(1,n.size,1)))
    )/(n.size))
    
    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,tpr_node=tpr.n,fnr_node1=fnr.n1))}
#Explore the TPR and FNR for 4 isolated clusters in random graph with different cluster size and connection probability
resultab22.sim11<-sim.test.pos22.cluster4(network=network,n.size=50,m=24,power=3,n.sim=100)
resultab22.sim12<-sim.test.pos22.cluster4(network=network,n.size=50,m=24,power=1,n.sim=100)
resultab22.sim21<-sim.test.pos22.cluster4(network=network,n.size=50,m=77,power=3,n.sim=100)
resultab22.sim22<-sim.test.pos22.cluster4(network=network,n.size=50,m=77,power=1,n.sim=100)

resultab22.sim13<-sim.test.pos22.cluster4(network=network,n.size=50,m=77,power=2,n.sim=100)
resultab22.sim23<-sim.test.pos22.cluster4(network=network,n.size=50,m=24,power=2,n.sim=100)
resultab22.sim223<-sim.test.pos22.cluster4(network=network,n.size=20,m=77,power=2,n.sim=100)
resultab22.sim213<-sim.test.pos22.cluster4(network=network,n.size=20,m=24,power=2,n.sim=100)
resultab22.sim214<-sim.test.pos22.cluster4(network=network,n.size=20,m=24,power=1.5,n.sim=100)
resultab22.sim211<-sim.test.pos22.cluster4(network=network,n.size=20,m=24,power=3,n.sim=100)
resultab22.sim212<-sim.test.pos22.cluster4(network=network,n.size=20,m=24,power=1,n.sim=100)
resultab22.sim221<-sim.test.pos22.cluster4(network=network,n.size=20,m=77,power=3,n.sim=100)
resultab22.sim222<-sim.test.pos22.cluster4(network=network,n.size=20,m=77,power=1,n.sim=100)

sim.resultab22.summary.223<-apply(resultab22.sim223,2,mean)
sim.resultab22.summary.23<-apply(resultab22.sim23,2,mean)
sim.resultab22.summary.13<-apply(resultab22.sim13,2,mean)

sim.resultab22.summary.21<-apply(resultab22.sim21,2,mean)
sim.resultab22.summary.22<-apply(resultab22.sim22,2,mean)
sim.resultab22.summary.12<-apply(resultab22.sim12,2,mean)
sim.resultab22.summary.11<-apply(resultab22.sim11,2,mean)

sim.resultab22.summary.211<-apply(resultab22.sim211,2,mean)
sim.resultab22.summary.212<-apply(resultab22.sim212,2,mean)
sim.resultab22.summary.221<-apply(resultab22.sim221,2,mean)
sim.resultab22.summary.222<-apply(resultab22.sim222,2,mean)



sim.resultab22.summary.13<-apply(resultab22.sim13,2,mean)
sim.resultab22.summary.23<-apply(resultab22.sim23,2,mean)

sim.resultab22.summary.223<-apply(resultab22.sim223,2,mean)


#####################simulation scenario 7,2 positive scoring clusters############################
#generate the Barabasi-Albert scale free random graph using the number of nodes in the network
#p:the probability of connecting edges
#n.size-the size of positive scoring cluster3,node scores in positive scoring clusters-uniform(0.1,1) or uniform(2,4)
#total 2 positive scoring clusters-within first high edge scores all positive-uniform(2,4) 
#with low node scores-uniform(0.1,1)
#within last positive scoring cluster-edge scores-low edge scores-uniform(0.1,1) and high node scores-uniform(2,4) 
#node scores in other nodes in the output random graph-uniform(-4,-0.1)
#edge scores in other edges in the output random graph-uniform(-4,-0.1)
network=network
m=24
power=3
n.size=10
randba.pos22.network5<-function(network,m,power,n.size){
  library(igraph)
  library(BioNet)
  g <- barabasi.game(vcount(network),m=m,power=power,directed=F)
  V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
  V(g)$score<-runif(vcount(g),min=-4,max=-0.1)
  E(g)$score<-runif(ecount(g),min=-4,max=-0.1)
  E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
  names(E(g)$score)<-E(g)$name
  from<-get.edgelist(g,names=T)[,1]
  to<-get.edgelist(g,names=T)[,2]
  #construct the positive scoring cluster using multivariate normal distribution
  #n.size, the size of the positive scoring cluter-the number of the nodes in the cluster
  #p, the percentage of connecting edges 
  library(MASS)
  ind.cluster1<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq(1,n.size,1))))==2)))
  if(sum(ind.cluster1)>=1){
    new.edgelist.cluster1<-get.edgelist(g,names=T)[ind.cluster1==1,]
    edge.name.cluster1<-paste(new.edgelist.cluster1[,1],new.edgelist.cluster1[,2],sep="_")
    n.edge1=length(edge.name.cluster1)
    score.edge1<-runif(n.edge1,min=2,max=4)
    names(score.edge1)<-edge.name.cluster1
    E(g)$score[ind.cluster1==1]<-score.edge1}
  
  score.node.cluster1<-runif(n.size,min=0.1,max=1)
  names(score.node.cluster1)<-V(g)$name[1:n.size]
  V(g)$score[1:n.size]<-score.node.cluster1
  
  ##################################################################################
  
  
  #moderate node score,moderate edge score cluster
  ind.cluster4<-unlist(apply(get.edgelist(g,names=F),1,function(x) return(sum(length(intersect(x,seq((vcount(g)-n.size+1),vcount(g),1))))==2)))
  if(sum(ind.cluster4)>=1){
    
    new.edgelist.cluster4<-get.edgelist(g,names=T)[ind.cluster4==1,]
    edge.name.cluster4<-paste(new.edgelist.cluster4[,1],new.edgelist.cluster4[,2],sep="_")
    n.edge4=length(edge.name.cluster4)
    score.edge4<-runif(n.edge4,min=0.1,max=1)
    names(score.edge4)<-edge.name.cluster4
    E(g)$score[ind.cluster4==1]<-score.edge4
  }
  score.node.cluster4<-runif(n.size,min=2,max=4)
  names(score.node.cluster4)<-V(g)$name[seq((vcount(g)-n.size+1),vcount(g),1)]
  V(g)$score[seq((vcount(g)-n.size+1),vcount(g),1)]<-score.node.cluster4
  
  return(g)}

length(intersect(c(1,3),1:3))


sim.test.pos22.cluster5<-function(network,n.size,m,power,n.sim){
  library(igraph)
  library(BioNet)
  tpr.e<-c()
  fnr.e1<-c()
  tpr.n<-c()
  fnr.n1<-c()
  
  for (i in 1:n.sim){
    network.test<-randba.pos22.network5(network=network,m=m,power=power,n.size=n.size)
    node.score<-V(network.test)$score
    names(node.score)<-V(network.test)$name
    
    edge.score<-E(network.test)$score
    names(edge.score)<-E(network.test)$name
    
    module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
    module.n<-runFastHeinz(network=network.test, scores=node.score)
    
    node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
    
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.e<-c(tpr.e,(length(intersect(node.num.e,seq(1,n.size,1)))
    )/(n.size))
    fnr.e1<-c(fnr.e1,(length(intersect(node.num.e,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
    
    node.num.n<-as.numeric(unlist(strsplit(V(module.n)$name,"r"))[seq(2,2*length(V(module.n)$name),2)])
    
    n.start.cluster2<-floor(vcount(network.test)/2)
    tpr.n<-c(tpr.n,(length(intersect(node.num.n,seq(1,n.size,1)))
    )/(n.size))
    
    fnr.n1<-c(fnr.n1,(length(intersect(node.num.n,seq((vcount(network.test)-n.size+1),vcount(network.test),1)))/n.size))
    
  }
  return(data.frame(tpr_edge=tpr.e,fnr_edge1=fnr.e1,tpr_node=tpr.n,fnr_node1=fnr.n1))}
#Explore the TPR and FNR for 4 isolated clusters in random graph with different cluster size and connection probability
resultab23.sim11<-sim.test.pos22.cluster5(network=network,n.size=50,m=24,power=3,n.sim=100)
resultab23.sim12<-sim.test.pos22.cluster5(network=network,n.size=50,m=24,power=1,n.sim=100)
resultab23.sim21<-sim.test.pos22.cluster5(network=network,n.size=50,m=77,power=3,n.sim=100)
resultab23.sim22<-sim.test.pos22.cluster5(network=network,n.size=50,m=77,power=1,n.sim=100)

resultab23.sim13<-sim.test.pos22.cluster5(network=network,n.size=50,m=77,power=2,n.sim=100)
resultab23.sim23<-sim.test.pos22.cluster5(network=network,n.size=50,m=24,power=2,n.sim=100)
resultab23.sim223<-sim.test.pos22.cluster5(network=network,n.size=20,m=77,power=2,n.sim=100)
resultab23.sim213<-sim.test.pos22.cluster5(network=network,n.size=20,m=24,power=2,n.sim=100)
resultab23.sim214<-sim.test.pos22.cluster5(network=network,n.size=20,m=24,power=1.5,n.sim=100)
resultab23.sim211<-sim.test.pos22.cluster5(network=network,n.size=20,m=24,power=3,n.sim=100)
resultab23.sim212<-sim.test.pos22.cluster5(network=network,n.size=20,m=24,power=1,n.sim=100)
resultab23.sim221<-sim.test.pos22.cluster5(network=network,n.size=20,m=77,power=3,n.sim=100)
resultab23.sim222<-sim.test.pos22.cluster5(network=network,n.size=20,m=77,power=1,n.sim=100)

sim.resultab23.summary.223<-apply(resultab22.sim223,2,mean)
sim.resultab23.summary.23<-apply(resultab22.sim23,2,mean)
sim.resultab23.summary.13<-apply(resultab22.sim13,2,mean)

sim.resultab23.summary.21<-apply(resultab22.sim21,2,mean)
sim.resultab23.summary.22<-apply(resultab22.sim22,2,mean)
sim.resultab23.summary.12<-apply(resultab22.sim12,2,mean)
sim.resultab23.summary.11<-apply(resultab22.sim11,2,mean)

sim.resultab23.summary.211<-apply(resultab22.sim211,2,mean)
sim.resultab23.summary.212<-apply(resultab22.sim212,2,mean)
sim.resultab23.summary.221<-apply(resultab22.sim221,2,mean)
sim.resultab23.summary.222<-apply(resultab22.sim222,2,mean)



sim.resultab23.summary.13<-apply(resultab22.sim13,2,mean)
sim.resultab23.summary.23<-apply(resultab22.sim23,2,mean)

sim.resultab23.summary.223<-apply(resultab22.sim223,2,mean)











































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
vcount(network.test)
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
  
  g.test<-randba.pos.network2(network=network,m=24,power=1,n.size=50)
  randba.pos.network2<-function(network,m,power,n.size){
    
    library(igraph)
    library(BioNet)
    g <- barabasi.game(vcount(network.test),m=m,power=power,directed=F)
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  #generate the random graph using the number of nodes in the network
  #p:the probability of connecting edges
  #node scores in the nodes in the output random graph-uniform(-2,2)
  #edge scores in the edges in the output random graph-uniform(-4,4)
  rand.neg.network<-function(network,p){
    library(igraph)
    library(BioNet)
    g <- erdos.renyi.game(vcount(network), p)
    V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
    V(g)$score<-runif(vcount(g),min=-2,max=2)
    E(g)$score<-runif(ecount(g),min=-4,max=4)
    E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
    names(E(g)$score)<-E(g)$name
    return(g)}
  
  #generate the random scale-free graphs according to the Barabasi-Albert model 
  #using the number of nodes in the network
  #m:Numeric constant, the number of edges to add in each time step 
  #m=24,approximately genearte the same number of edge connections as expected of connecting probability 0.1 in erdos renyi graph model
  #m=77,,approximately genearte the same number of edge connections as expected of connecting probability 0.3 in erdos renyi graph model
  #power: The power of the preferential attachment, distribution of degrees of the nodes in Barabasi-Albert model 
  #node scores in the nodes in the output random graph-uniform(-2,2)
  #edge scores in the edges in the output random graph-uniform(-4,4)
  
  g.test<-randba.neg.network2(network=network,m=24,power=1,n.size=50)
  randba.neg.network2<-function(network,m,power){
    
    library(igraph)
    library(BioNet)
    g <- barabasi.game(vcount(network),m=m,power=power,directed=F)
    V(g)$name<-paste("Var",as.character(1:vcount(network)),sep="")
    pval.sim<-runif(vcount(g)) 
    names(pval.sim)<-V(g)$name
    node.scores<-node.score(da.igraph=g,pval=pval.sim,fdr=0.2)
    V(g)$score<-node.scores
    E(g)$score<-runif(ecount(g),min=-4,max=4)
    E(g)$name<-paste(get.edgelist(g,name=T)[,1],get.edgelist(g,name=T)[,2],sep="_")
    names(E(g)$score)<-E(g)$name
    return(g)}
  sim.test.neg.cluster3<-function(network,m,power,n.sim){
    
    library(igraph)
    library(BioNet)
    
    fnr.e<-c()
    
    fnr.n<-c()
    
    for (i in 1:n.sim){
      network.test<-randba.neg.network2(network=network,m=m,power=power)
      node.score<-V(network.test)$score
      names(node.score)<-V(network.test)$name
      
      edge.score<-E(network.test)$score
      names(edge.score)<-E(network.test)$name
      
      
      module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
      module.n<-runFastHeinz(network=network.test, scores=node.score)
      
      fnr.n<-c(fnr.n,vcount(module.n)/vcount(network.test))
      fnr.e<-c(fnr.e,vcount(module.e)/vcount(network.test))
      
    }
    return(data.frame(fnr_edge=fnr.e,fnr_node=fnr.n))}
  
  sim.test.neg.cluster2<-function(network,p,n.sim){
    
    library(igraph)
    library(BioNet)
    
    fnr.e<-c()
    
    fnr.n<-c()
    
    for (i in 1:n.sim){
      network.test<-rand.neg.network(network=network,p=p)
      node.score<-V(network.test)$score
      names(node.score)<-V(network.test)$name
      
      edge.score<-E(network.test)$score
      names(edge.score)<-E(network.test)$name
      
      
      module.e<-runFastHeinz.e(network=network.test, node.scores=node.score,edge.scores=edge.score)
      module.n<-runFastHeinz(network=network.test, scores=node.score)
      
      node.num.e<-as.numeric(unlist(strsplit(V(module.e)$name,"r"))[seq(2,2*length(V(module.e)$name),2)])
      
      fnr.n<-c(fnr.n,vcount(module.n)/vcount(network.test))
      fnr.e<-c(fnr.e,vcount(module.e)/vcount(network.test))
      
      
      
    }
    return(data.frame(fnr_edge=fnr.e,fnr_node=fnr.n))}
  
  
  resultneg.sim13<-sim.test.neg.cluster2(network,p=0.5,n.sim=10)
  
  
  resultneg.sim11<-sim.test.neg.cluster2(network=network,p=0.3,n.sim=100)
  resultneg.sim12<-sim.test.neg.cluster2(network,p=0.1,n.sim=100)
  
  simneg.result.summary.12<-apply(resultneg.sim12,2,mean)
  simneg.result.summary.11<-apply(resultneg.sim11,2,mean)
  simneg.result.summary.13<-apply(resultneg.sim13,2,mean)
  
  
  resultabneg.sim13<-sim.test.neg.cluster3(network,p=0.5,n.sim=10)
  
  
  resultabneg.sim11<-sim.test.neg.cluster3(network,m=24,power=1,n.sim=1000)
  resultabneg.sim12<-sim.test.neg.cluster3(network,m=24,power=2,n.sim=1000)
  resultabneg.sim22<-sim.test.neg.cluster3(network,m=77,power=2,n.sim=1000)
  resultabneg.sim21<-sim.test.neg.cluster3(network,m=77,power=1,n.sim=1000)
  
  
  simabneg.result.summary.12<-apply(resultabneg.sim12,2,mean)
  simabneg.result.summary.11<-apply(resultabneg.sim11,2,mean)
  simabneg.result.summary.21<-apply(resultabneg.sim21,2,mean)
  simabneg.result.summary.22<-apply(resultabneg.sim22,2,mean)
  