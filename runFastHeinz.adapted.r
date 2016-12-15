
## try http if https is not available
source("https://bioconductor.org/biocLite.R")
biocLite("BioNet")

.getPathScore <- function(path, graph1, graph2, node.score)
{
  
  sum(c(node.score[V(graph1)[path]$name], node.score[V(graph2)[unique(unlist(V(graph1)[path]$clusters))]$name]))
}

library(igraph)
library(graph)
library(BioNet)
network<-da.igraph.test
node.scores<-node.scores
edge.scores<-edge.score[edge.score>-1]
length(edge.scores)
length(node.scores)
node.scores[which(node.scores==NA)]
names(scores.bm)[1:10]
vcount(network)
ecount(network)
length(node.scores)
max(node.scores)

V(pos.subgraph)$weight
#toy example to incoporate the edge weight with correlation among the nodes
#correlation between -1 to 1-rho
#T^2 statistic-1/2log((1+rho)/(1-rho))
#$Fisher's Z test \frac{1}{2}\log((1+rho)/(1-rho)) approximately normally distrubuted 
#with mean 0 and variance-1/(N-3) 
#large correlation among the nodes in the network

set.seed(123)
V <- LETTERS[1:8]
edL <- vector("list", length=8)
names(edL) <- V
for(i in 1:8)
  edL[[i]] <- list(edges=9-i, weights=0.9)
gR <- new("graphNEL", nodes=V, edgeL=edL)

edges(gR)[[1]]
names(E(network)$weight[[length(V(network))]])
get.edgelist(network)[1,]
myscores <- c(-2, 1, 0.25, 1.25,-1,2,3,-1.34,1.14)
edL[[1]][2]

### Run the code from runFastHeinz.code.R

network <- gR
node.scores <- myscores
names(scores) <- LETTERS[1:8]

mat.corr<-matrix(0.9,nrow=vcount(network),ncol=vcount(network))
diag(mat.corr)<-0
colnames(mat.corr)<-names(V(network))
rownames(mat.corr)<-names(V(network))
mat.corr[get.edgelist(network)[1,1],get.edgelist(network)[1,2]]
E(network)$weight

runFastHeinz<-function (network, node.scores,edge.scores) 
{
    net.flag <- FALSE
    if (is.null(names(node.scores))) {
        warning("Unnamed node scores")
    }
    if (any(is.na(names(node.scores)))) {
        warning("NA in names of node scores, the node score will be removed")
        node.scores <- node.scores[!is.na(names(node.scores))]
    }
    
    if (is.null(names(edge.scores))) {
      warning("Unnamed edge scores")
    }
    if (any(is.na(names(edge.scores)))) {
      warning("NA in names of edge scores, the edge score will be removed")
     edge.scores <- edge.scores[!is.na(names(edge.scores))]
    }   
    
    
    
    
    if (is(network, "graphNEL")) {
        network <- igraph.from.graphNEL(network)
        net.flag <- TRUE
    }
    if (is.null(V(network)$name)) {
        V(network)$name <- as.character(V(network))
    }
    
    network<-subnetwork.e(network,vid=names(node.scores),eid=names(edge.scores))
    
    
    V(network)$score <- node.scores[V(network)$name]
    E(network)$score[1:20]

    E(network)$name<-names(edge.scores)
    E(network)$score<-edge.scores
    typeof(network)
    edge_attr(pos.subgraph)
    edge_attr
    vertex_attr(network)
#the edges with edge scores>0 

    pos.edges<- names(edge.scores[which(edge.scores > 0)])
 #identify the names of the nodes connecting the edges  with edge scores>0  
  
    pos.nodes <- names(node.scores[which(node.scores > 0)])
   
    if (length(pos.nodes) == 0) {
        warning("No positive nodes and positive edges")
        module <- graph.empty(n = 0, directed = FALSE)
        if (net.flag) {
            nE <- ecount(module)
            module <- simplify(module, remove.multiple = TRUE)
            if (nE != ecount(module)) {
                warning("Multiple edges between two nodes had to be removed for calculation")
            }
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
    if (length(names.pos.edges.nodes) == 1) {
        module <- .subNetwork0(pos.nodes, network)
        if (net.flag) {
            nE <- ecount(module)
            module <- simplify(module, remove.multiple = TRUE)
            if (nE != ecount(module)) {
                warning("Multiple edges between two nodes had to be removed for calculation")
            }
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
  
    
    
    length(pos.nodes)
    
    
	#find positive scoring nodes 
    pos.subgraph <- subNetwork(pos.nodes, network)
 from.name.edge.pos<-  E(pos.subgraph)$from[E(pos.subgraph)$score>0]
 to.name.edge.pos<-  E(pos.subgraph)$to[E(pos.subgraph)$score>0]
 edge.score.pos.subgraph<- E(pos.subgraph)$score[E(pos.subgraph)$score>0]
 
 .pos.subgraph<- E(pos.subgraph)$score[E(pos.subgraph)$score>0]
 pos.subgraph1<-subgraph.edges(pos.subgraph,name.edge.pos , delete.vertices = F)
 length(from.name.edge.pos)
 length(to.name.edge.pos)
 length(edge.score.pos.subgraph)
 length(E(network)$name)
 
 relations.pos.edges <- data.frame(from= from.name.edge.pos,
                         to= from.name.edge.pos,
                         score=edge.score.pos.subgraph
                         )
 
 
 g <- graph_from_data_frame(relations, directed=F, vertices=pos.nodes)
 
 
 
 
 
 #find positive scoring nodes with the positive edges
 pos.subgraph <- subnetwork.e(network,pos.nodes,pos.edges )
 plotModule(pos.subgraph)
 
 #identify components of that subgraph via some algorithm (not clear to me how it works)
 #I believe the idea is that there may be components of the positive subgraph that are not 
 #actually connected. These non connected clusters are the 'graph components'
 #ultimately yields list of clusters of graph nodes, with each list element a separate 'component'
 conn.comp.graph <- decompose.graph(pos.subgraph)
 
 #get the sum of the scores of all component elements
 node.score.comp <- unlist(lapply(lapply(conn.comp.graph, vertex_attr, 
                                         "score"), sum))
 
  
edge.score.comp <- unlist(lapply(lapply(conn.comp.graph, edge_attr, 
                                   "score"), sum))
#sum of total node scores and edge scores in each component
score.comp<-edge.score.comp+node.score.comp

	#re-order the positive scoring clusters, and the sum of the cluster scores, 
	#such that the highest scores appear first 
   comp.graph <- score.comp [order(score.comp, decreasing = TRUE)]
    score.comp <- sort(score.comp, TRUE)
	
	#assign cumulative scores to each cluster 
    for (i in 1:length(conn.comp.graph)) {
        conn.comp.graph[[i]]$score <- score.comp[i]
    }
	
	#now working with the network input to the function (not just the positive scoring nodes, as above)
	#get vertices in the entire network 
    v.id <- seq(1, vcount(network))
    names(v.id) <- V(network)$name
cluster.score<-unlist(lapply(conn.comp.graph,graph_attr,"score"))
    names(cluster.score)<-new.names 
    v.score<-c(V(network)$score,cluster.score)
    names(v.score)<-c(V(network)$name,names(cluster.score))
	#get edges in network
    edgelist <- get.edgelist(network, FALSE)
    edgelist1 <- edgelist[, 1]
    edgelist2 <- edgelist[, 2]
    length(E(network)$name)
    edgelist1[1:10]
	dim(edgelist)
	#now loop through and assign the same id to edge list elements
	#in the same cluster of positively scoring nodes 
    for (i in 1:length(conn.comp.graph)) {
        new.id <- length(V(network)) + i
        for (j in as.character(v.id[V(conn.comp.graph[[i]])$name])) {
            edgelist1[which(edgelist1 == j)] <- new.id
            edgelist2[which(edgelist2 == j)] <- new.id
        }
    }
	edgelist[1:500]
	edgelist1.network[which(edgelist1==473)]
	
#also give more specific names to those ids as cluster1, cluster2, etc.
new.ids <- seq(length(V(network)) + 1, length(V(network)) + 
        length(conn.comp.graph))
new.names <- paste("cluster", seq(1:length(conn.comp.graph)), 
        sep = "")
new.names[order(new.names)]
names(new.ids) <- new.names
	
#put together an object of all vertex ids, including cluster ids, and their names 
v.id <- c(v.id, new.ids)
v.name <- names(v.id)
names(v.name) <- v.id
	
#now put together a new edgelist, with the negatively scoring nodes given unique ids, 
#and the positively scoring nodes still clustered as above 
new.edgelist <- cbind(v.name[as.character(edgelist1)], v.name[as.character(edgelist2)])

edgelist1[1:100]


length(edgelist1)
length(E(network)$score)
	dim(new.edgelist)
	length(E(network)$name)
cluster.score[	v.name[edgelist1[which((edgelist1+edgelist2)  >(2*vcount(network)+1))]]]

v.name[edgelist1[which((edgelist1+edgelist2)  >(2*vcount(network)+1))]]
v.name[edgelist2[which((edgelist1+edgelist2)  >(2*vcount(network)+1))]]
edgelist1[edgelist1>477]

ind.cluster.pair<-list()
for (i in 1:length(new.names)){

  ind.cluster.pair[[i]]<-intersect(grep(paste("cluster",as.character(i),sep=""),v.name[edgelist1]),grep(paste("cluster",as.character(i),sep=""),v.name[edgelist2]))
}

ind.cluster.pair.diff<-list()
for (i in 1:(length(new.names)-1)){
  ind.cluster.pair.diff[[i]]<-list()
  for( j in ((i+1):length(new.names))){
  ind.cluster.pair.diff[[i]][[j]]<-union(intersect(grep(paste("cluster",as.character(i),sep=""),v.name[edgelist1]),grep(paste("cluster",as.character(j),sep=""),v.name[edgelist2])),
                                    intersect(grep(paste("cluster",as.character(j),sep=""),v.name[edgelist1]),grep(paste("cluster",as.character(i),sep=""),v.name[edgelist2])))
}}
new.edge.score.cluster.diff<-list()
for (i in 1:(length(new.names)-1)){
    new.edge.score.cluster.diff[[i]]<-list()
  for( j in ((i+1):length(new.names))){
    if(length(ind.cluster.pair.diff[[i]][[j]])>1){
new.edge.score.cluster.diff[[i]][[j]]<-sum(edge.scores[ind.cluster.pair.diff[[i]][[j]]])
new.edge.score.cluster.diff[[i]]$name<-names(edge.scores)[ind.cluster.pair.diff[[i]][[j]]]
    }
}
}







lapply(ind.cluster.pair.diff,function(x){ifelse(length(x)>1,)}


paste(v.name)


lapply(1:length(new.names),grep())

edge.scores[paste(v.name[edgelist1.network[ind.cluster.pair.diff[[4]][[5]]]],v.name[edgelist2.network[ind.cluster.pair.diff[[4]][[5]]]],sep="_")]
paste(v.name[edgelist1[ind.cluster.pair.diff[[4]][[5]]]],v.name[edgelist2[ind.cluster.pair.diff[[4]][[5]]]],sep="_")

#get edges in network
edgelist <- get.edgelist(network, FALSE)
edgelist1.network <- edgelist[, 1]
edgelist2.network <- edgelist[, 2]

c(1,2,3,4)[-2]

unique(from.name[edgelist1.network[which(edgelist1==473)]])
V(conn.comp.graph[[1]])$name

lapply()










ind.cluster.pair1<-unlist(ind.cluster.pair)

v.name[edgelist1[ind.cluster.pair1]]

v.name[edgelist2[ind.cluster.pair1]]















lapply(1:ecount(network),function(x){return(sum(edge.scores[x]))})
	edgelist2[which((edgelist1+edgelist2)  >(2*vcount(network)+1))]
	degree<-degree(network)
#using degree of nodes not new degree of nodes
cluster.degree<-unlist(lapply(lapply(conn.comp.graph,vertex_attr,"name"),function(x){return(sum(degree[x]))}))
names(cluster.degree)<-new.names
cluster.new.degree<-

	
interactome$new.degree<-new.degree.inter
new.degree.inter<-c()
for (i in 1:length(V(interactome2)$name)){
  if(grepl("cluster",V(interactome2)$name[i])==T)
}
mat.inter.name<-as.matrix(V(interactome2)$name,nrow=length(V(interactome2)$name),ncol=1)
degree.inter<-apply(mat.inter.name,1,function(x){return(ifelse(grepl("cluster",x),degree(interactome2)[which(V(interactome2)$name==x)],new.degree[x]))})
names(new.degree.inter)<-V(interactome2)$name
max(new.degree.inter)==max(cluster.degree)
max(new.degree.inter[which(grepl("cluster",V(interactome2)$name)==F)])==max(new.degree)
V(interactome2)$name[which(grepl("cluster",V(interactome2)$name)==T)]
	
	#convert the edgelist to a graph
    interactome2 <- graph.edgelist(new.edgelist, FALSE)
    E( interactome2)$score<-edge.scores
    E( interactome2)$name<-names(edge.scores)
    V(interactome2)$score<-v.score[V(interactome2)$name]
    
    
from.name.cluster=c()
to.name.cluster=c()
#identify the nodes and the edges connecting in each cluster
from.name.cluster=c()
for (j in 1:length(from.name)){
  from.name.cluster[j]<-from.name[j]
  for(i in 1:length(conn.comp.graph)){   
   if(length(intersect(from.name[j],V(conn.comp.graph[[i]])$name))==1){
    from.name.cluster[j]<-new.names[i]
  }
  }
}

to.name.cluster=c()
for (j in 1:length(to.name)){
  to.name.cluster[j]<-to.name[j]
  for(i in 1:length(conn.comp.graph)){
    if(length(intersect(from.name[j],V(conn.comp.graph[[i]])$name))==1){
      to.name.cluster[j]<-new.names[i]
    }}
  }

length(from.name.cluster)
length(to.name.cluster)

length(from.name)
length(to.name)

length(from.name.cluster.1)
length(to.name.cluster.1)

#remove the self-loops in the network
from.to.name.cluster<-cbind(from.name.cluster,to.name.cluster)
for(i in 1:length(from.name.cluster)){
  if(grepl("cluster",from.name.cluster[i]) && grepl("cluster",to.name.cluster[i]) && from.name.cluster[i]==to.name.cluster[i]){
from.to.name.cluster[i,]<-c(NA,NA)
  }   
}


from.name.cluster.1<-from.to.name.cluster[,1][!is.na(from.to.name.cluster[,1])]  
to.name.cluster.1<-from.to.name.cluster[,2][!is.na(from.to.name.cluster[,2])]  

  
v.name.not.cluster<-V(interactome2)$name[which(grepl("cluster",V(interactome2)$name)==F)]
ind.rep.cluster<-list()
for(i in 1:length(v.name.not.cluster)){ 
  ind.rep<-c()
 for (j in 1:length(new.names)){
ind.rep<-c(ind.rep,union(intersect(grep(v.name.not.cluster[i],from.name.cluster) ,
                                                                grep(new.names[j],to.name.cluster)),
                          intersect(grep(v.name.not.cluster[i],to.name.cluster) , 
                                    grep(new.names[j],from.name.cluster))))
if(j==length(new.names)){ind.rep.cluster[[i]]<-ind.rep}
  } 
}
edgelist.12<-cbind(edgelist1,edgelist2)
da<-c(2,3,4,5,64,3,2,2,3)
duplicated(da)

da[duplicated(da)]
edgelist.12.unique<-edgelist.12[!duplicated(edgelist.12),]
unlist(ind.rep.cluster)

edgelist.12[duplicated(edgelist.12),][1:5,]
length(unique(edgelist.12[duplicated(edgelist.12),]))




edge.scores[paste(v.name[edgelist1.network[duplicated(edgelist.12),][1:5,1]],v.name[edgelist2.network[duplicated(edgelist.12),][1:5,2]],sep="_")]
length(from.name)
length(edgelist1.network)


intersect(grep(v.name.not.cluster[1],from.name.cluster),grep(new.names[1],to.name.cluster))
length(ind.rep.cluster)
from.name.cluster[grep(new.names[1],to.name.cluster)]

grep("cluster1",from.name.cluster)


to.name.cluster[1:500]

length(from.name.cluster)

v.name.not.cluster[1:200]








    apply(mat.from.name,1,lapply(lapply(conn.comp.graph,vertex_attr,"name"),function(x){}
    
  from.name.cluster<-apply(mat.from.name,1,function(x){ifelse(length(intersect(x,V(conn.comp.graph[[i]])$name))==1,new.names[i],x)})  
  
}
   intersect(V(conn.comp.graph[[1]])$name,
    
from.name<-unlist(strsplit(names(edge.scores),"_"))[seq(1,2*length(names(edge.scores)),by=2)]
mat.from.name<-matrix(from.name,nrow=length(from.name),ncol=1)

to.name<-unlist(strsplit(names(edge.scores,"_"))[seq(2,2*length(names(edge.scores)),by=2)]
mat.to.name<-matrix(to.name,nrow=length(to.name),ncol=1)               

grep(from.name,V(conn.comp.graph[[i]])$name)


                
             edge.name.sub<-as.matrix(cbind(from.name.sub,to.name.sub))
             edge.name.sub.node<-matrix(NA,nrow=dim(edge.name.sub)[1],ncol=dim(edge.name.sub)[2])
             for(i in 1:dim( edge.name.sub)[1]){
               if(length(intersect(vid,edge.name.sub[i,]))==2){
                 edge.name.sub.node[i,] <-edge.name.sub[i,]
               }
               else
                 edge.name.sub.node[i,] <-c(NA,NA)
               
             }
             dim( edge.name.sub.node)  
             
             
             
             
             
             
             
 cbind(v.name[edgelist1[grep("cluster2",v.name[edgelist1])]], v.name[edgelist2[grep("cluster2",v.name[edgelist1])]])            
             
             
 cbind(v.name[edgelist1[grep("cluster4",v.name[edgelist1])]], v.name[edgelist2[grep("cluster4",v.name[edgelist1])]])           
             v.name
             
   V(conn.comp.graph[[4]])        
             
             
             
    
    length( E(interactome2)$name)
    grep("cluster",E(interactome2)$name)
    length(edgelist1)
    length(edgelist2)
    length(E( interactome2)$score)
    v.name[edgelist1[1:100]]
    E( interactome2)$name
    plotModule(interactome2)
    
    names(edge.scores)[1:100]
	
	#now initialize weights for each edge
    E(interactome2)$weight <- rep(0, length(E(interactome2)))
    ecount(interactome2)
	
#remove any edges that start and end at the same node, and remove duplicate edges
#This will leave only one member of each of the positively scoring clusters, with the cumulative score 
#for that cluster assigned to it, and the rest of the negatively scoring nodes 
#Because all members of a cluster were assigned the same id and score, these all appeared as self loops
    #convert the edgelist to a graph
    interactome2 <- graph.edgelist(new.edgelist, FALSE)
    E( interactome2)$score<-edge.scores
    E( interactome2)$name<-names(edge.scores)
    V(interactome2)$score<-v.score[V(interactome2)$name]   
    
interactome2.1 <- simplify(interactome2, remove.loops = T, remove.multiple = T, edge.attr.comb=list(score="sum"))
interactome2.1.2 <- simplify(interactome2.1, remove.loops = T, remove.multiple = T, edge.attr.comb="sum")

interactome2.2 <- simplify(interactome2, remove.loops = F, remove.multiple = T, edge.attr.comb="sum")  

interactome2.2.1 <- simplify(interactome2.2, remove.loops = T, remove.multiple = T, edge.attr.comb="sum")  

sort(table(E(interactome2.1)$score))
ecount(interactome2.1)
ecount(interactome2.2)
ecount(interactome2.2.1)

sort(table(E(interactome2.2)$score))


length(E(interactome2.1)$score)

length(unique(cbind(edgelist1,edgelist2))  ) 

E(interactome2)$score[as.vector(table(E(interactome2)$score))>1000]
    
class(table(E(interactome2)$score))
    
sort(as.matrix(table(E(interactome2)$score)),decreasing =T )
    
g <- graph( c(1,2, 1,2, 1,2, 2,3, 3,4))

E(g)$weight <- 1:5
    
    ## print attribute values with the graph
    igraph.options(print.graph.attributes=TRUE)
    igraph.options(print.vertex.attributes=TRUE)
    igraph.options(print.edge.attributes=TRUE)
    
    ## new attribute is the sum of the old ones
   g1<- simplify(g, edge.attr.comb="sum")
    
    ## collect attributes into a string
   g2<- simplify(g, edge.attr.comb=toString)
    E(g1)$weight
    edge_attr(g2)

#############################################################################################################
#define the degree of the nodes based on the network with the network with positive edge scoring
#reflect the strength of connectness of the two nodes by the edge score>0 based on adjusted likelihood ration adjusting for FDR
#in the multiple testing of all the correlations among the vertex in the network
length(names(node.scores))
names(edge.scores[edge.scores>0])    
degree.network<-subnetwork.e(network,vid=names(node.scores),eid=names(edge.scores[edge.scores>0]))    
#the new degree of the nodes, if node score positive, (maximum of the positive scoring nodes-degree of the nodes+1),
#if node score negative,degree of the nodes
max.degree.pos.nodes<-max(degree(subnetwork.e(network,vid=names(node.scores[node.scores>0]),

                                              
                                                                                            eid=names(edge.scores[edge.scores>0]))))
new.degree<-c()
for(i in 1:length(V(degree.network))){
new.degree<-c(new.degree, ifelse(V(degree.network)$score[i]>0,
                                 max.degree.pos.nodes-degree(degree.network)[i]+1,degree(degree.network)[i]))
}
names(new.degree)<-V(network)$name    
V(network)$new.degree<-new.degree
#the degree of network by the network with positive connecting edge scores 
V(network)$degree<-degree(degree.network)
degree<-degree(degree.network)
names(degree)<-V(network)$name   

#degree for the clusters after removing the internal number of connecting edges within each cluster
#using degree of nodes not new degree of nodes
cluster.degree<-unlist(lapply(lapply(conn.comp.graph,vertex_attr,"name"),function(x){return(sum(degree[x]))}))
names(cluster.degree)<-new.names
cluster.new.degree<-max(cluster.degree)-cluster.degree+1
names(cluster.new.degree)<-new.names


get.edgelist(interactome2, 
             FALSE)[1:100, 1]


#get a vector of scores for the nodes in the graph
node.score1 <- node.scores[V(interactome2)$name]
V(interactome2)$score<-v.score[V(interactome2)$name]
max(V(interactome2)$score)
node.score2<-v.score[V(interactome2)$name]
names(node.score1) <- V(interactome2)$name
names(node.score2) <- V(interactome2)$name

length(node.score2)


grepl("cluster",names(node.score2))
score.degree<-c()
as.numeric(unlist(strsplit(names(node.score2)[i],"r"))[2])
for(i in 1:length(node.score2)){
  id.cluster<-as.numeric(unlist(strsplit(names(node.score2)[i],"r"))[2])
score.degree <- c(score.degree,ifelse(grepl("cluster",names(node.score2)[i]),
cluster.score[id.cluster]/(cluster.new.degree[id.cluster]+1),node.score2[i]/(new.degree[i]+1)))
}
V(interactome2)$score.degree <- score.degree


E(interactome2)$weight <- -(V(interactome2)[get.edgelist(interactome2, 
                                                         FALSE)[, 1]]$score.degree + V(interactome2)[get.edgelist(interactome2, 
                                                                                                                  V(interactome2)                                                                                                                  FALSE)[, 2]]$score.degree)

V(interactome2)$name











new.degree[V(conn.comp.graph[[1]])$name]
    
    degree(network)
    max(degree(network))-degree(network)+1
    V(network)$score[i]>0
		degree()
	#get a vector of scores for the nodes in the graph
    node.score1 <- node.scores[V(interactome2)$name]
    V(interactome2)$score<-v.score[V(interactome2)$name]
    max(V(interactome2)$score)
    node.score2<-v.score[V(interactome2)$name]
    length(node.score1)
    max(node.score2)
    length(node.score2)
    names(node.score1) <- V(interactome2)$name
    names(node.score2) <- V(interactome2)$name
    
	#all of the positively scoring nodes have been assigned names starting with 'cluster'
	#so they don't appear in the scores object above (which is one of this function's required inputs)
	#therefore, all of those scores will be NA, and should be reset to zero 
    node.score1[which(is.na(node.score1))] <- 0
  
    
    V(interactome2)$score[V(interactome2)$score>0]
	
	#assign a degree to each node (i.e., the score divided by the number of edges coming out of the node +1)
    unlist(lapply(lapply(interactome2,vertex_attr,"score"),function(x){
      return(ifelse(x>0,)
    }
    
    vertex_attr(interactome2,"score")
    score.degree <- node.score1/(igraph::degree(interactome2) + 1)
    V(interactome2)$score.degree <- score.degree
    max(igraph::degree(interactome2))
    min(igraph::degree(interactome2))
    
    
    
	plotModule(interactome2)
	length(score.degree)
	#assign an edge weight equal to the negative sum of the two input node 'score degrees' 
	#plus the edge degree connecting the two input nodes calculated in the previous step 
   
     E(interactome2)$weight <- -(V(interactome2)[get.edgelist(interactome2, 
        FALSE)[, 1]]$score.degree + V(interactome2)[get.edgelist(interactome2, 
        FALSE)[, 2]]$score.degree)
		length( E(interactome2)$weight )
		
		
		
		
		
		

	plotModule(mst)
	#compute sum of the scores for the positively scoring elements of each 
	#of the sub-clusters of the graph 
    sum.pos <- lapply(decomp.graphs, function(x) {
        sum(node.score[names(which(node.score[V(x)$name] > 0))])
    })
	
	#now restrict the graph to the cluster of nodes with the highest scoring positive elements
	#note here that the which.max() command is used to return the index of the maximum value
    interactome2 <- decomp.graphs[[which.max(sum.pos)]]
    rm(decomp.graphs)
	
	#now find a minimum spanning tree for the restricted graph
	#this will likely produce a still very large graph
    mst <- minimum.spanning.tree(interactome2, weights = E(interactome2)$weight)
	
	#identify the indices of the positively scoring clusters in the mst
    mst.cluster.id <- grep("cluster", V(mst)$name)
	
	#assign names to those indices corresponding to cluster id 
    names(mst.cluster.id) <- V(mst)[mst.cluster.id]$name
	
	#re-order by cluster number 
    mst.cluster.id <- mst.cluster.id[order(as.numeric(matrix(unlist(strsplit(names(mst.cluster.id), 
        "cluster")), nrow = 2)[2, ]))]
	
	#initiate vector for unknown reason
    all.ids <- c()
    if (length(mst.cluster.id) == 1) {
        neg.node.ids.2 = c()
    }
    else {
		
		#not completely clear on what is happening here
		#this loops through the positively scoring nodes and finds shortest paths to all other
		#positively scoring nodes 
        for (j in 1:(length(mst.cluster.id) - 1)) {
            path <- unlist(get.all.shortest.paths(mst, from = mst.cluster.id[j], 
                to = mst.cluster.id[(j + 1):length(mst.cluster.id)]))
            all.ids <- c(all.ids, path)
        }
		
		#this somehow identifies all nodes related to shortest paths in the graph
        all.ids <- unique(all.ids)
		
		#now this subsets the graph to include only those ids involved in shortest paths 
        sub.interactome2 <- subNetwork(V(mst)[all.ids]$name, 
            interactome2)
			plotModule(sub.interactome2)
		#now identify negative scoring nodes in the graph 
        neg.node.ids <- which(node.score[V(sub.interactome2)$name] < 
            0)
		
		#for each of the negatively scoring nodes, determine if a positively scoring 
		#cluster lies adjacent, and store the cluster ids in a list 
        for (i in neg.node.ids) {
            V(sub.interactome2)[i]$clusters <- list(neighbors(sub.interactome2, 
                v = i)[grep("cluster", V(sub.interactome2)[neighbors(sub.interactome2, 
                v = i)]$name)])
        }
        V(sub.interactome2)$clusters
		#initiate vector for negative node scores 
        score.neg.nodes <- c()
		
		#for each negative scoring node, if a positive scoring cluster is adjacent, sum the score 
		#for that node and the positive cluster score 
		#Otherwise, just assign the score of the node 
        for (i in neg.node.ids) {
            if (length(V(sub.interactome2)[i]$clusters[[1]]) > 
                0) {
                score.neg.nodes <- c(score.neg.nodes, sum(node.score[V(sub.interactome2)[c(i, 
                  V(sub.interactome2)[i]$clusters[[1]])]$name]))
            }
            else {
                score.neg.nodes <- c(score.neg.nodes, node.score[V(sub.interactome2)[i]$name])
            }
        }
		
		#id any negative nodes attached to positive clusters, whose sum total score is positive 
        neg.node.ids.2 <- neg.node.ids[score.neg.nodes > 0]
    }
	
	#if there are no negative nodes that connect to positive clusters such that the sum total 
	#score is positive, then grab the highest scoring cluster and assign it to object 'module' 
    if (length(neg.node.ids.2) == 0) {
        module <- unlist(lapply(conn.comp.graph, get.vertex.attribute, 
			#the term in the square brackets is a very complicated way of getting the number of the cluster with 
			#the highest score 
            "name")[as.numeric(matrix(unlist(strsplit(names(node.score.cluster)[which.max(node.score.cluster)], 
            "cluster")), nrow = 2)[2, ])])
        module <- .subNetwork0(module, network)
        if (net.flag) {
            nE <- ecount(module)
            module <- simplify(module, remove.multiple = TRUE)
            if (nE != ecount(module)) {
                warning("Multiple edges between two nodes had to be removed for calculation")
            }
            module <- igraph.to.graphNEL(module)
        }
        return(module)
    }
	
	#Otherwise, if there are some negative nodes attached to positive scoring clusters, such that 
	#the total combined score is positive, identify the largest component of the graph containing
	#those negative nodes 
    subg <- largestComp(induced.subgraph(sub.interactome2, neg.node.ids.2))
	plotModule( mst.subg)
	#find a minimum spanning tree for that graph
    mst.subg <- minimum.spanning.tree(subg, E(subg)$weight)
	
	#initialize maximum score 
    max.score <- 0
	
	#initialize best path vector 
    best.path <- c()
	
	#loop through the nodes in the mst subgraph 
    for (i in 1:(length(V(mst.subg)))) {
	
		#identify shortest paths from all nodes
        path <- get.all.shortest.paths(mst.subg, from = V(mst.subg)[i])
		
		#sum the score of the negative nodes and the adjacent positive clusters 
        path.score <- unlist(lapply(path$res, .getPathScore, 
            graph1 = mst.subg, graph2 = sub.interactome2, node.score = node.score))
			plotModule(sub.interactome2)
		#get the index of the highest scoring path 
        best.pos <- which.max(path.score)
		
		#if the score associated with that path is higher than the current 
		#high score, re-assign the high score and best path 
        if (path.score[[best.pos]] > max.score) {
            best.path <- path$res[[best.pos]]
            max.score <- path.score[[best.pos]]
        }
    }
    ? intersect
    best.path<-path$res[[1]]
    if (length(best.path) != 1) {
	
		#for each of the elements in the highest scoring path, determine the adjacent clusters 
        cluster.list <- V(mst.subg)[best.path]$clusters
        names.list <- as.character(1:length(cluster.list))
        names(cluster.list) <- names.list
        names(best.path) <- names.list
		
		#convoluted way of checking for overlap among clusters and filtering some out for some reason (I believe having to do with repeat clusters)
        for (i in names.list) {
            res <- lapply(cluster.list, intersect, cluster.list[[i]])
            if (length(intersect(unlist(cluster.list[as.character(which(as.numeric(names.list) < 
                as.numeric(i)))]), unlist(cluster.list[as.character(which(as.numeric(names.list) > 
                as.numeric(i)))]))) > 0) {
                if (length(setdiff(res[[i]], unique(unlist(res[names(res) != 
                  i])))) == 0) {
                  cluster.list <- cluster.list[names(cluster.list) != 
                    i]
                  names.list <- names.list[names.list != i]
                }
            }
        }
		
		#id vertices numbers in final best path 
        best.path <- best.path[names.list]
    }
	
	#get names of negative scoring vertices in final best path 
    module <- V(mst.subg)[best.path]$name

	#get names of adjacent positive scoring vertex clusters 
    pos.cluster <- V(sub.interactome2)[unique(unlist(V(mst.subg)[best.path]$clusters))]$name
	
	#combine into one character vector of the final vertex names 
    module <- c(module, unlist(lapply(conn.comp.graph, get.vertex.attribute, 
        "name")[as.numeric(matrix(unlist(strsplit(pos.cluster, 
        "cluster")), nrow = 2)[2, ])]))
		
	#make final subset igraph object 
    module <- subNetwork(module, network)
    plotModule(module)
    if (net.flag) {
        nE <- ecount(module)
        module <- simplify(module, remove.multiple = TRUE)
        if (nE != ecount(module)) {
            warning("Multiple edges between two nodes had to be removed for calculation")
        }
		
		#covnert from igraph to graphNEL
        module <- igraph.to.graphNEL(module)
    }
    return(module)
}



mydf <- read.table('thefile.txt', header=TRUE, sep="\t", fileEncoding="windows-1252")
str(mydf)








## try http if https is not available
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
## try http if https is not available
source("https://bioconductor.org/biocLite.R")


## try http if https is not available
source("https://bioconductor.org/biocLite.R")
biocLite("DLBCL")
biocLite("RBGL")

options(width=60)
ps.options(family="sans")

library(BioNet)
library(DLBCL)
data(dataLym)
data(interactome)

pvals <- cbind(t=dataLym$t.pval, s=dataLym$s.pval)
rownames(pvals) <- dataLym$label
pval <- aggrPvals(pvals, order=2, plot=FALSE)

##check the final optimized module identified
E(module)$score
E(module)$name
V(module)$score
V(module)$name
dat.finalmodule.v<-data.frame(v.name=V(module)$name,v.score=V(module)$score)
dat.finalmodule.e<-data.frame(e.name=E(module)$name,e.score=E(module)$score)
lst.finalmodule<-list(vertex=dat.finalmodule.v,edge=dat.finalmodule.e)

write.table(dat.finalmodule.v,"finalmodule_adaptedrunFastHeitz_v.txt",row.names=F,col.names=T)
write.table(dat.finalmodule.e,"finalmodule_adaptedrunFastHeitz_e.txt",row.names=F,col.names=T)
dat.finalmodule<-append(dat.finalmodule.v,dat.finalmodule.e)
class(dat.finalmodule)