#some simple function used for p-values for permutation test sampling of correlations 
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
#the function-trans.mat.vec.cross,transform the p1Xp2 correlation matrix to the vector with the length-p1*p2
#input:x- correlations between p1 variables in row and p2 variables in column
#output:the vector of all correlations between p1 variables and p2 variables with length-p1*p2

trans.mat.vec.cross<-function(x){
  corr.pearson.vec<-c()
  for (i in 1:dim(x)[1]){
    corr.pearson.vec<-c(corr.pearson.vec,x[i,])
  }
  return(corr.pearson.vec)
}

#the function-pval.perm.corr.cross to calculate the P-values of the correlations=0 
#of the variables between the two input datasets-dat1 and dat2
#based on permuation test sampling assuming these variables are pairwisely independent
#equivalently the correlation=0
#input: dat1,dat2-the two datasets containing the samples of the variables in columns
#input:nsim, the number of sample size for permutation test sampling
#input:bvar-block paired id variable for paired study
#assuming the paired observation at each pairid group has the same size
#bvar can't have any nonmissing observations
#output:the vector of p-values based on permuation resampling samples with size-p1*p2
#where p1-the number of variables(columns) in dat1;p2-the number of variables(columns) in dat2;
pval.perm.corr.cross<-function(dat1,dat2,nsim,bvar=NA){
  n<-dim(dat1)[1]
  if(all(is.na(bvar))==F && any(is.na(bvar))==T)
  {warning("NA in paired group id")
    break}
   if(any(is.na(bvar))==F){
  size.b<-length(unique(bvar))
  }
  #the number of variales for correlatins
  p1<-dim(dat1)[2]
  p2<-dim(dat2)[2]
  #the number of correlations from p variables
  pp<-p1*p2
  
  mat.corr.perm<-matrix(NA,ncol=pp,nrow=nsim)
  dat<-cbind(dat1,dat2)
  for ( i in 1:nsim){
    #permuation random samples for further correlation
    if(all(is.na(bvar))==T){
    var.perm<-apply(dat,2,function(x){return(sample(x,size=n))})}
   if(any(is.na(bvar))==F){
      var.perm<-apply(dat,2,function(x){return(x[rep(sample(unique(bvar),size=size.b),each=n/size.b)])})
   }

    #correlaton from nsim permutation samples
    corr.perm<-cor(var.perm[,1:ncol(dat1)],var.perm[,(1+ncol(dat1)):(ncol(dat1)+ncol(dat2))], use="complete.obs") 
    mat.corr.perm[i,]<-trans.mat.vec.cross(corr.perm)}
  #Pearson's correlation from real data
  corr.true<-trans.mat.vec.cross(cor(dat1,dat2, use="complete.obs"))
  #calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
  p.val.perm<-c()
  for (j in 1:pp){
    p.val.perm<-c(p.val.perm,(sum(abs(mat.corr.perm[,j])>=abs(corr.true[j]))+1)/(nsim+1))}
  from<-rep(names(dat1),ncol(dat2))
  to<-rep(names(dat2),each=ncol(dat1))
  names(p.val.perm)<-paste(from,to,sep="_")
  return(p.val.perm)}

#calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
#inputs:dat-the dataset containing the variables for correlations,nsim: the number of permutation samples
#input:bvar-block paired id variable for paired study
#assuming the paired observation at each pairid group has the same size
#bvar can't have any nonmissing observations
#output:the vector of p-values based on permuation resampling samples with size-p(p-1)/2
#where p-the number of variables(columns) in dat
pval.perm.corr.submat<-function(dat,nsim,bvar=NA){
  n<-dim(dat)[1]
  #the number of variales for correlatins
  p<-dim(dat)[2]
  if( all(is.na(bvar))==F && any(is.na(bvar))==T)
  {warning("NA in paired group id")
    break}
  if(any(is.na(bvar))==F){
    size.b<-length(unique(bvar))
  }
  #the number of correlations from p variables
  pp<-p*(p-1)/2
  mat.corr.perm<-matrix(NA,ncol=pp,nrow=nsim)
  
  for ( ii in 1:nsim){
    
    #permuation random samples for further correlation
    #permuation random samples for further correlation
    if(all(is.na(bvar))==T){
      var.perm<-apply(dat,2,function(x){return(sample(x,size=n))})}
    if(any(is.na(bvar))==F){
      var.perm<-apply(dat,2,function(x){return(x[rep(sample(unique(bvar),size=size.b),each=n/size.b)])})
    }
  
    #correlaton from nsim permutation samples
    corr.perm<-cor(var.perm, use="complete.obs") 
    mat.corr.perm[ii,]<-trans.mat.vec(corr.perm)}
  
  #Pearson's correlation from real data
  corr.true<-trans.mat.vec(cor(dat, use="complete.obs"))
  #calculate the P-values from the permutation samples with the assumptions of indepedent variables in dataset-dat
  p.val.perm<-c()
  for (j in 1:pp){
    p.val.perm<-c(p.val.perm,(sum(abs(mat.corr.perm[,j])>=abs(corr.true[j]))+1)/(nsim+1))}
  from<-c()
  to<-c()
  for (i in 1:(dim(dat)[2]-1)){
    from<-c(from,rep(names(dat)[i],length((i+1):dim(dat)[2])))
    to<-c(to,names(dat)[-c(1:i)])}
  names(p.val.perm)<-paste(from,to,sep="_")
  return(p.val.perm)}

#the function to calculate the P-values of all the correlations of the input-variable datasets
#based on permuation test sampling assuming these variables are pairwisely independent 
#input: dat-the dataset containing the samples of the variables in columns
#input:nsim, the number of sample size for permutation test sampling
#input:bvar-block paired id variable for paired study
#assuming the paired observation at each pairid group has the same size
#bvar can't have any nonmissing observations
#output:the vector of p-values based on permuation resampling samples with size-p(p-1)/2
#output:the vector of p-values based on permuation resampling samples with size-p(p-1)/2
#where p-the number of variables(columns) in dat
pval.perm.corr<-function(dat,nsim,bvar=NA){
  #permutation test 
  #sample size(the number of observations-n
  n<-dim(dat)[1]
  #the number of variales for correlatins
  p<-dim(dat)[2]
#thresholding of the edge.scores for the massive number of edges in the network-sparsity

edge.score.thre<--3

  n.submat<-ceiling(p/20)
  if(n.submat==1){
    return(pval.perm.corr.submat(dat,nsim,bvar))}
  if(n.submat>1){
    pval.perm.corr.sub<-c()
    pval.perm.corr.cross<-c()
    for(i in 1:n.submat){
      pval.perm.corr.sub<-c(pval.perm.corr.sub,pval.perm.corr.submat(dat=dat[,c(((i-1)*20+1):min(p,20*i))],nsim,bvar))}
    for(i in 1:(n.submat-1)){
      for (j in (i+1):n.submat){
        pval.perm.corr.cross<-c(pval.perm.corr.cross,pval.perm.corr.cross(dat1=dat[,c(((i-1)*20+1):min(p,20*i))],dat2=dat[,c(((j-1)*20+1):min(p,20*j))],nsim,bvar))}}
    return(c(pval.perm.corr.sub,pval.perm.corr.cross))
  }}


