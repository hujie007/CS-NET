rm(list = ls())
library("igraph")
library("clue")

#初始embedding向量经过row-normalize

##########################平方误差惩罚项#################################
sqlike_hood<-function(alpha, C_hat, hg, A, mu, lambda){ 
 #alpha为n*p阶矩阵 
 #C_hat为初始分组中心变量 G*p阶矩阵
 #hg为每个节点的潜在变量和组中心最小化距离的组标签
 #A为邻接矩阵
 #mu为sparse parameter
 #lambda为penalty parameter
 index<-lower.tri(A)
 n<-nrow(A)
 C_hat1<-C_hat[hg,]#每个点对应的组中心(坐标)
 B<-alpha%*%t(alpha)
 B_d<-diag(B)
 B<-B+mu
 B1<-B*A
 P_hat<-1+exp(B) 
 diag(P_hat)<-0
 alpha_t<-rowSums((alpha-C_hat1)^2) #每个点与组中心的距离二范数的平方
 likely<-sum(log(P_hat[index])-B1[index])/(n*(n-1)/2)+lambda*sum(alpha_t)
 return(likely)
 
}

sql_partial<-function(alpha, C_hat, hg, A, mu, lambda){ #分别对mu和alpha求偏导，以便于之后的坐标梯度下降法对两变量迭代
 index<-lower.tri(A)
 C_hat1<-C_hat[hg,]
 B<-alpha%*%t(alpha)
 B_d<-diag(B)
 P_hat<-1/(1+exp(-B-mu)) 
 diag(P_hat)<-0
 P_hat<-P_hat-A
 partial_mu<-2*sum(P_hat[index])/(n*(n-1))
 partial_alpha<-2*P_hat%*%alpha/(n*(n-1))+lambda*(2*(alpha-C_hat1))
 return(list(p_mu=partial_mu, p_alpha=partial_alpha ))
 
}


#聚类得到n个节点的分组 及 组中心
sq_cluster<-function(alpha, C_hat, ep, iter_max){
 G<-nrow(C_hat)
 hg<- numeric(nrow(alpha))
 alpha<-t(apply(alpha,1,function(x){x/sqrt(sum(x^2))})) #标准化潜在向量alpha
 iter<-0 #迭代次数
 delta<-1 #前后迭代得到的似然函数之差
 
 while (iter<iter_max && delta>ep) {
  C_hat1<-C_hat
  
  # 分别更新每个点的组
  for (i in 1:n){
   dist<- rowSums((matrix(1,G,1)%*%alpha[i,]-C_hat1)^2)
   hg[i]<- which.min(dist)
  }
  
  # 更新组中心
  for (g in 1:G) {
   # y<-colSums(alpha[hg==g, ])
   # C_hat[g,]<-y/sqrt(sum(y^2))
   
   C_hat[g,]<- colSums(alpha[hg==g, ])/sum(hg==g)
  } 
  delta<-sqrt(sum((C_hat-C_hat1)^2)/G) #计算新旧中心的距离，再据此判断是否继续迭代(直至中心变化很小)
  iter<-iter+1
 }
 
 for (i in 1:n){
  dist<- rowSums((matrix(1,G,1)%*%alpha[i,]-C_hat1)^2)
  hg[i]<- which.min(dist)
 }
 return(list(group=hg, center=C_hat1, iters=iter))
 
}


#根据分组中心及分组情况利用坐标梯度下降法估计出 每个节点的潜在变量 alpha、分组结果g_hat、组中心C_hat
sq_embedding<-function(A, G, p, mu , lambda, gamma, ep, iter_max){
 
 # #initial values
 n<-nrow(A) 
 pho<-sum(A)/(n*(n-1))
 pAdj<-A+as.numeric(0.25*pho)*matrix(1,nrow=n,ncol=n) #purmutated Adj  
 #compute the Laplace Matrix
 D <- apply(pAdj,1,sum)^(-0.5)
 D[D==Inf] <- 0
 H <- diag(D)
 L <- H%*%pAdj%*%H 
 U<-eigen(L)
 alpha<-U$vectors[,1:p]
 
 #alpha<- alpha/sqrt(apply(alpha^2,1,sum)) # 对特征向量进行 row-normalize
 
 # 对alpha做kmeans
 fitk<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
 C_hatk<- fitk$centers #kmeans初始分组中心坐标
 g_k<-fitk$cluster #kmeans初始聚类分组
 
 
 #原程序是先对alpha进行kmeans,然后再利用kmeans的结果用sq_cluster函数得到的值作为初始值
 #现程序直接用kmeans的值作为初始值
 C_hat<- C_hatk
 g_hat<- g_k
 # fit<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
 # g1<-fit$cluster#初始聚类分组
 # C_hat<-fit$centers#初始聚类中心
 # 
 # fit<-sq_cluster(alpha, C_hat=C_hat, ep=ep, iter_max=100)
 # C_hat<-fit$center
 # g_hat<-fit$group
 iter=0
 delta=1
 while(iter<iter_max && delta>ep){
  
  alpha1<-alpha
  C_hat1<-C_hat
  mu1<-mu
  g_hat1<-g_hat
  like1<-sqlike_hood(alpha1, C_hat1, g_hat1, A, mu1, lambda=lambda)
  
  fit<-sql_partial(alpha1, C_hat1, g_hat1, A, mu, lambda=lambda) 
  mu<-mu1-gamma*fit$p_mu 
  alpha<-alpha1-gamma*fit$p_alpha 
  #alpha<-up_alpha(alpha1, C_hat1, g_hat1, A, mu, gamma, lambda)
  
  fitk<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
  C_hat<- fitk$centers #kmeans初始分组中心坐标
  g_hat<-fitk$cluster
  like0<-sqlike_hood(alpha, C_hat, g_hat, A, mu, lambda=lambda)
  
  #如果新的似然大于旧的，则减小步长(gamma)再计算
  if(like0>like1){
   gamma<-gamma/2
   alpha<-alpha1
   C_hat<-C_hat1
   mu<-mu1
   g_hat<-g_hat1
  } 
  #如果新的似然小于旧的,寻找负对数似然的最小值(对数似然的最大值)
  if(like0<=like1){
   delta<-(like1-like0)/abs(like1)
   iter<-iter+1
  }
  
 }
 g_hat<-unname(g_hat)
 embedding<- alpha
 mu_hat<- mu
 return(list(kgroup=g_k, group=g_hat, embedding= embedding, mu_hat=mu_hat, center=C_hat, likely=like0, Gamma=gamma, iters=iter))
}

#############################余弦惩罚项##################################

coslike_hood<-function(alpha, C_hat, hg, A, mu, lambda){ 
 #alpha为n*p阶矩阵 
 #C_hat为初始分组中心变量 G*p阶矩阵
 #hg为每个节点的潜在变量和组中心最大化夹角的组标签
 #A为邻接矩阵
 #mu为sparse parameter
 #lambda为penalty parameter
 index<-lower.tri(A)
 n<-nrow(A)
 C_hat1<-C_hat[hg,]
 B<-alpha%*%t(alpha)
 B_d<-diag(B)
 B<-B+mu
 B1<-B*A
 P_hat<-1+exp(B) 
 diag(P_hat)<-0
 alpha_t<-rowSums(alpha*C_hat1)
 likely<-sum(log(P_hat[index])-B1[index])/(n*(n-1)/2)+lambda*sum(1-alpha_t/sqrt(B_d))
 return(likely)
 
}

cosl_partial<-function(alpha, C_hat, hg, A, mu, lambda){ #分别对mu和alpha求偏导，以便于之后的坐标梯度下降法对两变量迭代
 index<-lower.tri(A)
 C_hat1<-C_hat[hg,]
 B<-alpha%*%t(alpha)
 B_d<-diag(B)
 P_hat<-1/(1+exp(-B-mu)) 
 diag(P_hat)<-0
 alpha_t<-rowSums(alpha*C_hat1)
 P_hat<-P_hat-A
 partial_mu<-2*sum(P_hat[index])/(n*(n-1))
 partial_alpha<-2*P_hat%*%alpha/(n*(n-1))-lambda*(C_hat1/sqrt(B_d)
                                                  -(alpha_t/((B_d)^(3/2)))*alpha )
 return(list(p_mu=partial_mu, p_alpha=partial_alpha ))
 
}

# up_alpha<-function(alpha, C_hat, hg, A, mu, gamma, lambda){
#   n<-nrow(alpha)
#   for (i in 1:n) {
#    x<-alpha[i,]
#    y<-C_hat[hg[i],]
#    B<-alpha%*%x
#    phat<-1/(1+exp(-B[-i]-mu))-A[i,-i]
#    partial_alpha<-2*(t(alpha[-i,])%*%phat)/(n*(n-1))-lambda*(y/sqrt(B[i])
#                                                     -( sum(x*y)/(2*B[i]^(3/2)))*x ) 
#    alpha[i,]<-x-gamma*partial_alpha
#   }
#   return(alpha)
# }


#聚类得到n个节点的分组 及 组中心
cos_cluster<-function(alpha, C_hat, ep, iter_max){
 G<-nrow(C_hat)
 #alpha<-t(apply(alpha,1,function(x){x/sqrt(sum(x^2))})) #标准化潜在向量alpha
 iter<-0 #迭代次数
 delta<-1 #前后迭代得到的似然函数之差
 
 while (iter<iter_max && delta>ep) {
  C_hat1<-C_hat
  
  t_hat<-alpha%*%t(C_hat1)
  hg<-apply(t_hat, 1, which.max)
  
  for (g in 1:G) {
    y<-colSums(alpha[hg==g, ])
    C_hat[g,]<-y/sqrt(sum(y^2))
   
    # y<- alpha[hg==g, ]
    # h<- apply(y/sqrt(diag(y %*% t(y))),2,sum)
    # C_hat[g,]<- h/sqrt(t(h)%*%h)
  } 
  delta<-sqrt(sum((C_hat-C_hat1)^2)/G)
  iter<-iter+1
 }
 t_hat<-alpha%*%t(C_hat)
 hg<-apply(t_hat, 1, which.max)
 return(list(group=hg, center=C_hat1, iters=iter))
 
}


#根据分组中心及分组情况利用坐标梯度下降法估计出 每个节点的潜在变量 alpha、分组结果g_hat、组中心C_hat
cos_embedding<-function(A, G, p, mu , lambda, gamma, ep, iter_max){
 
 # #initial values
 n<-nrow(A) 
 pho<-sum(A)/(n*(n-1))
 pAdj<-A+as.numeric(0.25*pho)*matrix(1,nrow=n,ncol=n) #purmutated Adj  
 #compute the Laplace Matrix
 D <- apply(pAdj,1,sum)^(-0.5)
 D[D==Inf] <- 0
 H <- diag(D)
 L <- H%*%pAdj%*%H 
 
 U<-eigen(L)
 alpha<-U$vectors[,1:p]
 
 #alpha<- alpha/sqrt(apply(alpha^2,1,sum))# 对特征向量进行 row-normalize
 
 # fit<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
 # g1<-fit$cluster#初始聚类分组
 # C_hat<-fit$centers#初始聚类中心
 # 
 # fit<-cos_cluster(alpha, C_hat=C_hat, ep=ep, iter_max=100)
 # C_hat<-fit$center
 # g_hat<-fit$group
 
 fitk<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
 C_hatk<- fitk$centers #kmeans初始分组中心坐标
 g_k<-fitk$cluster #kmeans初始聚类分组
 
 fitsk<-skmeans(alpha,G) #利用skmeans计算初始聚类分组及聚类中心
 C_hatsk<-fitsk$prototypes  #skmeans初始分组中心坐标
 g_sk<-fitsk$cluster #skmeans初始聚类分组
 
 C_hat<- C_hatk
 g_hat<- g_k
 
 
 iter=0
 delta=1
 while(iter<iter_max && delta>ep){
  
  alpha1<-alpha
  C_hat1<-C_hat
  mu1<-mu
  g_hat1<-g_hat
  like1<-coslike_hood(alpha1, C_hat1, g_hat1, A, mu, lambda=lambda)
  
  fit<-cosl_partial(alpha1, C_hat1, g_hat1, A, mu, lambda=lambda) 
  mu<-mu1-gamma*fit$p_mu 
  alpha<-alpha1-gamma*fit$p_alpha 
  #alpha<-up_alpha(alpha1, C_hat1, g_hat1, A, mu, gamma, lambda)
  
  fit<-cos_cluster(alpha, C_hat=C_hat1, ep=ep, iter_max=100)
  C_hat<-fit$center
  g_hat<-fit$group
  like0<-coslike_hood(alpha, C_hat, g_hat, A, mu, lambda=lambda)
  
  #如果新的似然大于旧的，则减小步长(gamma)再计算
  if(like0>like1){
   gamma<-gamma/2
   alpha<-alpha1
   C_hat<-C_hat1
   mu<-mu1
   g_hat<-g_hat1
  } 
  #如果新的似然小于旧的,寻找负对数似然的最小值(对数似然的最大值)
  if(like0<=like1){
   delta<-(like1-like0)/abs(like1)
   iter<-iter+1
  }
  
 }
 g_hat<-unname(g_hat)
 return(list(skgroup=g_sk,kgroup=g_k, group=g_hat, embedding= alpha, mu_hat=mu, center=C_hat, likely=like0, Gamma=gamma, iters=iter))
}

#cos方法不经过normalize并且用kmeans结果作为初值 miscluster=0.0596
#sq方法不经过normalize并且用kmeans结果作为初值 miscluster=0.0731

##### Data input and pre-processing #####
setwd("C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/BMP")
foo=read.table("BritishMPcomm.txt",header=T)
mat=read.table("politicsuk-retweets.mtx",header=F)
#mat[1,]  # no of nodes, and no of edges (all 5 communities)
mat = mat[-1,]
A = matrix(0,nrow(foo),nrow(foo))
for (k in 1:nrow(mat)){ #adjacency matrix with 5 parties
 edgei = mat[k,1]
 edgej = mat[k,2]
 i = which(foo[,1]==edgei)
 j = which(foo[,1]==edgej)
 A[i,j]<-1;A[j,i]<-1
}
b.true <- foo[,2]
table(b.true) # membership sizes of the 5 parties
nodes1 = which(b.true=="conservative"); nodes2 = which(b.true=="labour"); nodes3 = which(b.true=="libdem")
nodes4 = which(b.true=="other"); nodes5 = which(b.true=="snp")
nodes = c(nodes1,nodes2,nodes3)
A <- A[nodes,nodes] #
unconnected<- which(rowSums(A)==0)
A<- A[-unconnected,-unconnected]
g<- b.true[nodes]
g<- g[-unconnected]
g[g=="conservative"]<- 1 #152个 1-152
g[g=="labour"]<- 2 #178个 153-330
g[g=="libdem"]<- 3 #39个 331-369
g<- as.numeric(g)


n<- nrow(A)
pho<-sum(A)/(n*(n-1))
mu<-log(pho/(1-pho))



#最优lambda=0.01, gamma=1
BMP_cos_results<- cos_embedding(A, G=3, p=3, mu=mu, lambda=0.01, gamma=1, ep=1e-7, iter_max=100)
cosmis<- miscluster(BMP_cos_results$group,g,3)
miscluster(BMP_cos_results$kgroup,g,3)
cos<- BMP_cos_results$center %*% t(BMP_cos_results$center)

# sum(BMP_cos_results$group[1:152]!=1)
# sum(BMP_cos_results$group[153:330]!=2)
# sum(BMP_cos_results$group[331:369]!=3)

write.csv(BMP_cos_results$embedding,file = "cosembedding.csv")
write.csv(BMP_cos_results$center,file = "coscenter.csv")
write.csv(g,file = "truegroup.csv")
write.csv(BMP_cos_results$kgroup,file = "kgroup.csv")
write.csv(BMP_cos_results$group,file = "cosgroup.csv")




mycolors <- c('green', 'blue', 'red')
color<- matrix(0,n,1)
color[1:ng1]<- mycolors[1]
color[(ng1+1):(ng1+ng2)]<- mycolors[2]
color[(ng1+ng2+1):(ng1+ng2+ng3)]<- mycolors[3]
library(scatterplot3d)
library(rgl)
plot3d(BMP_cos_results$embedding,col=color)
plot3d(BMP_cos_results$center,size=6,col="blue",add = TRUE)

#最优lambda=0.008
BMP_sq_results<- sq_embedding(A, G=3, p=3, mu=mu, lambda=0.008, gamma=3, ep=1e-7, iter_max=100)
sqmis<- miscluster(BMP_sq_results$group,g,3)
miscluster(BMP_sq_results$kgroup,g,3)


# sum(BMP_sq_results$group[1:152]!=3)
# sum(BMP_sq_results$group[153:330]!=2)
# sum(BMP_sq_results$group[331:369]!=1)

write.csv(BMP_sq_results$embedding,file = "sqembedding.csv")
write.csv(BMP_sq_results$center,file = "sqcenter.csv")
write.csv(BMP_sq_results$group,file = "sqgroup.csv")

library(rgl)
plot3d(BMP_sq_results$embedding,col=color)
plot3d(BMP_sq_results$center,size=6,col="blue",add = TRUE)