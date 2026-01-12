rm(list = ls())
library(igraph)
library(skmeans)

#初始embedding向量未经过row-normalize

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




#####美国两个党派博客数据 
#1-586 liberal; 587-1222 conservative
library(igraph)
setwd("C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/USblog")
G = read.graph("polblogs.gml", format = "gml")  # read the graph
G <- simplify(G)  # get rid of multiple edges and loops
G <- as.undirected(G,mode="each") # make undirected
A = get.adjacency(G,sparse=FALSE) 
foo <- which(clusters(G)$membership==1) # node labels for the largest connected component
A <- A[foo,foo]
A[A>1]<-1
aff <- get.vertex.attribute(G, "value", index=V(G)) # 0-liberal/ 1-conservative
aff <- aff[foo]
b.true <- aff+1


###################3
n<- nrow(A)
pho<-sum(A)/(n*(n-1))
mu<-log(pho/(1-pho))
degree<- rowSums(A)
#lambda=0.02,gamma=2
set.seed(0)
USB_sq_results<- sq_embedding(A, G=2, p=2, mu=mu, lambda=0.02, gamma=2, ep=1e-7, iter_max=100)

miscluster(USB_sq_results$kgroup,b.true,2)
USB_sq_mis<-miscluster(USB_sq_results$group,b.true,2)
USB_sq_mis

sum(USB_sq_results$group[1:586]!=1)
sum(USB_sq_results$group[587:1222]!=2)

#lambda=0.001,gamma=2
USB_cos_results<- cos_embedding(A, G=2, p=2, mu=mu, lambda=0.001, gamma=2, ep=1e-7, iter_max=100)

miscluster(USB_cos_results$skgroup,b.true,2)
miscluster(USB_cos_results$kgroup,b.true,2)
USB_cos_mis<- miscluster(USB_cos_results$group,b.true,2)
USB_cos_mis

sum(USB_cos_results$group[1:586]!=1)
sum(USB_cos_results$group[587:1222]!=2)

color<- NULL
color[1:586]<- 'green'
color[587:1222]<- 'blue'

write.csv(USB_cos_results$center,file = "coscenter.csv")

plot(USB_cos_results$embedding,col=color,xlab ='alpha[,1]',ylab = 'alpha[,2]')

#不要坐标轴
plot(USB_cos_results$embedding,col=color,xlab ='alpha[,1]',ylab = 'alpha[,2]',xaxt='n',yaxt='n')
axis(side=1,at=c(-0.08,-0.06,-0.04,-0.02,0),labels=c(-0.08,-0.06,-0.04,-0.02,0))
axis(side=2,at=c(-0.10,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06),labels=c(-0.10,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06))

plot(USB_sq_results$embedding,col=color,xlab ='alpha[,1]',ylab = 'alpha[,2]')

plot(USB_cos_results$center,pch=19,col="black")