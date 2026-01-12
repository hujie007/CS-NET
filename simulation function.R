#############################余弦惩罚项##################################
library(skmeans)
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
 #alpha<-t(apply(alpha,1,function(x){x/sqrt(sum(x^2))})) #标准化潜在向量alpha (数值模拟不标准化会好很多)
 iter<-0 #迭代次数
 delta<-1 #前后迭代得到的似然函数之差

 while (iter<iter_max && delta>ep) {
  C_hat1<-C_hat

  t_hat<-alpha%*%t(C_hat1)
  hg<-apply(t_hat, 1, which.max)

  for (g in 1:G) {
    # y<-colSums(alpha[hg==g, ])
    # C_hat[g,]<-y/sqrt(sum(y^2))

    y<- alpha[hg==g, ]
    h<- colSums(y/sqrt(diag(y %*% t(y)))) #h是一个向量
    C_hat[g,]<- h/sqrt(sum(h^2))
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
 
 #alpha<- alpha/sqrt(apply(alpha^2,1,sum)) # 对特征向量进行 row-normalize
 
 # 分别对alpha做kmeans和skmeans
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
  like1<-coslike_hood(alpha1, C_hat1, g_hat1, A, mu1, lambda=lambda)
  
  fit<-cosl_partial(alpha1, C_hat1, g_hat1, A, mu1, lambda=lambda) 
  mu<-mu1-gamma*fit$p_mu #更新了mu
  alpha<-alpha1-gamma*fit$p_alpha  #更新了alpha
  #alpha<-up_alpha(alpha1, C_hat1, g_hat1, A, mu, gamma, lambda)
  
  fit<- cos_cluster(alpha, C_hat1, ep, iter_max)
  C_hat<-fit$center
  g_hat<-fit$group
  
  
  # fit<-skmeans(alpha,G) #更新过后的alpha来skmeans聚类
  # C_hat<-fit$prototypes
  # g_hat<-fit$cluster
  
  like0<-coslike_hood(alpha, C_hat, g_hat, A, mu, lambda=lambda) #更新过后的alpha等参数计算新的似然函数
  
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
 return(list(kgroup=g_k,skgroup=g_sk, group=g_hat, embedding= embedding, mu_hat=mu_hat, center=C_hat, likely=like0, Gamma=gamma, iters=iter))
}


#为了保证运算速度以及方法差异，这里利用skmeans对输入的embedding聚类作为初始值而非kmeans
simcos_embedding<-function(A, G, p, mu ,embedding, lambda, gamma, ep, iter_max){
 
 #initial values
 alpha<-embedding #直接将随机生成的embedding向量作为初始向量
 
 #alpha<- alpha/sqrt(apply(alpha^2,1,sum)) # 对特征向量进行 row-normalize
 
 # fitk<-kmeans(alpha,G) #先利用K-Means计算初始聚类分组及聚类中心
 # g_k<-fitk$cluster #利用kmeans得到的聚类分组
 # C_hatk<-fitk$centers #利用kmeans得到的聚类中心
 
 fitsk<-skmeans(alpha,G) #利用skmeans计算初始聚类分组及聚类中心
 C_hatsk<-fitsk$prototypes  #skmeans初始分组中心坐标
 g_sk<-fitsk$cluster #skmeans初始聚类分组
 
 C_hat<- C_hatsk#一次迭代得到的初始聚类中心
 g_hat<-g_sk #一次迭代得到的初始聚类分组
 
 #mu<-log(pho/(1-pho))
 
 iter=0
 delta=1
 while(iter<iter_max && delta>ep){
  
  alpha1<-alpha
  C_hat1<-C_hat
  mu1<-mu
  g_hat1<-g_hat
  like1<-coslike_hood(alpha1, C_hat1, g_hat1, A, mu1, lambda=lambda)
  
  fit<-cosl_partial(alpha1, C_hat1, g_hat1, A, mu1, lambda=lambda) 
  mu<-mu1-gamma*fit$p_mu 
  alpha<-alpha1-gamma*fit$p_alpha 
  #alpha<-up_alpha(alpha1, C_hat1, g_hat1, A, mu, gamma, lambda)
  
  # fit<-skmeans(alpha,G) #更新过后的alpha来skmeans聚类
  # C_hat<-fit$prototypes
  # g_hat<-fit$cluster
  
  fit<- cos_cluster(alpha, C_hat1, ep, iter_max)
  C_hat<-fit$center
  g_hat<-fit$group
  # 
  like0<-coslike_hood(alpha, C_hat, g_hat, A, mu, lambda=lambda) #更新过后的alpha等参数计算新的似然函数
  
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
 return(list(skgroup=g_sk, group=g_hat, embedding= embedding, mu_hat=mu_hat, center=C_hat, likely=like0, Gamma=gamma, iters=iter))
}



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
 n<- nrow(alpha)
 hg<- numeric(n)
 alpha<-t(apply(alpha,1,function(x){x/sqrt(sum(x^2))})) #标准化潜在向量alpha
 iter<-0 #迭代次数
 delta<-1 #前后迭代得到的似然函数之差
 
 while (iter<iter_max && delta>ep) {
  C_hat1<-C_hat
  # 分别更新每个点的组
  dist<- as.matrix(dist(rbind(alpha,C_hat1)))[(n+1):(n+G),1:n]
  hg<- apply(dist,2,which.min)
  # 更新组中心
  for (g in 1:G) {
   # y<-colSums(alpha[hg==g, ]) #属于g组的这些点的坐标中心作为新的组中心
   # C_hat[g,]<-y/sqrt(sum(y^2)) 
   C_hat[g,]<- colSums(alpha[hg==g, ])/sum(hg==g)
  } 
  delta<-sqrt(sum((C_hat-C_hat1)^2)/G) #计算新旧中心的距离，再据此判断是否继续迭代(直至中心变化很小)
  iter<-iter+1
 }
 # for (i in 1:n){
 #  dist<- rowSums((matrix(1,G,1)%*%alpha[i,]-C_hat1)^2)
 #  hg[i]<- which.min(dist)
 # }
 
 dist<- as.matrix(dist(rbind(alpha,C_hat1)))[(n+1):(n+G),1:n]
 hg<- apply(dist,2,which.min)
 
 return(list(group=hg, center=C_hat1, iters=iter))
}




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
 
 alpha<- alpha/sqrt(apply(alpha^2,1,sum)) # 对特征向量进行 row-normalize
 
 # 对alpha做kmeans
 fitk<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
 C_hatk<- fitk$centers #kmeans初始分组中心坐标
 g_k<-fitk$cluster #kmeans初始聚类分组
 
 
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





simsq_embedding<-function(A, G, p, mu, embedding, lambda, gamma, ep, iter_max){
 
 #initial values
 alpha<-embedding #直接将随机生成的embedding向量作为初始向量
 
 #alpha<- alpha/sqrt(apply(alpha^2,1,sum)) # 对特征向量进行 row-normalize
 
 fit<-kmeans(alpha,G) #先利用K-Means计算初始聚类分组及聚类中心
 g_k<-fit$cluster #利用kmeans得到的聚类分组
 C_hatk<-fit$centers #利用kmeans得到的聚类中心
 # fit<-sq_cluster(alpha, C_hat=C_hat, ep=ep, iter_max=100)#先用一次cos惩罚函数
 # C_hat<-fit$center #一次迭代得到的初始聚类中心
 # g_hat<-fit$group #一次迭代得到的初始聚类分组
 
 C_hat<- C_hatk
 g_hat<- g_k
 #mu<-log(pho/(1-pho))
 
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
  
  # fit<-sq_cluster(alpha, C_hat=C_hat1, ep=ep, iter_max=100)
  # C_hat<-fit$center
  # g_hat<-fit$group
  # like0<-sqlike_hood(alpha, C_hat, g_hat, A, mu, lambda=lambda)
  
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


#rho<-matrix(c(0.5,0.1,0.1,0.1,0.5,0.1,0.1,0.1,0.5),3,3)
#rho<-matrix(c(0.3,0.1,0.3,0.1),2,2)
SBM_A_G <- function(n,G,rho){  #grouping #将n个节点均分成G个组 pho是group之间的连接参数
 b=NULL
 
 for(i in 1:G)
  if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i}
 else
 {b[((i-1)*round(n/G)+1):n]=i}
 
# b<- sample(1:2,n,replace=T,prob=c(0.4,0.6))
 
 
 A=P=matrix(0,n,n)
 for(i in 2:n){
  for(j in 1:(i-1)){
   index1=b[i]
   index2=b[j]
   P[i,j]=rho[index1,index2]
   A[i,j]=rbinom(1,1,P[i,j]) #以rho矩阵的第index1行 和 index2列 元素的概率为0-1分布的参数
   A[j,i]=A[i,j]
   P[j,i]=P[i,j]
  }
 }
 Z<- which(rowSums(A)!=0)
 A<- A[Z,Z]
 P<- P[Z,Z]
 b<- b[Z]
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 re<- list(n,A,b,mu,P)
 names(re)<- c('n','A','b','mu','P')
 return(re)
}



# 由Example 1和Example 3生成的数据(根据SBM和DCBM的表达式生成)
tuninglambda<- function(A,G,p,g,mu,lambda_n,gamma,ep,iter_max){
 cosl<- lapply(lambda_n,cos_embedding,A=A,G=G,p=p,mu=mu,gamma=gamma,ep=ep,iter_max=iter_max)
 sql<- lapply(lambda_n,sq_embedding,A=A,G=G,p=p,mu=mu,gamma=gamma,ep=ep,iter_max=iter_max)
 crit<- matrix(0,length(lambda_n),2)
 for (i in 1:length(lambda_n)){
  crit[i,1]<- miscluster(cosl[[i]]$group,g,G)$misrate
  crit[i,2]<- miscluster(sql[[i]]$group,g,G)$misrate
 }
 findmin<- apply(crit,2,which.min)
 re<- list(cosl[[findmin[1]]],sql[[findmin[2]]],lambda_n[findmin[1]],lambda_n[findmin[2]])
 names(re)<- c('cosres','sqres','lambdaofcos','lambdaofsq')
 return(re)
}


# 由Example 2和Example 4生成的数据(由embedding向量生成的数据)
simtuninglambda<- function(A, G, p, g, mu, embedding, lambda_n, gamma, ep, iter_max){
 cosl<- lapply(lambda_n,simcos_embedding,A=A,G=G,p=p,mu=mu,embedding=embedding,gamma=gamma,ep=ep,iter_max=iter_max)
 sql<- lapply(lambda_n,simsq_embedding,A=A,G=G,p=p,mu=mu,embedding=embedding,gamma=gamma,ep=ep,iter_max=iter_max)
 crit<- matrix(0,length(lambda_n),2)
 for (i in 1:length(lambda_n)){
  crit[i,1]<- miscluster(cosl[[i]]$group,g,G)$misrate
  crit[i,2]<- miscluster(sql[[i]]$group,g,G)$misrate
 }
 findmin<- apply(crit,2,which.min)
 re<- list(cosl[[findmin[1]]],sql[[findmin[2]]],lambda_n[findmin[1]],lambda_n[findmin[2]])
 names(re)<- c('cosres','sqres','lambdaofcos','lambdaofsq')
 return(re)
}


#rho<-matrix(c(0.5,0.05,0.05,0.05,0.5,0.05,0.05,0.05,0.5),3,3)
#rho<-matrix(c(0.5,0.05,0.5,0.05),2,2)
# DCSBM_A_G<- function(n,G,rho){        #(DCSBM) #rho是分组参数，theta是节点度的参数
#  theta<- matrix(0,n,1)
#  degreelevel<- seq(1,1.8,0.2)#seq(0.5,1.5,0.3)
#  pro<- seq(1,0.1,length.out=length(degreelevel))
#  for (g in 1:G){
#   theta[((g-1)*floor(n/G)+1):(g*floor(n/G))]<- matrix(sample(degreelevel,size=floor(n/G),replace=T,prob=(1/degreelevel)^2/sum((1/degreelevel)^2)),floor(n/G),1)
#   }
#  b=NULL
#  for(i in 1:G)
#   if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i} 
#  else
#  {b[((i-1)*round(n/G)+1):n]=i}
#  A=P=matrix(0,n,n)
#  for(i in 2:n){
#   for(j in 1:(i-1)){
#    index1=b[i]
#    index2=b[j]
#    P[i,j]<- theta[i]*rho[index1,index2]*theta[j]
#    P[i,j]<- P[i,j]*(P[i,j]<=1)+1*(P[i,j]>1)
#    P[j,i]=P[i,j]
#    A[i,j]=rbinom(1,1,P[i,j]) 
#    A[j,i]=A[i,j]
#   }
#  }
#  Z<- which(rowSums(A)!=0)
#  A<- A[Z,Z]
#  P<- P[Z,Z]
#  b<- b[Z]
#  n<- nrow(A)
#  pho<-sum(A)/(n*(n-1))
#  mu<-log(pho/(1-pho))
#  re<- list(n,A,b,mu,P)
#  names(re)<- c('n','A','b','mu','P')
#  return(re)
# }


# rho <- (1*1e-2)*(matrix(1,G,G)+diag(c(3,3)))
DCSBM_A_G<- function(n,G,rho){        #(DCSBM) #rho是分组参数，theta是节点度的参数
 v <- 4
 theta<-NULL
 dl<- seq(0.3,8,0.1)
 for (g in 1:G){
  theta[((g-1)*floor(n/G)+1):(g*floor(n/G))]<- sample(dl,size = n/G,replace = T,prob = 1/dl^2/sum(1/dl^2))
 }
 
 # for (g in 1:G){
 #  theta[((g-1)*floor(n/G)+1):(g*floor(n/G))]<- sample(c((v*2)/(v+1),(v-1)/(v+1),2/(v+1)),size = n/G,replace = T,prob = c(0.1,0.2,0.7))
 # }
 
 b=NULL
 for(i in 1:G)
  if(i!=G){b[((i-1)*n/G+1):(i*n/G)]=i} 
 else
 {b[((i-1)*n/G+1):n]=i}
 c<- diag(G)
 C<- c[b,]
 Theta<- diag(theta)
 P<- Theta%*%(C%*%rho%*%t(C))%*%Theta
 RM <- matrix(runif(n*n),n,n)
 RM[lower.tri(RM)] <- t(RM)[lower.tri(RM)] #Symmetrization
 A <- (RM < P)+0
 Z<- which(rowSums(A)!=0)
 A<- A[Z,Z]
 P<- P[Z,Z]
 b<- b[Z]
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 re<- list(n,A,b,mu,P)
 names(re)<- c('n','A','b','mu','P')
 return(re)
}




Adj_Generating_DCSBM <- function(n,Clusters_distribution,
                                 Lambda,Theta,
                                 sp = FALSE,del_0d = FALSE){
 num_clusters <- length(Clusters_distribution)
 label_clusters <- 1:num_clusters
 clusters <- sample(label_clusters,size=n,replace=T,
                    prob=Clusters_distribution)
 
 K_I <- diag(num_clusters)
 C_Mat <- K_I[clusters, ]
 
 Theta_Mat <- diag(Theta)
 
 P<- Theta_Mat%*%(C_Mat%*%Lambda%*%t(C_Mat))%*%Theta_Mat #Theta_Mat是每个点的degree参数
 
 RM <- matrix(runif(n*n),n,n)
 RM[lower.tri(RM)] <- t(RM)[lower.tri(RM)] #Symmetrization
 
 A <- (RM < P)+0
 
 if(del_0d){
  Ind_0 <- apply(A,1,sum)!=0
  A <- A[Ind_0,Ind_0]
  clusters <- clusters[Ind_0] 
 }
 
 if(sp){
  library(Matrix)
  A <- Matrix(A,sparse=T)
 }
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 re<- list(n,A,b,mu,P)
 names(re)<- c('n','A','b','mu','P')
 return(re)
}





library(MASS)
Embedding_A_G<- function(n,G,p,sigma){# embedding文章中的example生成数据方式
 b<- NULL
 #b<- numeric(n)
 mu<- -log(n)/n
 embedding<- matrix(0,n,p)
 # taudis=0
 # taucos=1
 #随机生成模为1的若干个中心(随机生成每个分量服从N(0,1)的p维向量，并将它标准化即得到模为1的向量)
 #while( taudis < 0.5 ) {
  center<- matrix(mvrnorm(G*p,0,1),G,p)
  center1<- center/sqrt(apply(center^2,1,sum)) #标准化每一个中心坐标
 #  inpro<- center1 %*% t(center1)
 #  #taucos<- max(inpro[upper.tri(inpro)]) #各个中心角度的余弦最大值(夹角的最小值)
 #  dist<- matrix(0, G, G)
 #  for (i in 1:G){
 #   dist[i,]<- sqrt(rowSums((matrix(1,G,1) %*% center[i,] - center)^2))
 #  }
 #  taudis<- min(dist[upper.tri(dist)]) #各个中心距离的最小值
 # }
 center<- unname(center)
 for (i in 1:n){
  b[i]<- floor(G*(i-1)/n+1)
 }
 for (k in 1:G){
  for (s in 1:p){
   embedding[b==k,s]<- rnorm(sum(b==k),center[k,s],sigma)
  }
 }
 #embedding<- embedding/sqrt(apply(embedding^2,1,sum))
 #各组节点所对应的embedding向量的第k个元素 的分布 共用同一个均值
 P<- (1+exp(-(mu+embedding %*% t(embedding))))^(-1)
 A<- matrix(0,n,n)
 A1<- rbinom(n*(n-1)/2,1,P[upper.tri(P)])
 A[upper.tri(A)]<- A1
 A<- A+t(A)
 Z<- which(rowSums(A)!=0)
 A<- A[Z,Z]
 P<- P[Z,Z]
 b<- b[Z]
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 re<- list(n,A,b,embedding,P,mu)
 names(re)<- c('n','A','b','embedding','P','mu')
 return(re)
}



#随机生成若干个角度中心，以它们为中心在与其有一定的夹角范围内随机生成n个节点进行分组
library(MASS)
COS_A_G<- function(n,G,p,sigma){
 b=NULL
 tau=1
 #随机生成模为1的若干个中心(随机生成每个分量服从N(0,1)的p维向量，并将它标准化即得到模为1的向量)
 while( tau > sqrt(1)/2 ) { #中心向量角度设定是否还需要？？
  center<- matrix(mvrnorm(G*p,0,1),G,p)#随机生成G个p维中心
  center<- center/sqrt(apply(center^2,1,sum)) #标准化每一个中心坐标
  inpro<- center %*% t(center)
  tau<- max(inpro[upper.tri(inpro)])
 }
 center<- unname(center)
 embedding<- matrix(0,n,p)
 for (i in 1:n){
  b[i]<- floor(G*(i-1)/n+1)
 }
 for (g in 1:G){
  ng<- sum(b==g)
  x<- mvrnorm(ng,center[g,],sigma*diag(p)) 
  embedding[b==g,]<- x/sqrt(apply(x^2,1,sum))*(abs(rnorm(ng,0,1))+0.1)
 }
 mu<- 0#-sqrt(log(n))
 P<- (1+exp(-(mu+embedding %*% t(embedding))))^(-1)
 A<- matrix(0,n,n)
 A1<- rbinom(n*(n-1)/2,1,P[upper.tri(P)])
 A[upper.tri(A)]<- A1
 A<- A+t(A)
 Z<- which(rowSums(A)!=0)
 A<- A[Z,Z]
 P<- P[Z,Z]
 b<- b[Z]
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 re<- list(n,A,P,b,embedding,mu)
 names(re)<- c('n','A','P','b','embedding','mu')
 return(re)
}




#针对G=2，p=2的场景

G2p2COS<- function(n,sigma){
 b=NULL
 z11<- rnorm(1,1,0.01)
 z12<- rnorm(1,0,0.01)
 z21<- rnorm(1,0,0.01)
 z22<- rnorm(1,1,0.01)
 center<- matrix(c(z11,z21,z12,z22),2,2)
 center<- center/sqrt(apply(center^2,1,sum))
 embedding<- matrix(0,n,2)
 for (i in 1:n){
  b[i]<- floor(2*(i-1)/n+1)
 }
 dl<- seq(1,2,0.5)
 for (g in 1:2){
  ng<- sum(b==g)
  x<- mvrnorm(ng,center[g,],sigma*diag(2))
  embedding[b==g,]<- x/sqrt(apply(x^2,1,sum))*(sample(dl,ng,replace=T,prob=1/dl^2/sum(1/dl^2)))
 }#
 mu<- -sqrt(log(n))
 P<- (1+exp(-(mu+embedding %*% t(embedding))))^(-1)
 A<- matrix(0,n,n)
 A1<- rbinom(n*(n-1)/2,1,P[upper.tri(P)])
 A[upper.tri(A)]<- A1
 A<- A+t(A)
 Z<- which(rowSums(A)!=0)
 A<- A[Z,Z]
 P<- P[Z,Z]
 b<- b[Z]
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 re<- list(n,A,P,b,embedding,mu)
 names(re)<- c('n','A','P','b','embedding','mu')
 return(re)
}



G2p2COS<- function(n,sigma){
 b=NULL
 omg<- c(3*pi/4,pi/6)
 embedding<- matrix(0,n,2)
 for (i in 1:n){
  b[i]<- floor(2*(i-1)/n+1)
 }
 d1<- c(2,8)
 d2<- c(1,2)
 scale1<- sample(d1,n/2,replace=T,prob=c(0.9,0.1))
 scale2<- sample(d2,n/2,replace=T,prob=c(0.9,0.1))
 # for (g in 1:2){
 #  ng<- sum(b==g)
 #  embedding[b==g,1]<- cos(omg[g]+rnorm(ng,0,sigma))*scale
 #  embedding[b==g,2]<- sin(omg[g]+rnorm(ng,0,sigma))*scale
 # }
 
 embedding[1:(n/2),1]<- cos(rnorm(n/2,omg[1],sigma))*scale1
 embedding[1:(n/2),2]<- sin(rnorm(n/2,omg[1],sigma))*scale1
 
 embedding[(n/2+1):n,1]<- cos(rnorm(n/2,omg[2],sigma))*scale2
 embedding[(n/2+1):n,2]<- sin(rnorm(n/2,omg[2],sigma))*scale2
 
 # for (g in 1:2){
 #  ng<- sum(b==g)
 #  embedding[b==g,1]<- cos(rnorm(ng,omg[g],sigma))*scale
 #  embedding[b==g,2]<- sin(rnorm(ng,omg[g],sigma))*scale
 # }

 mu<- -sqrt(log(100))
 
 P<- (1+exp(-(mu+embedding %*% t(embedding))))^(-1)
 A<- matrix(0,n,n)
 A1<- rbinom(n*(n-1)/2,1,P[upper.tri(P)])
 A[upper.tri(A)]<- A1
 A<- A+t(A)
 Z<- which(rowSums(A)!=0)
 A<- A[Z,Z]
 P<- P[Z,Z]
 b<- b[Z]
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 
 embedding[1:(n/2),1]<- cos(rnorm(n/2,omg[1],1/log(sqrt(n))) )*scale1
 embedding[1:(n/2),2]<- sin(rnorm(n/2,omg[1],1/log(sqrt(n))) )*scale1
 
 embedding[(n/2+1):n,1]<- cos(rnorm(n/2,omg[2],1/log(sqrt(n))) )*scale2
 embedding[(n/2+1):n,2]<- sin(rnorm(n/2,omg[2],1/log(sqrt(n))) )*scale2
 
 re<- list(n,A,P,b,embedding,mu)
 names(re)<- c('n','A','P','b','embedding','mu')
 return(re)
}





G2p2SQ<- function(n,sigma){
 b=NULL
 z11<- rnorm(1,0,0.1)
 z12<- rnorm(1,1,0.1)
 z21<- rnorm(1,1,0.1)
 z22<- rnorm(1,0,0.1)
 center<- matrix(c(z11,z21,z12,z22),2,2)
 center/sqrt(apply(center^2,1,sum))
 embedding<- matrix(0,n,2)
 for (i in 1:n){
  b[i]<- floor(2*(i-1)/n+1)
 }
 for (k in 1:G){
  for (s in 1:p){
   embedding[b==k,s]<- rnorm(sum(b==k),center[k,s],sigma)
  }
 }
 mu<- -sqrt(log(n))
 P<- (1+exp(-(mu+embedding %*% t(embedding))))^(-1)
 A<- matrix(0,n,n)
 A1<- rbinom(n*(n-1)/2,1,P[upper.tri(P)])
 A[upper.tri(A)]<- A1
 A<- A+t(A)
 Z<- which(rowSums(A)!=0)
 A<- A[Z,Z]
 P<- P[Z,Z]
 b<- b[Z]
 n<- nrow(A)
 pho<-sum(A)/(n*(n-1))
 mu<-log(pho/(1-pho))
 re<- list(n,A,P,b,embedding,mu)
 names(re)<- c('n','A','P','b','embedding','mu')
 return(re)
}






#######################计算聚类错误比例
zs<- function(x){#寻找一列数中的众数
 return(as.numeric(names(table(x)))[table(x) == max(table(x))])
}

miscluster<- function(x,g,G){#输入的分别是计算分组、真实分组、组的个数
 mis= NULL
 for (k in 1:G){
  mis[k]<- sum(x[g==k]!= zs(x[g==k])) #将对应的计算分组的众数作为该组编号
 }
 mistot<-  sum(mis)
 misrate<- mistot/length(x) #错分个数的比例
 re<- list(mistot,misrate)
 names(re)<- c('mistot','misrate')
 return(re)
} 
#######################




SBM_cos_sp_sq<- function(n,G,p,rho,lambda_n,gamma,ep,iter_max){
 net<- SBM_A_G(n,G,rho)
 n<- net$n
 A<- net$A
 g<- net$b
 mu<- net$mu
 P<- net$P
 #embedding<- net$embedding
 mu<- net$mu
 re<- tuninglambda(A,G,p,g,mu,lambda_n,gamma,ep,iter_max)
 cosres<- re$cosres
 sqres<- re$sqres
 
 cosP_hat<- 1/(1+exp(-cosres$mu_hat-cosres$embedding %*% t(cosres$embedding)))
 coserr_P<- 2*sum((P[upper.tri(P)]-cosP_hat[upper.tri(cosP_hat)])^2)/(n*(n-1))
 
 
 sqP_hat<- 1/(1+exp(-sqres$mu_hat-sqres$embedding %*% t(sqres$embedding)))
 sqerr_P<- 2*sum((P[upper.tri(P)]-sqP_hat[upper.tri(sqP_hat)])^2)/(n*(n-1))
 
 ratio_cosmis<- miscluster(cosres$group,g,G)$misrate #计算cos模型聚类错误的比例
 ratio_spmis<- miscluster(cosres$kgroup,g,G)$misrate #计算谱聚类聚类错误的比例
 ratio_sqmis<- miscluster(sqres$group,g,G)$misrate #计算sq模型聚类错误的比例
 
 ratio_mis<- list(ratio_cosmis,ratio_sqmis,coserr_P,sqerr_P)
 names(ratio_mis)<- c('cosmisrate','sqmisrate','coserr_P','sqerr_P')
 return(ratio_mis)
}


SBM_cos_km_sq<- function(n,G,p,rho,lambda_n,gamma,ep,iter_max){
 net<- G2p2SQ(n,sigma)
 n<- net$n
 A<- net$A
 g<- net$b
 mu<- net$mu
 P<- net$P
 #embedding<- net$embedding
 mu<- net$mu
 re<- tuninglambda(A,G,p,g,mu,lambda_n,gamma,ep,iter_max)
 cosres<- re$cosres
 sqres<- re$sqres
 
 cosP_hat<- 1/(1+exp(-cosres$mu_hat-cosres$embedding %*% t(cosres$embedding)))
 coserr_P<- 2*sum((P[upper.tri(P)]-cosP_hat[upper.tri(cosP_hat)])^2)/(n*(n-1))
 
 
 sqP_hat<- 1/(1+exp(-sqres$mu_hat-sqres$embedding %*% t(sqres$embedding)))
 sqerr_P<- 2*sum((P[upper.tri(P)]-sqP_hat[upper.tri(sqP_hat)])^2)/(n*(n-1))
 
 ratio_cosmis<- miscluster(cosres$group,g,G)$misrate #计算cos模型聚类错误的比例
 ratio_spmis<- miscluster(cosres$kgroup,g,G)$misrate #计算谱聚类聚类错误的比例
 ratio_sqmis<- miscluster(sqres$group,g,G)$misrate #计算sq模型聚类错误的比例
 
 ratio_mis<- list(ratio_cosmis,ratio_sqmis,coserr_P,sqerr_P)
 names(ratio_mis)<- c('cosmisrate','sqmisrate','coserr_P','sqerr_P')
 return(ratio_mis)
}



SBM_cos_km_sq<- function(n,G,p,sigma,lambda_n,gamma,ep,iter_max){
 net<- G2p2SQ(n,sigma)
 n<- net$n
 A<- net$A
 g<- net$b
 embedding<- net$embedding
 mu<- net$mu
 P<- net$P
 re<- simtuninglambda(A,G,p,g,mu,embedding,lambda_n,gamma,ep,iter_max)
 cosres<- re$cosres
 sqres<- re$sqres

 # cosres<- simcos_embedding(A,G,p,mu,embedding,lambda_n,gamma,ep,iter_max)
 # sqres<- simsq_embedding(A,G,p,mu,embedding,lambda_n,gamma,ep,iter_max)

 cosP_hat<- 1/(1+exp(-cosres$mu_hat-cosres$embedding %*% t(cosres$embedding)))
 coserr_P<- 2*sum((P[upper.tri(P)]-cosP_hat[upper.tri(cosP_hat)])^2)/(n*(n-1))


 sqP_hat<- 1/(1+exp(-sqres$mu_hat-sqres$embedding %*% t(sqres$embedding)))
 sqerr_P<- 2*sum((P[upper.tri(P)]-sqP_hat[upper.tri(sqP_hat)])^2)/(n*(n-1))

 ratio_cosmis<- miscluster(cosres$group,g,G)$misrate #计算cos模型聚类错误的比例
 ratio_spmis<- miscluster(cosres$kgroup,g,G)$misrate #计算谱聚类聚类错误的比例
 ratio_sqmis<- miscluster(sqres$group,g,G)$misrate #计算sq模型聚类错误的比例

 ratio_mis<- list(ratio_cosmis,ratio_sqmis,coserr_P,sqerr_P)
 names(ratio_mis)<- c('cosmisrate','sqmisrate','coserr_P','sqerr_P')
 return(ratio_mis)
}



DCBM_cos_sp_sq<- function(n,G,p,rho,lambda_n,gamma,ep,iter_max){
 net<- DCSBM_A_G(n,G,rho)
 n<- net$n
 A<- net$A
 g<- net$b
# embedding<- net$embedding
 mu<- net$mu
 P<- net$P
 re<- tuninglambda(A,G,p,g,mu,lambda_n,gamma,ep,iter_max)
 cosres<- re$cosres
 sqres<- re$sqres

 cosP_hat<- 1/(1+exp(-cosres$mu_hat-cosres$embedding %*% t(cosres$embedding)))
 coserr_P<- 2*sum((P[upper.tri(P)]-cosP_hat[upper.tri(cosP_hat)])^2)/(n*(n-1))
 
 
 sqP_hat<- 1/(1+exp(-sqres$mu_hat-sqres$embedding %*% t(sqres$embedding)))
 sqerr_P<- 2*sum((P[upper.tri(P)]-sqP_hat[upper.tri(sqP_hat)])^2)/(n*(n-1))
 
 ratio_cosmis<- miscluster(cosres$group,g,G)$misrate #计算cos模型聚类错误的比例
 ratio_spmis<- miscluster(cosres$kgroup,g,G)$misrate #计算谱聚类聚类错误的比例
 ratio_sqmis<- miscluster(sqres$group,g,G)$misrate #计算sq模型聚类错误的比例
 
 ratio_mis<- list(ratio_cosmis,ratio_sqmis,coserr_P,sqerr_P)
 names(ratio_mis)<- c('cosmisrate','sqmisrate','coserr_P','sqerr_P')
 return(ratio_mis)
}




DCBM_cos_km_sq<- function(n,G,p,sigma,lambda_n,gamma,ep,iter_max){
 net<- G2p2COS(n,sigma)
 n<- net$n
 A<- net$A
 g<- net$b
 embedding<- net$embedding
 mu<- net$mu
 P<- net$P
 re<- simtuninglambda(A,G,p,g,mu,embedding,lambda_n,gamma,ep,iter_max)
 cosres<- re$cosres
 sqres<- re$sqres
 
 cosP_hat<- 1/(1+exp(-cosres$mu_hat-cosres$embedding %*% t(cosres$embedding)))
 coserr_P<- 2*sum((P[upper.tri(P)]-cosP_hat[upper.tri(cosP_hat)])^2)/(n*(n-1))
 
 
 sqP_hat<- 1/(1+exp(-sqres$mu_hat-sqres$embedding %*% t(sqres$embedding)))
 sqerr_P<- 2*sum((P[upper.tri(P)]-sqP_hat[upper.tri(sqP_hat)])^2)/(n*(n-1))
 
 ratio_cosmis<- miscluster(cosres$group,g,G)$misrate #计算cos模型聚类错误的比例
 ratio_spmis<- miscluster(cosres$kgroup,g,G)$misrate #计算谱聚类聚类错误的比例
 ratio_sqmis<- miscluster(sqres$group,g,G)$misrate #计算sq模型聚类错误的比例
 
 ratio_mis<- list(ratio_cosmis,ratio_sqmis,coserr_P,sqerr_P)
 names(ratio_mis)<- c('cosmisrate','sqmisrate','coserr_P','sqerr_P')
 return(ratio_mis)
}

