SBM_A_G <- function(n,G,rho){  #grouping #将n个节点均分成G个组
  b=NULL
  for(i in 1:G)
    if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i} 
  else
  {b[((i-1)*round(n/G)+1):n]=i}
  
  A=matrix(0,n,n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      index1=b[i]
      index2=b[j]
      A[i,j]=rbinom(1,1,rho[index1,index2]) #以rho矩阵的第index1行 和 index2列 元素的概率为0-1分布的参数
      A[j,i]=A[i,j]
    }
  }
  A<- A[which(rowSums(A)!=0),which(rowSums(A)!=0)]
  b<- b[which(rowSums(A)!=0)]
  re<- list(A,b)
  names(re)<- c('A','b')
  return(re)
}


DCSBM_A_G<- function(n,G,rho,theta){        #(DCSBM) #rho是分组参数，theta是节点度的参数
 
 b=NULL
 for(i in 1:G)
  if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i} 
 else
 {b[((i-1)*round(n/G)+1):n]=i}
 
 A=matrix(0,n,n)
 for(i in 2:n){
  for(j in 1:(i-1)){
   index1=b[i]
   index2=b[j]
   pij<- theta[i]*rho[index1,index2]*theta[j]
   pij<- pij*(pij<=1)+1*(pij>1)
   A[i,j]=rbinom(1,1,pij) 
   A[j,i]=A[i,j]
  }
 }
 A<- A[which(rowSums(A)!=0),which(rowSums(A)!=0)]
 b<- b[which(rowSums(A)!=0)]
 re<- list(A,b)
 names(re)<- c('A','b')
 return(re)
}


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
  alpha<-t(apply(alpha,1,function(x){x/sqrt(sum(x^2))})) #标准化潜在向量alpha
  iter<-0 #迭代次数
  delta<-1 #前后迭代得到的似然函数之差
  
  while (iter<iter_max && delta>ep) {
    C_hat1<-C_hat
    
    t_hat<-alpha%*%t(C_hat1)
    hg<-apply(t_hat, 1, which.max)
    
    for (g in 1:G) {
      y<-colSums(alpha[hg==g, ])
      C_hat[g,]<-y/sqrt(sum(y^2))
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

star_time<-Sys.time()

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
fit<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
g1<-fit$cluster#初始聚类分组
C_hat<-fit$centers#初始聚类中心
fit<-cos_cluster(alpha, C_hat=C_hat, ep=ep, iter_max=100)
C_hat<-fit$center
g_hat<-fit$group

#mu<-log(pho/(1-pho))

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
return(list(sgroup=g1, group=g_hat, embedding= alpha, mu_hat=mu, center=C_hat, likely=like0, Gamma=gamma, iters=iter))
}

##########simulation
n<-400
G<-3
rho<-matrix(c(0.2,0.05,0.05,0.05,0.2,0.05,0.05,0.05,0.2),3,3)
#rho<-matrix(c(0.2,0.05,0.05,0.2),2,2)



#######SBM######
net<-SBM_A_G(n,G,rho)


####DCSBM####
theta<- matrix(sample(c(0.2,0.6,1),size=n,replace=T,prob=c(0.5,0.3,0.2)),n,1)
#theta[1]<- theta[2]<-  theta[n]<- theta[n-1]<- 3
#theta[n/G+1]<- theta[n/G+2]<-
net<-DCSBM_A_G(n,G,rho,theta)


A<- net$A
g<- net$b
ng1<- sum(g==1)
ng2<- sum(g==2)
ng3<- sum(g==3)
n<- nrow(A)
pho<-sum(A)/(n*(n-1))
mu<-log(pho/(1-pho))
degree<- rowSums(A)


simres1<- cos_embedding(A, G=2, p=2, mu=mu, lambda=0.005, gamma=1, ep=1e-8, iter_max=100)

sum(simres1$sgroup[1:586]!=2)
sum(simres1$sgroup[587:1222]!=1)
sum(simres1$sgroup!=simres1$group)

g[1:ng1]<- simres1$group[1]
g[(ng1+1):(ng1+ng2)]<- simres1$group[ng1+ng2]
g[(ng1+ng2+1):(ng1+ng2+ng3)]<- simres1$group[ng1+ng2+ng3]
sum(simres1$group!=g)


mycolors <- c('palegreen', 'skyblue', 'orange','red')
#mycolors <- c('palegreen', 'skyblue','red')
color<- matrix(0,n,1)
color[which(g==1),]<- mycolors[1]
color[which(g==2),]<- mycolors[2]
color[which(g==3),]<- mycolors[3]
#color[which.max(degree),]<- mycolors[4]
library(rgl)
plot3d(simres1$embedding,size=4,col=color)
plot3d(simres1$center,size=6,col="black",add = TRUE)


