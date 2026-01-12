block_N <- function(n,G){  #grouping
  b=NULL
  for(i in 1:G)
    if(i!=G){b[((i-1)*round(n/G)+1):(i*round(n/G))]=i} 
  else
  {b[((i-1)*round(n/G)+1):n]=i}
  return(b)
}



#generate A and E using SBM G=2
get_A<- function(n,b,rho){        #(Stochastic Block Model)
  A=E=matrix(0,n,n)
  for(i in 2:n){
    for(j in 1:(i-1)){
      index1=b[i]
      index2=b[j]
      A[i,j]=rbinom(1,1,rho[index1,index2]) 
      A[j,i]=A[i,j]
    }
  }
  return(A)
}

like_hood<-function(alpha, C_hat, hg, A, mu, lambda){
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

l_partial<-function(alpha, C_hat, hg, A, mu, lambda){
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
                                                   -(alpha_t/(2*B_d^(3/2)))*alpha )
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


cos_cluster<-function(alpha, C_hat, ep, iter_max){
  G<-nrow(C_hat)
  alpha<-t(apply(alpha,1,function(x){x/sqrt(sum(x^2))}))
  iter<-0
  delta<-1
  
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
fit<-kmeans(alpha,G)
g1<-fit$cluster
C_hat<-fit$centers
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
like1<-like_hood(alpha1, C_hat1, g_hat1, A, mu, lambda=lambda)

fit<-l_partial(alpha1, C_hat1, g_hat1, A, mu, lambda=lambda) 
mu<-mu1-gamma*fit$p_mu 
alpha<-alpha1-gamma*fit$p_alpha 
#alpha<-up_alpha(alpha1, C_hat1, g_hat1, A, mu, gamma, lambda)

fit<-cos_cluster(alpha, C_hat=C_hat1, ep=ep, iter_max=100)
C_hat<-fit$center
g_hat<-fit$group
like0<-like_hood(alpha, C_hat, g_hat, A, mu, lambda=lambda)

if(like0>like1){
  gamma<-gamma/2
  alpha<-alpha1
  C_hat<-C_hat1
  mu<-mu1
  g_hat<-g_hat1
} 
if(like0<=like1){
  delta<-(like1-like0)/abs(like1)
  iter<-iter+1
}

}
g_hat<-unname(g_hat)
return(list(sgroup=g1, group=g_hat, embedding= alpha, mu_hat=mu, center=C_hat, likely=like0, Gamma=gamma, iters=iter))
}

##########simulation

n<-100
G<-2
p<-2
gamma<-0.1
rho<-matrix(c(0.3,0.1,0.1,0.3),2,2)
g<-block_N(n,G)
A<-get_A(n,g,rho)
pho<-sum(A)/(n*(n-1))
mu<-log(pho/(1-pho))

cos_embedding(A, G=2, p=2, mu=mu, lambda=0.005, gamma=1, ep=1e-10, iter_max=1000)





