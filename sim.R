
repmean1<- repmean2<- repmean3<- repmean4<- matrix(0,4,4)
repsd1<- repsd2<-repsd3<- repsd4<- matrix(0,4,4)

n<- c(100,400,800,1000)
G<- 2 #分组个数
p<- 2 #embedding的维数 (G<=p)
gamma<- 2
ep<- 1e-7
iter_max<- 100
lambda<- 0.0001
lambda_n<- c(0.0001,0.0005,0.001,0.005,0.01,0.02)
sigma<- 0.1
rho<-matrix(c(0.5,0.1,0.5,0.1),2,2)



########################
#用SBM模型模拟
netSBM<- SBM_A_G(n[2],G,rho)
n<- netSBM$n
A<- netSBM$A
g<- netSBM$b
mu<- netSBM$mu
res1<- cos_embedding(A, G, p, mu , lambda, gamma, ep, iter_max)
res2<- sq_embedding(A, G, p, mu , lambda, gamma, ep, iter_max)

miscluster(res1$skgroup,g,G)
miscluster(res1$kgroup,g,G)
miscluster(res1$group,g,G)

miscluster(res2$kgroup,g,G)
miscluster(res2$group,g,G)

########################
#用DCSBM模型模拟
G<-2 #分组个数
p<- 2 #embedding的维数 (G<=p)
lambda<- 0.0003
gamma<- 2
ep<- 1e-7
iter_max<- 100
lambda_n<- seq(0.001,0.005,0.003)
rho <- (1*1e-2)*(matrix(1,G,G)+diag(c(3,5))) 
n<- c(200,400,800,1000) #
netDCBM<- DCSBM_A_G(n[1],G,rho)
n<- netDCBM$n
A<- netDCBM$A
g<- netDCBM$b
mu<- netDCBM$mu
res1<- cos_embedding(A, G, p, mu ,lambda, gamma, ep, iter_max)
res2<- sq_embedding(A, G, p, mu ,lambda, gamma, ep, iter_max)

miscluster(res1$group,g,G)
miscluster(res2$group,g,G)

##################
#用COS生成方式跑
netCOS<- G2p2(n,sigma)
n<- netCOS$n
A<- netCOS$A
g<- netCOS$b
mu<- netCOS$mu
embedding<- netCOS$embedding
res1<- simcos_embedding(A, G, p, mu ,embedding, lambda, gamma, ep, iter_max)
res2<- simsq_embedding(A, G, p, mu ,embedding, lambda, gamma, ep, iter_max)

miscluster(res1$group,g,G)
miscluster(res2$group,g,G)
############################
netEMB<- Embedding_A_G(n,G,p,sigma)
n<- netEMB$n
A<- netEMB$A
g<- netEMB$b
mu<- netEMB$mu
embedding<- netEMB$embedding
res1<- simcos_embedding(A, G, p, mu ,embedding, lambda, gamma, ep, iter_max)
res2<- simsq_embedding(A, G, p, mu ,embedding, lambda, gamma, ep, iter_max)

miscluster(res1$group,g,G)
miscluster(res2$group,g,G)



############################
G<- 2 #分组个数
p<- 2 #embedding的维数 (G<=p)
gamma<- 2
ep<- 1e-7
iter_max<- 100
n<- c(100,400,800,1000)
repmean1<- repmean2<- repmean3<- repmean4<- matrix(0,4,4)
repsd1<- repsd2<-repsd3<- repsd4<- matrix(0,4,4)
#lambda_n<- c(0.0001,0.0005,0.001,0.003)
lambda_n<- 0.003
#n<- c(150,300,600,900,1200)
rep<- 1
rho<-matrix(c(0.4,0.15,0.4,0.15),2,2) #SBM模型模拟用的rho
for (m in 1:4){
 res1<- replicate(rep,SBM_cos_sp_sq(n=n[m],G,p,rho,lambda_n,gamma,ep,iter_max))
 res1<- matrix(as.numeric(res1),4,rep)
 repmean1[,m]<- apply(res1,1,mean)
 repsd1[,m]<- apply(res1,1,sd)
}
repmean1
repsd1

write.csv(repmean1,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repmean1.csv")
write.csv(repsd1,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repsd1.csv")


today<- DCSBM_A_G(n[1],2,rho)
A<- today$A
n<-nrow(A) 
pho<-sum(A)/(n*(n-1))
pAdj<-A+as.numeric(0.25*pho)*matrix(1,nrow=n,ncol=n) #purmutated Adj  
#compute the Laplace Matrix
D <- apply(pAdj,1,sum)^(-0.5)
D[D==Inf] <- 0
H <- diag(D)
L <- H%*%pAdj%*%H 
U<-eigen(L)
alpha<-U$vectors[,1:2]
plot(alpha)

cos_embedding(A, G, p, mu , lambda, gamma, ep, iter_max)$group
sq_embedding(A, G, p, mu , lambda, gamma, ep, iter_max)$group

SBM_cos_sp_sq(n[1],2,2,rho,0.003,gamma,ep,iter_max)

########################################
G<- 2 #分组个数
p<- 2 #embedding的维数 (G<=p)
gamma<- 3
ep<- 1e-7
iter_max<- 100
sigma<- 0.5
lambda_n<- 0.003
lambda_n<- seq(0.001,0.03,0.003)
rep<- 1
n<- c(100,400,800,1000)
#SBM_cos_km_sq(100,G,p,sigma,lambda,gamma,ep,iter_max)
repmean2<- repsd2<- matrix(0,4,4)
for (m in 1:4){
 res2<- replicate(rep,SBM_cos_km_sq(n=n[m],G,p,sigma,lambda_n,gamma,ep,iter_max))
 res2<- matrix(as.numeric(res2),4,rep)
 repmean2[,m]<- apply(res2,1,mean)
 repsd2[,m]<- apply(res2,1,sd)
}
repmean2
repsd2


write.csv(repmean2,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repmean2.csv")
write.csv(repsd2,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repsd2.csv")




###########################################
G<-2 #分组个数
p<- 2 #embedding的维数 (G<=p)
gamma<- 2
ep<- 1e-7
iter_max<- 100
lambda_n<- c(0.0001,0.0005,0.001,0.01)
#set.seed(5678)
n<- c(100,400,800,1000) #n=400的最佳lambda_n可能为0.02
#n<- c(150,300,600,900,1200)
repmean1<- repmean2<- repmean3<- repmean4<- matrix(0,4,4)
repsd1<- repsd2<-repsd3<- repsd4<- matrix(0,4,4)
rho <- (1*1e-2)*(matrix(1,G,G)+diag(c(1,1)))
rep<- 10
for (m in 1:4){
 res3<- replicate(rep,DCBM_cos_sp_sq(n=n[m],G,p,rho,lambda_n,gamma,ep,iter_max))
 res3<- matrix(as.numeric(res3),4,rep)
 repmean3[,m]<- apply(res3,1,mean)
 repsd3[,m]<- apply(res3,1,sd)
}
repmean3
repsd3

write.csv(repmean3,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repmean3.csv")
write.csv(repsd3,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repsd3.csv")



SBM_cos_sp_sq(100,G,p,rho,lambda_n,gamma,ep,iter_max)
DCBM_cos_sp_sq(100,G,p,rho,lambda_n,gamma,ep,iter_max)

############################################

# DCBM_cos_km_sq(n,G,p,sigma,lambda_n,gamma,ep,iter_max)

# SBM_cos_km_sq(n=n,G=G,p=p,sigma=sigma,lambda_n=lambda_n,gamma=gamma,ep=ep,iter_max=iter_max)


# DCBM_cos_sp_sq(n=n,G=G,p=p,rho=rho,lambda_n=lambda_n,gamma=gamma,ep=ep,iter_max=iter_max)
rep<- 10
res3<- replicate(rep,DCBM_cos_sp_sq(n=n,G=G,p=p,rho=rho,lambda_n=lambda_n,gamma=gamma,ep=ep,iter_max=iter_max))
res3<- matrix(as.numeric(res3),4,rep)
repmean3<- apply(res3,1,mean)
repsd3<- apply(res3,1,sd)




G<- 2 #分组个数
p<- 2 #embedding的维数 (G<=p)
gamma<- 3
ep<- 1e-7
iter_max<- 100
#lambda_n<- 0.03
lambda_n<- c(0.0001,0.0005,0.001,0.01)
sigma<- 0.3
rep<- 10
n<- c(100,400,800,1000)
for (m in 1:4){
 res4<- replicate(rep,DCBM_cos_km_sq(n[m],G,p,sigma,lambda_n,gamma,ep,iter_max))
 res4<- matrix(as.numeric(res4),4,rep)
 repmean4[,m]<- apply(res4,1,mean)
 repsd4[,m]<- apply(res4,1,sd)
}
repmean4 #前两行是COS的估计错误率和概率误差, 后两行是SQ的估计错误率和概率误差
repsd4

write.csv(repmean4,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repmean4.csv")
write.csv(repsd4,file = "C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/repsd4.csv")

plot(G2p2COS(1000,0.05)$embedding)

sigma<- 0.3
today<- G2p2COS(100,sigma)
A<- today$A
mu<- today$mu
lambda<- 0.002
cos_embedding(A, G, p, mu , lambda, gamma, ep, iter_max)$group
sq_embedding(A, G, p, mu , lambda, gamma, ep, iter_max)$group


# 画图
plot(n,repmean4[1,],type="b",pch=15,cex=1.5, lty=1, lwd=2, col="red", ylim=c(0, 0.1),
     main="",
     xlab="the number of nodes n", ylab="misclustering rate")

lines(n, repmean4[2,], type="b",pch=15,cex=1.5, lty=1, lwd=2, col="blue")
#legend(5, 0.02, legend = c("K-means", "GMM"),
#       col = c("red", "blue"),lwd=1.5:1.5, lty = 1:2, cex = 0.6)
