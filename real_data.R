#install.packages(c("Matrix", "igraph", "blockmodels", "kernlab", "cluster", "dplyr", "ggplot2"))

rm(list = ls())
library(igraph)
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
 
  #fit<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
  #g1<-fit$cluster#初始聚类分组
  #C_hat<-fit$centers#初始聚类中心
 # 
   #fit<-cos_cluster(alpha, C_hat=C_hat, ep=ep, iter_max=1000)
   #C_hat<-fit$center
   #g_hat<-fit$group
 
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
  
  fit<-cos_cluster(alpha, C_hat=C_hat1, ep=ep, iter_max=iter_max)
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



calculate_ari <- function(true_labels, pred_labels) {
 if (length(true_labels) != length(pred_labels)) {
  stop("真实标签和预测标签长度必须一致！")
 }
 n <- length(true_labels)
 if (n < 2) {
  return(1)
 }
 
 true_f <- as.factor(true_labels)
 pred_f <- as.factor(pred_labels)
 conf_mat <- table(true_f, pred_f)
 n_ij <- conf_mat
 
 a <- sum(choose(n_ij, 2))
 n_i_dot <- rowSums(n_ij)
 n_dot_j <- colSums(n_ij)
 total_pairs <- choose(n, 2)
 
 e_a <- sum(choose(n_i_dot, 2)) * sum(choose(n_dot_j, 2)) / total_pairs
 e_ri <- (e_a + (total_pairs - sum(choose(n_i_dot, 2)) - sum(choose(n_dot_j, 2)) + e_a)) / total_pairs
 
 b <- total_pairs - sum(choose(n_i_dot, 2)) - sum(choose(n_dot_j, 2)) + a
 ri <- (a + b) / total_pairs
 
 max_ri <- (min(sum(choose(n_i_dot, 2)), sum(choose(n_dot_j, 2))) + b) / total_pairs
 min_ri <- (max(0, sum(choose(n_i_dot, 2)) + sum(choose(n_dot_j, 2)) - total_pairs) + b) / total_pairs
 
 if (max_ri - e_ri == 0) {
  return(0)
 }
 ari <- (ri - e_ri) / (max_ri - e_ri)
 
 return(ari)
}




calculate_nmi <- function(true_labels, pred_labels) {
 # 步骤1：输入校验
 if (length(true_labels) != length(pred_labels)) {
  stop("真实标签和预测标签长度必须一致！")
 }
 n <- length(true_labels)
 if (n < 2) {
  return(1)  # 单个样本无信息差异，NMI默认为1
 }
 
 # 步骤2：统一转换为因子（兼容数值型、字符型、因子型标签）
 true_f <- as.factor(true_labels)
 pred_f <- as.factor(pred_labels)
 
 # 步骤3：构建混淆矩阵（真实社区×预测社区的节点计数）
 conf_mat <- table(true_f, pred_f)
 n_ij <- conf_mat  # 每个单元格的节点数（真实社区i，预测社区j）
 n_i_dot <- rowSums(n_ij)  # 真实社区i的总节点数
 n_dot_j <- colSums(n_ij)  # 预测社区j的总节点数
 
 # 步骤4：计算联合概率、真实概率、预测概率
 p_ij <- n_ij / n  # 联合概率 P(C=i, C'=j)
 p_i <- n_i_dot / n  # 真实分组边缘概率 P(C=i)
 p_j <- n_dot_j / n  # 预测分组边缘概率 P(C'=j)
 
 # 步骤5：计算互信息（MI）—— 忽略p_ij=0的项（0*log0无意义）
 mi <- 0
 for (i in 1:nrow(p_ij)) {
  for (j in 1:ncol(p_ij)) {
   if (p_ij[i, j] > 0) {
    mi <- mi + p_ij[i, j] * log2(p_ij[i, j] / (p_i[i] * p_j[j]))
   }
  }
 }
 
 # 步骤6：计算真实分组熵（H(C)）和预测分组熵（H(C')）
 h_true <- -sum(p_i[p_i > 0] * log2(p_i[p_i > 0]))
 h_pred <- -sum(p_j[p_j > 0] * log2(p_j[p_j > 0]))
 
 # 步骤7：计算NMI（几何平均归一化，顶刊常用）
 if (h_true + h_pred == 0) {
  return(1)
 }
 nmi <- (2 * mi) / (h_true + h_pred)
 
 # 确保结果在[0,1]区间（避免浮点运算误差导致微小超出）
 nmi <- max(0, min(1, nmi))
 return(nmi)
}








#####美国三个党派政治类图书数据 
#nodes-105 links-441 groups-3
library(igraph)
setwd("C:/Users/lenovo/Desktop/2026revise/polbooks")
G = read.graph("polbooks.gml", format = "gml")  # read the graph
G <- simplify(G)  # get rid of multiple edges and loops
G <- as.undirected(G,mode="each") # make undirected
A = get.adjacency(G,sparse=FALSE) 
foo <- which(clusters(G)$membership==1) # node labels for the largest connected component
A <- A[foo,foo]
A[A>1]<-1
aff <- get.vertex.attribute(G, "value", index=V(G)) # l-liberal/ c-conservative /n-netural
TrueComm<- ifelse(aff == "n", 1,
                   ifelse(aff == "l", 2,
                          ifelse(aff == "c", 3, NA)))  # 非目标值替换为NA

n<- nrow(A)
pho<-sum(A)/(n*(n-1))
mu<-log(pho/(1-pho))
degree<- rowSums(A)
k <- 3

#===========================================cos惩罚项=================================
# misrate1<- matrix(0,2,100)
# for (num in 1:100){
#  PolBooks_cos_results<- cos_embedding(A, G=3, p=3, mu=mu, lambda=0.0006, gamma=3, ep=1e-9, iter_max=1000)
#  misrate1[1,num]<- miscluster(PolBooks_cos_results$skgroup,TrueComm,3)$mistot
#  misrate1[2,num]<- miscluster(PolBooks_cos_results$skgroup,TrueComm,3)$misrate
# }


library(cluster)
library(skmeans)
PolBooks_cos_results<- cos_embedding(A, G=3, p=3, mu=mu, lambda=0.0001, gamma=2, ep=1e-9, iter_max=100)
CorrectCluster1 <- 1 - miscluster(PolBooks_cos_results$skgroup,TrueComm,3)$misrate
ARI1 <- calculate_ari(TrueComm, PolBooks_cos_results$skgroup)
NMI1 <- calculate_nmi(TrueComm, PolBooks_cos_results$skgroup)
#===========================================谱聚类方法=======================
D_sqrt_inv <- diag(1 / sqrt(degree + 1e-8))  # 度矩阵的逆平方根（加小值避免除0）
laplacian_matrix <- diag(n) - D_sqrt_inv %*% A %*% D_sqrt_inv

# ---------------------- 步骤2：谱聚类建模与训练 ----------------------
# 2.1 对拉普拉斯矩阵进行特征值分解，提取低维嵌入
eigen_result <- eigen(laplacian_matrix)
# 提取前k个特征向量（k=3，对应真实社区数量），按特征值升序排列（谱聚类标准）

eigen_values <- eigen_result$values
eigen_vectors <- eigen_result$vectors
# 按特征值从小到大筛选前k个特征向量（去除第一个全1向量，对应特征值0）
idx <- order(eigen_values)[1:k]  # 跳过第一个特征值（接近0）
low_dim_embedding <- eigen_vectors[, idx]

# 2.2 对低维嵌入进行k-means聚类（谱聚类核心步骤）
#set.seed(01)  # 固定随机种子，保证结果可复现
kmeans_result <- kmeans(low_dim_embedding, centers = k, nstart = 20)
PolBooks_spectral_results <- kmeans_result$cluster  # 谱聚类预测分组（1/2/3）
CorrectCluster2 <- 1 - miscluster(PolBooks_spectral_results,TrueComm,3)$misrate
ARI2 <- calculate_ari(TrueComm, PolBooks_spectral_results)
NMI2 <- calculate_nmi(TrueComm, PolBooks_spectral_results)

# 
# 
# misrate2<- matrix(0,2,100)
# for (num in 1:100)
# {
#  kmeans_result <- kmeans(low_dim_embedding, centers = k, nstart = 20)
#  PolBooks_spectral_results <- kmeans_result$cluster  # 谱聚类预测分组（1/2/3）
#  misrate2[1,num] <-miscluster(PolBooks_spectral_results,TrueComm,3)$mistot
#  misrate2[2,num] <-miscluster(PolBooks_spectral_results,TrueComm,3)$misrate
# }



#===========================================DCBM方法=========================================
degree_correction <- degree / mean(degree)  # 度校正因子（归一化到平均度）
degree_correction[degree_correction == 0] <- 1e-8  # 避免0值导致后续计算报错

# 构建度校正后的邻接矩阵（手动实现公式，无小众函数调用）
A_dc <- A / outer(sqrt(degree_correction), sqrt(degree_correction), FUN = "*")
A_dc[is.na(A_dc)] <- 0  # 处理异常值（避免NA影响迭代）
A_dc[is.infinite(A_dc)] <- 0  # 处理无穷大值（保证矩阵有效性）

# ---------------------- 步骤3：EM算法优化社区划分（手动实现，无estimateSBM） ----------------------
# 3.1 初始化社区分配（随机分配节点到3个社区，基础R函数实现）
membership <- sample(1:k, n, replace = TRUE)

# 3.2 EM算法迭代（手动实现参数估计与社区更新，彻底避开estimateSBM）
max_iter <- 1000  # 最大迭代次数（确保收敛）
tol <- 1e-8      # 收敛阈值（似然值变化小于该值则停止）
likelihood_old <- -Inf  # 初始化旧似然值


#set.seed(123)  # 固定随机种子，保证结果可复现，无小众函数依赖
#misrate3<- matrix(0,2,100)
#for (num in 1:100){

for (iter in 1:max_iter) {
 # -------- E步：手动估计块概率矩阵（不同社区间的连接概率） --------
 block_prob <- matrix(0, nrow = k, ncol = k)  # 初始化k×k块概率矩阵
 for (i in 1:k) {
  for (j in 1:k) {
   # 提取第i社区和第j社区的节点索引
   idx_i <- which(membership == i)
   idx_j <- which(membership == j)
   
   if (length(idx_i) > 0 & length(idx_j) > 0) {
    # 手动计算该块内的平均连接概率（度校正后）
    block_sub <- A_dc[idx_i, idx_j]
    block_prob[i, j] <- mean(block_sub)
   }
   # 无向图保证块概率矩阵对称（手动处理，保证模型合理性）
   block_prob[j, i] <- block_prob[i, j]
  }
 }
 
 # -------- M步：手动更新节点社区分配（最大似然准则） --------
 membership_new <- rep(1, n)  # 初始化新的社区分配
 for (node in 1:n) {
  # 手动计算该节点属于每个社区的似然值
  node_likelihood <- rep(0, k)
  for (c in 1:k) {
   # 提取当前属于社区c的节点索引
   idx_c <- which(membership == c)
   if (length(idx_c) > 0) {
    # 手动结合块概率计算似然值（加1e-8避免log(0)报错）
    node_likelihood[c] <- sum(A_dc[node, idx_c] * log(block_prob[membership[node], c] + 1e-8))
   }
  }
  # 手动选择似然值最大的社区作为新分配
  membership_new[node] <- which.max(node_likelihood)
 }
 
 # -------- 收敛判断：手动计算当前似然值，对比旧似然值 --------
 likelihood_new <- 0
 for (i in 1:n) {
  for (j in 1:n) {
   c_i <- membership_new[i]
   c_j <- membership_new[j]
   likelihood_new <- likelihood_new + A_dc[i, j] * log(block_prob[c_i, c_j] + 1e-8)
  }
 }
 
 # 手动判断是否收敛，收敛则停止迭代
 if (abs(likelihood_new - likelihood_old) < tol) {
  cat("EM算法在第", iter, "步收敛\n")
  break
 }
 
 # 更新似然值和社区分配，进入下一轮迭代
 likelihood_old <- likelihood_new
 membership <- membership_new
}

# 3.3 提取最终DCBM预测分组结果（无任何小众函数调用）
PolBooks_DCBM_results <- membership_new
CorrectCluster3 <- 1 - miscluster(PolBooks_DCBM_results,TrueComm,3)$misrate
ARI3 <- calculate_ari(TrueComm, PolBooks_DCBM_results)
NMI3 <- calculate_nmi(TrueComm, PolBooks_DCBM_results)


#}

#====================================SBM方法=============================
set.seed(1)
# sbm_community <- cluster_spinglass(
#  G,
#  spins = 3,  # 预设分组数2（对应保守派/自由派）
#  update.rule = "config",  # 符合SBM的配置模型假设
# )

# 备选：鲁汶算法（自动选择分组数，效果同样优异）
# misrate4<- matrix(0,2,100)
# for (num in 1:100){
#  sbm_community <- cluster_louvain(G)
#  PolBooks_SBM_results <- membership(sbm_community)  # 提取每个节点的分组标签（从1开始）
#  misrate4[1,num]<- miscluster(PolBooks_SBM_results,TrueComm,3)$mistot
#  misrate4[2,num]<- miscluster(PolBooks_SBM_results,TrueComm,3)$misrate
# } 

sbm_community <- cluster_louvain(G)
PolBooks_SBM_results <- membership(sbm_community)  # 提取每个节点的分组标签（从1开始）
CorrectCluster4 <- 1 - miscluster(PolBooks_SBM_results,TrueComm,3)$misrate
ARI4 <- calculate_ari(TrueComm, PolBooks_SBM_results)
NMI4 <- calculate_nmi(TrueComm, PolBooks_SBM_results)

