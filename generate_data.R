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


library(igraph)
library(cluster)
library(skmeans)
# ---------------------- 配置DCBM核心参数 ----------------------
#set.seed(123)  # 随机种子，保证结果可复现
CorrectCluster1 <- matrix(0,1,10)
CorrectCluster2 <- matrix(0,1,10)
CorrectCluster3 <- matrix(0,1,10)
CorrectCluster4 <- matrix(0,1,10)
ARI1 <- matrix(0,1,10)
ARI2 <- matrix(0,1,10)
ARI3 <- matrix(0,1,10)
ARI4 <- matrix(0,1,10)
NMI1 <- matrix(0,1,10)
NMI2 <- matrix(0,1,10)
NMI3 <- matrix(0,1,10)
NMI4 <- matrix(0,1,10)

for (tms in 1:10)
{
K <- 3  # 社区数量
community_sizes <- c(100, 120, 80)  # 各社区节点数量
N <- sum(community_sizes)  # 总节点数（300）

# 1. 社区连接概率矩阵B
B <- matrix(
 data = c(
  0.2, 0.03, 0.03,
  0.03, 0.2, 0.03,
  0.03, 0.03, 0.18
 ),
 nrow = K,
 ncol = K,
 byrow = TRUE
)

# 2. 节点度偏好参数θ
theta <- runif(n = N, min = 0.1, max = 2)

# 3. 节点社区标签
community_labels <- rep(1:K, times = community_sizes)

# ---------------------- 手动构建DCBM邻接矩阵 ----------------------
dcbm_adjmatrix <- matrix(0, nrow = N, ncol = N)

# 遍历上三角节点对，避免重复计算
for (i in 1:(N-1)) {
 for (j in (i+1):N) {
  r <- community_labels[i]
  s <- community_labels[j]
  
  # DCBM连接概率计算
  connect_prob <- theta[i] * theta[j] * B[r, s]
  connect_prob <- min(connect_prob, 1)  # 截断超1的概率
  
  # 抽样生成边
  edge <- rbinom(n = 1, size = 1, prob = connect_prob)
  if (edge == 1) {
   dcbm_adjmatrix[i, j] <- 1
   dcbm_adjmatrix[j, i] <- 1
  }
 }
}

dcbm_graph <- graph_from_adjacency_matrix(
 adjmatrix = dcbm_adjmatrix,
 mode = "undirected",  # 无向网络
 diag = FALSE,         # 无自环
 weighted = NULL      # 非加权网络
)



n<- nrow(dcbm_adjmatrix)
pho<-sum(dcbm_adjmatrix)/(n*(n-1))
mu<-log(pho/(1-pho))
degree<- rowSums(dcbm_adjmatrix)
k<- 3

#===========================================cos惩罚项=================================
# misrate1<- matrix(0,2,100)
# for (num in 1:100){
#  PolBooks_cos_results<- cos_embedding(A, G=3, p=3, mu=mu, lambda=0.0006, gamma=3, ep=1e-9, iter_max=1000)
#  misrate1[1,num]<- miscluster(PolBooks_cos_results$skgroup,TrueComm,3)$mistot
#  misrate1[2,num]<- miscluster(PolBooks_cos_results$skgroup,TrueComm,3)$misrate
# }
PolBooks_cos_results<- cos_embedding(dcbm_adjmatrix, G=3, p=3, mu=mu, lambda=0.0001, gamma=2, ep=1e-9, iter_max=100)
CorrectCluster1[tms] <- 1 - miscluster(PolBooks_cos_results$kgroup,community_labels,3)$misrate
ARI1[tms] <- calculate_ari(community_labels, PolBooks_cos_results$kgroup)
NMI1[tms] <- calculate_nmi(community_labels, PolBooks_cos_results$kgroup)
#===========================================谱聚类方法=======================
D_sqrt_inv <- diag(1 / sqrt(degree + 1e-8))  # 度矩阵的逆平方根（加小值避免除0）
laplacian_matrix <- diag(n) - D_sqrt_inv %*% dcbm_adjmatrix %*% D_sqrt_inv

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
CorrectCluster2[tms] <- 1 - miscluster(PolBooks_spectral_results,community_labels,3)$misrate
ARI2[tms] <- calculate_ari(community_labels, PolBooks_spectral_results)
NMI2[tms] <- calculate_nmi(community_labels, PolBooks_spectral_results)


# misrate2<- matrix(0,2,100)
# for (num in 1:100)
# {
#  kmeans_result <- kmeans(low_dim_embedding, centers = k, nstart = 20)
#  PolBooks_spectral_results <- kmeans_result$cluster  # 谱聚类预测分组（1/2/3）
#  misrate2[1,num] <-miscluster(PolBooks_spectral_results,TrueComm,3)$mistot
#  misrate2[2,num] <-miscluster(PolBooks_spectral_results,TrueComm,3)$misrate
# }

#===========================================DCBM方法=========================================
node_degree <- rowSums(dcbm_adjmatrix)
fit_dcbm_manual <- function(adj_mat, K, max_iter = 100, tol = 1e-4) {
 node_degree <- rowSums(dcbm_adjmatrix)
 # 步骤1：初始化
 # 随机分配社区标签
 community <- sample(1:K, N, replace = TRUE)
 # 初始化θ（节点度标准化）
 theta <- node_degree / mean(node_degree)
 # 初始化B（社区连接概率矩阵）
 B <- matrix(0, nrow = K, ncol = K)
 for (r in 1:K) {
  for (s in 1:K) {
   idx_r <- which(community == r)
   idx_s <- which(community == s)
   if (length(idx_r) > 0 & length(idx_s) > 0) {
    sub_adj <- adj_mat[idx_r, idx_s]
    B[r, s] <- sum(sub_adj) / (length(idx_r) * length(idx_s) * mean(theta[idx_r]) * mean(theta[idx_s]))
   }
   B[r, s] <- min(max(B[r, s], 1e-6), 1 - 1e-6)  # 截断避免0/1
  }
 }
 
 ll_old <- -Inf  # 初始化似然值
 for (iter in 1:max_iter) {
  # 步骤2：E步（计算后验概率，简化为硬分配：直接更新社区标签）
  new_community <- rep(1, N)
  for (i in 1:N) {
   # 计算节点i属于各社区的似然值
   ll_community <- rep(0, K)
   for (r in 1:K) {
    for (j in 1:N) {
     if (i != j) {
      s <- community[j]
      p_ij <- theta[i] * theta[j] * B[r, s]
      p_ij <- min(max(p_ij, 1e-6), 1 - 1e-6)
      ll_community[r] <- ll_community[r] + adj_mat[i, j] * log(p_ij) + (1 - adj_mat[i, j]) * log(1 - p_ij)
     }
    }
   }
   new_community[i] <- which.max(ll_community)  # 硬分配：选择似然值最大的社区
  }
  
  # 步骤3：M步（更新θ和B）
  # 更新θ
  for (i in 1:N) {
   r <- new_community[i]
   denom <- 0
   for (j in 1:N) {
    if (i != j) {
     s <- new_community[j]
     denom <- denom + theta[j] * B[r, new_community[j]]
    }
   }
   theta[i] <- ifelse(denom > 0, node_degree[i] / denom, 1)
  }
  
  # 更新B
  new_B <- matrix(0, nrow = K, ncol = K)
  for (r in 1:K) {
   for (s in 1:K) {
    idx_r <- which(new_community == r)
    idx_s <- which(new_community == s)
    if (length(idx_r) > 0 & length(idx_s) > 0) {
     numer <- sum(adj_mat[idx_r, idx_s])
     denom <- sum(outer(theta[idx_r], theta[idx_s]))
     new_B[r, s] <- ifelse(denom > 0, numer / denom, 1e-6)
    }
    new_B[r, s] <- min(max(new_B[r, s], 1e-6), 1 - 1e-6)
   }
  }
  
  # 步骤4：计算似然值并判断收敛
  ll_new <- 0
  for (i in 1:(N-1)) {
   for (j in (i+1):N) {
    r <- new_community[i]
    s <- new_community[j]
    p_ij <- theta[i] * theta[j] * new_B[r, s]
    p_ij <- min(max(p_ij, 1e-6), 1 - 1e-6)
    ll_new <- ll_new + adj_mat[i, j] * log(p_ij) + (1 - adj_mat[i, j]) * log(1 - p_ij)
   }
  }
  
  # 收敛判断：似然值变化小于阈值或社区标签无变化
  if (abs(ll_new - ll_old) < tol | all(new_community == community)) {
   cat(paste("迭代收敛，迭代次数：", iter, "\n", sep = ""))
   break
  }
  
  # 更新参数
  community <- new_community
  B <- new_B
  ll_old <- ll_new
  
  # 最大迭代次数判断
  if (iter == max_iter) {
   cat("达到最大迭代次数，未完全收敛\n")
  }
 }
 
 return(list(
  pred_community = community,
  theta = theta,
  B = B
 ))
}


dcbm_manual_result <- fit_dcbm_manual(
 adj_mat = dcbm_adjmatrix,
 K = K,
 max_iter = 100,
 tol = 1e-4
)


# 3.3 提取最终DCBM预测分组结果（无任何小众函数调用）
PolBooks_DCBM_results <-  dcbm_manual_result$pred_community
CorrectCluster3[tms] <- 1 - miscluster(PolBooks_DCBM_results,community_labels,3)$misrate
ARI3[tms] <- calculate_ari(community_labels, PolBooks_DCBM_results)
NMI3[tms] <- calculate_nmi(community_labels, PolBooks_DCBM_results)


#}

#====================================SBM方法=============================
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

sbm_community <- cluster_louvain(dcbm_graph)
PolBooks_SBM_results <- membership(sbm_community)  # 提取每个节点的分组标签（从1开始）
CorrectCluster4[tms] <- 1 - miscluster(PolBooks_SBM_results,community_labels,3)$misrate
ARI4[tms] <- calculate_ari(community_labels, PolBooks_SBM_results)
NMI4[tms] <- calculate_nmi(community_labels, PolBooks_SBM_results)

}

# 打开新绘图窗口
dev.new()

# 第一步：绘制第一条曲线，奠定画布（指定横纵坐标范围）
plot(
 x = 1:100,
 y = CorrectCluster1,
 type = "l",  # 绘图类型：曲线（必填）
 lwd = 2,     # 线宽
 col = "#2E86AB",  # 曲线颜色
 xlab = "横坐标（1-100）",
 ylab = "纵坐标（0-1）",
 main = "多条曲线图（基础plot()实现）",
 xlim = c(1, 100),  # 横坐标固定为1-100
 ylim = c(0.5, 1),    # 纵坐标固定为0-1
 xaxt = "n"  # 后续自定义x轴刻度（可选，优化美观度）
)

# 第二步：添加自定义x轴刻度（1-100，每10个单位一个刻度）
axis(side = 1, at = seq(0, 100, 10), labels = seq(0, 100, 10))

# 第三步：用lines()添加后续曲线（不可用plot()，否则覆盖画布）
# 添加曲线2：下降趋势
lines(
 x = 1:100,
 y = CorrectCluster2,
 lwd = 2,
 col = "#A23B72",
 lty = 2  # 线型：虚线
)

# 添加曲线3：波动趋势
lines(
 x = 1:100,
 y = CorrectCluster3,
 lwd = 2,
 col = "#F18F01",
 lty = 3  # 线型：点线
)

# 添加曲线4：阶梯趋势
lines(
 x = 1:100,
 y = CorrectCluster4,
 lwd = 2,
 col = "#C73E1D",
 lty = 4  # 线型：点虚线
)

# 第四步：添加网格线（提升可读性）
grid(nx = 10, ny = 5, col = "gray80", lty = 2)

# 第五步：添加图例（区分多条曲线）
legend(
 x = "bottomright",  # 图例位置：右下角
 legend = c("上升趋势", "下降趋势", "波动趋势", "阶梯趋势"),
 col = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"),
 lwd = 2,
 lty = c(1, 2, 3, 4),
 bty = "n"  # 隐藏图例边框
)

