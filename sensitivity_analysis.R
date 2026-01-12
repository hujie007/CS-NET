library(skmeans)
library(MASS)
library(igraph)

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
cos_cluster<-function(alpha, C_hat, ep = 1e-6, iter_max = 100) {
 # ========== 步骤1：参数前置校验（避免初始NA） ==========
 G <- nrow(C_hat)
 if (G == 0) stop("C_hat的行数（G）不能为0！")
 if (ncol(alpha) != ncol(C_hat)) stop("alpha的列数必须等于C_hat的列数！")
 if (is.na(iter_max) || iter_max <= 0) iter_max <- 100  # 兜底最大迭代次数
 if (is.na(ep) || ep <= 0) ep <- 1e-6                    # 兜底收敛阈值
 
 # ========== 步骤2：初始化迭代变量（确保无NA） ==========
 iter <- 0                   # 迭代次数初始化为0
 delta <- 1                  # 初始误差（大于ep，确保进入循环）
 C_hat <- as.matrix(C_hat)   # 强制转为矩阵，避免数据类型异常
 alpha <- as.matrix(alpha)   # 强制转为矩阵
 
 # ========== 步骤3：带鲁棒校验的while循环 ==========
 while (iter < iter_max && delta > ep) {
  # 子步骤3.1：循环条件变量NA校验（提前拦截）
  if (any(is.na(c(iter, iter_max, delta, ep)))) {
   warning(paste("迭代", iter, "次：关键变量含NA，终止迭代"))
   break
  }
  
  # 子步骤3.2：保存上一轮的C_hat
  C_hat1 <- C_hat
  
  # 子步骤3.3：计算分组标签（hg）
  t_hat <- alpha %*% t(C_hat1)
  hg <- apply(t_hat, 1, which.max)  # 每行最大值对应的列（分组）
  
  # 子步骤3.4：更新C_hat（核心修复：空分组+分母防0）
  for (g in 1:G) {
   # 提取当前组的alpha行（空分组兜底）
   alpha_g <- alpha[hg == g, , drop = FALSE]  # drop=FALSE避免降维为向量
   
   if (nrow(alpha_g) == 0) {
    # 空分组：保留原C_hat[g,]，避免NA
    C_hat[g, ] <- C_hat1[g, ]
   } else {
    # 非空分组：计算y并防分母为0
    y <- alpha_g
    # 计算diag(y %*% t(y))，并将0替换为极小值（1e-10）防除以0
    diag_yy <- diag(y %*% t(y))
    diag_yy[diag_yy < 1e-10] <- 1e-10  # 兜底分母，避免0
    sqrt_diag <- sqrt(diag_yy)
    
    # 计算h向量（防NA/Inf）
    h <- colSums(y / sqrt_diag)
    # 归一化：若h全为0，保留原C_hat[g,]
    if (all(h == 0)) {
     C_hat[g, ] <- C_hat1[g, ]
    } else {
     C_hat[g, ] <- h / sqrt(sum(h^2))  # 标准化
    }
   }
  }
  
  # 子步骤3.5：计算delta（防NA/Inf）
  delta <- sqrt(sum((C_hat - C_hat1)^2, na.rm = TRUE) / G)
  # 兜底delta：若delta为NA/Inf，设为0（触发收敛）
  if (is.na(delta) || is.infinite(delta)) delta <- 0
  
  # 子步骤3.6：迭代次数+1
  iter <- iter + 1
 }
 
 # ========== 步骤4：最终分组计算 ==========
 t_hat_final <- alpha %*% t(C_hat1)
 hg_final <- apply(t_hat_final, 1, which.max)
 
 # ========== 步骤5：返回结果 ==========
 return(list(
  group = hg_final,
  center = C_hat1,
  iters = iter,
  final_delta = delta  # 新增：返回最终误差，便于调试
 ))
}



#根据分组中心及分组情况利用坐标梯度下降法估计出 每个节点的潜在变量 alpha、分组结果g_hat、组中心C_hat
cos_embedding_sens_anal<-function(A, G, p, mu, mu_nosy, lambda, gamma, gamma_nosy, ep, iter_max, ranosyx, ranosyy){
 
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
 alpha<-U$vectors[,1:p] #初始基准embedding向量
 alpha_nosy<-U$vectors[,1:p] + mvrnorm(n, 0, ranosyx) %*% t(mvrnorm(p, 0, ranosyx)) #加扰动项
 
 #alpha<- alpha/sqrt(apply(alpha^2,1,sum)) # 对特征向量进行 row-normalize
 
 # 分别对alpha做kmeans和skmeans
 fitk<-kmeans(alpha,G) #利用K-Means计算初始聚类分组及聚类中心
 C_hatk<- fitk$centers #kmeans初始分组中心坐标
 g_k<-fitk$cluster #kmeans初始聚类分组
 
 fitsk<-skmeans(alpha,G) #利用skmeans计算初始聚类分组及聚类中心
 C_hatsk<-fitsk$prototypes  #skmeans初始分组中心坐标
 g_sk<-fitsk$cluster #skmeans初始聚类分组
 
 C_hat<- C_hatk #初始基准中心向量
 C_hat_nosy<- C_hatk + mvrnorm(G, 0, ranosyy) %*% t(mvrnorm(p, 0, ranosyy)) #加扰动项
 g_hat<- g_k
 g_nosy<- g_k
 
 #initial value
 embedding_initial <- alpha
 embedding_nosy_initial <- alpha_nosy
 center_initial <- C_hat
 center_nosy_initial <- C_hat_nosy
 
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
 
 ###############用加了扰动项的中心向量、embedding向量计算################
 iter_nosy=0
 delta_nosy=1
 while(iter_nosy<iter_max && delta_nosy>ep){
  
  alpha1_nosy<-alpha_nosy
  C_hat1_nosy<-C_hat_nosy
  mu1_nosy<-mu_nosy
  g_hat1_nosy<-g_hat
  like1_nosy<-coslike_hood(alpha1_nosy, C_hat1_nosy, g_hat1_nosy, A, mu1_nosy, lambda=lambda)
  
  fit_nosy<-cosl_partial(alpha1_nosy, C_hat1_nosy, g_hat1_nosy, A, mu1_nosy, lambda=lambda) 
  mu_nosy<-mu1_nosy-gamma_nosy*fit_nosy$p_mu #更新了mu
  alpha_nosy<-alpha1_nosy-gamma_nosy*fit_nosy$p_alpha  #更新了alpha
  #alpha<-up_alpha(alpha1, C_hat1, g_hat1, A, mu, gamma, lambda)
  
  fit_nosy<- cos_cluster(alpha_nosy, C_hat1_nosy, ep, iter_max)
  C_hat_nosy<-fit_nosy$center
  g_hat_nosy<-fit_nosy$group
  
  
  like0_nosy<-coslike_hood(alpha_nosy, C_hat_nosy, g_hat_nosy, A, mu_nosy, lambda=lambda) #更新过后的alpha等参数计算新的似然函数
  
  #如果新的似然大于旧的，则减小步长(gamma)再计算
  if(like0_nosy>like1_nosy){
   gamma_nosy<-gamma_nosy/2
   alpha_nosy<-alpha1_nosy
   C_hat_nosy<-C_hat1_nosy
   mu_nosy<-mu1_nosy
   g_hat_nosy<-g_hat1_nosy
  } 
  #如果新的似然小于旧的,寻找负对数似然的最小值(对数似然的最大值)
  if(like0<=like1){
   delta_nosy<-(like1_nosy-like0_nosy)/abs(like1_nosy)
   iter_nosy<-iter_nosy+1
  }
  
 }
 
 g_hat_nosy <- unname(g_hat_nosy)
 embedding_nosy <- alpha_nosy
 mu_hat_nosy <- mu_nosy 
  
 g_hat<-unname(g_hat)
 embedding<- alpha
 mu_hat<- mu
 return(list(kgroup=g_k, skgroup=g_sk, group=g_hat, group_nosy=g_hat_nosy,
              embedding= embedding, embedding_nosy= embedding_nosy,
             embedding_initial=embedding_initial,embedding_nosy_initial=embedding_nosy_initial,
             center_initial=center_initial,center_nosy_initial=center_nosy_initial,
             mu_hat=mu_hat, mu_hat_nosy = mu_hat_nosy,
             center=C_hat, center_nosy = C_hat_nosy,
             likely=like0, likely_nosy=like0_nosy,
             Gamma=gamma, Gamma_nosy=gamma_nosy,
             iters=iter))
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


tms = 100
baseline = numeric(tms)
nosy = numeric(tms)
embedding_nosy_distance = numeric(tms)
center_nosy_distance = numeric(tms)

ranosy <- seq(from = 0, to = 0.5, by = 0.05)
xybaseline <- matrix(0,length(ranosy),length(ranosy))
xynosy <- matrix(0,length(ranosy),length(ranosy))
xyabs <- matrix(0,length(ranosy),length(ranosy))


for (x in 1:length(ranosy)){
for (y in 1:length(ranosy)){

for (t in 1:tms)
 {
 K <- 2  # 块数
 n_per_block <- c(50, 50)  # 每个块的节点数
 total_nodes <- sum(n_per_block)
 # 块间概率矩阵（确保无NA，且值在0-1之间）
 P <- matrix(c(0.3, 0.05,
               0.05, 0.2), nrow = K, byrow = TRUE)
 # 生成度校正参数（强制去除NA，归一化）
 theta <- rlnorm(total_nodes, meanlog = 0, sdlog = 0.5)
 theta[is.na(theta)] <- 1  # 替换NA为1（默认度校正系数）
 theta <- theta / mean(theta)  # 归一化
 
 # 生成块标签（确保标签在1-K范围内）
 block_labels <- rep(1:K, times = n_per_block)
 if (any(!block_labels %in% 1:K)) {
  stop("块标签超出1-K范围，请检查n_per_block和K的匹配性！")
 }

 # ===================== 2. 生成DCBM邻接矩阵（防NA） =====================
 adj_matrix <- matrix(0, nrow = total_nodes, ncol = total_nodes)
 
 for (i in 1:(total_nodes - 1)) {
  for (j in (i + 1):total_nodes) {
   b_i <- block_labels[i]
   b_j <- block_labels[j]
   
   # 步骤1：校验索引和概率值，避免NA
   if (is.na(b_i) || is.na(b_j) || b_i > K || b_j > K) {
    p_ij <- 0  # 无效块标签，设为无边
   } else {
    p_ij <- theta[i] * theta[j] * P[b_i, b_j]
   }

   # 步骤2：限制概率在0-1之间（避免p_ij>1或p_ij<0）
   p_ij <- max(0, min(1, p_ij))
   
   # 步骤3：伯努利抽样（此时p_ij必为有效值）
   adj_matrix[i, j] <- rbinom(1, size = 1, prob = p_ij)
   adj_matrix[j, i] <- adj_matrix[i, j]  # 无向网络对称赋值
  }
 }
 
 n<- nrow(adj_matrix)
 pho<-sum(adj_matrix)/(n*(n-1))
 
 mu<-log(pho/(1-pho))
 mu_nosy<-log(pho/(1-pho))
 
 lambda=0.01 
 gamma=2 
 gamma_nosy=2
 ep=1e-9
 iter_max=100
 
 sensitivity_analysis <- 
 cos_embedding_sens_anal(A=adj_matrix, G=K, p=K, mu, mu_nosy, lambda, gamma, gamma_nosy, ep, iter_max,ranosyx=ranosy[x],ranosyy=ranosy[y])
 
 baseline[t] <- miscluster(sensitivity_analysis$group,block_labels,2)$misrate
 nosy[t] <- miscluster(sensitivity_analysis$group_nosy,block_labels,2)$misrate
 
 #embedding_nosy_distance[t] <- sqrt(sum((sensitivity_analysis$embedding_initial-sensitivity_analysis$embedding_nosy_initial)^2))
 #center_nosy_distance[t] <- sqrt(sum(sensitivity_analysis$center_initial-sensitivity_analysis$center_nosy_initial)^2)
 }
 
 xybaseline[x,y] <- mean(baseline)
 xynosy[x,y] <- mean(nosy)
 xyabs[x,y] <- sum(abs(nosy-baseline))
 }
 }





# 生成A和B的所有组合（笛卡尔积）
dim3_new <- expand.grid(X=ranosy, Y=ranosy)
# 将C矩阵展开为一维向量，匹配组合顺序
dim3_new$Z <- as.vector(xyabs/100)  # as.vector默认按列展开，若需按行展开用as.vector(t(C))

write.csv(dim3_new,file = "C:/Users/lenovo/Desktop/2026revise/dim3_new.csv")
