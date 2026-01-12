rm(list = ls())
library("igraph")
library("clue")


#####Email数据（embedding文章中使用）
#使用来自大型欧洲研究机构的电子邮件数据生成，机构成员所有传入和传出电子邮件，所有
#人都属于研究所的42个部门，该网络表示email-EuAll网络的核心(因为还包含机构成员和外部人员的链接)
setwd("C:/Users/lenovo/Desktop/cosembedding/Numerical Experiments/Email")
network<- read.table("email-Eu-core.txt",header=F)+1
label<- read.table("email-Eu-core-department-labels.txt",header=F)+1
comm<- table(label[,2])
G<-graph_from_data_frame(network,directed=F)
A = matrix(0,nrow(label),nrow(label))
for (k in 1:nrow(network)){
 egdei=network[k,1]
 egdej=network[k,2]
 A[egdei,egdej]<-1
 A[egdej,egdei]<-1
}
diag(A)<- 0
degree<- t(rbind(label[,2],rowSums(A))) #第一列对应点的组别，第二列对应点的degree
tlabel<- table(degree[,1])


# example 1: 9 10 12 22

# example 2: 1 7 8 17 18

# example 3: 6 14 16 20 21 23

#repermu<- cbind(t(which(label[,2]==9)),t(which(label[,2]==10)),t(which(label[,2]==12)),t(which(label[,2]==22)))

#repermu<- cbind(t(which(label[,2]==1)),t(which(label[,2]==7)),t(which(label[,2]==8)),t(which(label[,2]==17)),t(which(label[,2]==18)))
 
repermu<- cbind(t(which(label[,2]==6)),t(which(label[,2]==14)),t(which(label[,2]==16)),t(which(label[,2]==20)),t(which(label[,2]==21)),t(which(label[,2]==23)))

exdata<- label[repermu,]

#,t(which(label[,2]==21)),t(which(label[,2]==23))
#exdata<- label[which(label[,2]==9|label[,2]==10|label[,2]==12|label[,2]==22),]

g<- exdata[,2]
a<- A[exdata[,1],exdata[,1]]
zerodegree<- which(rowSums(a)==0)
g<- g[-zerodegree]
a<- a[-zerodegree,-zerodegree]

n<- nrow(a)
pho<-sum(a)/(n*(n-1))
mu<-log(pho/(1-pho))

for (i in 1:length(table(g))){
 g[g==as.numeric(names(table(g))[i])]<- i
}
######可视化选择的社区组成的网络
# x<-graph_from_adjacency_matrix(a,mode=c("undirected"))
# plot(x)
#rowSums(a)

# 4 group的lambda=0.001,gamma=2
# 5 group的lambda=0.001,gamma=2
# 6 group的lambda=0.001,gamma=2


Email_cos_results<- cos_embedding(a, G=6, p=6, mu=mu, lambda=0.0005, gamma=2, ep=1e-7, iter_max=100)
Email_cos_results$group
Email_cos_results$kgroup
Email_cos_results$skgroup


miscluster(Email_cos_results$group,g,6)
miscluster(Email_cos_results$kgroup,g,6)
miscluster(Email_cos_results$skgroup,g,6)


 write.csv(g,file = "truegroup.csv")
 write.csv(Email_cos_results$group,file = "cosgroup.csv")
 
 write.csv(Email_cos_results$center,file = "coscenter4.csv")
 write.csv(Email_cos_results$center,file = "coscenter5.csv")
 write.csv(Email_cos_results$center,file = "coscenter6.csv")


set.seed(12)
Email_sq_results<- sq_embedding(a, G=4, p=4, mu=mu, lambda=0.0005, gamma=2, ep=1e-7, iter_max=100)
Email_sq_results$group
Email_sq_results$kgroup
Email_sq_results$skgroup

miscluster(Email_sq_results$group,g,4)
miscluster(Email_sq_results$kgroup,g,4)


write.csv(Email_sq_results$group,file = "sqgroup.csv")

