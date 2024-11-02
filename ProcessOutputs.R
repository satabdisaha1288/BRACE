library(BhGLM)
library(salso)
beta_clust_labels_train <- lapply(results, function(x) apply(x$beta_clust_labels, 2, as.numeric))

# Function to recreate a beta vector for a single iteration
recreate_beta_iteration <- function(iter, labels, theta) {
  zero_indices <- labels == 0
  beta_vector <- numeric(length(labels))
  non_zero_indices <- !zero_indices
  beta_vector[non_zero_indices] <- theta[labels[non_zero_indices]]
  return(beta_vector)
}

# Recreate new_beta_train 
new_beta_train <- lapply(1:repl, function(repp) {
  sapply((nburn + 1):(iterr), function(iter) {
    recreate_beta_iteration(iter, beta_clust_labels_train[[repp]][, iter], results[[repp]]$theta[[iter]])
  })
})

new_beta_train<-lapply(new_beta_train, function(x) round(x,2))
#Create sum of the beta vector 
new_beta_sum<-sapply(new_beta_train, colSums)

#Plot histogram
hist(new_beta_sum, breaks = 20)

# Mean of the sum of the beta vector for all replicates
mean_sum_new_beta = sapply(lapply(new_beta_train, colSums),mean)

# SD of the sum of the beta vector for all replicates
sd_sum_new_beta = sapply(lapply(new_beta_train, colSums),sd)

colmeans_list<-lapply(new_beta_train, rowMeans)
sd_new_beta = lapply(new_beta_train, function(x) apply(x,1,sd))
quantile_new_beta = lapply(new_beta_train, function(x) apply(x, 1, 
                                                             function(x) quantile(x, probs = c (0.025, 0.5,0.975))))

#Visualize the results
df<-list()
for(repp in 1:repl)
{
  df[[repp]]<-data.frame(beta = simdata$beta, beta_hat = quantile_new_beta[[repp]][2,], 
                         CI_low = quantile_new_beta[[repp]][1,],
                         CI_high =quantile_new_beta[[repp]][3,])
}

#Calculate prediction error and L2 Loss (accuracy)
PE<-list()
accuracy<-list()
a4<-list()
for(repp in 1: repl)
{
  yhat<-as.matrix(train_test_split$XlogRelAbun_split[[repp]]$test_data) %*% colmeans_list[[repp]]
  test_n = length(train_test_split$y_split[[repp]]$test_data)
  num1<-as.matrix(train_test_split$y_split[[repp]]$test_data) - yhat
  (PE[[repp]]<- (1/test_n) * crossprod(num1))
  accuracy[[repp]] = sqrt(sum((simdata$beta - as.matrix(colmeans_list[[repp]]))^2))
}
for(repp in 1:repl){
  y<- as.matrix(train_test_split$y_split[[repp]]$test_data)
  yhat<-as.matrix(train_test_split$XlogRelAbun_split[[repp]]$test_data) %*% colmeans_list[[repp]]
  a4[[repp]] <- measure.glm(y, yhat , family="gaussian")
}
a4.df<-do.call(rbind, a4)
a4.mean = colMeans(a4.df)
a4.sd = apply(a4.df, 2, sd)
#Prediction error
PE
#l2 loss
accuracy
(average_PE = mean(unlist(PE)))
(sd_PE = sd(unlist(PE)))
(average_accuracy = mean(unlist(accuracy)))
(sd_accuracy = sd(unlist(accuracy)))

RP = simdata$beta!=0
#Calculate error rates
FP4 = TP4 = TN4 = FN4 <- vector()
for(repp in 1:repl)
{
  #which(RP == TRUE)
  a = ((round(quantile_new_beta[[repp]][1,],2) < -10^-6) & (round(quantile_new_beta[[repp]][3,],2) > 10^-6))
  b = (round(quantile_new_beta[[repp]][1,],2) == 0) & (round(quantile_new_beta[[repp]][3,],2) == 0)
  c = (round(quantile_new_beta[[repp]][2,],2) == 0)
  PP4 = a | b | c
  PP4 = !PP4
  which(PP4 == TRUE)
  FP4[repp] = sum(PP4=="TRUE"&RP=="FALSE")
  TP4[repp] = sum(PP4=="TRUE"&RP=="TRUE")
  FN4[repp] = sum(PP4=="FALSE"&RP=="TRUE")
  TN4[repp] = sum(PP4=="FALSE"&RP=="FALSE")
}
mean(FP4); sd(FP4)
mean(FN4); sd(FN4)


#Compute rand index
#true clustering
labels <- match(simdata$beta, unique(simdata$beta))
#Find a point estimate for final cluster labels
beta_salso <-lapply(beta_clust_labels_train, function(x) salso::salso(t(x), 
                                                                      maxNClusters = 10))

(arand.mean<-mean(sapply(beta_salso, function(x) pdfCluster::adj.rand.index(labels,x ))))
(arand.sd <- sd(sapply(beta_salso, function(x) pdfCluster::adj.rand.index(labels,x ))))
