library(haven)
library(caret)
library(MASS)
library(randomForest)
library(PRROC)
library(reticulate)

set.seed(2021)
stem <- read.csv('dt17c3max.csv', header = T)
stem <- stem[,-c(2,3,21)]
stem <- stem[,-which(colnames(stem)=="c3v5")]
colnames(stem)[1:2] <- c("Drug_Name", 'risk')
stem$risk <- ifelse(stem$risk=='low', 'L', stem$risk)
stem$risk <- ifelse(stem$risk=='int', 'M', stem$risk)
stem$risk <- ifelse(stem$risk=='high', 'H', stem$risk)
stem <- stem[stem$Drug_Name!='Placebo',]
stem$risk <- factor(stem$risk, levels=c('L','M','H'))

# missing value
sum(!complete.cases(stem))
sum(is.na(stem))
drug <- unique(stem$Drug_Name)

# impute missing values
preProcValues <- preProcess(stem[,3:17], method = c("bagImpute"))
stem[,3:17] <- predict(preProcValues, stem[,3:17])

# p(risk>low)
stem$risk1 <- factor(ifelse(stem$risk=='L', 0, 1))

# p(risk>middle)
stem$risk2 <- factor(ifelse(stem$risk=='H', 1, 0))

set.seed(2021)
perm.count <- 10
system.time({
  # leave one drug out
  for(d in drug){
    print(d)
    index.test <- which(stem$Drug_Name==d)
    train <- stem[-index.test,]; dim(train)
    test <- stem[index.test,]; dim(test)
    
    # train base binary classifiers on un-permuted training data
    # L vs M, H
    model1 <- randomForest(risk1~c3v1+c3v2+c3v3+c3v4+c3v6+c3v7+c3v8+c3v9+c3v10+c3v11+c3v12+c3v13+c3v14+c3v15+c3v16, 
                           data=train)
    # H vs L, M
    model2 <- randomForest(risk2~c3v1+c3v2+c3v3+c3v4+c3v6+c3v7+c3v8+c3v9+c3v10+c3v11+c3v12+c3v13+c3v14+c3v15+c3v16, 
                           data=train)
    # permute each feature on test
    feature <- c('c3v1','c3v2','c3v3','c3v4','c3v6','c3v7','c3v8','c3v9','c3v10','c3v11','c3v12','c3v13','c3v14','c3v15','c3v16')
    for(f in feature){
      # repeat permutation
      #f <- feature[1]; f
      for(i in 1:perm.count){
        #i=1
        test[,which(colnames(test)==f)] <- sample(x = test[,which(colnames(test)==f)], 
                                                  size = length(test[,which(colnames(test)==f)]), replace=F)
        # p.L
        pred1 <- predict(model1, test, type = 'prob')[,2]
        p.L <- 1 - pred1 
        stem[index.test, paste(f, 'p.L', i, sep = '_')] <- p.L
        # p.M, p.H
        pred2 <- predict(model2, test, type = 'prob')[,2]
        p.H <- pred2
        p.M <- pred1-pred2
        p.M <- ifelse(p.M<0, 0, p.M)
        stem[index.test, paste(f, 'p.M', i, sep = '_')] <- p.M
        stem[index.test, paste(f, 'p.H', i, sep = '_')] <- p.H
      }
    }
  }
})

# calculate prediction per permuted predictor (obs)
l <- c('L','M','H')
acc.feature <- data.frame(matrix(ncol = length(feature), nrow = perm.count))
colnames(acc.feature) <- feature
for(f in feature){
  #f <- feature[1]; f
  for(j in 1:perm.count){
    #j <- 1; j
    stem[,paste(f, 'pred', j, sep = '_')] <- 
      unlist(apply(stem[,c(paste(f, 'p.L', j, sep = '_'),paste(f, 'p.M', j, sep = '_'),
                           paste(f, 'p.H', j, sep = '_'))],1,function(x) which.max(x)))
    # acc by obs
    stem[,paste(f, 'pred', j, sep = '_')] <- factor(l[stem[,paste(f, 'pred', j, sep = '_')]], 
                                                    levels=c('L','M','H'))
    acc.feature[j, f] <- mean(stem[,paste(f, 'pred', j, sep = '_')] == stem$risk)
  }
}

acc <- 0.831
decrease <- acc - apply(acc.feature, 2, mean)
imp <- sort(decrease/max(decrease), decreasing = T)
imp <- imp[imp>0]; imp
barplot(height = rev(imp), horiz = T, main='Permutation Predictor Importance')
imp.frame <- data.frame(predictor=names(imp), imp)
rownames(imp.frame) <- NULL
imp.frame$dataset <- 'wedge'
saveRDS(imp.frame, 'result/importance_wedge_rf_obs.rds')



