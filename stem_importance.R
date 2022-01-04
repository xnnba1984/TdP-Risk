library(haven)
library(caret)
library(MASS)
library(randomForest)

set.seed(2021)
stem <- data.frame(read_sas('allpred2.sas7bdat'))

stem$risk <- factor(stem$risk, levels=c('L','M','H'))
drug <- unique(stem$Drug_Name)

# impute missing values
preProcValues <- preProcess(stem[,6:12], method = c("bagImpute"))
stem[,6:12] <- predict(preProcValues, stem[,6:12])

# p(risk>low)
stem$risk1 <- factor(ifelse(stem$risk=='L', 0, 1))

# p(risk>middle)
stem$risk2 <- factor(ifelse(stem$risk=='H', 1, 0))

set.seed(2021)
perm.count <- 100
system.time({
  # leave one drug out
  for(d in drug){
    print(d)
    index.test <- which(stem$Drug_Name==d)
    train <- stem[-index.test,]; dim(train)
    test <- stem[index.test,]; dim(test)
    
    # train base binary classifiers on un-permuted training data
    # L vs M, H
    model1 <- randomForest(risk1~pred7+arrhythmia+fst_pro+maxpro+
                             foldaym+foldprolong+AMtwooutoffive, data=train)
    # H vs L, M
    model2 <- randomForest(risk2~pred7+arrhythmia+fst_pro+maxpro+
                             foldaym+foldprolong+AMtwooutoffive, data=train)
    # permute each feature on test
    feature <- c('pred7','arrhythmia','fst_pro','maxpro','foldaym','foldprolong','AMtwooutoffive')
    for(f in feature){
      # repeat permutation
      for(i in 1:perm.count){
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
  for(j in 1:perm.count){
    stem[,paste(f, 'pred', j, sep = '_')] <- 
      unlist(apply(stem[,c(paste(f, 'p.L', j, sep = '_'),paste(f, 'p.M', j, sep = '_'),
                           paste(f, 'p.H', j, sep = '_'))],1,function(x) which.max(x)))
    # acc by obs
    stem[,paste(f, 'pred', j, sep = '_')] <- factor(l[stem[,paste(f, 'pred', j, sep = '_')]], 
                                                    levels=c('L','M','H'))
    acc.feature[j, f] <- mean(stem[,paste(f, 'pred', j, sep = '_')] == stem$risk)
  }
}

acc <- 0.619
decrease <- acc - apply(acc.feature, 2, mean)
imp <- sort(decrease/max(decrease), decreasing = T); imp
barplot(height = rev(imp), horiz = T, main='Permutation Predictor Importance')
imp.frame <- data.frame(predictor=names(imp), imp)
rownames(imp.frame) <- NULL
imp.frame$dataset <- 'stem'
saveRDS(imp.frame, 'result/importance_stem_rf_obs.rds')


