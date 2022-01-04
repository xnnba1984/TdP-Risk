library(haven)
library(caret)
library(MASS)
library(randomForest)
library(PRROC)
library(reticulate)

set.seed(2021)
stem <- data.frame(read_sas('allpred2.sas7bdat'))

# impute missing values
preProcValues <- preProcess(stem[,6:12], method = c("bagImpute"))
stem[,6:12] <- predict(preProcValues, stem[,6:12])

# p(risk>low)
stem$risk1 <- factor(ifelse(stem$risk=='L', 0, 1))

# p(risk>middle)
stem$risk2 <- factor(ifelse(stem$risk=='H', 1, 0))
drug <- unique(stem$Drug_Name)
stem$risk <- factor(stem$risk, levels=c('L','M','H')); str(stem)

set.seed(2021)
boost.count <- 1000
system.time({
  # leave one drug out
  for(d in drug){
    #d <- drug[1]
    print(d)
    index.test <- which(stem$Drug_Name==d)
    train <- stem[-index.test,]; dim(train)
    test <- stem[index.test,]; dim(test)
    
    # stratified bootstrap
    index.table <- data.frame()
    drug.train <- unique(train$Drug_Name); length(drug.train)
    for(d.train in drug.train){
      #d.train <- drug.train[20]; d.train
      drug.index <- which(train$Drug_Name==d.train); length(drug.index)
      temp.table <- replicate(boost.count, 
                              sample(x = drug.index, size = length(drug.index), replace = T)); dim(temp.table)
      index.table <- rbind(index.table, temp.table); dim(index.table)
    }; dim(index.table)
    
    # modeling on stratified bootstrap sample
    for(i in 1:ncol(index.table)){
      #print(paste('boostrap',i, sep = '_'))
      # draw a stratified random sample
      #i=1
      train.bootstrap <- train[index.table[,i],]; dim(train.bootstrap)
      table(train.bootstrap$Drug_Name)
      
      # L vs M, H
      model1 <- randomForest(risk1~pred7+arrhythmia+fst_pro+maxpro+
                               foldaym+foldprolong+AMtwooutoffive, data=train.bootstrap)
      pred1 <- predict(model1, test, type = 'prob')[,2]
      p.L <- 1 - pred1 
      stem[index.test, paste('p.L', i, sep = '_')] <- p.L
      
      # H vs L, M
      model2 <- randomForest(risk2~pred7+arrhythmia+fst_pro+maxpro+foldaym+foldprolong+AMtwooutoffive,
                             data=train.bootstrap)
      pred2 <- predict(model2, test, type = 'prob')[,2]
      p.H <- pred2
      p.M <- pred1-pred2
      p.M <- ifelse(p.M<0, 0, p.M)
      stem[index.test, paste('p.M', i, sep = '_')] <- p.M
      stem[index.test, paste('p.H', i, sep = '_')] <- p.H
    }
  }
})

# calculate prediction per bootstrap
l <- c('L','M','H')
for(j in 1:boost.count){
  stem[,paste('pred',j, sep = '_')] <- 
    unlist(apply(stem[,c(paste('p.L', j, sep = '_'),paste('p.M', j, sep = '_'),
                         paste('p.H', j, sep = '_'))],1,function(x) which.max(x)))
  # acc by obs
  stem[,paste('pred',j, sep = '_')] <- factor(l[stem[,paste('pred',j, sep = '_')]], 
                                              levels=c('L','M','H'))
}

# save bootstrap result
saveRDS(stem, file = "result/randomforest_stem_bootstrap.RDS") 

######################################################################
# acc by obs
######################################################################
acc.obs <- c()
# confusion matrix
for(j in 1:boost.count){
  cm <- confusionMatrix(stem[,paste('pred',j, sep = '_')], stem$risk, mode = 'everything')
  acc.obs[j] <- cm[["overall"]][["Accuracy"]]
}

# CI boundary
ci.low <- boost.count * 0.025 + 1; ci.low
ci.high <- boost.count * 0.975; ci.high

# acc
cat('low 2.5%',sort(acc.obs)[ci.low], '\n')
cat('high 97.5%',sort(acc.obs)[ci.high], '\n')
cat('mean',mean(acc.obs), '\n')

######################################################################
# acc by drug: average
######################################################################
# drug risk
drug.pred.list <- list()
# drug prediction for each bootstrap
for(j in 1:boost.count){
  drug.pred <- unique(stem[,c('Drug_Name','risk')])
  temp <- stem[, c('Drug_Name', paste('p.L', j, sep = '_'), paste('p.M', j, sep = '_'), 
                   paste('p.H', j, sep = '_'))]; dim(temp)
  colnames(temp)[2:4] <- c('p.L','p.M','p.H')
  
  p.L <- aggregate(p.L~Drug_Name, data = temp, mean)
  p.M <- aggregate(p.M~Drug_Name, data = temp, mean)
  p.H <- aggregate(p.H~Drug_Name, data = temp, mean)
  L.M <- merge(p.L, p.M, by='Drug_Name')
  L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
  L.M.H.mean$pred <- apply(L.M.H.mean[,2:4], 1, function(x){
    which.max(x)
  })
  drug.pred <- merge(drug.pred, L.M.H.mean, by = 'Drug_Name')
  drug.pred$pred <- ifelse(drug.pred$pred==1, 'L', drug.pred$pred)
  drug.pred$pred <- ifelse(drug.pred$pred==2, 'M', drug.pred$pred)
  drug.pred$pred <- ifelse(drug.pred$pred==3, 'H', drug.pred$pred)
  drug.pred$pred <- factor(drug.pred$pred, levels=c('L','M','H'))
  
  drug.pred.list[[j]] <- drug.pred
}

acc.drug <- unlist(lapply(drug.pred.list, function(x){
  cm <- confusionMatrix(x$pred, x$risk)
  cm[["overall"]][["Accuracy"]]
}))

# acc
cat('low 2.5%',sort(acc.drug)[ci.low], '\n')
cat('high 97.5%',sort(acc.drug)[ci.high], '\n')
cat('mean',mean(acc.drug), '\n')

##############################################################
# binary obs: H vs ML
##############################################################
auroc.obs.H.ML <- c()

for(j in 1:boost.count){
  
  stem$risk.H.ML <- ifelse(stem$risk=='M' | stem$risk=='L', 'ML', 'H')
  stem$risk.H.ML <- factor(stem$risk.H.ML, levels=c('H','ML'))
  
  stem[[paste('pred.H.ML', j, sep='_')]] <- 
    ifelse(stem[[paste('pred', j, sep='_')]]=='M' | stem[[paste('pred', j, sep='_')]]=='L', 'ML', 'H')
  stem[[paste('pred.H.ML', j, sep='_')]] <- factor(stem[[paste('pred.H.ML', j, sep='_')]], levels=c('H','ML'))
  
  stem[[paste('p.ML', j, sep="_")]] <- 1-stem[[paste('p.H', j, sep="_")]]
  
  cm <- confusionMatrix(stem$risk.H.ML, stem[[paste('pred.H.ML', j, sep='_')]], mode = 'everything'); cm
  
  # roc curve
  auroc.obs.H.ML[j] <- roc.curve(scores.class0 = stem[stem$risk.H.ML=='H',paste('p.H', j, sep="_")], 
                        scores.class1 = stem[stem$risk.H.ML=='ML',paste('p.H', j, sep="_")], curve=F)$auc
}

# acc
cat('low 2.5%',sort(auroc.obs.H.ML)[ci.low], '\n')
cat('high 97.5%',sort(auroc.obs.H.ML)[ci.high], '\n')
cat('mean',mean(auroc.obs.H.ML), '\n')

##############################################################
# binary drug: H vs ML
##############################################################

auroc.drug.H.ML <- unlist(lapply(drug.pred.list, function(drug.pred){
  
  drug.pred$risk.H.ML <- ifelse(drug.pred$risk=='M' | drug.pred$risk=='L', 'ML', 'H')
  drug.pred$risk.H.ML <- factor(drug.pred$risk.H.ML, levels=c('H','ML'))
  
  drug.pred$pred.H.ML <- ifelse(drug.pred$pred=='M' | drug.pred$pred=='L', 'ML', 'H')
  drug.pred$pred.H.ML <- factor(drug.pred$pred.H.ML, levels=c('H','ML'))
  
  drug.pred$p.ML <- 1-drug.pred$p.H
  
  # roc curve
  roc.curve(scores.class0 = drug.pred[drug.pred$risk.H.ML=='H','p.H'], 
            scores.class1 = drug.pred[drug.pred$risk.H.ML=='ML','p.H'], curve=F)$auc
}))

# acc
cat('low 2.5%',sort(auroc.drug.H.ML)[ci.low], '\n')
cat('high 97.5%',sort(auroc.drug.H.ML)[ci.high], '\n')
cat('mean',mean(auroc.drug.H.ML), '\n')

##############################################################
# binary obs: HM vs L
##############################################################
auroc.obs.HM.L <- c()

for(j in 1:boost.count){
  
  stem$risk.HM.L <- ifelse(stem$risk=='H' | stem$risk=='M', 'HM', 'L')
  stem$risk.HM.L <- factor(stem$risk.HM.L, levels=c('HM','L'))
  
  stem[[paste('pred.HM.L', j, sep='_')]] <- 
    ifelse(stem[[paste('pred', j, sep='_')]]=='H' | stem[[paste('pred', j, sep='_')]]=='M', 'HM', 'L')
  stem[[paste('pred.HM.L', j, sep='_')]] <- factor(stem[[paste('pred.HM.L', j, sep='_')]], levels=c('HM','L'))
  
  stem[[paste('p.HM', j, sep="_")]] <- 1-stem[[paste('p.L', j, sep="_")]]
  
  # roc curve
  auroc.obs.HM.L[j] <- roc.curve(scores.class0 = stem[stem$risk.HM.L=='HM',paste('p.HM', j, sep="_")], 
                        scores.class1 = stem[stem$risk.HM.L=='L',paste('p.HM', j, sep="_")], curve=F)$auc
}

# acc
cat('low 2.5%',sort(auroc.obs.HM.L)[ci.low], '\n')
cat('high 97.5%',sort(auroc.obs.HM.L)[ci.high], '\n')
cat('mean',mean(auroc.obs.HM.L), '\n')

##############################################################
# binary drug: HM vs L
##############################################################

auroc.drug.HM.L <- unlist(lapply(drug.pred.list, function(drug.pred){
  
  drug.pred$risk.HM.L <- ifelse(drug.pred$risk=='H' | drug.pred$risk=='M', 'HM', 'L')
  drug.pred$risk.HM.L <- factor(drug.pred$risk.HM.L, levels=c('HM','L'))
  
  drug.pred$pred.HM.L <- ifelse(drug.pred$pred=='H' | drug.pred$pred=='M', 'HM', 'L')
  drug.pred$pred.HM.L <- factor(drug.pred$pred.HM.L, levels=c('HM','L'))
  
  drug.pred$p.HM <- 1-drug.pred$p.L
  
  # roc curve
  roc.curve(scores.class0 = drug.pred[drug.pred$risk.HM.L=='HM','p.HM'], 
            scores.class1 = drug.pred[drug.pred$risk.HM.L=='L','p.HM'], curve=F)$auc
}))

# acc
cat('low 2.5%',sort(auroc.drug.HM.L)[ci.low], '\n')
cat('high 97.5%',sort(auroc.drug.HM.L)[ci.high], '\n')
cat('mean',mean(auroc.drug.HM.L), '\n')

######################################################################
# concordance index by obs
######################################################################
cindex.obs <- c()
concordance_index <- import('lifelines.utils')

# confusion matrix
for(j in 1:boost.count){
  cindex.obs[j] <- concordance_index$concordance_index(as.numeric(stem$risk), as.numeric(stem[[paste('pred',j, sep = '_')]]))
}

cat('low 2.5%',sort(cindex.obs)[ci.low], '\n')
cat('high 97.5%',sort(cindex.obs)[ci.high], '\n')
cat('mean',mean(cindex.obs), '\n')

######################################################################
# concordance index by drugs
######################################################################

cindex.drug <- unlist(lapply(drug.pred.list, function(drug.pred){
  concordance_index$concordance_index(as.numeric(drug.pred$risk), as.numeric(drug.pred$pred))
}))

cat('low 2.5%',sort(cindex.drug)[ci.low], '\n')
cat('high 97.5%',sort(cindex.drug)[ci.high], '\n')
cat('mean',mean(cindex.drug), '\n')

######################################################################
# construct data frame for visualization
######################################################################
dataset <- rep('stem', 8000)
model <- rep('randomforest', 8000)
by <- rep(c('obs', 'drug'), each=4000)
measurement <- rep(rep(c('acc', 'auroc.H.ML', 'auroc.HM.L', 'cindex'), each=1000), 2)
value <- c(acc.obs, auroc.obs.H.ML, auroc.obs.HM.L, cindex.obs, acc.drug, auroc.drug.H.ML, auroc.drug.HM.L, cindex.drug)
figure.stem.randomforest <- data.frame(dataset, model, by, measurement, value); dim(figure.stem.randomforest)
saveRDS(figure.stem.randomforest, 'result/figure_stem_randomforest.RDS')

######################################################################
# drug correct rate
######################################################################
rate <- lapply(drug.pred.list, function(drug.pred){
  drug.pred$risk==drug.pred$pred
})
rate.drug <- do.call(cbind, rate); dim(rate.drug)
rate.drug <- cbind(drug.pred.list[[1]][,c('Drug_Name', 'risk')], apply(rate.drug, 1, mean))
colnames(rate.drug)[3] <- 'rate.correct'
rate.drug$dataset <- 'stem'
rate.drug$model <- "randomforest"
rate.drug[rate.drug=='D,l Sotalol'] <- 'Sotalol'
saveRDS(rate.drug, file = "result/randomforest_stem_bootstrap_drug_rate.RDS") 

######################################################################
# drug confusion matrix
######################################################################

pred.drug <- lapply(drug.pred.list, function(drug.pred){
  drug.pred$pred
})
pred.drug <- do.call(cbind, pred.drug); dim(pred.drug)
pred.drug <- apply(pred.drug, 1, function(x){
  c(sum(x==1), sum(x==2), sum(x==3))
})
pred.drug <- t(pred.drug)
pred.drug <- data.frame(pred.drug)
colnames(pred.drug) <- c("L", "M", "H")
pred.drug <- cbind(rate.drug, pred.drug)

saveRDS(pred.drug, file = "result/randomforest_stem_bootstrap_drug_confusion.RDS") 


