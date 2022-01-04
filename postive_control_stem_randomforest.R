library(haven)
library(caret)
library(MASS)
library(PRROC)
library(reticulate)

set.seed(2021)
stem <- data.frame(read_sas('allpred2.sas7bdat'))
######################################################
# ordinal logistic
######################################################
# missing value
sum(!complete.cases(stem))
sum(is.na(stem))

stem$risk <- factor(stem$risk, levels=c('L','M','H'))
stem$Cell_Type <- as.numeric(as.factor(stem$Cell_Type))-1
drug <- unique(stem$Drug_Name)

# impute missing values
preProcValues <- preProcess(stem[,6:12], method = c("bagImpute"))
stem[,6:12] <- predict(preProcValues, stem[,6:12])

# p(risk>low)
stem$risk1 <- factor(ifelse(stem$risk=='L', 0, 1))

# p(risk>middle)
stem$risk2 <- factor(ifelse(stem$risk=='H', 1, 0))

#############################################
# positive control: H risk
#############################################
pc.acc <- c()
acc.obs.all <- c()
acc.obs.correct <- c()
acc.drug.all <- c()
acc.drug.correct <- c()
auroc.H.ML.obs.all <- c()
auroc.H.ML.drug.all <- c()
auroc.H.ML.obs.correct <- c()
auroc.H.ML.drug.correct <- c()
auroc.HM.L.obs.all <- c()
auroc.HM.L.drug.all <- c()
auroc.HM.L.obs.correct <- c()
auroc.HM.L.drug.correct <- c()
cindex.obs.all <- c()
cindex.obs.correct <- c()
cindex.drug.all <- c()
cindex.drug.correct <- c()
drug.pcs <- c('D,l Sotalol')
drug.risk <- unique(stem[,c('Drug_Name','risk')])

set.seed(2021)
for(i in 1:length(drug.pcs)){
  drug.pc <- drug.pcs[i]
  print(drug.pc)
  stem$risk.pc <- drug.risk[drug.risk$Drug_Name==drug.pc, 'risk']
  index.pc <- which(stem$Drug_Name==drug.pc)
  test.pc <- stem[index.pc,]; dim(test.pc)
  drug <- drug[which(drug!=drug.pc)]
  
  # leave one drug out
  # predict test + pc
  for(d in drug){
    #d <- drug[1]
    #print(d)
    index.test <- which(stem$Drug_Name==d)
    train <- stem[-c(index.test, index.pc),]; dim(train)
    test <- stem[index.test,]; dim(test)
    
    # L vs M, H
    model1 <- randomForest(risk1~pred7+arrhythmia+fst_pro+maxpro+
                                  foldaym+foldprolong+AMtwooutoffive, data=train)
    pred1 <- predict(model1, test, type = 'prob')[,2]
    p.L <- 1 - pred1 
    stem[index.test, 'p.L'] <- p.L
    
    # predict on pc
    pred1.pc <- predict(model1, test.pc, type = 'prob')[,2]
    p.L.pc <- 1 - pred1.pc 
    stem[index.test, 'p.L.pc'] <- p.L.pc
    
    # H vs L, M
    model2 <- randomForest(risk2~pred7+arrhythmia+fst_pro+maxpro+
                                  foldaym+foldprolong+AMtwooutoffive, data=train)
    pred2 <- predict(model2, test, type = 'prob')[,2]
    p.H <- pred2
    p.M <- pred1-pred2
    p.M <- ifelse(p.M<0, 0, p.M)
    stem[index.test, 'p.M'] <- p.M
    stem[index.test, 'p.H'] <- p.H
    
    # predict on pc
    pred2.pc <- predict(model2, test.pc, type = 'prob')[,2]
    p.H.pc <- pred2.pc
    p.M.pc <- pred1.pc - pred2.pc
    p.M.pc <- ifelse(p.M.pc<0, 0, p.M.pc)
    stem[index.test, 'p.M.pc'] <- p.M.pc
    stem[index.test, 'p.H.pc'] <- p.H.pc
  } 
  stem <- stem[complete.cases(stem),]; dim(stem)
  stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')], 1, function(x) which.max(x)))
  stem$pred.pc <- unlist(apply(stem[,c('p.L.pc','p.M.pc','p.H.pc')], 1, function(x) which.max(x)))
  l <- c('L','M','H')
  stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))
  stem$pred.pc <- factor(l[stem$pred.pc], levels=c('L','M','H'))
  
  ##############################################################
  # acc obs all
  ##############################################################
  cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything')
  acc.obs.all[i] <- cm[["overall"]][["Accuracy"]]
  
  # acc on pc (obs)
  pc.acc[i] <- mean(stem$pred.pc=='H')
  
  ##############################################################
  # acc obs correct
  ##############################################################
  stem.correct <- stem[stem$pred.pc=='H',]; dim(stem.correct)
  cm <- confusionMatrix(stem.correct$pred, stem.correct$risk, mode = 'everything')
  acc.obs.correct[i] <- cm[["overall"]][["Accuracy"]]
  
  ##############################################################
  #acc drug all
  ##############################################################
  drug.pred <- unique(stem[,c('Drug_Name','risk')]); dim(drug.pred)
  p.L <- aggregate(p.L~Drug_Name, data = stem, mean)
  p.M <- aggregate(p.M~Drug_Name, data = stem, mean)
  p.H <- aggregate(p.H~Drug_Name, data = stem, mean)
  L.M <- merge(p.L, p.M, by='Drug_Name')
  L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
  drug.pred <- merge(drug.pred, L.M.H.mean, by = 'Drug_Name')
  drug.pred$pred <- apply(drug.pred[,3:5], 1, function(x){
    which.max(x)
  })
  drug.pred$pred <- ifelse(drug.pred$pred==1, 'L', drug.pred$pred)
  drug.pred$pred <- ifelse(drug.pred$pred==2, 'M', drug.pred$pred)
  drug.pred$pred <- ifelse(drug.pred$pred==3, 'H', drug.pred$pred)
  drug.pred$pred <- factor(drug.pred$pred, levels=c('L','M','H'))
  
  # cm on drugs
  cm <- confusionMatrix(drug.pred$pred, drug.pred$risk)
  acc.drug.all[i] <- cm[["overall"]][["Accuracy"]]
  
  ##############################################################
  #acc drug correct
  ##############################################################
  drug.pred.correct <- unique(stem[,c('Drug_Name','risk')]); dim(drug.pred.correct)
  p.L <- aggregate(p.L~Drug_Name, data = stem.correct, mean)
  p.M <- aggregate(p.M~Drug_Name, data = stem.correct, mean)
  p.H <- aggregate(p.H~Drug_Name, data = stem.correct, mean)
  L.M <- merge(p.L, p.M, by='Drug_Name')
  L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
  drug.pred.correct <- merge(drug.pred.correct, L.M.H.mean, by = 'Drug_Name')
  drug.pred.correct$pred <- apply(drug.pred.correct[,3:5], 1, function(x){
    which.max(x)
  })
  drug.pred.correct$pred <- ifelse(drug.pred.correct$pred==1, 'L', drug.pred.correct$pred)
  drug.pred.correct$pred <- ifelse(drug.pred.correct$pred==2, 'M', drug.pred.correct$pred)
  drug.pred.correct$pred <- ifelse(drug.pred.correct$pred==3, 'H', drug.pred.correct$pred)
  drug.pred.correct$pred <- factor(drug.pred.correct$pred, levels=c('L','M','H'))
  
  # cm on drugs
  cm <- confusionMatrix(drug.pred.correct$pred, drug.pred.correct$risk)
  acc.drug.correct[i] <- cm[["overall"]][["Accuracy"]]
  
  ##############################################################
  # binary obs: H vs ML all
  ##############################################################
  stem$risk.H.ML <- ifelse(stem$risk=='M' | stem$risk=='L', 'ML', 'H')
  stem$risk.H.ML <- factor(stem$risk.H.ML, levels=c('H','ML'))
  
  stem$pred.H.ML <- ifelse(stem$pred=='M' | stem$pred=='L', 'ML', 'H')
  stem$pred.H.ML <- factor(stem$pred.H.ML, levels=c('H','ML'))
  
  stem$p.ML <- 1-stem$p.H
  
  # roc curve
  r <- roc.curve(scores.class0 = stem[stem$risk.H.ML=='H','p.H'], scores.class1 = stem[stem$risk.H.ML=='ML','p.H'], curve=F)
  auroc.H.ML.obs.all[i] <- r$auc
  
  ##############################################################
  # binary obs: H vs ML correct
  ##############################################################
  stem.correct$risk.H.ML <- ifelse(stem.correct$risk=='M' | stem.correct$risk=='L', 'ML', 'H')
  stem.correct$risk.H.ML <- factor(stem.correct$risk.H.ML, levels=c('H','ML'))
  
  stem.correct$pred.H.ML <- ifelse(stem.correct$pred=='M' | stem.correct$pred=='L', 'ML', 'H')
  stem.correct$pred.H.ML <- factor(stem.correct$pred.H.ML, levels=c('H','ML'))
  
  stem.correct$p.ML <- 1-stem.correct$p.H
  
  # roc curve
  r <- roc.curve(scores.class0 = stem.correct[stem.correct$risk.H.ML=='H','p.H'], 
                 scores.class1 = stem.correct[stem.correct$risk.H.ML=='ML','p.H'], curve=F)
  auroc.H.ML.obs.correct[i] <- r$auc
  
  ##############################################################
  # binary drug: H vs ML all
  ##############################################################
  drug.pred$risk.H.ML <- ifelse(drug.pred$risk=='M' | drug.pred$risk=='L', 'ML', 'H')
  drug.pred$risk.H.ML <- factor(drug.pred$risk.H.ML, levels=c('H','ML'))
  
  drug.pred$pred.H.ML <- ifelse(drug.pred$pred=='M' | drug.pred$pred=='L', 'ML', 'H')
  drug.pred$pred.H.ML <- factor(drug.pred$pred.H.ML, levels=c('H','ML'))
  
  drug.pred$p.ML <- 1-drug.pred$p.H
  
  # roc curve
  r <- roc.curve(scores.class0 = drug.pred[drug.pred$risk.H.ML=='H','p.H'], scores.class1 = drug.pred[drug.pred$risk.H.ML=='ML','p.H'],
                 curve=F)
  auroc.H.ML.drug.all[i] <- r$auc
  
  ##############################################################
  # binary drug: H vs ML correct
  ##############################################################
  drug.pred.correct$risk.H.ML <- ifelse(drug.pred.correct$risk=='M' | drug.pred.correct$risk=='L', 'ML', 'H')
  drug.pred.correct$risk.H.ML <- factor(drug.pred.correct$risk.H.ML, levels=c('H','ML'))
  
  drug.pred.correct$pred.H.ML <- ifelse(drug.pred.correct$pred=='M' | drug.pred.correct$pred=='L', 'ML', 'H')
  drug.pred.correct$pred.H.ML <- factor(drug.pred.correct$pred.H.ML, levels=c('H','ML'))
  
  drug.pred.correct$p.ML <- 1-drug.pred.correct$p.H
  
  # roc curve
  r <- roc.curve(scores.class0 = drug.pred.correct[drug.pred.correct$risk.H.ML=='H','p.H'], scores.class1 = drug.pred.correct[drug.pred.correct$risk.H.ML=='ML','p.H'],
                 curve=F)
  auroc.H.ML.drug.correct[i] <- r$auc
  
  ##############################################################
  # binary obs: HM vs L all
  ##############################################################
  stem$risk.HM.L <- ifelse(stem$risk=='H' | stem$risk=='M', 'HM', 'L')
  stem$risk.HM.L <- factor(stem$risk.HM.L, levels=c('HM','L'))
  
  stem$pred.HM.L <- ifelse(stem$pred=='H' | stem$pred=='M', 'HM', 'L')
  stem$pred.HM.L <- factor(stem$pred.HM.L, levels=c('HM','L'))
  
  stem$p.HM <- 1-stem$p.L
  
  # roc curve
  r <- roc.curve(scores.class0 = stem[stem$risk.HM.L=='HM','p.HM'], scores.class1 = stem[stem$risk.HM.L=='L','p.HM'],
                 curve=F)
  auroc.HM.L.obs.all[i] <- r$auc
  
  ##############################################################
  # binary obs: HM vs L correct
  ##############################################################
  stem.correct$risk.HM.L <- ifelse(stem.correct$risk=='H' | stem.correct$risk=='M', 'HM', 'L')
  stem.correct$risk.HM.L <- factor(stem.correct$risk.HM.L, levels=c('HM','L'))
  
  stem.correct$pred.HM.L <- ifelse(stem.correct$pred=='H' | stem.correct$pred=='M', 'HM', 'L')
  stem.correct$pred.HM.L <- factor(stem.correct$pred.HM.L, levels=c('HM','L'))
  
  stem.correct$p.HM <- 1-stem.correct$p.L
  
  # roc curve
  r <- roc.curve(scores.class0 = stem.correct[stem.correct$risk.HM.L=='HM','p.HM'], scores.class1 = stem.correct[stem.correct$risk.HM.L=='L','p.HM'],
                 curve=F)
  auroc.HM.L.obs.correct[i] <- r$auc
  
  ##############################################################
  # binary drug: HM vs L all
  ##############################################################
  drug.pred$risk.HM.L <- ifelse(drug.pred$risk=='H' | drug.pred$risk=='M', 'HM', 'L')
  drug.pred$risk.HM.L <- factor(drug.pred$risk.HM.L, levels=c('HM','L'))
  
  drug.pred$pred.HM.L <- ifelse(drug.pred$pred=='H' | drug.pred$pred=='M', 'HM', 'L')
  drug.pred$pred.HM.L <- factor(drug.pred$pred.HM.L, levels=c('HM','L'))
  
  drug.pred$p.HM <- 1-drug.pred$p.L
  
  cm <- confusionMatrix(drug.pred$risk.HM.L, drug.pred$pred.HM.L, mode = 'everything'); cm
  
  # roc curve
  r <- roc.curve(scores.class0 = drug.pred[drug.pred$risk.HM.L=='HM','p.HM'], scores.class1 = drug.pred[drug.pred$risk.HM.L=='L','p.HM'],
                 curve=TRUE)
  auroc.HM.L.drug.all[i] <- r$auc
  
  ##############################################################
  # binary drug: HM vs L correct
  ##############################################################
  drug.pred.correct$risk.HM.L <- ifelse(drug.pred.correct$risk=='H' | drug.pred.correct$risk=='M', 'HM', 'L')
  drug.pred.correct$risk.HM.L <- factor(drug.pred.correct$risk.HM.L, levels=c('HM','L'))
  
  drug.pred.correct$pred.HM.L <- ifelse(drug.pred.correct$pred=='H' | drug.pred.correct$pred=='M', 'HM', 'L')
  drug.pred.correct$pred.HM.L <- factor(drug.pred.correct$pred.HM.L, levels=c('HM','L'))
  
  drug.pred.correct$p.HM <- 1-drug.pred.correct$p.L
  
  # roc curve
  r <- roc.curve(scores.class0 = drug.pred.correct[drug.pred.correct$risk.HM.L=='HM','p.HM'], scores.class1 = drug.pred.correct[drug.pred.correct$risk.HM.L=='L','p.HM'],
                 curve=TRUE)
  auroc.HM.L.drug.correct[i] <- r$auc
  
  ##############################################################
  # concordance index
  ##############################################################
  concordance_index <- import('lifelines.utils')
  # by obs all
  cindex.obs.all[i] <- concordance_index$concordance_index(as.numeric(stem$risk), as.numeric(stem$pred))
  # by obs correct
  cindex.obs.correct[i] <- concordance_index$concordance_index(as.numeric(stem.correct$risk), as.numeric(stem.correct$pred))
  # by drugs
  cindex.drug.all[i] <- concordance_index$concordance_index(as.numeric(drug.pred$risk), as.numeric(drug.pred$pred))
  # by drugs
  cindex.drug.correct[i] <- concordance_index$concordance_index(as.numeric(drug.pred.correct$risk), as.numeric(drug.pred.correct$pred))
}

result <- cbind(acc.obs.all, acc.obs.correct,acc.drug.all,acc.drug.correct, auroc.H.ML.obs.all, auroc.H.ML.drug.all, auroc.H.ML.obs.correct,
                auroc.H.ML.drug.correct, auroc.HM.L.obs.all,auroc.HM.L.drug.all, auroc.HM.L.obs.correct, auroc.HM.L.drug.correct,cindex.obs.all, 
                cindex.obs.correct, cindex.drug.all,cindex.drug.correct, pc.acc)
result <- data.frame(result)
result$postive_control <- drug.pcs


