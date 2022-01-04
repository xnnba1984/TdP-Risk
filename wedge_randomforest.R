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

# potential outlier drugs
#drug.outlier <- c('Domperidone', 'Cisapride')
#stem <- stem[!stem$Drug_Name%in%drug.outlier,]; dim(stem)

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

set.seed(38)
# leave one drug out
for(d in drug){
  index.test <- which(stem$Drug_Name==d)
  train <- stem[-index.test,]; dim(train)
  test <- stem[index.test,]; dim(test)
  
  # L vs M, H
  model1 <- randomForest(risk1~c3v1+c3v2+c3v3+c3v4+c3v6+c3v7+c3v8+c3v9+c3v10+c3v11+c3v12+c3v13+c3v14+c3v15+c3v16, 
                        data=train)
  #model1 <- randomForest(risk1~c3v9, data=train)
  pred1 <- predict(model1, test, type = 'prob')[,2]
  p.L <- 1 - pred1 
  stem[index.test, 'p.L'] <- p.L
  
  # H vs L, M
  model2 <- randomForest(risk2~c3v1+c3v2+c3v3+c3v4+c3v6+c3v7+c3v8+c3v9+c3v10+c3v11+c3v12+c3v13+c3v14+c3v15+c3v16, 
                        data=train)
  #model2 <- randomForest(risk2~c3v9, data=train)
  pred2 <- predict(model2, test, type = 'prob')[,2]
  p.H <- pred2
  p.M <- pred1-pred2
  p.M <- ifelse(p.M<0, 0, p.M)
  stem[index.test, 'p.M'] <- p.M
  stem[index.test, 'p.H'] <- p.H
} 

stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')],1,function(x) which.max(x)))

# acc by obs
l <- c('L','M','H')
stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))

# confusion matrix
cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything'); cm

#####################################################################
# acc by drug: average
#####################################################################
drug.pred <- unique(stem[,c('Drug_Name','risk')])
p.L <- aggregate(p.L~Drug_Name, data = stem, mean)
p.M <- aggregate(p.M~Drug_Name, data = stem, mean)
p.H <- aggregate(p.H~Drug_Name, data = stem, mean)
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

# cm on drugs
cm <- confusionMatrix(drug.pred$pred, drug.pred$risk, mode = 'everything'); cm

##############################################################
# binary obs: H vs ML
##############################################################
stem$risk.H.ML <- ifelse(stem$risk=='M' | stem$risk=='L', 'ML', 'H')
stem$risk.H.ML <- factor(stem$risk.H.ML, levels=c('H','ML'))

stem$pred.H.ML <- ifelse(stem$pred=='M' | stem$pred=='L', 'ML', 'H')
stem$pred.H.ML <- factor(stem$pred.H.ML, levels=c('H','ML'))

stem$p.ML <- 1-stem$p.H

cm <- confusionMatrix(stem$risk.H.ML, stem$pred.H.ML, mode = 'everything'); cm

# roc curve
plot(roc.curve(scores.class0 = stem[stem$risk.H.ML=='H','p.H'], scores.class1 = stem[stem$risk.H.ML=='ML','p.H'],
               curve=TRUE), xlab = 'False positive rate (1-specificity)', 
     ylab = 'True positive rate (sensitivity)')

##############################################################
# binary drug: H vs ML
##############################################################
drug.pred$risk.H.ML <- ifelse(drug.pred$risk=='M' | drug.pred$risk=='L', 'ML', 'H')
drug.pred$risk.H.ML <- factor(drug.pred$risk.H.ML, levels=c('H','ML'))

drug.pred$pred.H.ML <- ifelse(drug.pred$pred=='M' | drug.pred$pred=='L', 'ML', 'H')
drug.pred$pred.H.ML <- factor(drug.pred$pred.H.ML, levels=c('H','ML'))

drug.pred$p.ML <- 1-drug.pred$p.H

cm <- confusionMatrix(drug.pred$risk.H.ML, drug.pred$pred.H.ML, mode = 'everything'); cm

# roc curve
plot(roc.curve(scores.class0 = drug.pred[drug.pred$risk.H.ML=='H','p.H'], scores.class1 = drug.pred[drug.pred$risk.H.ML=='ML','p.H'],
               curve=TRUE), xlab = 'False positive rate (1-specificity)', 
     ylab = 'True positive rate (sensitivity)')

##############################################################
# binary obs: HM vs L
##############################################################
stem$risk.HM.L <- ifelse(stem$risk=='H' | stem$risk=='M', 'HM', 'L')
stem$risk.HM.L <- factor(stem$risk.HM.L, levels=c('HM','L'))

stem$pred.HM.L <- ifelse(stem$pred=='H' | stem$pred=='M', 'HM', 'L')
stem$pred.HM.L <- factor(stem$pred.HM.L, levels=c('HM','L'))

stem$p.HM <- 1-stem$p.L

cm <- confusionMatrix(stem$risk.HM.L, stem$pred.HM.L, mode = 'everything'); cm

# roc curve
plot(roc.curve(scores.class0 = stem[stem$risk.HM.L=='HM','p.HM'], scores.class1 = stem[stem$risk.HM.L=='L','p.HM'],
               curve=TRUE), xlab = 'False positive rate (1-specificity)', 
     ylab = 'True positive rate (sensitivity)')

##############################################################
# binary drug: HM vs L
##############################################################
drug.pred$risk.HM.L <- ifelse(drug.pred$risk=='H' | drug.pred$risk=='M', 'HM', 'L')
drug.pred$risk.HM.L <- factor(drug.pred$risk.HM.L, levels=c('HM','L'))

drug.pred$pred.HM.L <- ifelse(drug.pred$pred=='H' | drug.pred$pred=='M', 'HM', 'L')
drug.pred$pred.HM.L <- factor(drug.pred$pred.HM.L, levels=c('HM','L'))

drug.pred$p.HM <- 1-drug.pred$p.L

cm <- confusionMatrix(drug.pred$risk.HM.L, drug.pred$pred.HM.L, mode = 'everything'); cm

# roc curve
plot(roc.curve(scores.class0 = drug.pred[drug.pred$risk.HM.L=='HM','p.HM'], scores.class1 = drug.pred[drug.pred$risk.HM.L=='L','p.HM'],
               curve=TRUE), xlab = 'False positive rate (1-specificity)', 
     ylab = 'True positive rate (sensitivity)')

##############################################################
# concordance index
##############################################################
concordance_index <- import('lifelines.utils')
# obs
concordance_index$concordance_index(as.numeric(stem$risk), as.numeric(stem$pred))
# drugs
concordance_index$concordance_index(as.numeric(drug.pred$risk), as.numeric(drug.pred$pred))

