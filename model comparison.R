library(have)
library(caret)
library(MASS)
library(PRROC)
library(reticulate)
library(nnet)
library(randomForest)

##################################################################################################
# stem cell
##################################################################################################

set.seed(2021)
stem <- data.frame(read_sas('allpred2.sas7bdat'))

# missing value
sum(is.na(stem))

stem$risk <- factor(stem$risk, levels=c('L','M','H')); str(stem)
drug <- unique(stem$Drug_Name)

# impute missing values
preProcValues <- preProcess(stem[,6:12], method = c("bagImpute"))
stem[,6:12] <- predict(preProcValues, stem[,6:12])
str(stem)

##################################################################################################
# mutinomial logistic regression
##################################################################################################
# leave one drug out
for(d in drug){
  index.test <- which(stem$Drug_Name==d)
  train <- stem[-index.test,]; dim(train)
  test <- stem[index.test,]; dim(test)
  
  model <- multinom(risk ~ pred7+arrhythmia+fst_pro+maxpro+
                      foldaym+foldprolong+AMtwooutoffive, data = train, trace=F)
  
  # Predicting the values for test dataset
  pred <- predict(model, newdata = test, "probs")
  p.L <- pred[,1]
  p.M <- pred[,2]
  p.H <- pred[,3]
  stem[index.test, 'p.L'] <- p.L
  stem[index.test, 'p.M'] <- p.M
  stem[index.test, 'p.H'] <- p.H
} 
stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')],1,function(x) which.max(x)))

# acc by obs
l <- c('L','M','H')
stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))

# confusion matrix
cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything'); cm

# acc per drug
stem$correct <- stem$risk==stem$pred
drug.pred <- aggregate(correct~Drug_Name, data = stem, mean)

# drug risk
drug.risk <- unique(stem[,c('Drug_Name','risk')])
drug.pred <- merge(drug.pred, drug.risk, by = 'Drug_Name')

#acc by drug: average
p.L <- aggregate(p.L~Drug_Name, data = stem, mean)
p.M <- aggregate(p.M~Drug_Name, data = stem, mean)
p.H <- aggregate(p.H~Drug_Name, data = stem, mean)
L.M <- merge(p.L, p.M, by='Drug_Name')
L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
drug.pred$pred <- apply(L.M.H.mean[,2:4], 1, function(x){
  which.max(x)
})
drug.pred$pred <- ifelse(drug.pred$pred==1, 'L', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==2, 'M', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==3, 'H', drug.pred$pred)
drug.pred$pred <- factor(drug.pred$pred, levels=c('L','M','H'))
drug.pred <- merge(drug.pred, L.M.H.mean, by = 'Drug_Name')

# cm on drugs
cm <- confusionMatrix(drug.pred$pred, drug.pred$risk); cm

##################################################################################################
# classical ordinal logistic regression
##################################################################################################
set.seed(2021)
stem <- data.frame(read_sas('allpred2.sas7bdat'))

# missing value
sum(is.na(stem))

stem$risk <- factor(stem$risk, levels=c('L','M','H'), ordered = T); str(stem)
drug <- unique(stem$Drug_Name)

# impute missing values
preProcValues <- preProcess(stem[,6:12], method = c("bagImpute"))
stem[,6:12] <- predict(preProcValues, stem[,6:12])
stem$risk <-  factor(stem$risk, levels = c("L", "M", "H"), ordered = TRUE)
str(stem)

for(d in drug){
  index.test <- which(stem$Drug_Name==d)
  train <- stem[-index.test,]; dim(train)
  test <- stem[index.test,]; dim(test)
  
  model <- polr(risk ~ pred7+arrhythmia+fst_pro+maxpro+
                  foldaym+foldprolong+AMtwooutoffive, data = train)
  
  # Predicting the values for test dataset
  pred <- predict(model, newdata = test, "probs")
  p.L <- pred[,1]
  p.M <- pred[,2]
  p.H <- pred[,3]
  stem[index.test, 'p.L'] <- p.L
  stem[index.test, 'p.M'] <- p.M
  stem[index.test, 'p.H'] <- p.H
} 
stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')],1,function(x) which.max(x)))

# acc by obs
l <- c('L','M','H')
stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))

# confusion matrix
cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything'); cm

# acc per drug
stem$risk <- factor(stem$risk, ordered=F)
#str(stem)
stem$correct <- stem$risk==stem$pred
drug.pred <- aggregate(correct~Drug_Name, data = stem, mean)

# drug risk
drug.risk <- unique(stem[,c('Drug_Name','risk')])
drug.pred <- merge(drug.pred, drug.risk, by = 'Drug_Name')

#acc by drug: average
p.L <- aggregate(p.L~Drug_Name, data = stem, mean)
p.M <- aggregate(p.M~Drug_Name, data = stem, mean)
p.H <- aggregate(p.H~Drug_Name, data = stem, mean)
L.M <- merge(p.L, p.M, by='Drug_Name')
L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
drug.pred$pred <- apply(L.M.H.mean[,2:4], 1, function(x){
  which.max(x)
})
drug.pred$pred <- ifelse(drug.pred$pred==1, 'L', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==2, 'M', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==3, 'H', drug.pred$pred)
drug.pred$pred <- factor(drug.pred$pred, levels=c('L','M','H'))
drug.pred <- merge(drug.pred, L.M.H.mean, by = 'Drug_Name')

# cm on drugs
cm <- confusionMatrix(drug.pred$pred, drug.pred$risk); cm

##################################################################################################
# multinomial random forest
##################################################################################################
set.seed(2021)
stem <- data.frame(read_sas('allpred2.sas7bdat'))

# impute missing values
preProcValues <- preProcess(stem[,6:12], method = c("bagImpute"))
stem[,6:12] <- predict(preProcValues, stem[,6:12])
drug <- unique(stem$Drug_Name)
stem$risk <- factor(stem$risk, levels=c('L','M','H')); str(stem)

set.seed(8)
for(d in drug){
  index.test <- which(stem$Drug_Name==d)
  train <- stem[-index.test,]; dim(train)
  test <- stem[index.test,]; dim(test)
  
  model <- randomForest(risk~pred7+arrhythmia+fst_pro+maxpro+
                           foldaym+foldprolong+AMtwooutoffive, data=train)
  pred <- predict(model, test, type = 'prob')
  p.L <- pred[,1]
  p.M <- pred[,2]
  p.H <- pred[,3]
  stem[index.test, 'p.L'] <- p.L
  stem[index.test, 'p.M'] <- p.M
  stem[index.test, 'p.H'] <- p.H
} 
stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')],1,function(x) which.max(x)))

# acc by obs
l <- c('L','M','H')
stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))

# confusion matrix
cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything'); cm

# acc per drug
stem$correct <- stem$risk==stem$pred
drug.pred <- aggregate(correct~Drug_Name, data = stem, mean)

# drug risk
drug.risk <- unique(stem[,c('Drug_Name','risk')])
drug.pred <- merge(drug.pred, drug.risk, by = 'Drug_Name')

#acc by drug: average
p.L <- aggregate(p.L~Drug_Name, data = stem, mean)
p.M <- aggregate(p.M~Drug_Name, data = stem, mean)
p.H <- aggregate(p.H~Drug_Name, data = stem, mean)
L.M <- merge(p.L, p.M, by='Drug_Name')
L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
drug.pred$pred <- apply(L.M.H.mean[,2:4], 1, function(x){
  which.max(x)
})
drug.pred$pred <- ifelse(drug.pred$pred==1, 'L', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==2, 'M', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==3, 'H', drug.pred$pred)
drug.pred$pred <- factor(drug.pred$pred, levels=c('L','M','H'))
drug.pred <- merge(drug.pred, L.M.H.mean, by = 'Drug_Name')

# cm on drugs
cm <- confusionMatrix(drug.pred$pred, drug.pred$risk); cm

##################################################################################################
# wedge
##################################################################################################

##################################################################################################
# mutinomial logistic regression
##################################################################################################
set.seed(37)
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

# impute missing values
preProcValues <- preProcess(stem[,3:17], method = c("bagImpute"))
stem[,3:17] <- predict(preProcValues, stem[,3:17])
drug <- unique(stem$Drug_Name)

for(d in drug){
  index.test <- which(stem$Drug_Name==d)
  train <- stem[-index.test,]; dim(train)
  test <- stem[index.test,]; dim(test)
  model <- multinom(risk ~ c3v1+c3v2+c3v3+c3v4+c3v6+c3v7+c3v8+c3v9+c3v10+c3v11+c3v12+c3v13+c3v14+c3v15+c3v16,
                    data = train, trace=F)
  # Predicting the values for test dataset
  pred <- predict(model, newdata = test, "probs")
  p.L <- pred[,1]
  p.M <- pred[,2]
  p.H <- pred[,3]
  stem[index.test, 'p.L'] <- p.L
  stem[index.test, 'p.M'] <- p.M
  stem[index.test, 'p.H'] <- p.H
} 

stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')],1,function(x) which.max(x)))

# acc by obs
l <- c('L','M','H')
stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))

# confusion matrix
cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything'); cm

# acc by drug: average
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

##################################################################################################
# classical ordinal logistic regression
##################################################################################################
library(rms)
acc.obs <- c()
acc.drug <- c()
set.seed(56)
stem <- read.csv('dt17c3max.csv', header = T)
stem <- stem[,-c(2,3,21)]
stem <- stem[,-which(colnames(stem)=="c3v5")]
colnames(stem)[1:2] <- c("Drug_Name", 'risk')

stem$risk <- ifelse(stem$risk=='low', 'L', stem$risk)
stem$risk <- ifelse(stem$risk=='int', 'M', stem$risk)
stem$risk <- ifelse(stem$risk=='high', 'H', stem$risk)
stem <- stem[stem$Drug_Name!='Placebo',]
stem$risk <- factor(stem$risk, levels=c('L','M','H'), ordered = T)

# missing value
sum(!complete.cases(stem))
sum(is.na(stem))

# impute missing values
preProcValues <- preProcess(stem[,3:17], method = c("bagImpute"))
stem[,3:17] <- predict(preProcValues, stem[,3:17])
drug <- unique(stem$Drug_Name)

for(d in drug){
  index.test <- which(stem$Drug_Name==d)
  train <- stem[-index.test,]; dim(train)
  test <- stem[index.test,]; dim(test)
  
  model <- orm(risk ~ c3v1+c3v2+c3v3+c3v4+c3v6+c3v7+c3v8+c3v9+c3v10+c3v11+c3v12+c3v13+c3v14+c3v15+c3v16, 
                data = train)

  # Predicting the values for test dataset
  pred <- predict(model, newdata = test, type = 'fitted')
  p.L <- 1 - pred[,1]
  p.M <- pred[,1] - pred[,2]
  p.H <- pred[,2]
  stem[index.test, 'p.L'] <- p.L
  stem[index.test, 'p.M'] <- p.M
  stem[index.test, 'p.H'] <- p.H
} 
stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')],1,function(x) which.max(x)))

# acc by obs
l <- c('L','M','H')
stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))

# confusion matrix
cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything'); cm

# acc per drug
stem$risk <- factor(stem$risk, ordered=F)
#str(stem)
stem$correct <- stem$risk==stem$pred
drug.pred <- aggregate(correct~Drug_Name, data = stem, mean)

# drug risk
drug.risk <- unique(stem[,c('Drug_Name','risk')])
drug.pred <- merge(drug.pred, drug.risk, by = 'Drug_Name')

#acc by drug: average
p.L <- aggregate(p.L~Drug_Name, data = stem, mean)
p.M <- aggregate(p.M~Drug_Name, data = stem, mean)
p.H <- aggregate(p.H~Drug_Name, data = stem, mean)
L.M <- merge(p.L, p.M, by='Drug_Name')
L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
drug.pred$pred <- apply(L.M.H.mean[,2:4], 1, function(x){
  which.max(x)
})
drug.pred$pred <- ifelse(drug.pred$pred==1, 'L', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==2, 'M', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==3, 'H', drug.pred$pred)
drug.pred$pred <- factor(drug.pred$pred, levels=c('L','M','H'))
drug.pred <- merge(drug.pred, L.M.H.mean, by = 'Drug_Name')

# cm on drugs
cm <- confusionMatrix(drug.pred$pred, drug.pred$risk); cm

##################################################################################################
# multinomial random forest
##################################################################################################
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

# impute missing values
preProcValues <- preProcess(stem[,3:17], method = c("bagImpute"))
stem[,3:17] <- predict(preProcValues, stem[,3:17])
drug <- unique(stem$Drug_Name)

set.seed(11)
for(d in drug){
  index.test <- which(stem$Drug_Name==d)
  train <- stem[-index.test,]; dim(train)
  test <- stem[index.test,]; dim(test)
    
  model <- randomForest(risk~c3v1+c3v2+c3v3+c3v4+c3v6+c3v7+c3v8+c3v9+c3v10+c3v11+c3v12+c3v13+c3v14+c3v15+c3v16, 
                        data=train)
  pred <- predict(model, test, type = 'prob')
  p.L <- pred[,1]
  p.M <- pred[,2]
  p.H <- pred[,3]
  stem[index.test, 'p.L'] <- p.L
  stem[index.test, 'p.M'] <- p.M
  stem[index.test, 'p.H'] <- p.H
} 
stem$pred <- unlist(apply(stem[,c('p.L','p.M','p.H')],1,function(x) which.max(x)))

# acc by obs
l <- c('L','M','H')
stem$pred <- factor(l[stem$pred], levels=c('L','M','H'))

# confusion matrix
cm <- confusionMatrix(stem$pred, stem$risk, mode = 'everything'); cm
acc.obs[seed] <- cm$overall[1]
  
# acc per drug
stem$correct <- stem$risk==stem$pred
drug.pred <- aggregate(correct~Drug_Name, data = stem, mean)
  
# drug risk
drug.risk <- unique(stem[,c('Drug_Name','risk')])
drug.pred <- merge(drug.pred, drug.risk, by = 'Drug_Name')
  
#acc by drug: average
p.L <- aggregate(p.L~Drug_Name, data = stem, mean)
p.M <- aggregate(p.M~Drug_Name, data = stem, mean)
p.H <- aggregate(p.H~Drug_Name, data = stem, mean)
L.M <- merge(p.L, p.M, by='Drug_Name')
L.M.H.mean <- merge(L.M, p.H, by='Drug_Name')
drug.pred$pred <- apply(L.M.H.mean[,2:4], 1, function(x){
  which.max(x)
})
drug.pred$pred <- ifelse(drug.pred$pred==1, 'L', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==2, 'M', drug.pred$pred)
drug.pred$pred <- ifelse(drug.pred$pred==3, 'H', drug.pred$pred)
drug.pred$pred <- factor(drug.pred$pred, levels=c('L','M','H'))
drug.pred <- merge(drug.pred, L.M.H.mean, by = 'Drug_Name')
  
# cm on drugs
cm <- confusionMatrix(drug.pred$pred, drug.pred$risk); cm