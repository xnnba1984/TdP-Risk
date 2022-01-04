library(ggplot2)
library(tidyverse)
library(ggh4x)

######################################################################################
# model uncertainty
######################################################################################
stem.logistic <- readRDS('result/figure_stem_logistic.RDS')
stem.randomforest <- readRDS('result/figure_stem_randomforest.RDS')
wedge.logistic <- readRDS('result/figure_wedge_logistic.RDS')
wedge.randomforest <- readRDS('result/figure_wedge_randomforest.RDS')

result <- rbind(stem.logistic, stem.randomforest, wedge.logistic, wedge.randomforest); dim(result)
result[result=='stem'] <- 'Stem cell'
result[result=='wedge'] <- 'Wedge'
result[result=='obs'] <- 'By observations'
result[result=='drug'] <- 'By drugs'
result[result=='acc'] <- 'Accuracy'
result[result=='auroc.H.ML'] <- 'AUROC\n(H vs ML)'
result[result=='auroc.HM.L'] <- 'AUROC\n(HM vs L)'
result[result=='cindex'] <- 'Concordance\nindex'
result[result=='logistic'] <- 'Ordinal logistic regression'
result[result=='randomforest'] <- 'Ordinal random forest'
result$by = factor(result$by, levels=c('By observations','By drugs'))

# only plot by drugs
ggplot(result[result$by=='By drugs',], aes(x=measurement, y=value, fill=model)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(dataset)) +
  theme_bw() +
  theme(text=element_text(size = 25), axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size=15))

# only plot by observations
ggplot(result[result$by=='By observations',], aes(x=measurement, y=value, fill=model)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(dataset)) +
  theme_bw() +
  theme(text=element_text(size = 25), axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size=15))

######################################################################################
# drug correct rate
######################################################################################
rate.stem.logistic <- readRDS("result/logistic_stem_bootstrap_drug_rate.RDS"); head(rate.stem.logistic)
rate.stem.rf <- readRDS("result/randomforest_stem_bootstrap_drug_rate.RDS"); head(rate.stem.rf)
rate.wedge.logistic <- readRDS("result/logistic_wedge_bootstrap_drug_rate.RDS"); head(rate.wedge.logistic)
rate.wedge.rf <- readRDS("result/randomforest_wedge_bootstrap_drug_rate.RDS"); head(rate.wedge.rf)

result.rate <- rbind(rate.stem.logistic, rate.stem.rf, rate.wedge.logistic, rate.wedge.rf); dim(result.rate)
result.rate$Drug_Name <- as.factor(result.rate$Drug_Name)
result.rate$risk = factor(result.rate$risk, levels=c('H', 'M', 'L'))
result.rate[result.rate=='stem'] <- 'Stem cell'
result.rate[result.rate=='wedge'] <- 'Wedge'
result.rate[result.rate=='logistic'] <- 'Ordinal logistic regression'
result.rate[result.rate=='randomforest'] <- 'Ordinal random forest'

ggplot(result.rate, aes(x=rate.correct, y=reorder(Drug_Name, rate.correct))) +
  geom_line(aes(group = Drug_Name), linetype = "dashed")+
  geom_point(aes(color=model, shape=model), size=4.5, alpha=.7) +
  facet_grid(risk~dataset, scales = 'free', space = 'free')+
  xlab('\nCorrect rate') +
  scale_x_continuous(label=scales::percent) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(), text=element_text(size = 25),
        axis.title.y=element_blank(), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13),
        legend.title=element_blank(), legend.position="bottom", axis.title.x = element_text(size = 15),
        legend.text = element_text(size=15))

######################################################################################
# permutation predictor importance
######################################################################################
imp.stem <- readRDS('result/importance_stem_rf_obs.rds')
imp.wedge <- readRDS('result/importance_wedge_rf_obs.rds')
result.imp <- rbind(imp.stem, imp.wedge)
result.imp[result.imp=='stem'] <- 'Stem cell'
result.imp[result.imp=='wedge'] <- 'Wedge'
result.imp <- result.imp[result.imp$imp>0.1,]

ggplot(result.imp, aes(x=reorder(predictor, imp), y=imp)) +
  geom_bar(stat="identity", fill='steelblue', width=0.7) +
  coord_flip() +
  ylab('Normalized permutation predictor importance') +
  geom_text(aes(label=scales::percent(imp, accuracy = 1L)), hjust=1.1, col='white', size=6) +
  theme(axis.ticks.x = element_blank(),  axis.ticks.y = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(), text=element_text(size = 25), axis.text.y = element_text(size = 15)) +
  facet_grid(rows = vars(dataset), scales = 'free', space = 'free') 
  
######################################################################################
# outlier drug
######################################################################################
result.outlier <- read.csv('result/result.outlier.csv')
result.outlier$by = factor(result.outlier$by, levels=c('By observations','By drugs'))

# only plot by drugs
ggplot(result.outlier[result.outlier$by=='By drugs',], aes(x=measurement, y=value, fill=drug.remove)) +
  geom_bar(stat="identity", position=position_dodge(width=0.55), width=0.5) +
  facet_grid(rows = vars(model), cols = vars(dataset)) +
  coord_flip(ylim=c(.5, 1)) +
  theme_bw() +
  scale_x_discrete(limits = rev) +
  theme(legend.position="bottom", legend.title=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), text=element_text(size = 23), axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), legend.text = element_text(size=15),
        panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Dark2") 

# only plot by observations
ggplot(result.outlier[result.outlier$by=='By observations',], aes(x=measurement, y=value, fill=drug.remove)) +
  geom_bar(stat="identity", position=position_dodge(width=0.55), width=0.5) +
  facet_grid(rows = vars(model), cols = vars(dataset)) +
  coord_flip(ylim=c(.5, 1)) +
  theme_bw() +
  scale_x_discrete(limits = rev) +
  theme(legend.position="bottom", legend.title=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), text=element_text(size = 23), axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), legend.text = element_text(size=15),
        panel.grid.major = element_blank()) +
  scale_fill_brewer(palette="Dark2") 

######################################################################################
# Univariate vs Multivariate
######################################################################################
result.variable <- read.csv('result/mutivariate.csv')
result.variable[result.variable$model=='Ordinal logistic regression', 'model'] <- 'Ordinal\nLogistic regression'
result.variable[result.variable$model=='Ordinal random forest', 'model'] <- 'Ordinal\nRandom forest'

ggplot(result.variable, aes(x=model, y=acc, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge(width=0.55), width=0.5) +
  coord_flip(ylim=c(.1, 0.85)) +
  facet_grid(rows = vars(dataset), cols = vars(by)) +
  theme_bw() + 
  scale_x_discrete(limits = rev) +
  theme(legend.position="bottom", legend.title=element_blank(), axis.title.x = element_text(size = 15),
        axis.title.y=element_blank(), text=element_text(size = 23), axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13), legend.text = element_text(size=15),
        panel.grid.major = element_blank()) +
  ylab('Accuracy') +
  scale_fill_brewer(palette="Set1") 
