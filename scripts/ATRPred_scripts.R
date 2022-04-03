rm(list=ls())

ra_meta=read.table('ra_final.tsv', sep='\t', header = T)
ra_npx=read.table('ra_npx.tsv', sep='\t', header = T)

library(tidyverse)
ra=left_join(ra_meta,ra_npx, by='StudyID')
ra_resp=ra[!is.na(ra$nice),]
ra_nice=ra_resp[c(10,11,4,7,5,12:363)]

hpc <- caret::preProcess(ra_nice[,-c(1:3,5)], method = c("knnImpute","center", "scale"))
transformed <- predict(hpc, newdata = ra_nice[,-c(1:3,5)])
pca <- prcomp(transformed[,-c(1,354:356)], rank=352)#, center = F,scale. = F)
summary(pca)
p=ggplot(transformed, aes(x=pca$x[,2], y=pca$x[,3],size=DAS)) +
  geom_point(aes(fill=Gender), shape=21) + #21 to 25 have colour and fill
  #scale_color_manual(values=c("black", "white")) +
  geom_hline(yintercept=0) + 
  theme_bw()+
  #scale_fill_discrete(name = "Genders", labels = c("Male", "Female")) +
  scale_fill_manual(values=c("white", "black"), name= "Genders", guide = guide_legend(reverse = TRUE), labels = c("Female", "Male")) + 
  #labs(fill = "Gender") +
  #ggtitle("PCA")+theme(plot.title = element_text(hjust = 0.5))+
  xlab("PC2") + ylab("PC3")

p

## Machine Learning
library(glmnet)
library(caret)
library(pROC)
library(beanplot)
library(RANN)

# Feature Importance
getFI=function(data,feature,outcome,n_run=100, div=0.8, frac=0.8) {
  df=data[c(feature,outcome)]
  run=list()
  for(i in 1:n_run) {
    trainIndex <- createDataPartition(df[[outcome]], p = div, list = FALSE)
    Train <- df[trainIndex,]
    Test  <- df[-trainIndex,]
    x_train=as.matrix(Train[,!names(Train) %in% outcome])
    y_train=Train[,outcome]
    if(div*nrow(data)-sum(y_train) < 8) next
    x_test=as.matrix(Test[,!names(Train) %in% outcome])
    y_test=Test[,outcome]
    if(length(y_test)-sum(y_test) == 0) next
    cvfit = cv.glmnet(x_train, y_train, type.measure = "class", alpha=0.9, nfolds = 10, family="binomial")
    y_test_pred=predict(cvfit, newx = x_test, s = "lambda.min")
    #if(div*nrow(data)-sum(y_test) == 0) next
    roc_test <- roc(as.numeric(y_test),as.numeric(y_test_pred))
    auc_test=pROC::auc(roc_test)
    if(auc_test>=0.5) {
      beta=as.matrix(coef(cvfit, s = "lambda.min"))
      prot=rownames(beta)[beta[,1]!=0]
      run=append(run,prot)
    }
  }
  freq=data.frame(table(as.factor(apply(cbind(run), 1, unlist))))[-1,]
  #return(freq)
  return(freq[order(freq["Freq"], decreasing = T),])
}
set.seed(200)
fi_lasso=getFI(transformed,features_prot, outcome,500,0.8,1)
h=fi_lasso$Freq[1:30]
names(h)=fi_lasso$Var1[1:30]
barplot(h,las=2,ylab="No. of times protein appeared in the 500 models")

mcc <- function (conf_matrix) {
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  mcc_num <- (TP*TN - FP*FN)
  mcc_den <- 
    as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))
  
  mcc_final <- mcc_num/sqrt(mcc_den)
  return(mcc_final)
}

acc <- function (conf_matrix) {
  TP <- conf_matrix$table[1,1]
  TN <- conf_matrix$table[2,2]
  FP <- conf_matrix$table[1,2]
  FN <- conf_matrix$table[2,1]
  
  acc_final <- (TP+TN)/(TN+TP+FN+FP)
  return(acc_final)
}

library(stringr)
calcPerf=function(data,feature,outcome, k=5) {
  new=data[c(feature,outcome)]
  flds <- createFolds(new[[outcome]], k = 5, list = TRUE, returnTrain = FALSE)
  auc_train_5f=0
  auc_test_5f=0
  se_train_5f=0
  sp_train_5f=0
  se_test_5f=0
  sp_test_5f=0
  mcc_train_5f=0
  mcc_test_5f=0
  acc_train_5f=0
  acc_test_5f=0
  
  #Fold1
  Train <- new[-flds$Fold1,]
  Test  <- new[flds$Fold1,]
  x_train=as.matrix(Train[,!names(Train) %in% outcome])
  y_train=Train[,outcome]
  x_test=as.matrix(Test[,!names(Train) %in% outcome])
  y_test=Test[,outcome]
  cvfit = cv.glmnet(x_train, y_train, type.measure = "class", alpha=0.9, nfolds = 10, family="binomial")
  y_train_pred=predict(cvfit, newx = x_train, s = "lambda.min")
  y_test_pred=predict(cvfit, newx = x_test, s = "lambda.min")
  roc_train <- roc(as.numeric(y_train),as.numeric(y_train_pred))
  roc_test <- roc(as.numeric(y_test),as.numeric(y_test_pred))
  delong=roc.test(roc_train, roc_test, method = "delong", paired = F)
  auc_train1=auc(roc_train)
  auc_test1=auc(roc_test)
  th_train=coords(roc_train, "b", ret = "t", best.method = "youden", transpose = T)
  y_train_pred_b <- as.factor(ifelse(y_train_pred > th_train,1,0))
  y_test_pred_b <- as.factor(ifelse(y_test_pred > th_train,1,0))
  cm_train=caret::confusionMatrix(data = y_train_pred_b, reference = as.factor(y_train))
  cm_test=caret::confusionMatrix(data = y_test_pred_b, reference = as.factor(y_test))
  mcc_train1=mcc(cm_train)
  mcc_test1=mcc(cm_test)
  acc_train1=acc(cm_train)
  acc_test1=acc(cm_test)
  se_train1=cm_train$byClass['Sensitivity']
  sp_train1=cm_train$byClass['Specificity']
  se_test1=cm_test$byClass['Sensitivity']
  sp_test1=cm_test$byClass['Specificity']
  auc_train_5f=auc_train_5f + auc_train1
  auc_test_5f=auc_test_5f +auc_test1
  se_train_5f=se_train_5f +se_train1
  sp_train_5f=sp_train_5f+ sp_train1
  se_test_5f=se_test_5f+se_test1
  sp_test_5f=sp_test_5f+sp_test1
  mcc_train_5f=mcc_train_5f + mcc_train1
  mcc_test_5f=mcc_test_5f +mcc_test1
  acc_train_5f=acc_train_5f + acc_train1
  acc_test_5f=acc_test_5f +acc_test1
  
  #Fold2
  Train <- new[-flds$Fold2,]
  Test  <- new[flds$Fold2,]
  x_train=as.matrix(Train[,!names(Train) %in% outcome])
  y_train=Train[,outcome]
  x_test=as.matrix(Test[,!names(Train) %in% outcome])
  y_test=Test[,outcome]
  cvfit = cv.glmnet(x_train, y_train, type.measure = "class", alpha=0.9, nfolds = 10, family="binomial")
  y_train_pred=predict(cvfit, newx = x_train, s = "lambda.min")
  y_test_pred=predict(cvfit, newx = x_test, s = "lambda.min")
  roc_train <- roc(as.numeric(y_train),as.numeric(y_train_pred))
  roc_test <- roc(as.numeric(y_test),as.numeric(y_test_pred))
  delong=roc.test(roc_train, roc_test, method = "delong", paired = F)
  auc_train2=auc(roc_train)
  auc_test2=auc(roc_test)
  th_train=coords(roc_train, "b", ret = "t", best.method = "youden", transpose = T)
  y_train_pred_b <- as.factor(ifelse(y_train_pred > th_train,1,0))
  y_test_pred_b <- as.factor(ifelse(y_test_pred > th_train,1,0))
  cm_train=caret::confusionMatrix(data = y_train_pred_b, reference = as.factor(y_train))
  cm_test=caret::confusionMatrix(data = y_test_pred_b, reference = as.factor(y_test))
  mcc_train2=mcc(cm_train)
  mcc_test2=mcc(cm_test)
  acc_train2=acc(cm_train)
  acc_test2=acc(cm_test)
  se_train2=cm_train$byClass['Sensitivity']
  sp_train2=cm_train$byClass['Specificity']
  se_test2=cm_test$byClass['Sensitivity']
  sp_test2=cm_test$byClass['Specificity']
  auc_train_5f=auc_train_5f + auc_train2
  auc_test_5f=auc_test_5f +auc_test2
  se_train_5f=se_train_5f +se_train2
  sp_train_5f=sp_train_5f+ sp_train2
  se_test_5f=se_test_5f+se_test2
  sp_test_5f=sp_test_5f+sp_test2
  mcc_train_5f=mcc_train_5f + mcc_train2
  mcc_test_5f=mcc_test_5f +mcc_test2
  acc_train_5f=acc_train_5f + acc_train2
  acc_test_5f=acc_test_5f +acc_test2
  
  #Fold3
  Train <- new[-flds$Fold3,]
  Test  <- new[flds$Fold3,]
  x_train=as.matrix(Train[,!names(Train) %in% outcome])
  y_train=Train[,outcome]
  x_test=as.matrix(Test[,!names(Train) %in% outcome])
  y_test=Test[,outcome]
  cvfit = cv.glmnet(x_train, y_train, type.measure = "class", alpha=0.9, nfolds = 10, family="binomial")
  y_train_pred=predict(cvfit, newx = x_train, s = "lambda.min")
  y_test_pred=predict(cvfit, newx = x_test, s = "lambda.min")
  roc_train <- roc(as.numeric(y_train),as.numeric(y_train_pred))
  roc_test <- roc(as.numeric(y_test),as.numeric(y_test_pred))
  delong=roc.test(roc_train, roc_test, method = "delong", paired = F)
  auc_train3=auc(roc_train)
  auc_test3=auc(roc_test)
  th_train=coords(roc_train, "b", ret = "t", best.method = "youden", transpose = T)
  y_train_pred_b <- as.factor(ifelse(y_train_pred > th_train,1,0))
  y_test_pred_b <- as.factor(ifelse(y_test_pred > th_train,1,0))
  cm_train=caret::confusionMatrix(data = y_train_pred_b, reference = as.factor(y_train))
  cm_test=caret::confusionMatrix(data = y_test_pred_b, reference = as.factor(y_test))
  mcc_train3=mcc(cm_train)
  mcc_test3=mcc(cm_test)
  acc_train3=acc(cm_train)
  acc_test3=acc(cm_test)
  se_train3=cm_train$byClass['Sensitivity']
  sp_train3=cm_train$byClass['Specificity']
  se_test3=cm_test$byClass['Sensitivity']
  sp_test3=cm_test$byClass['Specificity']
  auc_train_5f=auc_train_5f + auc_train3
  auc_test_5f=auc_test_5f +auc_test3
  se_train_5f=se_train_5f +se_train3
  sp_train_5f=sp_train_5f+ sp_train3
  se_test_5f=se_test_5f+se_test3
  sp_test_5f=sp_test_5f+sp_test3
  mcc_train_5f=mcc_train_5f + mcc_train3
  mcc_test_5f=mcc_test_5f +mcc_test3
  acc_train_5f=acc_train_5f + acc_train3
  acc_test_5f=acc_test_5f +acc_test3
  
  #Fold4
  Train <- new[-flds$Fold4,]
  Test  <- new[flds$Fold4,]
  x_train=as.matrix(Train[,!names(Train) %in% outcome])
  y_train=Train[,outcome]
  x_test=as.matrix(Test[,!names(Train) %in% outcome])
  y_test=Test[,outcome]
  cvfit = cv.glmnet(x_train, y_train, type.measure = "class", alpha=0.9, nfolds = 10, family="binomial")
  y_train_pred=predict(cvfit, newx = x_train, s = "lambda.min")
  y_test_pred=predict(cvfit, newx = x_test, s = "lambda.min")
  roc_train <- roc(as.numeric(y_train),as.numeric(y_train_pred))
  roc_test <- roc(as.numeric(y_test),as.numeric(y_test_pred))
  delong=roc.test(roc_train, roc_test, method = "delong", paired = F)
  auc_train4=auc(roc_train)
  auc_test4=auc(roc_test)
  th_train=coords(roc_train, "b", ret = "t", best.method = "youden", transpose = T)
  y_train_pred_b <- as.factor(ifelse(y_train_pred > th_train,1,0))
  y_test_pred_b <- as.factor(ifelse(y_test_pred > th_train,1,0))
  cm_train=caret::confusionMatrix(data = y_train_pred_b, reference = as.factor(y_train))
  cm_test=caret::confusionMatrix(data = y_test_pred_b, reference = as.factor(y_test))
  mcc_train4=mcc(cm_train)
  mcc_test4=mcc(cm_test)
  acc_train4=acc(cm_train)
  acc_test4=acc(cm_test)
  se_train4=cm_train$byClass['Sensitivity']
  sp_train4=cm_train$byClass['Specificity']
  se_test4=cm_test$byClass['Sensitivity']
  sp_test4=cm_test$byClass['Specificity']
  auc_train_5f=auc_train_5f + auc_train4
  auc_test_5f=auc_test_5f +auc_test4
  se_train_5f=se_train_5f +se_train4
  sp_train_5f=sp_train_5f+ sp_train4
  se_test_5f=se_test_5f+se_test4
  sp_test_5f=sp_test_5f+sp_test4
  mcc_train_5f=mcc_train_5f + mcc_train4
  mcc_test_5f=mcc_test_5f +mcc_test4
  acc_train_5f=acc_train_5f + acc_train4
  acc_test_5f=acc_test_5f +acc_test4
  
  #Fold5
  Train <- new[-flds$Fold5,]
  Test  <- new[flds$Fold5,]
  x_train=as.matrix(Train[,!names(Train) %in% outcome])
  y_train=Train[,outcome]
  x_test=as.matrix(Test[,!names(Train) %in% outcome])
  y_test=Test[,outcome]
  cvfit = cv.glmnet(x_train, y_train, type.measure = "class", alpha=0.9, nfolds = 10, family="binomial")
  y_train_pred=predict(cvfit, newx = x_train, s = "lambda.min")
  y_test_pred=predict(cvfit, newx = x_test, s = "lambda.min")
  roc_train <- roc(as.numeric(y_train),as.numeric(y_train_pred))
  roc_test <- roc(as.numeric(y_test),as.numeric(y_test_pred))
  delong=roc.test(roc_train, roc_test, method = "delong", paired = F)
  auc_train5=auc(roc_train)
  auc_test5=auc(roc_test)
  th_train=coords(roc_train, "b", ret = "t", best.method = "youden", transpose = T)
  y_train_pred_b <- as.factor(ifelse(y_train_pred > th_train,1,0))
  y_test_pred_b <- as.factor(ifelse(y_test_pred > th_train,1,0))
  cm_train=caret::confusionMatrix(data = y_train_pred_b, reference = as.factor(y_train))
  cm_test=caret::confusionMatrix(data = y_test_pred_b, reference = as.factor(y_test))
  mcc_train5=mcc(cm_train)
  mcc_test5=mcc(cm_test)
  acc_train5=acc(cm_train)
  acc_test5=acc(cm_test)
  se_train5=cm_train$byClass['Sensitivity']
  sp_train5=cm_train$byClass['Specificity']
  se_test5=cm_test$byClass['Sensitivity']
  sp_test5=cm_test$byClass['Specificity']
  auc_train_5f=auc_train_5f + auc_train5
  auc_test_5f=auc_test_5f +auc_test5
  se_train_5f=se_train_5f +se_train5
  sp_train_5f=sp_train_5f+ sp_train5
  se_test_5f=se_test_5f+se_test5
  sp_test_5f=sp_test_5f+sp_test5
  mcc_train_5f=mcc_train_5f + mcc_train5
  mcc_test_5f=mcc_test_5f +mcc_test5
  acc_train_5f=acc_train_5f + acc_train5
  acc_test_5f=acc_test_5f +acc_test5
  
  #Final
  auc_train_5f=auc_train_5f/5
  auc_test_5f=auc_test_5f/5
  se_train_5f=se_train_5f/5
  sp_train_5f=sp_train_5f/5
  se_test_5f=se_test_5f/5
  sp_test_5f=sp_test_5f/5
  mcc_train_5f=mcc_train_5f/5
  mcc_test_5f=mcc_test_5f/5
  acc_train_5f=acc_train_5f/5
  acc_test_5f=acc_test_5f/5
  
  auc_train_5f=mean(c(auc_train1,auc_train2,auc_train3,auc_train4,auc_train5), na.rm=T)
  auc_test_5f=mean(c(auc_test1,auc_test2,auc_test3,auc_test4,auc_test5), na.rm=T)
  se_train_5f=mean(c(se_train1,se_train2,se_train3,se_train4,se_train5), na.rm=T)
  sp_train_5f=mean(c(sp_train1,sp_train2,sp_train3,sp_train4,sp_train5), na.rm=T)
  se_test_5f=mean(c(se_test1,se_test2,se_test3,se_test4,se_test5), na.rm=T)
  sp_test_5f=mean(c(sp_test1,sp_test2,sp_test3,sp_test4,sp_test5), na.rm=T)
  mcc_train_5f=mean(c(mcc_train1,mcc_train2,mcc_train3,mcc_train4,mcc_train5), na.rm=T)
  mcc_test_5f=mean(c(mcc_test1,mcc_test2,mcc_test3,mcc_test4,mcc_test5), na.rm=T)
  acc_train_5f=mean(c(acc_train1,acc_train2,acc_train3,acc_train4,acc_train5), na.rm=T)
  acc_test_5f=mean(c(acc_test1,acc_test2,acc_test3,acc_test4,acc_test5), na.rm=T)
  
  c=str_c(names(new),collapse = "; ")
  
  return(data.frame(c,auc_train1,auc_test1,acc_train1,acc_test1,se_train1,se_test1,sp_train1,sp_test1,mcc_train1,mcc_test1,auc_train2,auc_test2,acc_train2,acc_test2,se_train2,se_test2,sp_train2,sp_test2,mcc_train2,mcc_test2,auc_train3,auc_test3,acc_train3,acc_test3,se_train3,se_test3,sp_train3,sp_test3,mcc_train3,mcc_test3,auc_train4,auc_test4,acc_train4,acc_test4,se_train4,se_test4,sp_train1,sp_test4,mcc_train4,mcc_test4,auc_train5,auc_test5,acc_train5,acc_test5,se_train5,se_test5,sp_train5,sp_test5,mcc_train5,mcc_test5,auc_train_5f,auc_test_5f,acc_train_5f,acc_test_5f,se_train_5f,se_test_5f,sp_train_5f,sp_test_5f,mcc_train_5f,mcc_test_5f))
}
#Running Signature loop
c=data.frame()
n=32 #should be greater than 2
set.seed(200)
for(i in 2:n){
  a=calcPerf(transformed,feat[1:i],outcome)
  c=rbind(c,a)
}

xdata=0:(nrow(c)-1)

# plot the first curve by calling plot() function
# First curve is plotted
plot(xdata, c$auc_train_5f, type="l", col="black", pch="o", lty=1 , ylim=c(0,1), xlab="No. of Proteins", ylab="AUC")

# Add second curve to the same plot by calling points() and lines()
# Use symbol '*' for points.
# points(xdata, c$auc_test_5f, col="black", pch="*")
lines(xdata, c$auc_test_5f, col="black",lty=2)

# Add Third curve to the same plot by calling points() and lines()
# Use symbol '+' for points.
# points(xdata, c$sp_test_5f, col="dark red",pch="+")
# lines(xdata, c$sp_test_5f, col="dark red", lty=3)

legend(20, 0.2,col=c("black", "black"), legend=c("Training set", "Test set"), lty=1:2, cex=0.8)

png("es.png")
p=coef(cvfit, s = "lambda.min")
plot(p[-1], xlab="Decreasing Importance", ylab="Effect Size", xlim = c(0,21))
abline(h=0)
text(p, row.names(p), cex=0.6, pos=4, col="black", offset=-0.8) 
dev.off()

library(beeswarm)
beeswarm(fin ~ nice, data = transformed,
         method = 'swarm',
         pch = 16, #pwcol = eular+1,
         xlab = '', ylab = 'Score',
         labels = c('NR', 'R'))

boxplot(fin ~ nice, data = transformed, add = T, col="#00000000", names=c("NR","R"))#, names = c("",""), col="#0000ff22")
abline(h=th)


## Cluster Validation
var_explained_df <- data.frame(PC=1:89, var_explained=100*(pca$sdev)^2/sum((pca$sdev)^2))

var_explained_df %>%
  ggplot(aes(x=PC,y=var_explained))+
  geom_col()+
  ylab("Percentage of explained variances")+
  xlab("Principal Components")+
  geom_line(aes(y=1),linetype="dotted") +
  scale_x_continuous(limits = c(0, 31)) +
  theme_minimal()

factoextra::fviz_eig(pca, ncp=352)

par(cex=.5)
factoextra::fviz_eig(pca, ncp=50, ggtheme = theme(axis.text = element_text(size = 7)))

screeplot(pca, type = "l", npcs = 50, main = "Screeplot of the first 50 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

cumpro <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
plot(cumpro[0:50], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
pc=4
abline(v = pc, col="blue", lty=5)
abline(h = cumpro[pc], col="blue", lty=5)
legend("topleft", legend=paste0(c("Cut-off @ PC"),pc),
       col=c("blue"), lty=5, cex=0.6)

X <- as.matrix(transformed[,-c(1,354:356)])
res <- sinkr::pca_loocv(X[,1:20])
res2 <- lapply(res, colSums)

COL <- 2:4
LTY <- 1:3
op <- par(mar=c(4,4,2,1), tcl=-0.25, mgp=c(2.5,0.5,0))
for(i in seq(res)){
  if(i==1) {
    plot(res2[[i]], t="n", ylim=range(unlist(res2)), 
         xlab="Number of Principal Components", ylab="PRESS")
    grid()
  } 
  lines(res2[[i]], t="b", bg=c(NaN,COL[i])[(res2[[i]]==min(res2[[i]])) + 1],
        col=COL[i], lty=LTY[i], pch=21)
}
legend("topright", legend=c("naive", "approximate", "pseudoinverse"),
       col=COL, lty=LTY, pch=21, bty="n")
par(op)

