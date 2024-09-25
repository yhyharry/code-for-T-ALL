rm(list = ls()) 
options(stringsAsFactors = F)

#### load data ####
library(tidyverse)
library(broom)
library(glmnet)
k<-read.csv("T-ALL 24sample fpkm.csv",row.names = 1)

library(readxl)
gene<-read_xlsx("genelist.xlsx")
k<-k[rownames(k)%in%gene$gene,]
boxplot(k)
k1 <- log2(k)
boxplot(k1)
k<-as.data.frame(t(k))

fen<-read_xlsx("grouplist.xlsx")
fen <- as.data.frame(fen)
rownames(fen) <- fen[,1]
group_list <- fen$`ETP-like status`
lasso = factor(group_list)
fen <- data.frame(row.names = rownames(k), lasso) 

k<-k[rownames(fen),]
k<-cbind(fen,k)
write.csv(k,file = 'group and gene expression.csv')

states<-as.matrix(k)
x<-states[,-1]
y<-states[,1]
y <- as.numeric(y)
x <- glmnet::makeX(train = k[, !names(k) == "mdv"])
x <- x[,-c(1,2)]

#### lasso regression ####
set.seed(2024)
cvfit=cv.glmnet(data.matrix(x),y,type.measure = "mse",nfolds = 5,alpha=1)
plot(cvfit)
cvfit$lambda.min
c(cvfit$lambda.min, cvfit$lambda.1se)
lasso<-glmnet(x,y,family="binomial",alpha=1,nlambda = 100)
coef(lasso, s=c(0.08362403,0.16802016))
score <- coef(lasso, s=c(0.08362403,0.16802016))

print(lasso) 
plot(lasso,label=T) 
plot(lasso,xvar="lambda",label=T) 
lasso.coef<-predict(lasso,s=0.4,type="coefficients") 
plot(lasso,xvar="dev",label=T) 
lasso.y<-predict(lasso,newx=x,type="response",s=0.4) 
plot(lasso.y,y,xlab="Predicted",ylab="Actual",main="Lasso Regression")

coef <- coef(lasso, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
geneCoef     
write.csv(geneCoef,file = 'lasso gene.csv')

lassoGene <- as.matrix(lassoGene)
lassoGene <- lassoGene[-1,]
actCoef <- as.matrix(actCoef)
actCoef <- actCoef[-1,]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
geneCoef


#### validation ####
library(ggpubr) 
library(ROCR)   
library(glmnet)
library(caret)
myFun <-  function(x){crossprod(as.numeric(x),actCoef)}

###### validation1 Target
test1 <- read_xlsx('Test1_Target fpkm.xlsx')
test1 <- as.data.frame(test1)
rownames(test1) <- test1[,1]
test1 <- test1[,-1]
test1 <- log2(test1+1)
test1 <- test1[c("JUP","MGLL","PRR5L","SH3RF1","SLC9A7","SLC45A3","ERMP1"),]
test1<-as.data.frame(t(test1))
grouptest1<-read_xlsx("Test1_Target grouplist.xlsx")
grouptest1 <- as.data.frame(grouptest1)
rownames(grouptest1) <- grouptest1[,1]
group_list <- grouptest1$`ETP status`
lasso = factor(group_list)
group_list <- data.frame(row.names = rownames(test1), lasso) 

diagnosticScoretest1 <-  apply(test1,1,myFun)
outCol <-  c('lasso', lassoGene)
diagnostictest1 <-  as.vector(ifelse(diagnosticScoretest1 > median(diagnosticScoretest1), "ETP", "nonETP"))
dat_test1 <- cbind(test1, diagnosticScoretest1=as.vector(diagnosticScoretest1), diagnostictest1) 
dat_test1 <- cbind(dat_test1,group_list)

predtest1 <- prediction(dat_test1$diagnosticScoretest1, dat_test1$lasso)
ROCtest1 <- performance(predtest1,"tpr","fpr")
AUCtest1 <- performance(predtest1,"auc")
AUCtest1@y.values

plot(ROCtest1,colorize=FALSE, col="red", print.auc =TRUE,xlab="1-Specificity", ylab="Sensitivity") +  
  lines(c(0,1),c(0,1),col = "black", lty = 10 ) + 
  legend(0.5,0.2,'AUC=0.97')


###### validation2 NODE
test2 <- read.csv('Test2_NODE TPM.csv')
rownames(test2) <- test2[,1]
test2 <- test2[,-1]
test2 <- log2(test2+1)
test2 <- test2[c("JUP","MGLL","PRR5L","SH3RF1","SLC9A7","SLC45A3","ERMP1"),]
test2<-as.data.frame(t(test2))
grouptest2<-read_xlsx("Test2_NODE grouplist.xlsx")
grouptest2 <- as.data.frame(grouptest2)
rownames(grouptest2) <- grouptest2[,1]
group_list <- grouptest2$`ETP status`
lasso = factor(group_list)
group_list <- data.frame(row.names = rownames(test2), lasso) 

diagnosticScoretest2 <-  apply(test2,1,myFun)
outCol <-  c('lasso', lassoGene)
diagnostictest2 <-  as.vector(ifelse(diagnosticScoretest2 > median(diagnosticScoretest2), "ETP", "nonETP"))
dat_test2 <- cbind(test2, diagnosticScoretest2=as.vector(diagnosticScoretest2), diagnostictest2) 
dat_test2 <- cbind(dat_test2,group_list)

predtest2 <- prediction(dat_test2$diagnosticScoretest2, dat_test2$lasso)
ROCtest2 <- performance(predtest2,"tpr","fpr")
AUCtest2 <- performance(predtest2,"auc")   
AUCtest2@y.values 

plot(ROCtest2,colorize=FALSE, col="red", print.auc =TRUE,xlab="1-Specificity", ylab="Sensitivity") + 
  lines(c(0,1),c(0,1),col = "black", lty = 10 ) +  
  legend(0.5,0.2,'AUC=0.8125') 

###### validation3 GSE146901
test3 <- read_xlsx('Test3_GSE146901 fpkm.xlsx')
test3 <- as.data.frame(test3)
rownames(test3) <- test3[,1]
test3 <- test3[,-1]
test3 <- log2(test3+1)
test3 <- test3[c("JUP","MGLL","PRR5L","SH3RF1","SLC9A7","SLC45A3","ERMP1"),]
test3<-as.data.frame(t(test3))
grouptest3<-read_xlsx("Test3_GSE146901 group list.xlsx")
grouptest3 <- as.data.frame(grouptest3)
rownames(grouptest3) <- grouptest3[,1]
group_list <- grouptest3$`ETP status`
lasso = factor(group_list)
group_list <- data.frame(row.names = rownames(test3), lasso) 

diagnosticScoretest3 <-  apply(test3,1,myFun)
outCol <-  c('lasso', lassoGene)
diagnostictest3 <-  as.vector(ifelse(diagnosticScoretest3 > median(diagnosticScoretest3), "ETP", "nonETP"))
dat_test3 <- cbind(test3, diagnosticScoretest3=as.vector(diagnosticScoretest3), diagnostictest3) 
dat_test3 <- cbind(dat_test3,group_list)

predtest3 <- prediction(dat_test3$diagnosticScoretest3, dat_test3$lasso)
ROCtest3 <- performance(predtest3,"tpr","fpr")
AUCtest3 <- performance(predtest3,"auc")  
AUCtest3@y.values 

plot(ROCtest3,colorize=FALSE, col="red", print.auc =TRUE,xlab="1-Specificity", ylab="Sensitivity") +  
  lines(c(0,1),c(0,1),col = "black", lty = 10 ) +  
  legend(0.5,0.4,'AUC=0.925') 

###### validation4 GSE141140
test4 <- read_xlsx('Test4_GSE141140 fpkm.xlsx')
test4 <- as.data.frame(test4)
rownames(test4) <- test4[,1]
test4 <- test4[,-1]
test4 <- log2(test4+1)
test4 <- test4[c("JUP","MGLL","PRR5L","SH3RF1","SLC9A7","SLC45A3","ERMP1"),]
test4<-as.data.frame(t(test4))
grouptest4<-read_xlsx("Test4_GSE141140 group list.xlsx")
grouptest4 <- as.data.frame(grouptest4)
rownames(grouptest4) <- grouptest4[,1]
group_list <- grouptest4$`ETP status`
lasso = factor(group_list)
group_list <- data.frame(row.names = rownames(test4), lasso) 

diagnosticScoretest4 <-  apply(test4,1,myFun)
outCol <-  c('lasso', lassoGene)
diagnostictest4 <-  as.vector(ifelse(diagnosticScoretest4 > median(diagnosticScoretest4), "ETP", "nonETP"))
dat_test4 <- cbind(test4, diagnosticScoretest4=as.vector(diagnosticScoretest4), diagnostictest4) 
dat_test4 <- cbind(dat_test4,group_list)

predtest4 <- prediction(dat_test4$diagnosticScoretest4, dat_test4$lasso)
ROCtest4 <- performance(predtest4,"tpr","fpr")
AUCtest4 <- performance(predtest4,"auc")  
AUCtest4@y.values

plot(ROCtest4,colorize=FALSE, col="red", print.auc =TRUE,xlab="1-Specificity", ylab="Sensitivity") +  
  lines(c(0,1),c(0,1),col = "black", lty = 10 ) +  
  legend(0.5,0.2,'AUC=0.84') 

