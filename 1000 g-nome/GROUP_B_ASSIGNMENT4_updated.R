#1000 Genome Project
#Group B
#Members: ADEDAMOLA OGUNBONA
#         BRITTANY MAC INTYRE
#         HOC TRAN
#         ZHIKANG CHEN
          
# Load Necessary packages 
library(vcfR)
library(genetics)
library(glmnet)
library(ggplot2)
library(GGally)
library(pROC)
library(MASS)
library(lars)
library(caret)
library(class)
library(gtools)
library(ROCR)
library(randomForest)

### Load the vcf files for the BRCA1 gene, the vcf files have been filtered
BRCAEAS <- read.vcfR("filt_ba_freq_BRCA1EAS.vcf.gz")
BRCAEASgt <- extract.gt(BRCAEAS, as.numeric = T)
BRCAEASgt2 <- t(BRCAEASgt)
BRCAEASgt2<-cbind(BRCAEASgt2,0)
colnames(BRCAEASgt2) <- c(colnames(BRCAEASgt2)[1:661],"Population")

BRCAEUR <- read.vcfR("filt_ba_freq_BRCA1EUR.vcf.gz")
BRCAEURgt <- extract.gt(BRCAEUR, as.numeric=T)
BRCAEURgt2<-t(BRCAEURgt)
BRCAEURgt2<-cbind(BRCAEURgt2,1)
colnames(BRCAEURgt2) <- c(colnames(BRCAEURgt2)[1:737],"Population")

### Load the vcf file for the IGF1 gene, this file has also been filtered
IGF1EAS <- read.vcfR("filt_ba_freq_EAS_IGF1.vcf.gz")
IGF1EASgt <- extract.gt(IGF1EAS, as.numeric=TRUE)
IGF1EASgt2<-t(IGF1EASgt)
IGF1EASgt2<-cbind(IGF1EASgt2,0)
colnames(IGF1EASgt2) <- c(colnames(IGF1EASgt2)[1:549],"Population")

IGF1EUR <- read.vcfR("filt_ba_freq_EUR_IGF1.vcf.gz")
IGF1EURgt <- extract.gt(IGF1EUR, as.numeric=TRUE)
IGF1EURgt2<-t(IGF1EURgt)
IGF1EURgt2<-cbind(IGF1EURgt2,1)
colnames(IGF1EURgt2) <- c(colnames(IGF1EURgt2)[1:528],"Population")


### Load the vcf file for the HMGA2 gene, this file has also been filtered
HMGA2EAS <- read.vcfR("filt_ba_freq_EAS_HMGA2.vcf.gz")
HMGA2EASgt <- extract.gt(HMGA2EAS, as.numeric=TRUE)
HMGA2EASgt2<-t(HMGA2EASgt)
HMGA2EASgt2<-cbind(HMGA2EASgt2,0)
colnames(HMGA2EASgt2) <- c(colnames(HMGA2EASgt2)[1:999],"Population")

HMGA2EUR <- read.vcfR("filt_ba_freq_EUR_HMGA2.vcf.gz")
HMGA2EUR@fix[726,3] <- "rs377603591.1" #Duplicated variant ID, modifying one of them
HMGA2EURgt <- extract.gt(HMGA2EUR, as.numeric=TRUE)
HMGA2EURgt2<-t(HMGA2EURgt)
HMGA2EURgt2<-cbind(HMGA2EURgt2,1)
colnames(HMGA2EURgt2) <- c(colnames(HMGA2EURgt2)[1:1007],"Population")


# RANDOM FOREST   ############################################################################################
### Combining both EUR and EAS matrices for the gene and convert all NAs to 0
BRCA1gt2 <- smartbind(BRCAEASgt2, BRCAEURgt2, fill=0)
### Delete the row with the same values
BRCA1gt2 <- Filter(var, BRCA1gt2)
which(colnames(BRCA1gt2) == "Population")
Population <- BRCA1gt2[,503]
BRCA1gt2 <- BRCA1gt2[,-503]
BRCA1gt2 <- cbind(BRCA1gt2, Population)
class(BRCA1gt2)

### Repeat the same for IGF1
### Combining both EUR and EAS matrices for the gene and convert all NAs to 0
IGF1gt2 <- smartbind(IGF1EASgt2, IGF1EURgt2, fill = 0)
IGF1gt2 <- Filter(var, IGF1gt2)
which(colnames(IGF1gt2) == "Population")
Population <- IGF1gt2[,394]
IGF1gt2 <- IGF1gt2[,-394]
IGF1gt2 <- cbind(IGF1gt2, Population)

### Repeat for HMga2
HMGA2gt2 <- smartbind(HMGA2EASgt2, HMGA2EURgt2, fill=0)
HMGA2gt2 <- Filter(var, HMGA2gt2)
which(colnames(HMGA2gt2) == "Population")
Population <- HMGA2gt2[,730]
HMGA2gt2 <- HMGA2gt2[,-730]
HMGA2gt2 <- cbind(HMGA2gt2, Population)

#Creating a function that will create a classifier and generate ROC curves based on predictive performance of the classifiers when trying to classify into either EUR or EAS populations
Classifier <- function (data){
  set.seed(44) #Ensure results are repeatable
  train.index <- sample(1:nrow(data), round(0.70*nrow(data),0)) #70% of the data will be training, the rest are for testing
  training <- data[train.index,]
  testing <- data [-train.index,] 
  x = which(colnames(data) == "Population")
  #Creating the classifier, x has been set to include all columns except the last one (population) as predictor variables, mtry set to the default value. Importance set to true to later examine variable importance
  rf_pop_classifier <- randomForest(x = training[, 1:(length(data)-1)], y = as.factor(training$Population), ntree = 500, importance = TRUE) 
  print(rf_pop_classifier) #Check OOB Error Rate, class error rate, and the confusion matrix
  plot(rf_pop_classifier) 
  
  varImpPlot(rf_pop_classifier) #Check variable importance
  
  #cross.validate <- rfcv(trainx = training[, 1:(length(data)-1)], trainy = as.factor(training$Population), cv.fold = 10)
  #cross.validate
  #with(cross.validate, plot(n.var, error.cv, log="x", type="o", lwd=2))
  
  #Generating ROC Curves
  classes <- levels(as.factor(training$Population)) #Generating classes to be used in ROC Curve
  colours <- c("RED","BLACK")#Colours to be used in the ROC Curves
  
  
  predict_for_rf_pop_classifier <- predict(rf_pop_classifier, testing[, 1:(length(data)-1)], type = "prob")
  predict_for_rf_pop_classifier1 <- predict(rf_pop_classifier, testing[, 1:(length(data)-1)])
  
  print(table(observed=testing[,length(data)],predicted=predict_for_rf_pop_classifier1))
  
  
  
  for (i in 1:2){ #For loop set to i in 1:2 as there are two classes present, EAS and EUR
    true_values <- ifelse(testing$Population==classes[i],1,0) #Using an ifelse statement to correctly attribute samples to the right class
    pred <- prediction(predict_for_rf_pop_classifier[,i],true_values) #Calculating probability values for samples in each class using the inputted predict_for_rf_pop_classifier taken in by the function
    perf <- performance(pred, "tpr", "fpr") #Using pred as the prediction object, "fpr" (false positive rate) is the x axis for the ROC Curve and "tpr" is the y axis 
    if (i==1) #Plotting for EUR
    {
      plot(perf,main="ROC Curve",col=colours[i]) #Using red as colour choice
    }
    else #Plotting for EAS
    {
      plot(perf,main="ROC Curve",col=colours[i],add=TRUE)  #Using black as colour choice
    }
    
    auc.perf <- performance(pred, measure = "auc") #Calculates the area under the curve value for each class
    print(auc.perf@y.values) #Prints the AUC value
  }
}

Classifier(BRCA1gt2) #OOB Error Rate 9.91%,  AUC 0.967 

Classifier(IGF1gt2) #OOB Error Rate 15.41%,  AUC 0.9029

Classifier(HMGA2gt2) #OOB estimate of  error rate: 7.7%   AUC 0.982



# ELASTIC NET    #############################################################################################
### Combining both EUR and EAS matrixes for the gene
BRCA1gt2 <- smartbind(BRCAEASgt2, BRCAEURgt2, fill=0) #Repeat the smartbind to include the data we need for analysis, without this, the transposed gene GT gets filtered incorrectly and continuously in every analysis
rownames(BRCA1gt2)<-c(1:908)

### Delete the row with all values same
BRCA1gt2 <- Filter(var, BRCA1gt2)
which(colnames(BRCA1gt2) == "Population")

IGF1gt2 <- smartbind(IGF1EASgt2, IGF1EURgt2, fill=0)
rownames(IGF1gt2)<-c(1:908)
IGF1gt2 <- Filter(var, IGF1gt2)
which(colnames(IGF1gt2) == "Population")

HMGA2gt2 <- smartbind(HMGA2EASgt2, HMGA2EURgt2, fill=0)
rownames(HMGA2gt2)<-c(1:908)
HMGA2gt2 <- Filter(var, HMGA2gt2)
which(colnames(HMGA2gt2) == "Population")

### Train and test sets for the BRCA1 gene
set.seed(44)
train.index.brca1 <- sample(1:nrow(BRCA1gt2), round(0.70*nrow(BRCA1gt2),0))
brca1train<-BRCA1gt2[train.index.brca1,]
brca1test<-BRCA1gt2[-train.index.brca1,]
which(names(BRCA1gt2)=="Population")
brca1trainx<-data.matrix(brca1train[,-(which(names(BRCA1gt2)=="Population"))])
brca1trainx<-data.matrix(brca1train[,-503])
brca1trainy<-as.factor(brca1train[,503])
brca1testx<-data.matrix(brca1test[,-503])
brca1testy<-as.factor(brca1test[,503])

### Repeat for the IGF1 gene
set.seed(44)
train.index.igf1 <- sample(1:nrow(IGF1gt2), round(0.70*nrow(IGF1gt2),0))
igf1train<-IGF1gt2[train.index.igf1,]
igf1test<-IGF1gt2[-train.index.igf1,]
igf1trainx<-data.matrix(igf1train[,-394])
igf1trainy<-as.factor(igf1train[,394])
igf1testx<-data.matrix(igf1test[,-394])
igf1testy<-as.factor(igf1test[,394])

### And the HMGA2 gene
set.seed(44)
train.index.hmga2 <- sample(1:nrow(HMGA2gt2), round(0.70*nrow(HMGA2gt2),0))
hmga2train<-HMGA2gt2[train.index.hmga2,]
hmga2test<-HMGA2gt2[-train.index.hmga2,]
hmga2trainx<-data.matrix(hmga2train[,-730])
hmga2trainy<-as.factor(hmga2train[,730])
hmga2testx<-data.matrix(hmga2test[,-730])
hmga2testy<-as.factor(hmga2test[,730])

### Perform 10-fold CV Elastic Net and tune lambda
set.seed(44)
enet <- trainControl(method = "cv", number = 10)
tunegrid=expand.grid(.alpha = seq(0.1,0.9,0.01),.lambda = seq(0,0.5,0.025))

### Perform elastic net on the BRCA1 gene
set.seed(44)
system.time(brca1_net <- train(brca1trainx,brca1trainy,method = "glmnet", trControl = enet,tuneGrid = tunegrid))

### Perform elastic net on the IGF1 gene
set.seed(44)
system.time(igf1_net <- train(igf1trainx,igf1trainy,method = "glmnet", trControl = enet,tuneGrid = tunegrid))

### And on the HMGA2 gene
set.seed(44)
system.time(hmga2_net <- train(hmga2trainx,hmga2trainy,method = "glmnet", trControl = enet,tuneGrid = tunegrid))

### Best model
set.seed(44)
w1 <- which.max(brca1_net$results[,3])
cvbrca1 <- cv.glmnet(brca1trainx,brca1trainy, nfolds = 10, family="binomial", alpha=brca1_net$results[w1,1], type.measure = "auc")

set.seed(44)
w2 <- which.max(igf1_net$results[,3])
cvigf1 <- cv.glmnet(igf1trainx,igf1trainy, nfolds = 10, family="binomial", alpha=igf1_net$results[w2,1], type.measure = "auc")

set.seed(44)
w3 <- which.max(hmga2_net$results[,3])
cvhmga2 <- cv.glmnet(hmga2trainx,hmga2trainy, nfolds = 10, family="binomial", alpha=hmga2_net$results[w3,1], type.measure = "auc")

### Prediction by elastic net for each gene
prds.train.brca1 <- predict(cvbrca1,newx = brca1trainx, type = "response", s=cvbrca1$lambda.min)[,1]
prds.test.brca1 <- predict(cvbrca1,newx = brca1testx, type = "response", s=cvbrca1$lambda.min)[,1]

prds.train.igf1 <- predict(cvigf1,newx = igf1trainx, type = "response", s=cvigf1$lambda.min)[,1]
prds.test.igf1 <- predict(cvigf1,newx = igf1testx, type = "response", s=cvigf1$lambda.min)[,1]

prds.train.hmga2 <- predict(cvhmga2,newx = hmga2trainx, type = "response", s=cvhmga2$lambda.min)[,1]
prds.test.hmga2 <- predict(cvhmga2,newx = hmga2testx, type = "response", s=cvhmga2$lambda.min)[,1]

### Threshold for elastic net in each gene
auc.train.brca1 <- roc(brca1trainy,prds.train.brca1)
auc.test.brca1 <- roc(brca1testy,prds.test.brca1)
snsp.train.brca1 <- cbind(auc.train.brca1$sensitivities,auc.train.brca1$specificities)
snsp.test.brca1 <- cbind(auc.test.brca1$sensitivities,auc.test.brca1$specificities)
indx.brca1 <- which.max(apply(snsp.train.brca1,1,min))
cutoff.brca1 <- auc.train.brca1$thresholds[indx.brca1]
cutoff.brca1

auc.train.igf1 <- roc(igf1trainy,prds.train.igf1)
auc.test.igf1 <- roc(igf1testy,prds.test.igf1)
snsp.train.igf1 <- cbind(auc.train.igf1$sensitivities,auc.train.igf1$specificities)
snsp.test.igf1 <- cbind(auc.test.igf1$sensitivities,auc.test.igf1$specificities)
indx.igf1 <- which.max(apply(snsp.train.igf1,1,min))
cutoff.igf1 <- auc.train.igf1$thresholds[indx.igf1]
cutoff.igf1

auc.train.hmga2 <- roc(hmga2trainy,prds.train.hmga2)
auc.test.hmga2 <- roc(hmga2testy,prds.test.hmga2)
snsp.train.hmga2 <- cbind(auc.train.hmga2$sensitivities,auc.train.hmga2$specificities)
snsp.test.hmga2 <- cbind(auc.test.hmga2$sensitivities,auc.test.hmga2$specificities)
indx.hmga2 <- which.max(apply(snsp.train.hmga2,1,min))
cutoff.hmga2 <- auc.train.hmga2$thresholds[indx.hmga2]
cutoff.hmga2

### Sensitivity and specificity, and AUC
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

sn.sp(table(brca1testy, as.numeric(prds.test.brca1>cutoff.brca1)))
auc.test.brca1<-roc(igf1testy,prds.test.brca1)
auc.test.brca1

sn.sp(table(igf1testy, as.numeric(prds.test.igf1>cutoff.igf1)))
auc.test.igf1<-roc(igf1testy,prds.test.igf1)
auc.test.igf1

sn.sp(table(hmga2testy, as.numeric(prds.test.hmga2>cutoff.hmga2)))
auc.test.hmga2<-roc(hmga2testy,prds.test.hmga2)
auc.test.hmga2

### ROC Curve for all 3 classification
snsp.test.brca1<-cbind(auc.test.brca1$sensitivities,auc.test.brca1$specificities)[,-3][,-3]
snsp.test.igf1<-cbind(auc.test.igf1$sensitivities,auc.test.igf1$specificities)[,-3][,-3]
snsp.test.hmga2<-cbind(auc.test.hmga2$sensitivities,auc.test.hmga2$specificities)[,-3][,-3]

colnames(snsp.test.brca1) <- c("sn","sp")
colnames(snsp.test.igf1) <- c("sn","sp")
colnames(snsp.test.hmga2) <-c("sn","sp")

snsp.test.brca1<-as.data.frame(snsp.test.brca1)
snsp.test.igf1<-as.data.frame(snsp.test.igf1)
snsp.test.hmga2<-as.data.frame(snsp.test.hmga2)

snsp.test.brca1$model <- factor("BRCA1",c("BRCA1","IGF1","HMGA2"))
snsp.test.igf1$model <- factor("IGF1",c("BRCA1","IGF1","HMGA2"))
snsp.test.hmga2$model <- factor("HMGA2",c("BRCA1","IGF1","HMGA2"))

roc_dat <- rbind(snsp.test.brca1,snsp.test.igf1,snsp.test.hmga2)


ggplot(data=roc_dat, aes(x=1-sp, y=sn, color=model)) +
  ggtitle("Elastic Net ROC Curves")+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  geom_point()


# RANDOM FOREST FOR BETTER COMPARISON OF ROC CURVE VIA GGPLOT ################################################
### Combining both EUR and EAS matrices for the gene and convert all NAs to 0
BRCA1gt2 <- smartbind(BRCAEASgt2, BRCAEURgt2, fill=0) #Repeat the smartbind to include the data we need for analysis, without this, the transposed gene GT gets filtered incorrectly and continuously in every analysis
### Delete the row with the same values
BRCA1gt2 <- Filter(var, BRCA1gt2)
which(colnames(BRCA1gt2) == "Population")
Population <- BRCA1gt2[,503]
BRCA1gt2 <- BRCA1gt2[,-503]
BRCA1gt2 <- cbind(BRCA1gt2, Population)
class(BRCA1gt2)

### Combining both EUR and EAS matrices for the gene and convert all NAs to 0
IGF1gt2 <- smartbind(IGF1EASgt2, IGF1EURgt2, fill = 0)
### Repeat the same for IGF1
IGF1gt2 <- Filter(var, IGF1gt2)
which(colnames(IGF1gt2) == "Population")
Population <- IGF1gt2[,394]
IGF1gt2 <- IGF1gt2[,-394]
IGF1gt2 <- cbind(IGF1gt2, Population)

HMGA2gt2 <- smartbind(HMGA2EASgt2, HMGA2EURgt2, fill=0)
### Repeat for HMga2
HMGA2gt2 <- Filter(var, HMGA2gt2)
which(colnames(HMGA2gt2) == "Population")
Population <- HMGA2gt2[,730]
HMGA2gt2 <- HMGA2gt2[,-730]
HMGA2gt2 <- cbind(HMGA2gt2, Population)

set.seed(44) #Ensure results are repeatable
train.index <- sample(1:nrow(BRCA1gt2), round(0.70*nrow(BRCA1gt2),0)) #70% of the data will be training, the rest are for testing
training <- BRCA1gt2[train.index,]
testing <- BRCA1gt2 [-train.index,] 
x = which(colnames(BRCA1gt2) == "Population")
#Creating the classifier, x has been set to include all columns except the last one (population) as predictor variables, mtry set to the default value. Importance set to true to later examine variable importance
rf_pop_classifier <- randomForest(x = training[, 1:(length(BRCA1gt2)-1)], y = as.factor(training$Population), ntree = 500, importance = TRUE) 
print(rf_pop_classifier) #Check OOB Error Rate, class error rate, and the confusion matrix
plot(rf_pop_classifier) 

varImpPlot(rf_pop_classifier) #Check variable importance

#Generating ROC Curves
classes <- levels(as.factor(training$Population)) #Generating classes to be used in ROC Curve
#colours <- c("RED","BLACK")#Colours to be used in the ROC Curves


predict_for_rf_pop_classifier <- predict(rf_pop_classifier, testing[, 1:(length(BRCA1gt2)-1)], type = "prob")
predict_for_rf_pop_classifier1 <- predict(rf_pop_classifier, testing[, 1:(length(BRCA1gt2)-1)])

print(table(observed=testing[,length(BRCA1gt2)],predicted=predict_for_rf_pop_classifier1))

for (i in 1:2){ #For loop set to i in 1:2 as there are two classes present, EAS and EUR
  true_values <- ifelse(testing$Population==classes[i],1,0) #Using an ifelse statement to correctly attribute samples to the right class
  pred <- prediction(predict_for_rf_pop_classifier[,i],true_values) #Calculating probability values for samples in each class using the inputted predict_for_rf_pop_classifier taken in by the function
  perf <- performance(pred, "tpr", "fpr") #Using pred as the prediction object, "fpr" (false positive rate) is the x axis for the ROC Curve and "tpr" is the y axis 
}

# Create a dataframe to plot the ROC curve
newdf <- as.data.frame(perf@x.values)
colnames(newdf)<- "x.values"
newdf2 <- as.data.frame(perf@y.values)
colnames(newdf2) <- "y.values"
brcaplot <- cbind(newdf,newdf2)
ggplot(data=brcaplot, aes(x=x.values, y=y.values))+
  geom_point(colour = "blue")

auc.perf <- performance(pred, measure = "auc") #Calculates the area under the curve value for each class
print(auc.perf@y.values) 


set.seed(44) #Ensure results are repeatable
train.index <- sample(1:nrow(HMGA2gt2), round(0.70*nrow(HMGA2gt2),0)) #70% of the data will be training, the rest are for testing
training <- HMGA2gt2[train.index,]
testing <- HMGA2gt2 [-train.index,] 
x = which(colnames(HMGA2gt2) == "Population")
#Creating the classifier, x has been set to include all columns except the last one (population) as predictor variables, mtry set to the default value (could tune this parameter like the example?). Importance set to true to later examine variable importance
rf_pop_classifier <- randomForest(x = training[, 1:(length(HMGA2gt2)-1)], y = as.factor(training$Population), ntree = 500, importance = TRUE) 
print(rf_pop_classifier) #Check OOB Error Rate, class error rate, and the confusion matrix
plot(rf_pop_classifier) 

varImpPlot(rf_pop_classifier) #Check variable importance

#Generating ROC Curves
classes <- levels(as.factor(training$Population)) #Generating classes to be used in ROC Curve
#colours <- c("RED","BLACK")#Colours to be used in the ROC Curves


predict_for_rf_pop_classifier <- predict(rf_pop_classifier, testing[, 1:(length(HMGA2gt2)-1)], type = "prob")
predict_for_rf_pop_classifier1 <- predict(rf_pop_classifier, testing[, 1:(length(HMGA2gt2)-1)])

print(table(observed=testing[,length(HMGA2gt2)],predicted=predict_for_rf_pop_classifier1))

for (i in 1:2){ #For loop set to i in 1:2 as there are two classes present, EAS and EUR
  true_values <- ifelse(testing$Population==classes[i],1,0) #Using an ifelse statement to correctly attribute samples to the right class
  pred <- prediction(predict_for_rf_pop_classifier[,i],true_values) #Calculating probability values for samples in each class using the inputted predict_for_rf_pop_classifier taken in by the function
  perf <- performance(pred, "tpr", "fpr") #Using pred as the prediction object, "fpr" (false positive rate) is the x axis for the ROC Curve and "tpr" is the y axis 
}

# Create a dataframe to plot the ROC curve
newdf <- as.data.frame(perf@x.values)
colnames(newdf)<- "x.values"
newdf2 <- as.data.frame(perf@y.values)
colnames(newdf2) <- "y.values"
hmgaplot <- cbind(newdf,newdf2)
ggplot(data=hmgaplot, aes(x=x.values, y=y.values))+
  geom_point(colour = "blue")

auc.perf <- performance(pred, measure = "auc") #Calculates the area under the curve value for each class
print(auc.perf@y.values) 



set.seed(44) #Ensure results are repeatable
train.index <- sample(1:nrow(IGF1gt2), round(0.70*nrow(IGF1gt2),0)) #70% of the data will be training, the rest are for testing
training <- IGF1gt2[train.index,]
testing <- IGF1gt2 [-train.index,] 
x = which(colnames(IGF1gt2) == "Population")
#Creating the classifier, x has been set to include all columns except the last one (population) as predictor variables, mtry set to the default value (could tune this parameter like the example?). Importance set to true to later examine variable importance
rf_pop_classifier <- randomForest(x = training[, 1:(length(IGF1gt2)-1)], y = as.factor(training$Population), ntree = 500, importance = TRUE) 
print(rf_pop_classifier) #Check OOB Error Rate, class error rate, and the confusion matrix
plot(rf_pop_classifier) 

varImpPlot(rf_pop_classifier) #Check variable importance

#Generating ROC Curves
classes <- levels(as.factor(training$Population)) #Generating classes to be used in ROC Curve
#colours <- c("RED","BLACK")#Colours to be used in the ROC Curves


predict_for_rf_pop_classifier <- predict(rf_pop_classifier, testing[, 1:(length(IGF1gt2)-1)], type = "prob")
predict_for_rf_pop_classifier1 <- predict(rf_pop_classifier, testing[, 1:(length(IGF1gt2)-1)])

print(table(observed=testing[,length(IGF1gt2)],predicted=predict_for_rf_pop_classifier1))

for (i in 1:2){ #For loop set to i in 1:2 as there are two classes present, EAS and EUR
  true_values <- ifelse(testing$Population==classes[i],1,0) #Using an ifelse statement to correctly attribute samples to the right class
  pred <- prediction(predict_for_rf_pop_classifier[,i],true_values) #Calculating probability values for samples in each class using the inputted predict_for_rf_pop_classifier taken in by the function
  perf <- performance(pred, "tpr", "fpr") #Using pred as the prediction object, "fpr" (false positive rate) is the x axis for the ROC Curve and "tpr" is the y axis 
}

# Create a dataframe to plot the ROC curve
newdf <- as.data.frame(perf@x.values)
colnames(newdf)<- "x.values"
newdf2 <- as.data.frame(perf@y.values)
colnames(newdf2) <- "y.values"
igf1plot <- cbind(newdf,newdf2)
ggplot(data=igf1plot, aes(x=x.values, y=y.values))+
  geom_point(colour = "blue")

auc.perf <- performance(pred, measure = "auc") #Calculates the area under the curve value for each class
print(auc.perf@y.values) 

# Enhance the dataframe for ROC curve for the 3 genes
brcaplot$model <- factor("BRCA1",c("BRCA1", "IGF1","HMGA2"))
igf1plot$model <- factor("IGF1",c("BRCA1","IGF1","HMGA2"))
hmgaplot$model <- factor("HMGA2",c("BRCA1","IGF1","HMGA2"))

roc_data <- rbind(brcaplot,igf1plot,hmgaplot)

ggplot(data = roc_data, aes(x=x.values, y=y.values, color=model)) +
  ggtitle("Random Forest ROC Curves") +
  xlab("Specificity") +
  ylab("Sensitivity") +
  geom_point()