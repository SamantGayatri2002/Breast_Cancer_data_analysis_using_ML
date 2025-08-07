#-----------------------#
# DEG Analysis - Part 1
#-----------------------#
# 1. Reading CEL files
# 2. Expression Normalization (RMA)
# 3. Pre and Post Normalization Expression Visualization


# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("limma")
BiocManager::install("hgu133plus2.db")
BiocManager::install("hgu133plus2cdf")


# Load packages
library(limma)
library(affy)
library(affyPLM)
library(hgu133plus2.db)
library(hgu133plus2cdf)
library(IRanges)
library(RColorBrewer)

getwd()
setwd("E:/Nyberman_internship/GSE21422_RAW_CELfiles")
getwd()

targets <- readTargets("target_H_DCIS.txt") 
targets

#Read CEL Files
data <- ReadAffy(filenames = targets$FileName)
data


#=====================#
# RMA Normalization
#=====================#
eset <- rma(data)
normset <- exprs(eset)

write.csv(normset, "ExpSet_PostNorm.csv", quote = F)


#=====================#
#   Box Plot
#=====================#
par(mfrow=c(1,2))

#Boxplot Before Normalization
tiff(file="Boxplot_Pre-Normalization.tiff", bg="transparent", width=400, height=500)
par(mar = c(12, 4, 6, 2) + 0.1); # This sets the plot margins
boxplot(data,col="red", main="Boxplot Pre-Normalization", las=2, cex.axis=0.74, ylab="Intensities" )
title(xlab = "Sample Array", line = 8); # Add x axis title
dev.off()

#Boxplot After Normalization
tiff(file="Boxplot_Post-Normalization.tiff", bg="transparent", width=400, height=500)
par(mar = c(12, 4, 6, 2) + 0.1); # This sets the plot margins
boxplot(normset,col="blue", main="Boxplot Post-Normalization", las=2, cex.axis=0.74, ylab="Intensities") #, col=colors 
title(xlab = "Sample Array", line = 8); # Add x axis title
dev.off()



#-----------------------#
# DEG Analysis - Part 2
#-----------------------#
# 1. Data Cleaning - Dimensionality Reduction (PCA)
# 2. Differentially Expressed Genes (DEGs) Identification 
# 3. DEGs Visualization (Volcano plot & Heatmap)


#=====================#
#   PCA Plot
#=====================#
install.packages("factoextra")
library(tidyverse)
library(factoextra)
library(limma)


data <- read.csv("ExpSet_PostNorm.csv")

nrow(data)  # Check the number of rows in your data
length(c(rep("DCIS", 9), rep("Healthy", 5)))  # Check the length of your group labels

# Adjust Group Labels
data_t <- t(data[,-1 ])  # Exclude the first col (gene names) during transpose
data_t <- as.data.frame(data_t)
data_t$Group <- c(rep("DCIS", 9), rep("Healthy", 5))

# Perform PCA
pca_res <- prcomp(data_t[, -ncol(data_t)])

# view all PC scores 
head(pca_res$x)
library(factoextra)


png("pca.png")
# Generate PCA plot
fviz_pca_ind(pca_res,
             geom.ind = c("point", "text"),
             col.ind = data_t$Group,
             palette = c("blue", "red"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Group",
             labelsize = 2
)
dev.off()

# Remove outliers
normset <- read.csv("ExpSet_PostNorm.csv", h=TRUE)
head(normset)
#normset <- normset[1:(ncol(normset)-1)]# this code is used if the oulier is present in last column and we need to remove it
normset <- normset[, -10]# i have used this code because my outlier was is 10th colum, i.e DCIS_9
head(normset)

#=======================#
# DEGs Identification
#=======================#
# Model Matrix Design
# Let's create a model matrix using the factor() function to represent the condition labels ("Healthy" or "DCIS")
design <- model.matrix(~factor(c("DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "DCIS", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy")))
# Now assign names ("Healthy" and "DCIS") to the columns of the model matrix
colnames(design) <- c("DCIS", "Healthy")

# Fits a linear model for each gene based on the given series of arrays
# It estimates the relationship between gene expression and conditions

fit <- lmFit(normset[,1:ncol(normset)], design) 
fit

# Contrast Matrix Design
# Define the specific comparison between conditions you want to analyze.
cont.matrix = makeContrasts(DCIS-Healthy, levels=design)
cont.matrix

# Fitting model with Contrasts(2 Groups), so apply the defined contrast to the previously fitted model (fit)
fit2 <- contrasts.fit(fit, cont.matrix)

# Model optimization / Empirical Bayes Moderation
# Improves the estimation of variances for genes with low expression
# Computes moderated t-statistics and log-odds (B-stats) of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
# contrast-specific information from fit2 is incorporated into the fit object, which is then passed to eBayes()
fit2 <- eBayes(fit)  
fit2

# Result Top Table
topTable(fit2, coef = 2, adjust.method = "BH") 

DEGS <- topTable(fit2, coef=2, adjust="BH", sort.by="logFC", number=100000); #inf
DEGS
write.csv(DEGS, "Result_Table_logFCsorted.csv", quote = F, row.names = FALSE)


#After this the the result table log file is manually sorted in excel file and annotated using biodbnet


# Filter & Save DEGs
#logFC_1 <- DEGS[DEGS$P.Value < 0.05 & (DEGS$logFC > 2 | DEGS$logFC < -2), ]
#write.csv(logFC_1,"pval_0.05_logFC_2.csv", quote = F)




#==========================#
# Annotate filtered DEGs
#==========================#
#BiocManager::install("hgu133plus2.db")
#library("hgu133plus2.db")

#probes=row.names(logFC_1)
#Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))

# Combine gene annotations with raw data
#deg_anno = cbind(probes,Symbols, logFC_1)
#write.csv(deg_anno, "DEGs_Annotated.csv", quote = F, row.names = F)


#-----------------------#
# DEG Analysis - Part 3
#-----------------------#
# DEGs Visualization 
#   - Volcano plot 
#   - Heatmap


#==================================================#
# Filter & Save final DEGs based on Pvalue & logFC
#==================================================#
# Read data of topTable
DEGs <- read.csv("Result_Table_logFCsorted.csv", header = TRUE)

# Filter & Save DEGs
final_DEGs <- DEGs[DEGs$P.Value < 0.05 & (DEGs$logFC > 2 | DEGs$logFC < -2), ]
write.csv(final_DEGs,"finalDEGs.csv", quote = F, row.names = F)



#=============================================#
# Annotate(getting Gene Symbols) filtered DEGs
#=============================================#
#BiocManager::install("hgu133plus2.db")
library("hgu133plus2.db")

DEGs <- read.csv("finalDEGs.csv", header = TRUE)
head(DEGs)

#open the the finalDEGs file manually and write the heading of affimetric IDs as Probe_IDs
probes=DEGs$Probe_ID
head(probes)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
head(Symbols)

# Combine gene annotations with raw data
deg_anno = cbind(probes,Symbols, DEGs)
write.csv(deg_anno, "DEGs_Annotated.csv", quote = F, row.names = F)



#================================================#
# DEG Viz (Volcano plot)
#================================================#
install.packages("gdata")
install.packages("gplots")
library(gdata)
library(gplots)

DEGs <- read.csv("Result_Table_logFCsorted.csv", h=T)
head(DEGs)

png(filename = "VolcanoPlot_FC_2.png")
with(DEGs, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot"))
with(subset(DEGs, P.Value < 0.05 & logFC > 2 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(DEGs, P.Value < 0.05 & logFC < -2), points(logFC, -log10(P.Value), pch=20, col="green"))
dev.off()



#========================#
# Probe/Gene Annotation
#========================#
BiocManager::install("hgu133plus2.db")
library("hgu133plus2.db")

normset <- read.csv("ExpSet_PostNorm.csv", h=TRUE)
head(normset)

# Match probe IDs and retrieve Gene SYMBOLS
probes=normset$X
head(probes)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))

# Combine gene annotations with raw data
normset_anno = cbind(probes,Symbols,normset)
write.csv(normset_anno, "ExpSet_PostNorm_Annotated.csv", quote = F, row.names = F)




#====================================#
# DEG Viz (Heatmap DEG Expression)
#====================================#
data <- read.csv(file = "heatmap_expdata.csv", h=T)
head(data)

rnames <- data[,1]
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- rnames
my_palette <- colorRampPalette(c("red", "blue"))(n = 299)
col_breaks = c(seq(0,5,length=100), # for red
               seq(5.1,10,length=100), # combo of red & blue
               seq(10.1,15,length=100)) # for blue
tiff("heatmap_exp_deg_cluster.tiff",     
     width = 6*300,        # 5 x 300 pixels
     height = 6*300,
     res = 300,            # 300 pixels per inch
     pointsize = 8)        # smaller font size

heatmap.2(mat_data,
          main = "Heatmap", # heat map title
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks, 
          dendrogram="both",     # only draw a row dendrogram
          Colv="T" ,         # turn off column clustering
          lhei = c(1,7)         # Key size width adjustment
)            
dev.off()


# ML (Machine Learning) 
# Model Development for Gene Expression Analysis


# Load the required packages
install.packages("xlsx")
install.packages("caret")
install.packages("glmnet")

# Importing required libraries
library("readxl")  
library("tidyverse") # data manipulation and visualization package
library("caret")  # machine learning library package
library("glmnet") # implement computing penalized regression - elasticNet

# Prerequisites 
graphics.off()  # clear all graphs from RStudio
rm(list = ls()) # remove all files from your workspace
# Ctrl+L -> Clear Console


#======================#
## load Data
#======================#
df <- read_excel("shortTable.xlsx", sheet = 1)
head(df)

#============================================#
# Data Pre-processing |   EDA  
#============================================#

{  # Check the structure of the dataset
  
  str(df)
  
  dim(df)
  
  
  # remove rows that contain NA values
  df <- df[complete.cases(df), ]
  head(df)
  dim(df)
  
  
  #Calculate Mean of duplicate genes
  x <- df
  x <- data.frame(x)
  x <- do.call(rbind,lapply(lapply(split(x,x$Symbols),`[`,2:ncol(x)),colMeans))
  dim(x)
  
  
  #Convert rownames as a 1st column with header Symbols -> which became rownames after previous operation
  library(tibble) # from tidyverse
  x <- data.frame(x)
  x <- tibble::rownames_to_column(x, var="Symbols")
  head(x)
  dim(x) 
  
  df <- x
  
  # Transpose table 
  install.packages("sjmisc")
  library(sjmisc)  
  df_t <- rotate_df(df, cn=T)
  Symbols <- colnames(df[-1])
  df_t <- cbind(Symbols, df_t)
  write.csv(df_t, "transposed_table.csv", row.names=F)
  
  df_t <- read.csv("transposed_table.csv", h=T)
  dim(df_t)
  df_t[1]
  
  # Healthy_1.CEL to N/T
  df_t[,1] <- gsub("_.*$", "", df_t[,1])
  df_t[1]
  
  # convert Healthy/DCIS from char to factor
  df_t[1] <- factor(df_t$Symbols)
  str(df_t)
  
  # df_t -> df
  df <- df_t # df is df_NT
  
}


# view transformed data
str(df)
df <- data.frame(df)
head(df)

#===============================================#
# Step 3.  Visualize Dataset - Figures - Plots 
#===============================================#
#####  Box and Whisker Plots  ##### 
# Given that the input variables are numeric, we can create box and whisker plots of each
png("box_and_whisker_plots.png")
par(mfrow=c(2,4))
for(i in 2:9) {
  boxplot(x[,i], main=names(df)[i], col="blue")
}
dev.off()


#####  Sample matrix  ##### 
library(ggplot2)
library(caret) 
# split input and output
x <- df[,2:ncol(df)]  # x -  inputs attributes
y <- df[,1]     # y -  outputs attributes
y<- as.factor(y)
plot(y, col="blue")
#featurePlot(x=x, y=y, plot="ellipse") # time intensive


#======================#
# Step 2. Data Splitting
#======================#
# create 70%/30% for training and testing dataset
library(caret)
set.seed(101)
split <- createDataPartition(df$Symbols, p=0.70, list=FALSE)  # Return the row indices as a matrix/vector, not as a list.
train <- df[split,]
test <- df[-split,]

# dimensions of dataset, train, test
dim(df)
dim(train)
dim(test)


#set cross-validation control for training
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

#{
 # library(yardstick)  # for metric_set()
  #metric1 <- metric_set(rmse, rsq, mae, ccc) #need to check
 # metric2 <- metric_set(accuracy, kap, sens, spec, roc_auc, pr_auc) #need to check
#}

#=========================#
# Step 3. Build ML Models   
#=========================#
{ #extra start.....
  ##This part of code has not been executed in class
  
  head(df_t)
  
  # kNN (k Nearest Neighbour)
  set.seed(101)
  fit.knn <- train(Symbols~., data=train,
                   method="knn",  metric = metric ,trControl=control)
  fit.knn
  summary(fit.knn)
  
  # Feature Importance
  varImp(fit.knn)
  plot(varImp(fit.knn, scale=FALSE), top=20) 
  
  # export plot  
  tiff("kNN_varimp.tiff")
  plot(varImp(fit.knn, main="k-Nearest Neighbors", scale=TRUE), top=30)
  dev.off()
  
  # Make Predictions # Confusion Matrix #model evaluation
  predictions <- predict(fit.knn, test)
  conf_matrix <- confusionMatrix(predictions, as.factor(test$Symbols))
  conf_matrix
  
  
  
  # SVM (Support Vector Machine) - Radial
  install.packages("kernlab")
  library(kernlab)
  set.seed(101)
  fit.svm <- train(Symbols~., data=train, method="svmRadial", trControl=control)
  fit.svm
  
  # Feature Importance
  varImp(fit.svm)
  plot(varImp(fit.svm, scale=FALSE), top=20)
  # export plot
  png("svm_varimp.png")
  plot(varImp(fit.svm, main="Support Vector Machines with Radial Basis", scale=TRUE), top=30)
  dev.off()
  
  # Make Predictions # Confusion Matrix
  predictions <- predict(fit.svm, test)
  conf_matrix <- confusionMatrix(predictions, as.factor(test$Symbols))
  conf_matrix
  
  # extract the actual class labels from the test data
  # actual_classes <- test$Symbols
  # convert the class probabilities to predicted class labels
  # predicted_classes <- ifelse(predicted_probs > 0.5, "versicolor", "setosa")
  
  {
    ## Installing the package 
    install.packages("caTools")    # For Logistic regression 
    install.packages('pROC')       # For ROC curve to evaluate model 
    library(caTools)
    library(pROC)  
    library(ggplot2)
    
    # Logistic regression 
    #  E-net  glmnet  
    set.seed(7)
    fit.enet <- train(Symbols ~ ., data = train, method ='glmnet', #metric="Accuracy", 
                      type.measure="deviation", family="binomial",
                      tuneGrid = expand.grid(alpha = seq(0,1,length=10), lambda = seq(0.0001,0.2,length=5)),
                      trControl = control)
    print(fit.enet)    
    #print(fit.enet$finalModel)
    varImp(fit.enet)
    plot(varImp(fit.enet), top=20, main="ENET") #, col="red")
    plot(fit.enet, main="ENET")
    
    pred.enet = predict(fit.enet, newdata = test)#newdata=test[-3], type='response')
    
    # Model Evaluation
    confusionMatrix(pred.enet, test$Symbols, positive = "T")
    c5 <- confusionMatrix(pred.enet, test$Symbols, positive = "T")$table
  }
  
  
  # Naive Bayes  # package - klaR 
  library(klaR)
  set.seed(7)
  fit.nb <- train(Symbols~., data=train, 'nb', trControl=control)# , usekernel = T    # trControl=trainControl(method='cv', number=10)
  fit.nb
  varImp(fit.nb)
  #Plot Variable performance
  plot(varImp(fit.nb), top=20, main="nvBayes") 
  
  # Model Evaluation (Confusion Matrix)
  confusionMatrix(pred.nb, test$Symbols, positive = "T")
  cm_nv <- confusionMatrix(pred.nb, test$Symbols, positive = "T")$table
  
  
} # extra end
## The below codes can be executed for ML analysis


head(df_t)


#=========================#
# Build ML Models   
#=========================#

# ML (Machine Learning) 
# Model Development for Gene Expression Analysis - [Part 2]


# Install and Load the required packages
install.packages("caret")
library("caret")


# 1... kNN(k-Nearest Neighbor) - [Model 1]
#-------------------------------------------
set.seed(7)
fit.knn <- train(Symbols~., 
                 data=train, 
                 method="knn", 
                 metric=metric, 
                 trControl=control)
fit.knn


install.packages("cowplot")     # Only once
library(cowplot)                # Load every time you use plot_grid()

# check important variables
varImp(fit.knn)
p1 <- plot(varImp(fit.knn), top = 20, main="kNN")
p2 <- plot(fit.knn, main="kNN")
plot_grid(p1, p2)

# make predictions using trained model on new/test
pred.knn <- predict(fit.knn, newdata = test)

# Model Evaluation
confusionMatrix(pred.knn, test$Symbols, positive = "DCIS")
c1 <- confusionMatrix(pred.knn, test$Symbols, positive = "DCIS")$table

#========================#
# plot Confusion Matrix
#========================#
library(ggplot2)
library(dplyr)
table <- data.frame(c1)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = Freq)) + # alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(high = "#009194", low="#FF9966")) +
  #scale_fill_gradient(low="white", high="#009194") +
  theme_bw() +
  xlim(rev(levels(table$Reference)))



# 2.SVM model     [model-2]
#-----------------------------------
# For bioinformatics tasks like gene expression classification:
# Use svmLinear when classes are clearly separable (e.g., PCA shows clusters),
# Use svmRadial for complex, non-linear patterns (common in omics data).

install.packages("kernlab")
library(kernlab)
set.seed(101)
fit.svm <- train(Symbols~., data=train, method="svmRadial", metric=metric, trControl=control)
fit.svm

# Feature Importance
varImp(fit.svm)
plot(varImp(fit.svm, scale=FALSE), top=20)
# export plot
png("svm_varimp.png")
plot(varImp(fit.svm, main="Support Vector Machines with Radial Basis", scale=TRUE), top=30)
dev.off()

# Make Predictions # Confusion Matrix
predictions <- predict(fit.svm, test)
conf_matrix <- confusionMatrix(predictions, as.factor(test$Symbols))
conf_matrix




# 3 .Random Forest Model - [Model 3]
#----------------------------------------
install.packages("randomForest")
library(randomForest)

set.seed(123)
fit.rf <- train(Symbols~.,
                data=train,
                method="rf",
                mtry=2,# this line is optional, it tells how may decision trees are needed, in RF model. without this also code runs.
                metric=metric,
                trControl=control)
fit.rf

# view important genes
varImp(fit.rf)


# save the list of important genes in a file
var <- varImp(fit.rf)
var_df <- as.data.frame(var$importance)
var_df_sorted <- var_df[order(-var_df$Overall), , drop = FALSE]
write.csv(var_df_sorted, "random_forest_variable_importance.csv", row.names = TRUE)

# visualize the important genes
plot(varImp(fit.rf), top = 30)

# make predictions using trained model on new/test
pred.rf <- predict(fit.rf, newdata = test)

# Model Evaluation
confusionMatrix(pred.rf, test$Symbols, positive = "DCIS")
c1 <- confusionMatrix(pred.rf, test$Symbols, positive = "DCIS")$table



#========================#
# plot Confusion Matrix
#========================#
library(ggplot2)
library(dplyr)
table <- data.frame(c1)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = Freq)) + # alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(high = "#009194", low="#FF9966")) +
  #scale_fill_gradient(low="white", high="#009194") +
  theme_bw() +
  xlim(rev(levels(table$Reference)))



