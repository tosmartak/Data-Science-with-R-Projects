if(!require(pacman)) install.packages(pacman)
pacman::p_load(  # Use p_load function from pacman
  readr,
  caTools,
  FNN,
  tree,
  kernlab,
  neuralnet,
  caret,
  inspectdf,
  plotly,
  dplyr,
  tidyr,
  broom,
  ggplot2,
  car,
  lmtest,
  lessR,
  psych,
  janitor,
  visdat,
  gapminder,
  esquisse,
  readxl
)

# Load Thrive data
thrive_vulnerability_data <- read_excel("thrive_vulnerability data.xlsx")
View(thrive_vulnerability_data)

#Using K-means for two clusters----
thrive <- thrive_vulnerability_data
head(thrive)
thrive_work <- thrive[c(5,11,13,15)]
head(thrive_work)
as.double(thrive_work$endline_rsci_score)
thrive_work[, c(1,2,3,4)] <- lapply(thrive_work[,c(1,2,3,4)], as.numeric)
str(thrive_work[, c(1,2,3,4)])

thriveD <- data.frame(scale(thrive_work[,c(1,2,3,4)]))
#create a numeric vector, wss, to keep track of the within sum of squares values
wss<-numeric(10)

# use the following line of code to assign the withinss for 1 through 10 clusters

for(i in 1:10){wss[i]<-sum(kmeans(thriveD,centers=i)$withinss)}
# plot wss
plot(wss,type="b") #the plot didn't help much because the intervals appear equally distributed 
##but since the desire is generate two clusters i used "2" centres" in the codes below

thriveKM<-kmeans(thriveD,centers=2)
#assign the cluster numbers to a new variable in the original data


thrive_work$KMSegment<-thriveKM$cluster
thrive_work


#calculate the *unscaled* means for your clusters with the aggregate() function: 

aggregate(thrive_work[c("hh_size","endline_fcs_score", "endline_rsci_score", "endline_hhs_score")],by=list(KMSegment=thrive_work$KMSegment), mean)
#Based on the aggregates, the interpretation of the cluster is "1" means recovered, "2" means unrecoverd

# Examining the Clusters on 2-dimension
# Perform PCA on the scaled data (thriveD)
pca_result <- prcomp(thriveD)
# Extract the first two principal components
pc1 <- pca_result$x[, 1]
pc2 <- pca_result$x[, 2]
# Create a scatter plot using the first two principal components with cluster coloring
plot(pc1, pc2, col = thriveKM$cluster, pch = 19, 
     main = "PCA with K-Means Clustering - 2D", xlab = "Principal Component 1", ylab = "Principal Component 2")
# Calculate cluster centers manually
cluster_centers <- aggregate(thriveD, by=list(thriveKM$cluster), FUN=mean)
# Add cluster centers to the plot
points(cluster_centers[, 1], cluster_centers[, 2], col = 1:2, pch = 8, cex = 2)
# Add a legend
legend("topright", legend = unique(thriveKM$cluster), col = 1:2, pch = 19,
       title = "Cluster")

#I also tried out hierrachical clusters below to see whether it generates better cluster output

#Using Hierrachichal clusters----
library(readr)
thriveH <- thrive_vulnerability_data
# check variable types, make sure abv and ibu are doubles, recast as necessary
head(thriveH)

# I used the code below to set up new scaled dataframe with only the variables we need
thriveDH<-data.frame(scale(thrive_work[,c(1,2,3,4)]))

# I used the dist() function to calculate euclidean distances between each household, and
# assign the result to the object thriveDist.
thriveDist<-dist(thriveDH,method = "euclidean", diag = FALSE)
#I used hclust to perform clustering on thriveDist,and assigned it to thriveHC
thriveHC <-hclust(thriveDist,method="complete")

#I cut the tree using cutree() into 2 to be able to compare with the output of Kmeans output
# i.e. use thrive$HCSegment<-cutree()
CutTree<-cutree(thriveHC,2)
thrive_work$HCSegment<-CutTree

#I then generated cluster means for key variables using the aggregate() function
aggregate(thrive_work[c("hh_size","endline_fcs_score", "endline_rsci_score", "endline_hhs_score")],by=list(HCSegment=thrive_work$HCSegment), mean)

# Create a scatter plot using the first two principal components with hcluster
plot(pc1, pc2, col = thrive_work$HCSegment, pch = 19, 
     main = "PCA with Hierarchical Clustering - 2D", xlab = "Principal Component 1", ylab = "Principal Component 2")
# Calculate cluster centers manually
cluster_centers <- aggregate(thriveD, by=list(thrive_work$HCSegment), FUN=mean)
# Add cluster centers to the plot
points(cluster_centers[, 1], cluster_centers[, 2], col = 1:2, pch = 8, cex = 2)
# Add a legend
legend("topright", legend = unique(thrive_work$HCSegment), col = 1:2, pch = 19,
       title = "Cluster")


##Comment: The output of the hierrachical clustering is similar to that of Kmeans but 
###didn't quite separate the HHsize like the Kmeans model. So,on this basis, I will go with the output 
####of the later.


# To further validate how distinct the two clusters are, i test statistical significance of each of 
##the variables in the model below

# Testing statistical significance betw for the 2 unsupervised learning model outputs----

## for FCS score
t_test_result <- t.test(thrive_work$endline_fcs_score ~ thrive_work$KMSegment)
t_test_result

t_test_result <- t.test(thrive_work$endline_fcs_score ~ thrive_work$HCSegment)
t_test_result


## for rCSI

t_test_result <- t.test(thrive_work$endline_rsci_score ~ thrive_work$KMSegment)
t_test_result

t_test_result <- t.test(thrive_work$endline_rsci_score ~ thrive_work$HCSegment)
t_test_result


## for HHS
t_test_result <- t.test(thrive_work$endline_hhs_score ~ thrive_work$KMSegment)
t_test_result

t_test_result <- t.test(thrive_work$endline_hhs_score ~ thrive_work$HCSegment)
t_test_result

## for HH size
t_test_result <- t.test(thrive_work$hh_size ~ thrive_work$KMSegment)
t_test_result

t_test_result <- t.test(thrive_work$hh_size ~ thrive_work$HCSegment)
t_test_result

# Getting out an output file in CSV----

write.csv(thrive_work, "output_file.csv", row.names = FALSE)

output_with_id <- cbind(thrive$hh_id,thrive_work)
head(output_with_id)
write.csv(thrive_work, "output_with_id_file.csv", row.names = FALSE)

#Comment on Model fit----
#I used Kmeans because the two models separate the FCS, HHS and rCSI in similar way 
## but kmeans seems to produce higher significant values compared to Hierrarchical clusters



## Supervised Learning: Classification----
#K_Nearest Neighbours
library(readr) # load the required packages into R
library(class)
library(caTools)

# Load Thrive data with label
output_with_id_file <- read_csv("output_with_id_file.csv", show_col_types = FALSE)


thrivecl <- output_with_id_file
head(thrivecl)
thrive_trim <- thrivecl[,c(1:5)] 
head(thrive_trim)

# Convert the class labels to factors with levels
thrive_trim$KMSegment <- factor(thrive_trim$KMSegment)

set.seed(1)
sample<-sample.split(thrive_trim$KMSegment, SplitRatio = .80)
train<-subset(thrive_trim, sample == TRUE)
test<-subset(thrive_trim, sample == FALSE)

# Prepare data for all models
X_train <- train[, -which(names(train) == "KMSegment")]
y_train <- train$KMSegment
X_test <- test[, -which(names(test) == "KMSegment")]
y_test <- test$KMSegment


# Train knn model
knn1<-knn(X_train, X_test, y_train, k=1) 
knn1

#Accuracy of the K-Nearest Neighbours Model----

#use the table function to calculate the precision of your predictions
CF<-table(knn1,y_test)
CF
Precision_k1 <-CF[2,2]/(CF[2,1]+CF[2,2])
Precision_k1

confusion <- confusionMatrix(knn1, y_test)
confusion

# Calculate the F1-score
f1_score <- confusionMatrix(knn1, y_test)$byClass["F1"]
f1_score

## One nearest neighbour gives the most precise model, The model precision decreases with increasing number of K



# Using Logistic Regression Method ----
thrive_trim$KMSegment<-as.factor(thrive_trim$KMSegment)
thrive_trim

sum(is.na(thrive_trim[, c("endline_fcs_score", "endline_rsci_score", "endline_hhs_score")]))


#Excluding household size from the model
thriveLM <- glm(KMSegment ~ endline_fcs_score + endline_rsci_score + endline_hhs_score, family = binomial(), data = thrive_trim)

summary(thriveLM)
exp(coef(thriveLM))

test$recover <- predict(thriveLM,test,type= "response" )
test$recover.d <-round(test$recover) #this is to round up the numbers
test

##Accuracy and Precision of the Logistics Regression Model----
unique(test$recover.d)
test$recover.d <- factor(test$recover.d, levels = c(0,1), labels = c("unrecovered", "recovered"))
test

unique(test$KMSegment)
test$KMSegment <- factor(test$KMSegment, levels = c(1,2), labels = c("unrecovered", "recovered"))
test

AccuracyofLM <- confusionMatrix(test$recover.d,test$KMSegment)
AccuracyofLM

## Naive Bayes----
library(readr)
library(caTools)

thrivecl <- output_with_id_file
head(thrivecl)
thrive_trim <- thrivecl[,c(1:5)]
head(thrive_trim)
thrivecl$hh_size <- as.factor(thrivecl$hh_size)
thrive_trim$endline_fcs_score <- as.factor(thrive_trim$endline_fcs_score)
thrive_trim$endline_rsci_score <- as.factor(thrive_trim$endline_rsci_score)
thrive_trim$endline_hhs_score <- as.factor(thrive_trim$endline_hhs_score)
thrive_trim$KMSegment <-  as.factor(thrive_trim$KMSegment)
set.seed(1)
sample<-sample.split(thrive_trim$KMSegment, SplitRatio = .80)
train<-subset(thrive_trim, sample == TRUE)
train.att <- train[-5] #Unlabeled training data
train.res <- train$KMSegment
test<-subset(thrive_trim, sample == FALSE)
test.att <- test[-5] #Unlabeled testing data
test.res <- test$KMSegment

library(e1071)
NBayes<-naiveBayes(train.att,train.res)
NBayesPredict<-predict(NBayes,test.att,type="class",laplace=0)

#Accuracy and Precision of Naive Bayes ----
library(gmodels)
CrossTable(NBayesPredict,test.res,prop.chisq = FALSE,prop.t =FALSE,prop.c = FALSE,prop.r = FALSE)
PerCorrect<-sum(NBayesPredict==test.res)/nrow(test)
PerCorrect
#Bayes's accuracy is comparable with Logistics regression(93%) but and KNN yields 95%

#### Using Tree Method----

library(readr)
library(class)
library(caTools)
library(tree)
thrivecl <- output_with_id_file
head(thrivecl)
thrive_trim <- thrivecl[,c(1:5)]
head(thrive_trim)
thrive_trim$KMSegment <- as.factor(thrive_trim$KMSegment)

set.seed(123) 
sample<-sample.split(thrive_trim$KMSegment, SplitRatio = .80)
train<-subset(thrive_trim, sample == TRUE)
test<-subset(thrive_trim, sample == FALSE)

thriveTree <- tree(KMSegment ~ hh_size + endline_hhs_score+endline_rsci_score+ endline_fcs_score, data = train)
#plot and visualize your tree
plot(thriveTree)
text(thriveTree, cex=.5)
summary(thriveTree)
#use the predict function to see how well your model performs on the test subset
Predrecovery <- predict(thriveTree, test, type="class")
 ## Accuracy/Precision of the Tree Method Model----
CF<-table(test$KMSegment,Predrecovery)
CF
Precision<-CF[2,2]/(CF[2,1]+CF[2,2])
Precision
###Tree method is less accurate compared to KNN





#________________________________________________________________________________________________________________
#install.packages("DALEX")
library(DALEX)
library(caret)
library(rpart)
library(rpart.plot)

# Try KNN with additional cv method
# Define the KNN model
knn_model <- train(x = X_train, y = y_train, method = "knn", trControl = trainControl(method = "cv", number = 5))

# Make predictions with the KNN model
knn_predictions <- predict(knn_model, newdata = X_test)

# Create a confusion matrix
knn_confusion_matrix <- confusionMatrix(knn_predictions, y_test)
knn_confusion_matrix

# Extract F1 score from the confusion matrix
knn_f1_score <- knn_confusion_matrix$byClass["F1"]
print(paste("KNN Model F1 Score:", knn_f1_score))




# Recreate a new logistic regression model
# Define the logistic regression model
lr_model <- train(x = X_train, y = y_train, method = "glmnet", trControl = trainControl(method = "cv", number = 5))

# Make predictions with the LR model
lr_predictions <- predict(lr_model, newdata = X_test)

# Create a confusion matrix
lr_confusion_matrix <- confusionMatrix(lr_predictions, y_test)
lr_confusion_matrix

# Extract F1 score from the confusion matrix
lr_f1_score <- lr_confusion_matrix$byClass["F1"]
print(paste("Logistic Regression Model F1 Score:", lr_f1_score))

# Get variable importance scores (absolute values of coefficients)
lr_importance_scores <- varImp(lr_model, scale = FALSE)$importance

# Create a data frame for plotting
lr_importance_df <- data.frame(
  Feature = rownames(lr_importance_scores),
  Importance = lr_importance_scores
)
# Create a bar plot of feature importance
ggplot(data = lr_importance_df, aes(x = reorder(Feature, -Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Logistic Regression Model Feature Importance", x = "Feature", y = "Importance")




# Train Decision Tree
# Examine the tree plot
dtree_model <- rpart(KMSegment ~ ., data = train, method = "class") # decision tree explain

# Make Predictions with the dtree Model
dtree_predictions <- predict(dtree_model, test, type = "class")

# Create a confusion matrix
dtree_confusion_matrix <- confusionMatrix(dtree_predictions, y_test,)
dtree_confusion_matrix

# Extract F1 score from the confusion matrix
dtree_f1_score <- dtree_confusion_matrix$byClass["F1"]
print(paste("Decision Tree Model F1 Score:", dtree_f1_score))

# Plot the decision tree
rpart.plot(dtree_model, extra = 4)




# Train Random Forest
# Define the Random Forest model
rf_model <- train(x = X_train, y = y_train, method = "rf",
                  trControl = trainControl(method = "cv", number = 5, search = "grid"))

# Make Predictions with the LR Model
rf_predictions <- predict(rf_model, newdata = X_test)

# Create a confusion matrix
rf_confusion_matrix <- confusionMatrix(rf_predictions, y_test)
rf_confusion_matrix

# Extract F1 score from the confusion matrix
rf_f1_score <- rf_confusion_matrix$byClass["F1"]
print(paste("Random Forest Model F1 Score:", rf_f1_score))

# Get variable importance scores (absolute values of coefficients)
rf_importance_scores <- varImp(rf_model, scale = FALSE)$importance

# Create a data frame for plotting
rf_importance_df <- data.frame(
  Feature = rownames(rf_importance_scores),
  Importance = rf_importance_scores
)
# Create a bar plot of feature importance
ggplot(data = rf_importance_df, aes(x = reorder(Feature, -Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Random Forest Model Feature Importance", x = "Feature", y = "Importance")




# Train Bayesian Model
# Define the Bayesian model
baye_model <- train(x = X_train, y = y_train, method = "bayesglm", trControl = trainControl(method = "cv", number = 5))

# Make predictions with the LR model
baye_predictions <- predict(baye_model, newdata = X_test)

# Create a confusion matrix
baye_confusion_matrix <- confusionMatrix(baye_predictions, y_test)
baye_confusion_matrix

# Extract F1 score from the confusion matrix
baye_f1_score <- baye_confusion_matrix$byClass["F1"]
print(paste("Bayesian Model F1 Score:", baye_f1_score))

summary(baye_model)




library(pROC)
# Plot ROC Curves for all the Models
# Calculate ROC curves and AUC for Random Forest, Decision Tree, KNN, and Logistic Regression
roc_rf <- roc(y_test, as.numeric(rf_predictions))
roc_dtree <- roc(y_test, as.numeric(dtree_predictions))
roc_knn <- roc(y_test, as.numeric(knn_predictions))
roc_lr <- roc(y_test, as.numeric(lr_predictions))
roc_baye <- roc(y_test, as.numeric(baye_predictions))

# Calculate AUC values for all three models
auc_rf <- auc(roc_rf)
auc_dtree <- auc(roc_dtree)
auc_knn <- auc(roc_knn)
auc_lr <- auc(roc_lr)
auc_baye <- auc(roc_baye)

# Create a plot with all ROC curves and AUC values in the legend labels
plot(roc_rf, col = "blue", main = "ROC Curves for the Different Models")
lines(roc_dtree, col = "green")
lines(roc_knn, col = "yellow")
lines(roc_lr, col = "red")
lines(roc_baye, col = "orange")
# Create legend labels with AUC values for all models
legend_labels <- c(
  paste("Random Forest (AUC = ", round(auc_rf, 2), ")", sep = ""),
  paste("Decision Tree (AUC = ", round(auc_dtree, 2), ")", sep = ""),
  paste("KNN (AUC = ", round(auc_knn, 2), ")", sep = ""),
  paste("Logistic Reg. (AUC = ", round(auc_lr, 2), ")", sep = ""),
  paste("Bayesian (AUC = ", round(auc_baye, 2), ")", sep = "")
)
# Add the legend with modified labels
legend("bottomright", legend = legend_labels, col = c("blue", "green", "yellow", "red", "orange"), lty = 1)



# Plot the F1 Scores
# Create a vector of F1 scores
f1_scores <- c(knn_f1_score, lr_f1_score, dtree_f1_score, rf_f1_score, baye_f1_score)
# Names of the models
model_names <- c("KNN", "Logistic Reg.", "Decision Tree", "Random Forest", "Bayesian")
# Create a data frame for the F1 scores and model names
pltdata <- data.frame(Model = model_names, F1_Score = f1_scores)
# Set y-axis limits
y_limits <- c(0, 1.1)
# Create a bar chart with ggplot2
p <- ggplot(pltdata, aes(x = Model, y = F1_Score, fill = Model)) +
  geom_bar(stat = "identity") +
  ylim(y_limits) +
  labs(title = "F1 Scores for Different Models", x = "Models", y = "F1 Score") +
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = round(F1_Score, 2)), vjust = -0.5, size = 3, color = "black")
# Print the plot
print(p)

#________________________________________________________________________________________________________________


