### Testing for phenotypic limits
### Authors: Angelo Soto-Centeno & Pedro Ivo Monico
### Data: 16 linear phenotypic skull and dentary characters
# date: 14 of September, 2023
################################################


# Packages
library(mice)
library(ade4)
library(tidyverse)
library(caret)
library(MASS)
library(car)
library(gridExtra)
library(hrbrthemes)
library(gridExtra)
library(ggplot2)
library(tidyverse)

#for more information about packages
#?package.name

###########################################################################################

#Set working directory with the downloaded data 
setwd("/Users/Pedro/Desktop/Eptesicus_fuscus/morphometrics")

#Upload dataset
morpho_final<-read.csv(file = "./morphodata_COMPLETE_3GROUPS.csv", header = T)
dim(morpho_final)
View(morpho_final)


##########################################################################################

# Extract the group names
sp.names <- morpho_final[19]
head(sp.names)

# Log transform to "normalize" the values
# Displacing the data into normal distribution (one of the model assumptions);
sp.ln <- log(morpho_final[, c(3:18)])
sp.ln 

# Combine names + ln data
sp <- cbind(sp.names, sp.ln)
head(sp)

##########################################################################################

# Create training/testing splits
set.seed(123) #this is a random number; analysis will start at "'123" position 
training <- sp$Locality %>% createDataPartition(p = 0.75, list = F) #3/4 of total data used to training purposes 
# training set
train.sp <- sp[training, ]
dim(train.sp)
head(train.sp)
# testing set
test.sp <- sp[-training, ]

# Visualize the data splits; Scatterplots show if data is Gausian; 
# In this case, the arrangements are done to allow visualization
scatterplotMatrix(train.sp[2:8]) 
scatterplotMatrix(train.sp[9:17])

# Create boxplots help you to see the distribution of means
boxplot(train.sp[, c(2:17)], main = "Raw Data")

# Preprocess to center and scale the data (i.e. normalize)
prep.train.sp <- preProcess(train.sp[, c(2:17)], method = c("center", "scale"))
prep.train.sp.d <- predict(prep.train.sp, train.sp[, c(2:17)])

# Boxplot to verify normalization
boxplot(prep.train.sp.d, main = "Normalized Data")

# Fit the LDA model to training dataset
model.fit <- train(Locality~., data = train.sp, preProcess = c("center", "scale"), method = "lda")

# Check results
model.fit$finalModel

##########################################################################################

# Make the model predictions to testing dataset
predictions <- predict(model.fit, newdata = test.sp)
# Make a table to ensure both data & reference are factors
t <- table(factor(predictions), factor(test.sp$Locality))
t

# calculate the Confusion Matrix
confusionMatrix(t)

###########################################################################################
# CROSS VALIDATION 
# The number of reps depends on the size of your data (for more info, read: https://towardsdatascience.com/5-reasons-why-you-should-use-cross-validation-in-your-data-science-project-8163311a1e79)

kfoldcv <- trainControl(method = "cv", number = 5) 
performance.metric <- "Accuracy"


##########################################################################################
#FINAL LDA MODELS
set.seed(123)
# run Linear Discriminant Analysis LDA
brownbat.lda <- lda(Locality ~ ., data = sp)
brownbat.lda

## Display the results of LDA. Note the proportion of trace
## The values presented by the LDA function is the % separation achieved by each discriminant function


# Confusion matrix of LDA machine learning classification; Assess the prediction accuracy of the LDA from the model to the actual data
brownbat.lda.predict <- train(Locality ~ .,  data = sp, method = "lda", metric = performance.metric,
                              trControl = kfoldcv, preProcess = c("center", "scale"))

# calculate a confusion matrix to see the accuracy of classification
confusionMatrix(as.factor(sp$Locality), predict(brownbat.lda.predict, sp))  

## To interpret Kappa, see here: https://stats.stackexchange.com/questions/82162/cohens-kappa-in-plain-english

##########################################################################################

# create stacked histograms of the LDA values for csv subspecies
brownbat.lda.values <- predict(brownbat.lda)
ldahist(brownbat.lda.values$x[ ,1], g = sp$Locality, 
        type = "histogram", col = "light gray") 

## note: if "margins are too large", expand the plot viewer & do this:
## par(mar = c(4,3,1,1)) # tweak each side of the plot to make it view able

# create a scatterplot of LD to visualize discrimination; convert LD data into a data.frame for the analyzed dataset
brownbat_LD <- data.frame(type = sp[,1], lda = brownbat.lda.values$x)
head(brownbat_LD)


##########################################################################################
#VISUALIZING DATA
## step-by-step ggplot for LDA with all features

# set theme
theme_set(theme_light())
colrs <- c("#9e9ac8", "#807dba", "#6a51a3", "#4a1486")

F1 <- ggplot(brownbat_LD, aes(lda.LD1, lda.LD2, color = type)) + 
  geom_point(size = 2) + stat_ellipse(level = 0.68)

# ADJUST ARGUMENT "labs" (x="LD1" TO THE PRODUCED VALUES)
F1a <- F1 + scale_color_manual(values = colrs,
                               aesthetics = c("colour", "fill"),
                               labels = c("Caribbean", "Mesoamerica", "Southeastern North America")) +
  geom_rug() +
  labs(x = "LD1 (98.53%)", y = "LD2 (1.47%)") +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.position = "none")

#VISUALIZE
F1a


##########################################################################################
# MARGINAL DENSITY PLOT 
# For more info, read:  http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
brownbatFig1density <- ggplot(brownbat_LD, aes(lda.LD1, fill = type)) + 
  geom_density(alpha = 0.75) + 
  scale_fill_manual(values = c("#9e9ac8", "#807dba", "#6a51a3", "#4a1486"),
                    labels = c("Caribbean", "Mesoamerica", "Southeastern North America")) + 
  theme(legend.position = "right", axis.title.x = element_blank(), axis.title.y = element_text(size = 18))

grid.arrange(F1a, brownbatFig1density,
             ncol = 1, nrow = 2, heights = c(4, 1.4))

##########################################################################################
#PCA
#Examine PC contributions & rotation
vil_pca <- prcomp(sp.ln, center = T, scale. = T)

summary(vil_pca)

vil_pca$rotation

#principle component contribution
qplot(c(1:16), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

pca.table<-as.data.frame(results$x)
caribbean.names.pca<-morpho_final$Locality
new.names.caribbean<-cbind(caribbean.names.pca,pca.table)

theme_set(theme_light())
colrs <- c("#9e9ac8", "#807dba", "#6a51a3", "#4a1486")


F2 <- ggplot(new.names.caribbean, aes(PC1, PC2, color = caribbean.names.pca)) + 
  geom_point(size = 3) + stat_ellipse(level = 0.68)  
#why using 68%?? https://www.jstatsoft.org/article/view/v017i06

F2

# ADJUST ARGUMENT "labs" (x="LD1" TO THE PRODUCED VALUES)
F2a <- F2 + scale_color_manual(values = colrs,
                               aesthetics = c("colour", "fill"),
                               labels = c("Caribbean", "Mesoamerica", "Southeastern North America")) +
  geom_rug() +
  labs(x = "PC1 (63.3%)", y = "PC2 (8.94%)") +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        legend.position = "none")

#VISUALIZE
F2a

##########################################################################################


# FOR THE BOXPLOTS  
# create the boxplots to visualize characters variance

#upload data
morpho_final<-read.csv(file = "./morphodata_COMPLETE_3GROUPS.csv", header = T)

#make a colors vector
clrs = c("#9e9ac8", "#807dba", "#6a51a3")

#substitute "MASTOID_34" for the specific trait abbreviation contained in the Supplementary Table S3
MASTOID_34 <- ggplot(morpho_final, aes(x =Locality, y =MASTOID_34)) +
  geom_boxplot(outlier.alpha = 1, coef = 1, color = "black", fill = clrs, width = 0.5, alpha = 0.75) +
  geom_jitter(color = "black", size = 2, alpha = 0.4) +
  labs(x = "Locality", y = "insert.character.name") +
  theme_light() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12)) + scale_x_discrete(labels=c("Caribbean", "Mesoamerica", "S.E_North.America"))
MASTOID_34


##########################################################################################
#SUMMARY STATISTICS FOR CHARACTERS
#substitute "MAX.TOOTH_h" for the specific trait abbreviation contained in the Supplementary Table S3
morpho_final %>% group_by(Locality) %>% 
  summarise(Mean = mean(MAX.TOOTH_h), Median = median(MAX.TOOTH_h), Std = sd(MAX.TOOTH_h),
            Min = min(MAX.TOOTH_h), Max = max(MAX.TOOTH_h))

##########################################################################################
# END #