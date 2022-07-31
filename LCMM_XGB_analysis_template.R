## R code for LCMM and XGB analysis
## NOTE: code removed for parsimony where repetitive 

library(lcmm)
library(reshape)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(mice)
library(caret)
library(xgboost)
library(pROC)
##
##
##
## LCMM ANALYSIS
## CODE ONLY INCLUDED FOR SAMPLE A (IDENTICAL ANALYSIS FOR SAMPLE B)
##
##
##

# ONLY USING PCL SCORES FOR LCMM SO PULLING OUT ONLY THOSE COLUMNS FROM BASELINE DATA
set.seed(42)
sampleA_PCL_timepoints<-dplyr::select(sampleA_goodvars_goodsubs, ID, PCL_total_T1, PCL_total_T2, PCL_total_T3)

# RENAME COLUMNS TO SOMETHING READABLE
colnames(sampleA_PCL_timepoints)<-c("ID", "T1", "T2", "T3")

# RESHAPE DATA TO LONG FORMAT (BY PCL AT 3 TIMEPOINTS)
sampleA_reshape_PCL_timepoints<- melt(data=sampleA_PCL_timepoints, id=c("ID"), measures=c(colnames(sampleA_PCL_timepoints[,2:ncol(sampleA_PCL_timepoints)])))

# RENAME COLUMNS TO SOMETHING READABLE
colnames(sampleA_reshape_PCL_timepoints)<-c("ID", "TIME", "PCL")
sampleA_reshape_PCL_timepoints$TIME<-as.numeric(sampleA_reshape_PCL_timepoints$TIME)

# CALCULATE LCMM MODELS FOR 1-6 CLASS SOLUTIONS
d1<-lcmm(PCL~TIME, random=~TIME, subject="ID", ng=1, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear")

# GRID SEARCH EMPLOYED FOR ng=2-6, sampleATING FROM INITIAL VALUES FROM ng=1 SOLUTION ABOVE
# 50 REPETITIONS (DEPARTURES FROM RANDOM INTIAL VALUES)
# 100 ITERATIONS TO OPTIMIZE ALGORITHM
d2<-gridsearch(minit=d1, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=2, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d3<-gridsearch(minit=d1, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=3, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d4<-gridsearch(minit=d1, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=4, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d5<-gridsearch(minit=d1, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=5, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d6<-gridsearch(minit=d1, maxiter=100, rep=50, lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=6, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))

# COMPARE MODEL FIT PARAMETERS
summarytable(d1, d2, d3, d4, d5, d6, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

#
# REPEAT ABOVE ANALYSIS WITH INCLUSION OF A QUADRATIC TERM
#
# CALCULATE QUADRATIC TIME SLOPE
sampleA_reshape_PCL_timepoints$TIME2 <- as.numeric(sampleA_reshape_PCL_timepoints$TIME)^2

# FIT LCMM MODEL
d1_2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", ng=1, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear")
# AGAIN, USE INITIAL VALUES FROM ng=1 SOLUTION FOR GRID SEARCH IMPLEMENTED FOR ng=2-6 SOLUTIONS
d2_2<-gridsearch(minit=d1_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=2, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d3_2<-gridsearch(minit=d1_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=3, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d4_2<-gridsearch(minit=d1_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=4, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d5_2<-gridsearch(minit=d1_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=5, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))
d6_2<-gridsearch(minit=d1_2, maxiter=100, rep=50, lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=6, idiag=TRUE, data=sampleA_reshape_PCL_timepoints, link="linear"))

# COMPARE MODEL FIT PARAMETERS
summarytable(d1_2, d2_2, d3_2, d4_2, d5_2, d6_2, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

##
##
##
## SPLIT HALF VALIDATION OF LCMM 
##
##
##
# ASSESS SAME MODEL FITS (AS ABOVE) USING SPLIT-HALF CROSS VALIDATION (AS IN GALATZER-LEVY ET AL 2017 PAPER)
# SPLIT HALF CROSS VALIDATION

# TO APPLY PARAMETERS FROM PREVIOUS MODEL SOLUTION USE B="MODEL$BEST"
# SEE: https://github.com/CecileProust-Lima/lcmm/issues/59

# RANDOMLY SPLIT ORIGINAL DATASET INTO 2 EQUAL HALVES
set.seed(10)
half1rownums<-sample(1:nrow(sampleA_PCL_timepoints), nrow(sampleA_PCL_timepoints)/2, replace=FALSE)

# FIRST HALF
half1_temp<-sampleA_PCL_timepoints[half1rownums,]
# ORGANIZE half1 TO BE ORDERED BY ID
half1<-half1_temp[order(half1_temp$ID),]
# SECOND HALF
half2<-sampleA_PCL_timepoints[-half1rownums,]


## CALCULATE MODELS FOR FIRST HALF OF DATA USING ESTIMATES FROM FULL SAMPLE ABOVE
reshape_half1<- melt(data=half1, id=c("ID"), measures=c(colnames(half1[,2:ncol(half1)])))
colnames(reshape_half1)<-c("ID", "TIME", "PCL")
reshape_half1$TIME <- as.numeric(reshape_half1$TIME)
reshape_half1$TIME2 <- as.numeric(reshape_half1$TIME)^2

# LINEAR TERM ONLY
d1_half1<-lcmm(PCL~TIME, random=~TIME, subject="ID", ng=1, idiag=TRUE, data=reshape_half1, link="linear")
d2_half1<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=2, idiag=TRUE, data=reshape_half1, link="linear")
d3_half1<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=3, idiag=TRUE, data=reshape_half1, link="linear")
d4_half1<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=4, idiag=TRUE, data=reshape_half1, link="linear")
d5_half1<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=5, idiag=TRUE, data=reshape_half1, link="linear")
d6_half1<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=6, idiag=TRUE, data=reshape_half1, link="linear")

# LINEAR AND QUADRATIC TERMS
d1_2_half1<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", ng=1, idiag=TRUE, data=reshape_half1, link="linear", B=d1_2$best)
d2_2_half1<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=2, idiag=TRUE, data=reshape_half1, link="linear", B=d2_2$best)
d3_2_half1<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=3, idiag=TRUE, data=reshape_half1, link="linear", B=d3_2$best)
d4_2_half1<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=4, idiag=TRUE, data=reshape_half1, link="linear", B=d4_2$best)
d5_2_half1<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=5, idiag=TRUE, data=reshape_half1, link="linear", B=d5_2$best)
d6_2_half1<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=6, idiag=TRUE, data=reshape_half1, link="linear", B=d6_2$best)


## CALCULATE MODELS FOR SECOND HALF OF DATA
reshape_half2<- melt(data=half2, id=c("ID"), measures=c(colnames(half2[,2:ncol(half2)])))
colnames(reshape_half2)<-c("ID", "TIME", "PCL")
reshape_half2$TIME <- as.numeric(reshape_half2$TIME)
reshape_half2$TIME2 <- as.numeric(reshape_half2$TIME)^2

# LINEAR TERM ONLY
d1_half2<-lcmm(PCL~TIME, random=~TIME, subject="ID", ng=1, idiag=TRUE, data=reshape_half2, link="linear")
d2_half2<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=2, idiag=TRUE, data=reshape_half2, link="linear")
d3_half2<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=3, idiag=TRUE, data=reshape_half2, link="linear")
d4_half2<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=4, idiag=TRUE, data=reshape_half2, link="linear")
d5_half2<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=5, idiag=TRUE, data=reshape_half2, link="linear")
d6_half2<-lcmm(PCL~TIME, random=~TIME, subject="ID", mixture=~TIME, ng=6, idiag=TRUE, data=reshape_half2, link="linear")

# LINEAR AND QUADRATIC TERMS
d1_2_half2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", ng=1, idiag=TRUE, data=reshape_half2, link="linear", B=d1_2$best)
d2_2_half2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=2, idiag=TRUE, data=reshape_half2, link="linear", B=d2_2$best)
d3_2_half2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=3, idiag=TRUE, data=reshape_half2, link="linear", B=d3_2$best)
d4_2_half2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=4, idiag=TRUE, data=reshape_half2, link="linear", B=d4_2$best)
d5_2_half2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=5, idiag=TRUE, data=reshape_half2, link="linear", B=d5_2$best)
d6_2_half2<-lcmm(PCL~TIME + TIME2, random=~TIME + TIME2, subject="ID", mixture=~TIME + TIME2, ng=6, idiag=TRUE, data=reshape_half2, link="linear", B=d6_2$best)


## MODEL FIT COMPARISON

# LINEAR TERM ONLY
summarytable(d1_half1, d2_half1, d3_half1, d4_half1, d5_half1, d6_half1, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

summarytable(d1_half2, d2_half2, d3_half2, d4_half2, d5_half2, d6_half2, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

# LINEAR AND QUADRATIC TERM
summarytable(d1_2_half1, d2_2_half1, d3_2_half1, d4_2_half1, d5_2_half1, d6_2_half1, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))
summarytable(d1_2_half2, d2_2_half2, d3_2_half2, d4_2_half2, d5_2_half2, d6_2_half2, which=c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "%class"))

# PLOT ALL CLASS SOLUTIONS
# FOR PARSIMONY CODE BELOW ONLY DEPICTS PLOTTING K=1 SOLUTION
# PROCESS WAS REPEATED FOR K=2-6 SOLUTIONS
# DESIGNATE COLOR PALATTE FOR GGPLOT
colpal="Oranges"

# ENSURE ID ISN'T TREATED AS NUMERIC
sampleA_reshape_PCL_timepoints$ID <- as.character(sampleA_reshape_PCL_timepoints$ID)

# FOR EACH MODEL GRAB THE POSTERIOR PROBABILITIES
people1 <- as.data.frame(d1_2$pprob[,1:2])

# PUT CLASS ASSIGNMENTS FROM POSTERIOR PROBABILITIES INTO OUR DATAFRAME WITH THE PCL SCORES SO WE CAN PLOT CLASSES
sampleA_reshape_PCL_timepoints$group1 <- factor(people1$class[sapply(sampleA_reshape_PCL_timepoints$ID, function(x) which(people1$ID==x))])

# PLOT LATENT CLASS ASSIGNMENTS WITH INDIVIDUAL (group=ID) AND MEAN GROUP LINES (colour=group#)
p1<-ggplot(sampleA_reshape_PCL_timepoints, aes(TIME, PCL, group=ID, colour=group1)) + geom_line(data=na.omit(sampleA_reshape_PCL_timepoints)) + geom_smooth(aes(group=group1), method="loess", size=2, se=T) + labs(title="", x="Time",y="PCL Total Severity",colour="Class") + theme_classic() + theme(plot.title = element_text(hjust = 0.5, size=14,face="bold"), legend.title=element_text(size=10, face="bold"), legend.text=element_text(size=10), axis.title = element_text(size=12, face="bold"), axis.text=element_text(size=12)) + scale_color_brewer(palette = colpal) +
  geom_hline(yintercept=32, linetype="dashed", size=0.5) + scale_x_continuous(breaks=c(1,2,3))


# FINAL CLASS SOLUTION SELECTED
# REPORT POSTERIOR PROBABILITIES
postprob(d4_2)

# REPORT PREDICTED PCL-5 VALUES FROM FINAL LCMM MODEL
# CREATE DUMMY DATAFRAME
x=c(1,1)
y=c(2,4)
z=c(3,9)
df_temp<-as.data.frame(rbind(x,y,z))
colnames(df_temp)<-c("TIME", "TIME2")
# PREDICT PCL-5 FROM DUMMY DATAFRAME
predictY(d4_2, df_temp)

# ATTACH CLASSES TO BASELINE DATAFRAME
sampleA_goodvars_goodsubs_classes<-cbind(sampleA_goodvars_goodsubs, people4$class)
colnames(sampleA_goodvars_goodsubs_classes)[colnames(sampleA_goodvars_goodsubs_classes) == 'people4$class'] <- 'Latent_Class'
sampleA_goodvars_goodsubs=within(sampleA_goodvars_goodsubs, rm("useablesubs"))

# REPORT NUMBER OF SUBS PER CLASS
table(sampleA_goodvars_goodsubs_classes$Latent_Class)

##
##
##
##
##
##
## TAKE LCMM SOLUTIONS TO XGB ANALYSIS
## 1) SPLIT TRAINING TEST DATA SETS
## 2) IMPUTE MISSING DATA
## 3) DUMMY CODE CATEGORICAL DATA
## 4) XGB WITH RECURSIVE FEATURE ELIMINATION 
## 5) XGB WITH ALL VARIABLES
##
##

# EVALUATE MISSINGNESS WITH VIM
aggr(sampleA_goodvars_goodsubs_classes, numbers=TRUE, prop=FALSE, labels=names(trainSet_sampleA_M1), cex.axis=.55, gap=1, combined=TRUE)
aggr(sampleB_goodvars_goodsubs_classes, numbers=TRUE, prop=FALSE, labels=names(trainSet_sampleA_M1), cex.axis=.55, gap=1, combined=TRUE)

# FOR REPLICABILITY OF RESULTS, SET SEED
set.seed(11)

# SPLIT DATA INTO TRAINING AND TEST SETS (75/25 SPLIT)
splits<-createDataPartition(sampleA_goodvars_goodsubs_classes$Latent_Class_nonRemitVSall, p=0.75, list=FALSE)

# APPLY SPLITS TO FIRST MODEL DATASET
trainSet_sampleA_M1<-sampleA_goodvars_goodsubs_classes[splits,]
testSet_sampleA_M1<-sampleA_goodvars_goodsubs_classes[-splits,]

# MAKE SURE CATEGORICAL VARIABLES ARE FACTORS
# IMPUTE MISSING CASES USING MULTIPLE IMPUTATION
# M = 20, NUMBER OF MULTIPLE IMPUTATIONS
# MAXIT = 20 MEANS 20 ITERATIONS PER IMPUTATION
# PMM = PREDICTIVE MEAN MATCHING
imp_sampleA_train<-mice::mice(trainSet_sampleA_M1, m=20, maxit=20, method="pmm")
imp_sampleA_test<-mice(testSet_sampleA_M1, m=20, maxit=20, method="pmm")
imp_sampleB<-mice(sampleB_goodvars_goodsubs_classes, m=20, maxit=20, method="pmm")

# VERIFY CONVERGENCE, LOOKING FOR NO REAL PATTERN IN PLOTS ACROSS IMPUTATIONS
plot(imp_sampleA_train)
plot(imp_sampleA_test)
plot(imp_sampleB)

# INSERT IMPUTED CASES TO DATASET
# RANDOMLY SELECT IMPUTATION TO CARRY FORWARD
sampleA_train_imputed<-complete(imp_sampleA_train,round(runif(1, min=1, max=20), 0))
sampleA_test_imputed<-complete(imp_sampleA_test,round(runif(1, min=1, max=20), 0))
sampleB_imputed<-complete(imp_sampleB,round(runif(1, min=1, max=20), 0))

# REMOVE EXTRANEOUS ID VARIABLES
sampleA_train_imputed=within(sampleA_train_imputed, rm("X", "ID"))
sampleA_test_imputed=within(sampleA_test_imputed, rm("X", "ID"))
sampleB_imputed=within(sampleB_imputed, rm("ID"))

# DUMMY CODE CATEGORICAL VARIABLES
# CREATED LOOP TO DO THIS ACROSS TRAINING AND TEST DATASETS
MODEL_NAMES=c("sampleA_train_imputed", "sampleA_test_imputed", "sampleB_imputed")

# dummyVars ONLY IDENTIFIES WHICH VARIABLES SHOULD BE DUMMY CODED
# fullRank=T WILL CREATE ONLY (N-1) COLUMNS FOR CATEGORY WITH N LEVELS
for (MODEL in 1:length(MODEL_NAMES)) {
  dmy<-dummyVars("~.", data=eval(parse(text= MODEL_NAMES[MODEL])), fullRank=T)
  temp_df<-data.frame(predict(dmy, newdata=eval(parse(text= MODEL_NAMES[MODEL]))))
  # CONVERT OUR OUTCOME VARIABLE BACK TO FACTOR AHEAD OF SVM
  temp_df$Latent_Class<-factor(temp_df$Latent_Class)
  # SAVE OUT DUMMY CODED DF WITH CORRESPONDING NAME FROM ORIGINAL DF
  assign(strcat(c(MODEL_NAMES[MODEL], "_wDummies")), temp_df)
}

# NOW THAT OUR MODEL VARIABLE SETS HAVE BEEN CLEANED AND PREPPED, WE BEGIN OUR XGB CLASSIFICATION

##
##
## TRAIN XGB ON sampleA, TEST on sampleB
##
##

# REMOVE THE PCL VARIABLES FROM EACH TIME POINT AND THE LATENT CLASS VARIABLE FROM THE PREDICTOR DF IN EACH OF THE TRAINING AND TEST DATASETS

sampleA_train_imputed_wDummies=within(sampleA_train_imputed_wDummies, rm("PCL_total_T1", "PCL_total_T2", "PCL_total_T3", "Latent_Class"))
sampleA_test_imputed_wDummies=within(sampleA_test_imputed_wDummies, rm("PCL_total_T1", "PCL_total_T2", "PCL_total_T3", "Latent_Class"))
sampleB_imputed_wDummies=within(sampleB_imputed_wDummies, rm("PCL_total_T1", "PCL_total_T2", "PCL_total_T3", "Latent_Class"))

####
####
####
#### XTREME GRADIENT BOOSTING 
#### WITH RECURSIVE FEATURE ELIMINATION (RFE)
####
####
####

# FOR REPLICABILITY OF RESULTS, AND DUE TO RANDOM SAMPLING NATURE OF CV PROCEDURES, SET SEED AHEAD OF RFE
set.seed(1)

# SET UP THE CONTROL PARAMETERS FOR THE RFE FUNCTION BEFORE EXECUTION
ctrl_xgb<-rfeControl(functions=caretFuncs, 
                     method="repeatedcv", 
                     number=5, repeats=10,  # FOR CV: 5 X 10 FOLD. 10 FOLDS, REPEATED 5 TIMES.
                     index = createMultiFolds(sampleA_train_imputed_wDummies[, "Latent_Class_nonRemitVSall"], times = 5), # POINT TO THE TARGET VARIABLE, CREATING FOLDS BASED ON DISTRIBUTION OF LATENT CLASS VARIABLE
                     verbose = TRUE, 
                     allowParallel = TRUE)

# EXECUTION OF THE RFE METHODS USING THE ABOVE CONTROL PARAMETERS
# THIS TAKES SOME TIME (10+ MIN) DUE TO THE MANY FOLDS/ITERATIONS
# search = "random" INSTATES RANDOM HYPERPARAMETER GRID SEARCH
# PREPROCESS DATA BY CENTERING, SCALING, AND REMOVING VARS WITH NON-ZERO VARIANCES
xgb_rfe_M1<-rfe(Latent_Class_nonRemitVSall~., 
                data=sampleA_train_imputed_wDummies, 
                method="xgbLinear", 
                rfeControl=ctrl_xgb, 
                preProc=c("center", "scale", "nzv"))

# SUMMARY OF RFE SOLUTION
xgb_rfe_M1

# PULL BEST PREDICTOR SET FROM xgb-RFE
preds_M1<- predictors(xgb_rfe_M1)
# WHAT WERE THE BEST RFE PREDICTORS?
preds_M1

# TEST AND CALCULATE CLASSIFICATION ACCURACY ON TEST DATA SET
# PREDICT TRAJECTORIES USING RFE VARS
xgb_M1_bestvars_internal=predict(xgb_rfe_M1, newdata=sampleA_test_imputed_wDummies[,-35])
xgb_M1_bestvars_external=predict(xgb_rfe_M1, newdata=sampleB_imputed_wDummies[,-35])

# EVALUATE ACCURACY OF PREDICTIONS
confusionMatrix(xgb_M1_bestvars_internal, sampleA_test_imputed_wDummies$Latent_Class_nonRemitVSall, mode="everything")
confusionMatrix(xgb_M1_bestvars_external, sampleB_imputed_wDummies$Latent_Class_nonRemitVSall, mode="everything")

# CALCULATE AUC CONFIDENCE INTERVAL
# FIRST CALCULATE FOR INTERNAL VALIDATION
roc_rfe_sampleA_test<-roc(as.numeric(sampleA_test_imputed_wDummies$Latent_Class_nonRemitVSall), as.numeric(xgb_M1_bestvars_internal))
# WHAT IS AUC?
pROC::auc(roc_rfe_sampleA_test)
# WHAT IS AUC 95% CONFIDENCE INTERVAL? USE BOOTSTRAP SAMPLING
pROC::ci.auc(roc_rfe_sampleA_test, method="bootstrap")

# REPEAT AUC CALCULATION FOR EXTERNAL VALIDATION SET
roc_rfe_sampleB_test<-roc(as.numeric(sampleB_imputed_wDummies$Latent_Class_nonRemitVSall), as.numeric(xgb_M1_bestvars_external))
pROC::auc(roc_rfe_sampleB_test)
pROC::ci.auc(roc_rfe_sampleB_test, method="bootstrap")

####
####
####
#### REPEAT XGB ANALYSIS FOR FULL VARIABLE SET
####
####
####

# SET SEED FOR REPLICABILITY OF CV PROCESS
set.seed(110)

# AGAIN USING 5 X 10 FOLD CV METHOD
# INSTANTIATE RANDOM GRID SEARCH PROCESS FOR HYPERPARAMETER TUNING
fitControl<-trainControl(method = "repeatedcv", 
                         number = 5, 
                         repeats = 10, 
                         search = "random")

# RUN XGB ANALYSIS
# USE SAME PREPROCESSING STEPS AS IN RFE ANALYSIS
xgb_M1<-train(Latent_Class_nonRemitVSall ~ . , 
              data=sampleA_train_imputed_wDummies, 
              method="xgbLinear", 
              trControl=fitControl, 
              preProc = c("center", "scale", "nzv"))

# TEST AND CALCULATE CLASSIFICATION ACCURACY ON TEST DATA SET
# CALCULATE TRAJECTORY PREDICTIONS USING ALL VARS
xgb_M1_internal=predict(xgb_M1, newdata=sampleA_test_imputed_wDummies[,-35])
xgb_M1_external=predict(xgb_M1, newdata=sampleB_imputed_wDummies[,-35])

# EVALUATE PREDICTION PERFORMANCE ON INTERNAL AND EXTERNAL VALIDATION SETS
confusionMatrix(xgb_M1_internal, sampleA_test_imputed_wDummies$Latent_Class_nonRemitVSall, mode="everything")
confusionMatrix(xgb_M1_external, sampleB_imputed_wDummies$Latent_Class_nonRemitVSall, mode="everything")

# CALCULATE AUC CONFIDENCE INTERVAL
# FIRST CALCULATE FOR INTERNAL VALIDATION
roc_sampleA_test<-roc(as.numeric(sampleA_test_imputed_wDummies$Latent_Class_nonRemitVSall), as.numeric(xgb_M1_internal))
pROC::auc(roc_sampleA_test)
pROC::ci.auc(roc_sampleA_test, method="bootstrap")

# REPEAT FOR EXTERNAL VALIDATION
roc_sampleB_test<-roc(as.numeric(sampleB_imputed_wDummies$Latent_Class_nonRemitVSall), as.numeric(xgb_M1_external))
pROC::auc(roc_sampleB_test)
pROC::ci.auc(roc_sampleB_test, method="bootstrap")


# PLOT SHAPLEY VALUES FROM FULL VAR SET
# SHAP SUMMARY PLOT REQUIRES DF IS IN MATRIX FORM
test_mat_m1<-as.matrix(sampleB_imputed_wDummies[,-35])

shap_plot_M1<-
  xgb.ggplot.shap.summary(data=test_mat_m1, model = xgb_M1$finalModel, top_n=10) + 
  theme_bw() + 
  ggtitle("Model 1: Non-remitting vs. all") + 
  labs(y="SHAP Value", x="", colour="Feature Value")





## THE ABOVE XGB WITH RFE AND XGB WITH FULL VARIABLE SET WAS REPEATED FOR MODELS 2 AND 3 (I.E. NON-REMITTING VS. RESILIENT, AND NON-REMITTING VS. REMITTING)
## CODE EDITED OUT FOR PARSIMONY

