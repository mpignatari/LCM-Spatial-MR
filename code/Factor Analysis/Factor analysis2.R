df_reduced1<-read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/New Data/factor.csv")
anyNA(df_reduced1)
nrow(df_reduced1)

#install.packages("psych")
library(psych)
library(dplyr)      # data handling
library(missMDA)    # handy EM helpers (optional)
library(lavaan)
library(mice)

# Select the relevant variables
vars <- c("Clim_1","Clim_2","Clim_3","Clim_5","Clim_6","TR_1","TR_2","TR_3","TR_4","TR_5","Rel_1","Rel_2","Rel_3","Rel_4","Rel_5","Rel_6","Pol_1","Pol_2","Pol_3","Pol_4","Pol_5",    "Pol_6","Schw_1","Schw_2","Schw_3","Schw_4","Schw_5","Schw_6","Schw_7","Schw_8","Schw_9","Schw_10","Schw_11")
sapply(df_reduced1[, vars], is.ordered)


## 1.  Tell mice that *every* variable is ordinal
df_reduced1 <- df_reduced %>%                                     # keep other cols intact if you want
  mutate(across(all_of(vars), ~ ordered(.)))                  # 1 < 2 < … < 7

sapply(df_reduced1[, vars], is.ordered)

## 1.  One-shot EM imputation - impute EM on Nas.
df_reduced1 <- complete(mice(df_reduced1[ , vars], method = "polr", m = 1, maxit = 10))
## 2.  Polychoric correlation on the completed data
R_imp <- polychoric(df_reduced)$rho
tail(R_imp)






R_imp <- polychoric(df_reduced[vars])$rho
tail(R_imp, round)
round(R_imp,2)




## come back to numeric
df_reduced <- as.data.frame(lapply(df_reduced, function(x) as.numeric(as.character(x))))

R_imp <- polychoric(df_reduced)$rho
print(R_imp, row.names=FALSE) ######### proposal
anyNA(df_reduced)


# Run PCA using principal components (no rotation)
pca_result <- principal(df_reduced, nfactors = ncol(df_reduced), 
                        rotate = "none", scores = TRUE, method = "principal")

# Print summary
print(pca_result)

# Scree plot
scree(pca_result)

fa.parallel(df_reduced, fa = "pc", n.iter = 1000, quant = 0.5)


## Factor analysis without Nas.
fa_result <- fa(df_reduced, nfactors = 5, rotate = "promax", fm="ml")



# Step 1: Extract factors using Maximum Likelihood Estimation
fa_result <- fa(r = df_reduced,
                nfactors = 5,              # Choose number of factors (adjust as needed)
                fm = "ml",                  # MLE method
                rotate = "promax",         # Varimax rotation
                scores = TRUE)
print(fa_result, cut = 0.3)
# Extract fit statistics
fa_result$STATISTIC      # chi-square value
fa_result$dof            # degrees of freedom
fa_result$PVAL           # p-value

# Residual matrix
residuals <- fa_result$residual
round(residuals, 4)

# RMSR
fa_result$rms
# Partial correlation matrix
partial_corr <- residuals(fa_result, type = "partial")
round(partial_corr, 4)

# Assuming your dataset is called `frailty`
alpha(df_reduced[, c("Clim_1","Clim_2","Clim_3","Clim_5","Clim_6")],check.keys = TRUE)
alpha(df_reduced[, c("TR_1","TR_2","TR_3","TR_4","TR_5")])
alpha(df_reduced[, c("Rel_1","Rel_2","Rel_3","Rel_5","Rel_6","Schw_3")])
## 2nd round drop REl_6?
alpha(df_reduced[, c("Pol_2","Pol_3","Pol_4","Pol_5", "Pol_6")])
#Drop Pol_6
alpha(df_reduced[, c("Schw_2","Schw_4","Schw_5","Schw_6","Schw_7","Schw_8","Schw_9","Schw_10","Schw_11")])
# Sch2, Drop 10 and 11

alpha(df_reduced[, c("Clim_1","Clim_2","Clim_3","Clim_5","Clim_6")],check.keys = TRUE)
alpha(df_reduced[, c("TR_1","TR_2","TR_3","TR_4","TR_5")])
alpha(df_reduced[, c("Rel_1","Rel_2","Rel_3","Rel_5")])
## 2nd round drop REl_6? Then drop Scw_3
alpha(df_reduced[, c("Pol_2","Pol_3","Pol_4","Pol_5")])
#Drop Pol_6
alpha(df_reduced[, c("Schw_4","Schw_5","Schw_6","Schw_7","Schw_8","Schw_9")])
# Sch2, Drop 10 and 11

# Visualize item-total correlations and alpha if item dropped
plot(alpha(df_reduced[, c("TR_1","TR_2","TR_3","TR_4","TR_5")]))

#Vars without 
vars <- c("Clim_1","Clim_2","Clim_3","Clim_5","Clim_6","TR_1","TR_2","TR_3","TR_4","TR_5","Rel_1","Rel_2","Rel_3","Rel_5","Pol_1","Pol_2","Pol_3","Pol_4","Pol_5","Schw_4","Schw_5","Schw_6","Schw_7","Schw_8","Schw_9")




fa_result <- fa(r = df_reduced[ , vars],
                nfactors = 5,              # Choose number of factors (adjust as needed)
                fm = "ml",                  # MLE method
                rotate = "promax",         # Varimax rotation
                scores = TRUE)
print(fa_result, cut = 0.3)


# Step 2: View factor loadings
print(fa_result$loadings, cutoff = 0.3)

df_reduced1<-read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/New Data/factor.csv")

# Step 3: Calculate Factor Scores (like PROC SCORE)
# If you want to compute factor scores for individuals:
scores <- factor.scores(df_reduced[,vars], fa_result, method = "regression")$scores
nrow(scores)

# Optional: Combine with original data
df_scores2 <- cbind(df_reduced1, scores)
colnames(df_scores2)
head(df_scores2)
setwd("C:/Users/A02369659/Documents/12-2020/Marine Project/Chapter 3")
write.csv(df_scores2,"df_scores2.csv")

# ML 5 Equalitarians/social justice
# ML 2 Conservatives/risk averse
# ML 4 Anthropocenters
# ML 3 Ecocentrists
# ML 1 Climate concerners

#Now the CFA!

library(lavaan)

model <- '
  # Latent factors
  ML1 =~ Clim_1 + Clim_2 + Clim_3 + Clim_5 + Clim_6
  ML3 =~ TR_1 + TR_2 + TR_3 + TR_4 + TR_5
  ML4 =~ Rel_1 + Rel_2 + Rel_3 + Rel_5
  ML2 =~ Pol_5 + Schw_4 + Schw_5 + Schw_6 + Schw_7 + Schw_8 + Schw_9 
  ML5 =~ Pol_1 + Pol_2 + Pol_3 + Pol_4  

  # Optional: fix latent variances to 1
  ML1 ~~ 1*ML1
  ML2 ~~ 1*ML2
  ML3 ~~ 1*ML3
  ML4 ~~ 1*ML4
  ML5 ~~ 1*ML5
  
  # Optional: fix correlation between Body and Speed to 0 if matching SAS
  ML4 ~~ 0*ML2
'

fit <- cfa(model, data = df_reduced[,vars])

summary(fit, fit.measures = TRUE, standardized = TRUE)
modindices(fit, sort.=TRUE)[1:10, ]

model_refined <- '
  # Latent constructs
  ML1 =~ Clim_1 + Clim_2 + Clim_5 + Clim_6
  ML2 =~ Pol_5 + Schw_4 + Schw_5 + Schw_6 + Schw_7 + Schw_8 + Schw_9
  ML3 =~ TR_1 + TR_2 + TR_3 + TR_4 + TR_5
  ML4 =~ Rel_1 + Rel_2 + Rel_3 + Rel_5
  ML5 =~ Pol_1 + Pol_2 + Pol_3 + Pol_4

  # Fix non-significant covariance
  ML4 ~~ 0*ML2

  # New correlated residuals
  Clim_1 ~~ Clim_2
'
fit_refined <- cfa(model_refined, data = df_reduced[,vars])

summary(fit_refined, fit.measures = TRUE, standardized = TRUE)

anova(fit, fit_refined)

#install.packages("semTools")
library(semTools)
reliability(fit_refined)

# ML1: Omega and alpha very low → this factor might be misspecified or contains contradictory items. Check again whether Clim_3 and Clim_6 should stay.

#ML5: Borderline reliability → keep if theoretically important, but acknowledge as a limitation.

## FInal Factor

vars <- c("Clim_1","Clim_2","Clim_5","Clim_6","TR_1","TR_2","TR_3","TR_4","TR_5","Rel_1","Rel_2","Rel_3","Rel_5","Pol_1","Pol_2","Pol_3","Pol_4","Pol_5","Schw_4","Schw_5","Schw_6","Schw_7","Schw_8","Schw_9")

vars <- c("Clim_1","Clim_2","Clim_5","Clim_6","TR_1","TR_2","TR_3","TR_4","TR_5","Rel_1","Rel_2","Rel_3","Rel_5","Pol_2","Pol_3","Pol_4","Schw_4","Schw_5","Schw_9") ## variables used

fa_result <- fa(r = df_reduced[ , vars],
                nfactors = 5,              # Choose number of factors (adjust as needed)
                fm = "minres",                  # MLE method
                rotate = "promax",         # Varimax rotation
                scores = TRUE)
print(fa_result, cut = 0.3)


# Step 2: View factor loadings
print(fa_result$loadings, cutoff = 0.5)

df_reduced1<-read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/New Data/factor.csv")

# Step 3: Calculate Factor Scores (like PROC SCORE)
# If you want to compute factor scores for individuals:
scores <- factor.scores(df_reduced[,vars], fa_result, method = "Bartlett")$scores
nrow(scores)

# Optional: Combine with original data
df_scores2 <- cbind(df_reduced1, scores)
colnames(df_scores2)
head(df_scores2)
setwd("C:/Users/A02369659/Documents/12-2020/Marine Project/Chapter 3")
#write.csv(df_scores2,"df_scoresMR_50.csv")


set.seed(123)
n <- nrow(df_reduced[,vars])
idx <- sample(seq_len(n), size = n/2)

your_data_frame<-df_reduced[,vars]
train_data <- your_data_frame[idx, ]
test_data  <- your_data_frame[-idx, ]
anyNA(your_data_frame)
anyNA(train_data)

anyNA(train_data)
model <- '
  # Latent factors
  ML1 =~ Clim_1 + Clim_2 + Clim_5 + Clim_6
  ML3 =~ TR_1 + TR_2 + TR_3 + TR_4 + TR_5
  ML4 =~ Rel_1 + Rel_2 + Rel_3 + Rel_5
  ML2 =~ Schw_4 + Schw_5 + Schw_9 
  ML5 =~ Pol_2 + Pol_3 + Pol_4  

  # Optional: fix latent variances to 1
  ML1 ~~ 1*ML1
  ML2 ~~ 1*ML2
  ML3 ~~ 1*ML3
  ML4 ~~ 1*ML4
  ML5 ~~ 1*ML5
  
  # Optional: fix correlation between Body and Speed to 0 if matching SAS
  ML4 ~~ 0*ML2
  
  # New correlated residuals
  Clim_1 ~~ Clim_2
  '



# Fit CFA in training sample
fit_train <- cfa(model, data = train_data)
names(train_data)

# Re-fit in test sample (fixed structure)
fit_test <- cfa(model, data = test_data)

# Compare model fit or test invariance
anova(fit_train, fit_test)

your_data_frame$group <- ifelse(seq_len(nrow(your_data_frame)) %in% idx, "train", "test")

fit_configural <- cfa(model, data = your_data_frame, group = "group")
fit_metric     <- cfa(model, data = your_data_frame, group = "group", group.equal = "loadings")
fit_scalar     <- cfa(model, data = your_data_frame, group = "group", group.equal = c("loadings", "intercepts"))

anova(fit_configural, fit_metric, fit_scalar)

## EFA with train data and CFA with test data

fa_result <- fa(r = train_data,
                nfactors = 5,              # Choose number of factors (adjust as needed)
                fm = "minres",                  # MLE method
                rotate = "promax",         # Varimax rotation
                scores = TRUE)
print(fa_result, cut = 0.3)

fa_result <- fa(r = train_data,
                nfactors = 5,              # Choose number of factors (adjust as needed)
                fm = "minres",                  # MLE method
                rotate = "promax",         # Varimax rotation
                scores = TRUE)
print(fa_result, cut = 0.3)

fit_test <- cfa(model, data = test_data)
summary(fit_test, fit.measures=TRUE)

modindices(fit_test, sort = TRUE)

model <- '
  ML1 =~ Clim_1 + Clim_2 + Clim_5 + Clim_6
  ML3 =~ TR_1 + TR_2 + TR_3 + TR_4 + TR_5
  ML4 =~ Rel_1 + Rel_2 + Rel_3 + Rel_5
  ML2 =~ Schw_4 + Schw_5 + Schw_9 
  ML5 =~ Pol_2 + Pol_3 + Pol_4  

  ML1 ~~ 1*ML1
  ML2 ~~ 1*ML2
  ML3 ~~ 1*ML3
  ML4 ~~ 1*ML4
  ML5 ~~ 1*ML5

  ML4 ~~ 0*ML2

  # Residual correlations
  Clim_1 ~~ Clim_2
  Clim_1 ~~ Clim_5
  TR_1 ~~ TR_2
  TR_4 ~~ TR_5
  Rel_2 ~~ Rel_3
  Schw_4 ~~ Schw_9
'

fit_test <- cfa(model, data = test_data)
summary(fit_test, fit.measures=TRUE)
fitMeasures(fit_test)

df_reduced1<-read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/New Data/factor.csv")
anyNA(df_reduced)
anyNA(df_reduced1)

fa_result <- fa(r = df_reduced[, vars], nfactors = 5, rotate = "promax", scores = TRUE,    fm = "minres")
print(fa_result, cut = 0.3)
fa_result <- fa(r = df_reduced1[, vars], nfactors = 5, rotate = "promax", scores = TRUE,    fm = "minres")
print(fa_result, cut = 0.3)

nrow(df_reduced[, vars])



##GPT
#⚙️ Why Your Factor Analysis Results Changed
#1. mice() uses random draws when imputing NAs
#When you ran:

#df_reduced <- complete(mice(df_reduced[, vars], method = "polr", m = 5))each run generated slightly different imputed values for your missing data, because:The mice() algorithm draws random samples from predictive distributions.Those random draws affect the correlation matrix used in fa().Even tiny differences in correlations will slightly change loadings, communalities (h2, u2), and fit indices (TLI, RMSEA, etc.).This is why you cannot reproduce the exact same factor loadings unless you control the random seed.

#2. Your models differ only due to random imputation variationCompare the two results:

#  Metric	Old	New	Difference
#Clim_1 loading	0.93	0.91	minimal
#TLI	0.929	0.93	minimal
#RMSEA	0.055	0.054	minimal
#BIC	-321	-325	minimal