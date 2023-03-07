#Load the data
library("TH.data")
library("car")
library("leaps")
data("bodyfat")
plot(bodyfat)
rownames(bodyfat) <- 1:nrow(bodyfat)

create_model <- function(dataset){
  model <- lm(DEXfat ~., data = dataset)
  return(model)
}

get_res <- function(model) {
  n = nrow(bodyfat)
  p = ncol(bodyfat)-1
  h_ii = lm.influence(model)$hat
  res = residuals(model)
  MS_res = sum(res^2)/(n-p) 
  res_std   = res / sqrt(MS_res)
  res_stud  = res / sqrt(MS_res*(1-h_ii))
  res_press = res / (1 - h_ii)
  S2 = ((n-p)*MS_res-res^2/(1-h_ii))/(n-p-1)
  res_rstud <- res/sqrt(S2*(1-h_ii))
  return(list("res"=res,
              "res_std"=res_std,
              "res_stud"=res_stud,
              "res_rstud"=res_rstud,
              "res_press"=res_press))}

plot_normal <- function(model){
  res = get_res(model)
  qqnorm(res$res_stud)
  qqline(res$res_stud)
}

plot_res_vs_fitted <- function(model){
  fit <- fitted(model)
  res = get_res(model) 
  plot(fit, res$res_stud)
}

transform_boxcox <- function(dataset){
  bf = dataset
  bc <- boxcox(DEXfat ~., data=bf)
  (lambda <- bc$x[which.max(bc$y)])
  bf$DEXfat = bf$DEXfat^lambda
  return(bf)
}

transform_boxTidwell <- function(dataset){
  bf = dataset
  DEXfat <- bf$DEXfat
  dt <- bf[,1,drop=FALSE]
  dt <- cbind(dt, DEXfat)
  
  #Transformation of regressors
  (bt <- boxTidwell(DEXfat ~., data=dt)) #Vary to get lambda for all regressors. 
  
  bf$age = bf$age^-1.6591
  bf$waistcirc = bf$waistcirc^0.60716
  bf$hipcirc = bf$hipcirc^-0.77347
  bf$elbowbreadth = bf$elbowbreadth^1.8084
  bf$kneebreadth = bf$kneebreadth^2.908
  bf$anthro3a = bf$anthro3a^4.7541
  bf$anthro3b = bf$anthro3b^5.6786
  bf$anthro3c = bf$anthro3c^3.3807
  bf$anthro4 = bf$anthro4^5.2287
  return(bf)
}

cross_validate <- function(dataset, k=10){
  bf = dataset
  MSE <- numeric(511)
  ADJR2 <- numeric(511)
  R2 <- numeric(511)
  n <- round(nrow(bf)/k, digits = 0)
  for (i in 1:k){
    bf <- bf[sample(1:nrow(bf)),] #Shuffle the data
    train_ind <- sample(seq_len(nrow(bf)), size = k)
    from <- n*(i-1)
    to <- min(n*i, nrow(bf))
    test_set <- bf[from:to,]
    train_set <- bf[-(from:to),]
    y_train <- train_set$DEXfat
    y_test <- test_set$DEXfat
    subs <- regsubsets(DEXfat ~., data=train_set, method = "exhaustive", nbest = 126, nvmax=9, really.big = TRUE)
    (test_set <- test_set[,-2])
    (train_set <- train_set[,-2])
    for (idx in 1:511){
      (model <- summary(subs)$which[idx,])
      (model <- model[-1])
      (train_set_model <- train_set[,model,drop=FALSE])
      (test_set_model <- test_set[,model,drop=FALSE])
      (train_set_model <- cbind(train_set_model, y_train))
      lm_model = lm(y_train ~., data=train_set_model)
      
      (pred <- predict(lm_model,newdata = test_set_model))
      
      (SS_res <- sum((y_test-pred)^2))
      (SS_R <- sum((pred-mean(y_test))^2))
      (SS_T <- SS_R + SS_res)
      (MSE[idx] <- MSE[idx] + mean((y_test-pred)^2))
      (summary(lm_model))
      (R2[idx]<- R2[idx] + 1-SS_res/SS_T)
      (ADJR2[idx] <- ADJR2[idx] + 1-(SS_res/(nrow(train_set)-ncol(test_set)))/(SS_T/(nrow((train_set)-1))))
    }
  }
  (MSE <- MSE/k)
  (ADJR2 <- ADJR2/k)
  (R2 <- R2/k)
  
  subs <- regsubsets(DEXfat ~., data=bf, method = "exhaustive", nbest = 126, nvmax=9, really.big = TRUE)
  (best_model_idx_adjr2 <- which.max(ADJR2))
  (best_model_idx_mse <- which.min(MSE))
  
  (best_model_cols_adjr2 <- summary(subs)$which[best_model_idx_adjr2,])
  (best_model_cols_mse <- summary(subs)$which[best_model_idx_mse,])
  
  (best_model_cols_adjr2 <- best_model_cols_adjr2[-1])
  (best_model_cols_mse <- best_model_cols_mse[-1])
  
  (DEXfat <- bf$DEXfat)
  bf <- bf[-2]
  
  (dataset_adjr2 <- bf[,best_model_cols_adjr2,drop=FALSE])
  (dataset_mse <- bf[,best_model_cols_mse,drop=FALSE])
  
  (dataset_adjr2 <- cbind(dataset_adjr2,DEXfat))
  (dataset_mse <- cbind(dataset_mse,DEXfat))
  
  model_adjr2 <- create_model(dataset_adjr2)
  model_mse <- create_model(dataset_mse)
  
  return(list("adjr2" = model_adjr2, "mse" = model_mse))
}

#Initial model
full_model <- create_model(bodyfat)
(mse_full_model<-mean((summary(full_model)$residuals)^2))
summary(full_model)
plot_normal(full_model) #Normal probability plot
plot_res_vs_fitted(full_model) #Plot of Residuals against the Fitted Values
avPlots(full_model) #Partial regression plots

#Box Cox
data_boxcox <- transform_boxcox(bodyfat)
boxcox_model <- create_model(data_boxcox)
summary(boxcox_model)
plot_normal(boxcox_model) #Normal probability plot
plot_res_vs_fitted(boxcox_model) #Plot of Residuals against the Fitted Values
avPlots(boxcox_model) #Partial regression plots

#Box Tidwell
data_boxtidwell <- transform_boxTidwell(bodyfat)
boxtidwell_model <- create_model(data_boxtidwell)
summary(boxtidwell_model)
plot_normal(boxtidwell_model) #Normal probability plot
plot_res_vs_fitted(boxtidwell_model) #Plot of Residuals against the Fitted Values
avPlots(boxtidwell_model) #Partial regression plots

#Cross validation
boxcox_models <- cross_validate(data_boxcox)
summary(boxcox_models$adjr2)
summary(boxcox_models$mse)


boxtidwell_models <- cross_validate(data_boxtidwell)
summary(boxtidwell_models$adjr2)
summary(boxtidwell_models$mse)


##We should move this \/

# Prediction variability power
(press = sum(res2$res_press^2))
SS_T = sum((bodyfat$DEXfat - mean(bodyfat$DEXfat))^2)
R2_pred = 1 - press / SS_T

# Examination of residuals - Seemingly ~3 outliers 73, 87, 94
plot(res2$res_stud)
plot(res2$res_press)
plot(res2$res_rstud)
View(res2$res_stud)

# Covratio with cutoff lines - A fair amount of seemingly influential points for precision
# Influential points - bad for precision 73 (0.31), 87 (0.12), 92 (0.56) and 94 (0.40)
covrat = covratio(model2)
plot(covrat)
abline(h=1 - 3*p/n)
abline(h=1 + 3*p/n)
View(covrat)

# Cooks distance - No major influential points displacing the model parameters
# The most influential points 71 (0.12), 73 (0.10), 87 (0.24) and 94 (0.17)
cooksD = cooks.distance(model2)
plot(cooksD)
cutoff = qf(0.5, p, n-p)
abline(h=cutoff)
View(cooksD)

# Leverage points - About 4 points far enough away to be considered leverage points
# Leverage points include 71 (0.31), 81 (0.26), 112 (0.26) and 113 (0.29)
H_diag = lm.influence(model2)$hat
plot(H_diag)
abline(h=2*p/n)
View(H_diag)

# Checking for normality - 73, 87 and 94 seems to violate the normality condition
qqinfo = qqPlot(res2$res_stud)
View(qqinfo)


