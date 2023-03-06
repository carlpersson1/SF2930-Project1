#Load the data
library("TH.data")
library("car")

data("bodyfat")
head(bodyfat)
plot(bodyfat)
summary(bodyfat)

#Create the full model
model <- lm(DEXfat ~., data = bodyfat)
summary(model)

#Residual analysis
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

res <- get_res(model)

plot(res$res)
plot(res$res_stud)
plot(res$res_press)
plot(res$res_rstud)

#Normal probability plot
qqnorm(res$res_stud)
qqline(res$res_stud)

#(s_res_stud <- sort(res_stud))
#c_prob <- ((1:n)-1/2)/n
#plot(s_res_stud, c_prob)

#Plot of Residuals against the Fitted Values
fit <- fitted(model)
plot(fit, res$res_stud)

#Plot of Residuals against the Regressor
plot(bodyfat$age, res$res_stud)
plot(bodyfat$waistcirc, res$res_stud)
plot(bodyfat$hipcirc, res$res_stud)
plot(bodyfat$elbowbreadth, res$res_stud)
plot(bodyfat$kneebreadth, res$res_stud)
plot(bodyfat$anthro3a, res$res_stud)
plot(bodyfat$anthro3b, res$res_stud)
plot(bodyfat$anthro3c, res$res_stud)
plot(bodyfat$anthro4, res$res_stud)

#Partial regression plots
avPlots(model)

#Partial residual plots
crPlots(model)

#PRESS statistic
(press = sum(res$res_press^2))

#Transformation of regressors
#bodyfat$waistcirc = 1 / bodyfat$waistcirc
#bodyfat$hipcirc = 1 / bodyfat$hipcirc
bodyfat$anthro3a = 1 / bodyfat$anthro3a
bodyfat$anthro3b = 1 / bodyfat$anthro3b
bodyfat$anthro3c = 1 / bodyfat$anthro3c
bodyfat$anthro4 = 1 / bodyfat$anthro4


#Variance-stabilization of the model
bodyfat$DEXfat <- sqrt(bodyfat$DEXfat)
model2 <- lm(DEXfat ~., data = bodyfat)

plot(bodyfat)

#Plot of Residuals against the Fitted Values
res2 <- get_res(model2)
fit2 <- fitted(model2)
plot(fit2, res2$res_stud)

#Normal probability plot
qqnorm(res2$res_stud)
qqline(res2$res_stud)

#Plot of Residuals against the Regressor
plot(bodyfat$age, res2$res_stud)
plot(bodyfat$waistcirc, res2$res_stud)
plot(bodyfat$hipcirc, res2$res_stud)
plot(bodyfat$elbowbreadth, res2$res_stud)
plot(bodyfat$kneebreadth, res2$res_stud)
plot(bodyfat$anthro3a, res2$res_stud)
plot(bodyfat$anthro3b, res2$res_stud)
plot(bodyfat$anthro3c, res2$res_stud)
plot(bodyfat$anthro4, res2$res_stud)

#Partial regression plots
avPlots(model2)


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


