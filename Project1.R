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
n = nrow(bodyfat)
p = ncol(bodyfat)-1
res <- residuals(model)
SS_res = sum(res^2)
MS_res = SS_res/(n-p)
h_ii = lm.influence(model)$hat
S2 = ((n-p)*MS_res-res^2/(1-h_ii))/(n-p-1)

res_std <- res/sqrt(MS_res)
res_stud <- res/sqrt(MS_res*(1-h_ii))
res_press <- res/(1-h_ii)
res_rstud <- res/sqrt(S2*(1-h_ii))

plot(res)
plot(res_stud)
plot(res_press)
plot(res_rstud)

#Normal probability plot
qqnorm(res_stud)
qqline(res_stud)
#(s_res_stud <- sort(res_stud))
#c_prob <- ((1:n)-1/2)/n
#plot(s_res_stud, c_prob)

#Plot of Residuals against the Fitted Values
fit <- fitted(model)
plot(fit, res_stud)

#Plot of Residuals against the Regressor
plot(bodyfat$age, res_stud)
plot(bodyfat$waistcirc, res_stud)
plot(bodyfat$hipcirc, res_stud)
plot(bodyfat$elbowbreadth, res_stud)
plot(bodyfat$kneebreadth, res_stud)
plot(bodyfat$anthro3a, res_stud)
plot(bodyfat$anthro3b, res_stud)
plot(bodyfat$anthro3c, res_stud)
plot(bodyfat$anthro4, res_stud)

#Partial regression plots
avPlots(model)

#Partial residual plots
crPlots(model)


#Diagnostics and handling of outliers, leverage and influential observations 

#PRESS statistic (hör nog inte riktigt hemma här men kan vara användbar senare)
(press = sum(res_press^2))



