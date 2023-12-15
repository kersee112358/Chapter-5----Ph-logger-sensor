# data sets ---------------------------------------------------------------
library(readxl)
library(xlsx)
library(dplyr)
library(stringr)
library(MethodCompare)
library(tidyr)

setwd("C:/Users/kersee/Dropbox/Thesis/Chapter_5/R-Script")
# the following script comes from [28]
source("MAMMA.R")

# The data frame will be named as df and the format of the data will be long.
# In the current application we have the following columncs
# 1. the repeated measurements
# 2. The temperature of the solution
# 3. The subject ID
# 4. The quantity of the solution (in ml)
# 5. The method (open-source or reference)
# 6. The type of the fruit
# 7. The pH measurement

#The following line names the columns properly according to the rest of the script
colnames(df) <- c("rep", "Temp", "subject", "quantity", "method","turn","type", "meas")

df$subject<-factor(df$subject)
df$method<-factor(df$method)
df$type<-factor(df$type)

df$meas<-round(df$meas,2)

#converting from long to wide for exploratory analysis
os<-subset(df,df$method=="OS")
ref<-subset(df, df$method=="IND")
wide<-data.frame(os = os$meas, ref = ref$meas, id = os$subject)


# callibration curves -----------------------------------------------------

x=seq(0,3,0.1)
y1=-6.27615*x+16.4456
y2=-6.134969*x+16.25767
y3=-6.0241*x+16.09036
y4=-6.04351*x+16.1196615

data<-matrix(c(y1,y2,y3,y4), byrow=F, nrow=31)
l1<-0:2.1
ly<-c(4,4,4)
lx<-c(2.03,2.03,2.03,2.03)
ly2<-c(0,1,2,4)

plot(x,y1, type="o",col = "blue", lty=1, xlim=c(1.5,2.5), ylim = c(2,5),
     xaxt="n",yaxt="n",ylab="pH",cex=2,xlab="Volts",main="Calibration Line")
points(x, y2, col="red", pch="*")
lines(x, y2, col="red",lty=2)
points(x, y3, col="black",pch="+")
lines(x, y3, col="black", lty=3)
points(x, y4, col="orange",pch="$")
lines(x, y4, col="orange", lty=4)
axis(2,at=seq(0,14,0.5),cex.axis=1.2)
axis(1,at=seq(0,3,0.1),cex.axis=1.2)
legend(5,15,legend=c("y1","y2","y3","y4"), col=c("blue","red","black","orange"),
       pch=c("o","*","+","$"),lty=c(1,2,3,4), ncol=1)


# Exploratory Analysis ----------------------------------------------------
## Extended Bland-Altman plot via "MethodCompare" package ------------------

bland_altman_plot(wide, new ="os", Ref = "ref", ID ="id")

## Bland-Atman plot and scatter plot with id numbers based on [11]------------------
levels(df$method)<-c("Reference","Open-Source")

par(mfrow = c(1, 2), pty = "m", oma = c(0.5, 0, 0, 0),
    mar = c(4.5,4.5, 2.5, 0.5) + 0.1, 
    cex.lab = 2.5, cex.main = 2.5, cex.axis = 2.5, cex.sub=2.5)

# The following code is for Bland-Altman plot with no 95% limits of agreement
plot(1, type = "n", xlim = range(df$meas), ylim = range(df$meas), 
     xlab = "pH, Reference", ylab = "pH, Open-Source", 
     panel.first = grid(lty = "solid"))

# split data by subject id
df.split <- split(df, df$subject)
for (i in 1:length(df.split)) {
  
  # further split data by rep number
  yy <- split(df.split[[i]], df.split[[i]]$rep, drop = TRUE)
  yy1 <- numeric(length(yy))
  yy2 <- numeric(length(yy))
  for (j in 1:length(yy)) {
    
    # extract the measurements by the two methods
    yy1[j] <- yy[[j]][yy[[j]]$method == "Reference", "meas"]
    yy2[j] <- yy[[j]][yy[[j]]$method == "Open-Source", "meas"]
  }
  # finally plot
  text(yy1, yy2, as.character(i), cex = 3)
}
abline(a = 0, b = 1)
title(main = "(a)")

#The following code is for the scatter plot
plot(1, type = "n", xlim = range(df$meas), ylim = c(-0.4, 0.4), xlab = "average", 
     ylab = "difference, Open-Source (-) Reference",
     panel.first = grid(lty = "solid"))
for (i in 1:length(df.split)) {
  # further split data by rep number
  yy <- split(df.split[[i]], df.split[[i]]$rep, drop = TRUE)
  yy1 <- numeric(length(yy))
  yy2 <- numeric(length(yy))
  for (j in 1:length(yy)) {
    # extract the measurements by the two methods
    yy1[j] <- yy[[j]][yy[[j]]$method == "Reference", "meas"]
    yy2[j] <- yy[[j]][yy[[j]]$method == "Open-Source", "meas"]
  }
  
  # finally plot
  text((unlist(yy1) + unlist(yy2))/2, unlist(yy2) - unlist(yy1), as.character(i), cex = 3)
}
abline(h = 0)
title(main = "(b)")

## Trellis plot------------------
plot(groupedData(meas ~ method | subject, data = df, order.groups = TRUE), 
     ylim = c(-1, 31),xlab = "pH", ylab = "Sorted Subject ID")

## Interaction plots -------------------------------------------------------

par(mfrow = c(1, 2),cex.lab = 2, cex.main = 2, cex.axis = 1.5)
with(df, {
  interaction.plot(method, subject, meas, legend = FALSE, type = "o", 
                   lty = 1, pch = 19, cex =1, lwd = 0.8, xlab = "method", 
                   ylab = "", las = 1)
  title(ylab = "average measurement", mgp=c(2.7,1,0),main = "(a)", cex.main = 2)
  
  interaction.plot(rep, subject, meas, legend = FALSE, type = "o", 
                   lty = 1, pch = 19, cex = 1, lwd = 0.8, xlab = "time",ylab="", 
                    las = 1)
  title(ylab = "average measurement", mgp=c(3,1,0), main = "(b)", cex.main = 2)
})


# Statistical tools to assess agreement and Similarity -------------------------------------------------------------------
# Modeling  ---------------------------------------------------------------

# First the homoscedastic model is fit
# The following use file "MAMMA.R" found in [28]

levels(df$method)<-c("method1","method2")

##Fitting homoscedastic unlinked model ----------------------------------
#Function lme.unlinked.repmeas.fit, fits a homoscedastic 
#mixed-effect model with unlinked repeated measurements. 
#The function returns the model and its parameters as in table 5.

fit.results <- lme.unlinked.repmeas.fit(df)

#this function returns the parameters of the model, CCC, TDI, similarity measures,
#and repeatability CCC,TDI.
conf.results <- conf.measures.unlinked(df, fit.results$param.hat, 
                                       prob = 0.9)

#The following line extracts the homoscedastic fixed-effects model
fit.homosced<-fit.results$fit

##Fitting heteroscedastic unlinked model ----------------------------------
# v is defined as the variance covariate and is considered as an argument
# for the variance function. In this case the average of both device's measurements
#is used. In other cases we use the average of the reference device's measurements

v <- get.var.covariate(df, "avg.mean")

#vgrid is a grid of 20 values of the variance covariate
vgrid <- seq(min(df$meas), max(df$meas), l = 20)

#the following line of code fits the heteroscedastic mixed-effect model using
#the power as variance function
fit.results.power <- lme.repmeas.hetero.fit(df, v, "power", vgrid)

#the following line extracts the power heteroscedastic mixed-effects model
fit.power <- fit.results.power$fit


## Comparing the heteroscedastic and the homoscedastic model---------------------
anova(fit.homosced,fit.power)


## Diagnostics -------------------------------------------------------------
### Plot of residuals against fitted values ---------------------------------
orig.method <- df$method

hetero.stdres <- resid(fit.power, type = "pear")

# standardized residuals against fitted values
p1 <- xyplot(hetero.stdres ~ fitted(fit.power) | orig.method, 
             grid = T, data = df, xlab = "fitted values", 
             ylab = "standardized residuals", panel = function(...) {
               panel.abline(h = 0, lty = 1, col = "black")
               panel.xyplot(...)
             })

# absolute residuals against fitted values
p2 <- xyplot(abs(hetero.stdres) ~ fitted(fit.power) | orig.method, 
             grid = T, data = df, xlab = "fitted values",
             ylab = "| standardized residuals |")

print(p1, split = c(1, 1, 1, 2), more = TRUE) # ?print.trellis
print(p2, split = c(1, 2, 1, 2))
dev.off()

# Checking for normality of errors
qqnorm(fit.power, ~resid(., type = "pear") | method, abline = c(0,1))

## Model with covariates ---------------------------------------------------
fit5.sub <- lme(meas ~ method - 1 + type + Temp + turn+quantity, method = "ML", random = ~1 | subject, data = df, 
                weights = varComb(varIdent(form = ~1 | method), varPower(form = ~v | 
                                                                           method)))
# Compare the heteroscedastic models with and without covariates
anova(fit5.sub, fit.power)

# Agreement ---------------------------------------------------------------
#Function lme.unlinked.hetero.repmeas.fit, fits a heteroscedastic 
#mixed-effect model with unlinked repeated measurements. 
#The function returns the model and its parameters as in table 5.
conf.results <- conf.measures.hetero.repmeas(df, fit.results.power$param.hat, 
                                             v, "power", vgrid)
## Bland-Altman plot and 95% limits of agreement (inter/intra) ---------------------------

plot(vgrid, fit.results.power$loa$loa.inter[, 1], 
     ylim = c(-0.5,0.1), type = "l", xlab = "pH", ylab = "difference", 
     panel.first = grid(lty = "solid"))
lines(vgrid, fit.results.power$loa$loa.inter[, 2], lty = 1)
lines(vgrid, fit.results.power$loa$loa.intra.1[, 1], lty = 2)
lines(vgrid, fit.results.power$loa$loa.intra.1[, 2], lty = 2)
lines(vgrid, fit.results.power$loa$loa.intra.2[, 1], lty = 3)
lines(vgrid, fit.results.power$loa$loa.intra.2[, 2], lty = 3)
legend("bottomleft", legend = c("inter-method", "Hanna", "Open-Source"), 
       lty = 1:3, bty = "n", cex = 3)


## CCC ---------------------------------------------------------------------
# CCC
par(mfrow=c(1,2))
plot(vgrid, conf.results$ccc[, "lcb.ccc"], cex.axis = 2,yaxt="n", type = "l", 
     xlab = "pH", 
     ylab = "", ylim = c(0.3, 1), panel.first = grid(lty = "solid"),lwd= 2,
     cex.lab=2,cex.sub=2,cex.main=2)
mtext("CCC",side = 2,line=1.8,cex=2)
lines(vgrid, conf.results$rep.ccc[, "lcb.ccc.1"], lty = 2,lwd= 2)
lines(vgrid, conf.results$rep.ccc[, "lcb.ccc.2"], lty = 3,lwd= 2)
legend(3, 0.4, legend = c("inter-method", "Reference", "Open-Source"), 
       lty = 1:3, bty = "n",cex=3)
axis(2, at = seq(0,1,0.1),cex=3)
title(main = "(a)")
round(apply(conf.results$ccc, MAR=2, range)[,c(1,5,6)],4)

## TDI ---------------------------------------------------------------------
plot(vgrid, conf.results$tdi[, "ucb.tdi"],ylim = c(-0.5, 0.3), type = "l", 
     xlab = "pH", ylab = "", panel.first = grid(lty = "solid"),lwd= 2, cex.lab=2,
     yaxt="n"  ,cex.lab=2,cex.sub=2,cex.main=2)
mtext("TDI",side = 2,line=1.8,cex=2)
lines(vgrid, -conf.results$tdi[, "ucb.tdi"], lty = 1,lwd= 2)

lines(vgrid, conf.results$rep.tdi[, "ucb.tdi.1"], lty = 2,lwd= 2)
lines(vgrid, -conf.results$rep.tdi[, "ucb.tdi.1"], lty = 2,lwd= 2)
lines(vgrid, conf.results$rep.tdi[, "ucb.tdi.2"], lty = 3,lwd= 2)
lines(vgrid, -conf.results$rep.tdi[, "ucb.tdi.2"], lty = 3,lwd= 2)
legend(3, -0.4, legend = c("inter-method", "Reference", "Open-Source"), lty = 1:3, 
       bty = "n", cex = 3)
axis(2, at = seq(-0.4,0.3,0.1))
title(main = "(b)")
round(apply(conf.results$tdi, MAR=2, range)[,c(1,5,6)],4)


# Investigating possible sources of disagreenent --------------------------
## Similarity assessment ---------------------------------------------------
### Lamda -------------------------------------------------------------------

plot(vgrid, conf.results$lambda[, "est"], yaxt="n",xaxt="n",
     ylim = c(min(conf.results$lambda[,"lcl.par"]), 
              max(conf.results$lambda[, "ucl.par"])), 
     cex.main = 2.5, 
     cex.axis = 2.5,
     cex.lab = 2.5,
     type = "l", 
     xlab = "pH", ylab = "",
     panel.first = grid(lty = "solid"))
lines(vgrid, conf.results$lambda[, "lcl.par"], lty = 2)
lines(vgrid, conf.results$lambda[, "ucl.par"], lty = 2)
axis(2, at = seq(0,6,0.2), cex.axis = 2)
axis(1, at = seq(2.7,3.8,0.1), cex.axis =2)
title(ylab = "precision ratio", line = 2.6, cex.lab=2.5)  

paste0("precision ratio for 2.78 : ",round(conf.results$lambda[, "est"][1],4))
paste0("precision ratio for 3.7 : ",round(conf.results$lambda[, "est"][20],4))

# Fixed bias --------------------------------------------------------------
fit.results.power$param.hat[1]

## Assessing Repeatability -------------------------------------------------
### limits of agreement (intra) ---------------------------------------------
#### For the reference device ------------------------------------------------
summary(fit.results.power$loa$loa.intra.1)

#### For the open-source device ------------------------------------------------
summary(fit.results.power$loa$loa.intra.2)

### CCC (intra) ---------------------------------------------
#### For the reference device ------------------------------------------------
round(apply(conf.results$rep.ccc, MAR=2, range)[,c(1,5,6)],4)
#### For the open-source device ------------------------------------------------
round(apply(conf.results$rep.ccc, MAR=2, range)[,c(1,10,11)],4)
### TDI (intra) ---------------------------------------------
#### For the reference device ------------------------------------------------
round(apply(conf.results$rep.tdi, MAR=2, range)[,c(1,5,6)],4)
#### For the open-source device ------------------------------------------------
round(apply(conf.results$rep.tdi, MAR=2, range)[,c(1,10,11)],4)


# Recalibration -----------------------------------------------------------
#creating a new copy of the data frame

df_recal <- df
df_recal[df_recal$method == "method1", "meas"] <- df_recal[df_recal$method == "method1", 
                                               "meas"] + fit.results.power$param.hat["alpha"]
fit2.results <- lme.repmeas.hetero.fit(df_recal, v, "power", vgrid)
conf.results2 <- conf.measures.hetero.repmeas(df_recal, fit2.results$param.hat, 
                                             v, "power", vgrid)
# CCC after recalibration
round(conf.results2$ccc,4)[,c(1,5,6)]
# TDI after recalibration
round(conf.results2$tdi,4)[,c(1,5,6)]


# redefined functions for the MethodCompare package-----------------------------------------------------
## B_A plot source code ----------------------------------------------------
bland_altman_plot <- function(data,new="y1",Ref="y2",ID="id",fill=TRUE){
  data_sub <- (data[,c(ID,new,Ref)])
  colnames(data_sub) <- c("id","y1","y2")
  #### calculate the difference and average
  data_sub$Diff_M <- data_sub$y1-data_sub$y2
  data_sub$AVG_M <- (data_sub$y1+data_sub$y2)/2
  #### Calculate the harmonic means
  y1 <- as.numeric(data.frame(table(na.omit(data_sub)$id))$Freq)
  y2 <- as.numeric(data.frame(table(data_sub$id))$Freq)
  coef_1<- 1/mean(1/y1)
  coef_1 <- 1-1/coef_1
  coef_2<- 1/mean(1/y2)
  coef_2 <- 1- 1/coef_2
  #### With-object variance
  model_y2 <- lme(y2~1,data=data_sub,random = ~1|id, na.action = na.exclude)
  VIF<- model_y2$sigma^2*coef_2
  if (max(y1) > 1){
    model_y1 <- lme(y1~1,data=data_sub,random = ~1|id,na.action = na.exclude)
    VIF<- model_y2$sigma^2*coef_2+model_y1$sigma^2*coef_1
  }
  ####    Model on difference based on average
  data_sub_1 <- na.omit(data_sub)
  Model_1<- lm(Diff_M~AVG_M,data=data_sub_1)
  data_sub_1$fitted <- fitted(Model_1)
  data_sub_1$resid_fitted <- residuals(Model_1)
  data_sub_1$resid_fitted_abs <- abs(data_sub_1$resid_fitted)
  #####   Regression on absoulte residuals based on averge
  Model_2 <- lm(resid_fitted_abs~AVG_M,data=data_sub_1)
  data_sub_1$sig2_abs_res=fitted(Model_2)*sqrt(pi/2)
  data_sub_1$sig2_abs_res <- data_sub_1$sig2_abs_res^2
  data_sub_1$upper <- data_sub_1$fitted+1.96*sqrt(data_sub_1$sig2_abs_res+VIF)
  data_sub_1$lower <- data_sub_1$fitted-1.96*sqrt(data_sub_1$sig2_abs_res+VIF)
  ####
  if (fill) {
    mean_y1 <- aggregate(y1~id,data=data_sub,mean)
    colnames(mean_y1)[2] <- "y1.mean"
    data_sub <-merge(data_sub,mean_y1,by="id")
    data_sub$y1 <- ifelse(is.na(data_sub$y1),data_sub$y1.mean,data_sub$y1)
    data_sub$Diff_M <- data_sub$y1-data_sub$y2
    data_sub$AVG_M <- (data_sub$y1+data_sub$y2)/2
  }
  ### model for plot
  Model_3 <- lm(upper~AVG_M,data=data_sub_1)
  Model_4 <- lm(lower~AVG_M,data=data_sub_1)
  max <- max(abs(data_sub$Diff_M),na.rm=TRUE)

  ##### final plot
  par(mar=c(3.5,3.5,2,2)+0.1)
  plot(data_sub$AVG_M,data_sub$Diff_M,xlab = "",ylab = "",axes = F,
       col="grey",ylim=c(-0.45,0.05))
  title(main="Extended Bland-Altman limits of agreement (LoA) plot",cex.main=1.2)
  xlab="Average:(Open-Source (+) Reference)/2"
  ### Add the y axis
  axis(2,col="black",las=1)
  mtext("Difference:Open-Source (-) Reference",side = 2,line=2)
  box(col="black")
  ### Add the x axis
  axis(1)
  mtext("Average:(Open-Source (+) Reference)/2",side=1,col="black",line=2)
  abline(Model_1$coefficients,col="blue",lwd=2,lty=1)
  abline(h=0,col="black",lwd=2)
  abline(Model_3$coefficients,col="blue",lty=2,lwd=2)
  abline(Model_4$coefficients,col="blue",lty=2,lwd=2)
  legend("topleft",legend=c("Regression line","95% LoA"),col = c("blue","blue"),
         y.intersp = 0.7,yjust=0.5,lty=c(1,2),bty = "n",cex=0.8)
}

measure_compare <- function(data,new="y1",Ref="y2",ID="id"){
  #data_sub <- (data[,c("id","y1","y2")])
  #View(data_sub)
  data_sub <-(data[,c(ID,new,Ref)])
  colnames(data_sub) <- c("id","y1","y2")
  #### Calculate average value for reference method and categorize it into 10
  #### or 5 categories
  Mean_y2 <- aggregate(y2~id,data=data_sub,mean)
  N_cat <- ifelse(dim(Mean_y2)[1]>=100,11,6)
  Mean_y2$cat_y2_mean <- cut(Mean_y2$y2,
                             breaks=quantile(Mean_y2$y2, probs=seq(0,1,length.out = N_cat)),
                             include.lowest=TRUE)
  colnames(Mean_y2)[2] <- "y2_mean"
  data_sub <- merge(data_sub,Mean_y2,by="id")
  ### create a dataset for reference method by excluding missing value
  data_sub_y2 <- data_sub[ , c("y2")]
  data_sub_nomiss_y2 <- data_sub[complete.cases(data_sub_y2), ]
  ### Model 1: mixed model for BLUP estimate of x
  model_1 <- lme(y2 ~ 1, data = data_sub_nomiss_y2, random = ~ 1|id,
                 weights=varIdent(form = ~1|cat_y2_mean), na.action = na.exclude) # breaks down when n2 small...
  #model_1 <- lme(y2 ~ 1, data = data_sub_nomiss_y2, random = ~ 1|id, na.action = na.exclude)
  data_sub_nomiss_y2$y2_hat <- fitted(model_1)
  ### create a dataset with y2_hat
  data_sub_y2_hat <- data_sub_nomiss_y2[ , c("id","y2_hat")]
  data_sub_y2_hat <- aggregate(y2_hat~id,data=data_sub_y2_hat,mean)
  data_sub <- merge(data_sub,data_sub_y2_hat,by="id")
  ### create a dataset for new method by excluding missing value of new method y1
  data_sub_y1 <- data_sub[ , c("y1")]
  data_sub_nomiss_y1 <- data_sub[complete.cases(data_sub_y1), ]
  data_newmethod <- data_sub_nomiss_y1
  ########         Models on the reference method
  #### Model 2: regression of y2 based on BLUP of x
  model_2 <- lm(y2 ~ 0+y2_hat, data = data_sub, na.action = na.exclude) # not useful...
  #model_2 <- lm(y2 ~ y2_hat, data = data_sub, offset = rep(0, length(y)) # not useful...
  #model_2$coefficients
  data_sub$fitted_y2 <- data_sub$y2_hat
  data_sub$resid_y2 <- data_sub$y2-data_sub$y2_hat
  data_sub$resid_y2_abs <- abs(data_sub$resid_y2)
  #### Model 3: Estimation of variance function for y2
  model_3 <- lm(resid_y2_abs~y2_hat,data=data_sub, na.action = na.exclude)
  ### Smooth standard deviation estimate
  data_sub$sig_resid_y2 <- fitted(model_3)*sqrt(pi/2)
  #######         Models on the new method - original measured value
  #### Model 4: Regression of y1 based on BLUP of x
  model_4 <- lm(y1~y2_hat,data=data_newmethod)
  #### Differential and proportional bias of new method
  Bias <- cbind(model_4$coefficients,confint(model_4))
  rownames(Bias) <- c("Differential bias","Proportional bias")
  colnames(Bias)[1] <- "Estimate"
  #### Residuals and corrected y1
  data_newmethod$resid_y1 <- residuals(model_4)
  data_newmethod$resid_y1_abs <- abs(data_newmethod$resid_y1)
  data_newmethod$y1_corrected <- (data_newmethod$y1-Bias[1])/Bias[2]
  data_newmethod$Bias <- Bias[1]+data_newmethod$y2_hat*(Bias[2]-1)
  #### Model 5: Variance estimation for y1 based on BLUP of x
  model_5<- lm(resid_y1_abs ~ y2_hat,data=data_newmethod, na.action = na.exclude)
  data_newmethod$fitted_y1 <- fitted(model_5)
  ### Smooth standard deviation estimate
  data_newmethod$sig_resid_y1 <- fitted(model_5)*sqrt(pi/2)
  ########        Models on new method- corrected value
  #### Model 6: regression on corrected y1 using BLUP_y2
  model_6 <- lm(y1_corrected~y2_hat,data=data_newmethod, na.action = na.exclude)
  data_newmethod$y1_corrected_fitted <- fitted(model_6)
  data_newmethod$resid_y1_corrected <- residuals(model_6)
  data_newmethod$resid_y1_corrected_abs <- abs(data_newmethod$resid_y1_corrected)
  #### Model 7: variance function estimation for corrected y1
  model_7 <- lm(resid_y1_corrected_abs ~ y2_hat,data=data_newmethod, na.action = na.exclude)
  data_newmethod$sig_resid_y1_corrected <- fitted(model_7)*sqrt(pi/2)
  ######  Final output
  ## 1. Bias (differential and proportional) from New method
  ## 2. Output values including y2_hat, sig_resid_y2, sig_resid_y1_corr, Bias,
  ## 3. list of all models: mixed model
  Model_list <-  list(model_1,model_2, model_3,model_4,model_5,
                      model_6,model_7)
  data_sub <- data_sub[,c("id","y1","y2","y2_hat","sig_resid_y2")]
  data_newmethod <- data_newmethod[,c("id","y1","y2_hat","Bias","y1_corrected",
                                      "sig_resid_y1_corrected")]
  output <- list(Bias=Bias,Models=Model_list,Ref=data_sub,New=data_newmethod)
  return(output)
}

bias_plot <- function(object){
  ### extract the objects from the output
  Bias <- object$Bias
  data_old <- object$Ref
  data_new <- object$New
  Models <- object$Models
  subtitle <- paste("Differential bias:",round(Bias[1,1],3),"; ",
                    "Proportional bias:",round(Bias[2,1],3),sep="")
  min_y <- min(min(data_old$y2,na.rm=TRUE),min(data_new$y1,na.rm=TRUE))
  max_y <- max(max(data_old$y2,na.rm=TRUE),max(data_new$y1,na.rm=TRUE))
  min_y <- floor(min_y)
  range = max_y-min_y
  max_y <- ceiling(max_y+range*0.2)
  par(mar=c(3.5,3.5,3,4)+0.1)
  ### plot scatter for y1 and y2 with respect to BLUP of x
  plot(data_old$y2_hat,data_old$y2,axes=F,xlab="",ylab="",
       col="grey",ylim = c(min_y,max_y),cex=0.5)
  title(main="Bias",cex.main=0.9)
  points(data_new$y2_hat,data_new$y1,col="blue",pch=19,cex=0.5)
  
  ### Add the fitted line from model 2 and model 4,
  #abline(Models[[2]]$coefficients,lwd=2)
  abline(c(0,1),lwd=2)
  abline(Models[[4]]$coefficients,lwd=2,lty=2,col="blue")
  
  ### Add the subtitle
  mtext(subtitle,side=3,cex=0.8)
  ### Add the left y axis
  axis(2,col="black",las=1)
  mtext("Reference and Open-Source",side = 2,line=2.3)
  box(col="black")
  
  ### Add second plot: bias axis
  par(new=TRUE)
  ## Plot the bias plot and put axis scale on right
  plot(data_new$y2_hat,data_new$Bias,xlab="",ylab="",axes=FALSE,
       col="red",lty=1,type="l")
  abline(h=0,col="black",lwd=2)
  ## Add the right y axis and label
  mtext("Bias",side=4,col="black",line=2.5)
  axis(4,col="black",col.axis="black",las=1)
  ### Draw the x axis and add the label
  axis(1)
  mtext("BLUP of x",side=1,col="black",line=2)
  ### Add the legend
  lm_bias <- lm(data_new$Bias~data_new$y2_hat)
  if (coef(lm_bias)[2] <0){
    legend("top",legend=c("Reference","Open-Source","Bias"),
           pch = c(1,19),col = c("black","blue","red"),pt.cex=c(1,1,0.0001),
           y.intersp = 0.7,yjust=0.2,lty=c(1,2,1),bty = "n",cex=0.8)}
  else {legend("topleft",legend=c("Reference","Open-Source","Bias"),
               pch = c(1,19),col = c("black","blue","red"),pt.cex=c(1,1,0.0001),
               y.intersp = 0.7,yjust=0.2,lty=c(1,2,1),bty = "n",cex=0.8)}
}

compare_plot <- function(object){
  data_old <- object$Ref
  data_new <- object$New
  Models <- object$Models
  max <- max(data_old$y2,data_new$y1,data_new$y1_corrected,na.rm=TRUE)
  min <- min(data_old$y2,data_new$y1,data_new$y1_corrected,na.rm=TRUE)
  range <- max-min
  par(mar=c(3.5,3.5,2,2)+0.1)
  plot(data_old$y2_hat,data_old$y2,pch=1,cex=1,col="grey",axes = F,
       xlab="",ylab="",ylim=c(min-range*0.1,max+range*0.2))
  title(main="Comparison of the methods",cex.main=1.2)
  xlab="Average:(Reference (+) Open-Source)/2"
  ### Add the y axis
  axis(2,col="black",las=1)
  mtext("Measurement",side = 2,line=2)
  box(col="black")
  ### Add the x axis
  axis(1)
  mtext("BLUP of x",side=1,col="black",line=2)
  points(data_new$y2_hat,data_new$y1,pch=19,col="blue",cex=1)
  points(data_new$y2_hat,data_new$y1_corrected,pch=18,col="red",cex=1)
  #abline(Models[[2]]$coefficients,lwd=2)
  abline(c(0,1),lwd=2)
  abline(Models[[4]]$coefficients,lwd=2,lty=1,col="blue")
  abline(Models[[6]]$coefficients,lwd=2,lty=2,col="red")
  legend("topleft",legend=c("Reference","Open-Source",
                            "Open-Source (corrected)"),
         pch=c(1,19,18),lty = c(1,1,2),col=c("black","blue","red"),
         y.intersp = 0.7,yjust=0.2,bty = "n",cex=1)
}

precision_plot <- function(object){
  data_old <- object$Ref
  data_new <- object$New
  min_y <- min(min(data_old$sig_resid_y2,na.rm=TRUE),min(data_new$sig_resid_y1_corrected,na.rm=TRUE))
  max_y <- max(max(data_old$sig_resid_y2,na.rm=TRUE),max(data_new$sig_resid_y1_corrected,na.rm=TRUE))
  min_y <- floor(min_y)
  range = max_y - min_y
  max_y <- ceiling(max_y+range*0.2)
  par(mar=c(3.5,3.5,2,2)+0.1)
  plot(data_old$y2_hat,data_old$sig_resid_y2,xlab="",
       ylab="",axes = F,
       cex=0.8,ylim=c(min_y,max_y))
  title(main="Precision plot",cex.main=0.9)
  ### Add the y axis
  axis(2,col="black",las=1)
  mtext("Standard deviation of measurement errors",side = 2,line=2,cex=0.9)
  box(col="black")
  ### Add the x axis
  axis(1)
  mtext("BLUP of x",side=1,col="black",line=2,cex=0.9)
  points(data_new$y2_hat,data_new$sig_resid_y1_corrected,
         cex=0.8,pch=19,col="blue")
  legend("topleft",legend=c("Reference","Open-Source(corrected)"),
         pch=c(1,19),col=c("black","blue"),xpd=T,horiz = F,
         box.lwd=0,cex=0.7)
}
