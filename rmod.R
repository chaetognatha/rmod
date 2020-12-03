#!/usr/bin/Rscript
#USAGE
#------------------------------------------------------------------------------
#source('~/rmods/rmod/rmod.R')
#DEPENDENCIES
#------------------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(car)
library(lmodel2)
library(lmtest)
library(tidyverse)
library(gridExtra)
library(Hmisc)
library(MuMIn)
library(ppcor)
library(rgl)
library(knitr)
library(rmarkdown)
#SETTINGS
#------------------------------------------------------------------------------
my.par <- par(mfrow=c(2,2)) #show plots in a 2x2 window, reset with c(1,1)
my.palette <- terrain.colors(12)
# FUNCTIONS
#------------------------------------------------------------------------------
# Function to perform a kolmogorov test with simulated dataset v1
kolmogorov.test <- function(x){
  y <- pnorm(summary(x), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  return(ks.test(x,y))
}
#------------------------------------------------------------------------------
# Histogram with a normal distribution trendline v2
# Takes a dataset and a name of dataset, default to My Histogram
dist.line.hist <- function(x, my.name="My Histogram"){
  h <- hist(x, main = my.name, xlab = my.name)
  xfit <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length=100)
  yfit <- dnorm(xfit,mean=mean(x, na.rm = TRUE), sd=sd(x, na.rm = TRUE))
  yfit <- yfit*max(h$counts)/max(yfit)
  lines(xfit,yfit, col="red", lwd=2)
}
#------------------------------------------------------------------------------
# Function to perform common tests of normality and review normally distributed data v2
tests.of.normality <- function(x){
  summary(x)
  str(x)
  par(mfrow=c(2,3))
  # in the first case we have a column that is a factor
  if(is.factor(x[,1])){
for(i in 1:(length(x)-1)){
      cat("Test results of: ", colnames(x[i+1]), "\n ----------------- \n")
      print(shapiro.test(x[,i+1]))
      print(kolmogorov.test(x[,i+1]))
      qqnorm(x[,i+1])
      qqline(x[,i+1])
      hist(x[,i+1], col=my.palette, main = colnames(x[i+1]), xlab = colnames(x[i+1]))
      dist.line.hist(x[,i+1], colnames(x[i+1]))
    }
  }
  else{
    # in this case we only have columns that we want to analyse
    for(i in 1:(length(x))){
      cat("Test results of: ", colnames(x[i]), "\n ----------------- \n")
      print(shapiro.test(x[,i]))
      print(kolmogorov.test(x[,i]))
      qqnorm(x[,i])
      qqline(x[,i])
      hist(x[,i], col=my.palette, main = colnames(x[i]), xlab = colnames(x[i]))
      dist.line.hist(x[,i], colnames(x[i]))
    }
  }
  par(my.par) #reset
}
#---------------------------------------------------------------------------------------------------
#no function for getting standard error
# so I made one :)
std.err <- function(sel.col){
  sd(sel.col, na.rm=T)/sqrt(length(!is.na(sel.col)))
}
#-----------------------------------------------------------------------------------
# for testing a model used for anova
anova.model.analysis <- function(model){
  #Check assumptions
  #Independent measurements
  print("Checking for independent measurements:")
  print(dwtest(model))
  print("IF Not significant, means no autocorrelation")
  
  
  #Homogeneity of variances and normality of residuals
  par(mfrow=c(2,2))
  plot(model)
  par(mfrow=c(1,1))
  print("Checking for normality in residuals:")
  print(shapiro.test(model$residuals))
  print("If the test is not significant, the data is normally distributed!")
  x <- model$residuals
  y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
  print("And with the Kolmogorov test:")
  print(ks.test(x, y))
  print("If the P value is small, conclude that the two groups were sampled from populations with different distributions. The populations may differ in median, variability or the shape of the distribution.")
  print("Levene Test results:")
  print(leveneTest(model))
  print("non-significance means the two variances are approximately equal")
  
  #Get results of one-way anova
  print("Here's an F-test:")
  print(anova(model))
  print("And here you can find your estimates:")
  print(summary(model))
}
#--------------------------------------------------------------------------------
# arcsine squareroot transform function
# generally recommended for proportional data
asinTransform <- function(p) { asin(sqrt(p)) }
    