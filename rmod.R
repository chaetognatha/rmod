#!/usr/bin/Rscript
#USAGE
#------------------------------------------------------------------------------
#source('~/rmods/rmod.R')
#DEPENDENCIES
#------------------------------------------------------------------------------
library("car")
library("lmodel2")
library("lmtest")
library("tidyverse")
#SETTINGS
#------------------------------------------------------------------------------
my.par <- par(mfrow=c(2,2)) #show plots in a 2x2 window, reset with c(1,1)
my.palette <- terrain.colors(12)
# FUNCTIONS
#------------------------------------------------------------------------------
# Function to perform a komogorov test with simulated dataset v1
komogorov.test <- function(x){
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
      print(komogorov.test(x[,i+1]))
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
      print(komogorov.test(x[,i]))
      qqnorm(x[,i])
      qqline(x[,i])
      hist(x[,i], col=my.palette, main = colnames(x[i]), xlab = colnames(x[i]))
      dist.line.hist(x[,i], colnames(x[i]))
    }
  }
  par(my.par) #reset
}
#no function for getting standard error
# so I made one :)
std.err <- function(sel.col){
  sd(sel.col, na.rm=T)/sqrt(length(!is.na(sel.col)))
}
