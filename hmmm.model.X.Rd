\name{hmmm.model.X}
\alias{hmmm.model.X}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EMRS model as proposed by Colombi R., Giordano S., and Tutz G. (JEBS, 2021)
}
\description{
A function for constructing the marginal model of response-style latent 
variables, with one latent variable per response. 
It uses the hmmm package syntax.
}
\usage{
hmmm.model.X(marg, lev, names, Formula = NULL, strata = 1, fnames = NULL, cocacontr = NULL, ncocacontr = NULL, replace = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marg}{
A string specifying the marginals where the interactions are defined. See 
the hmmm package for reference. 
Both univariate and bivariate marginals must be included.
}
  \item{lev}{
Numeric vector indicating the levels of the latent variables. 
Usually c(2, 2, …), as response-style variables are dichotomous.
}
  \item{names}{
Usually names = c("LA", "LB", …). 
It specifies the names of the latent variables, one for each response.
}
  \item{Formula}{
A list of formulas describing how covariates affect interactions 
among latent variables. Refer to the hmmm package. 
Interactions of three or more variables must 
be set to zero (e.g., LA.LB.LC = "zero"). 
}
  \item{strata}{
Numeric vector indicating the number of levels for each covariate.
}
  \item{fnames}{
Names of the covariates, listed in the same order as their 
levels appear in the strata.
}
  \item{cocacontr}{
This option is currently not used.
}
  \item{ncocacontr}{
This option is currently not used.
}
  \item{replace}{
Logical value. If TRUE, the object describing the model is created; 
if FALSE, only partial information is reported. Default is TRUE.
}
}
\details{
The sequence of the estimation procedure is as follows:
a) create a data matrix where each column contains the joint response for 
a given combination of covariate levels (see examples). Continuous covariates
are not considered. The order in which the categories change depends on the order 
in which the variables are introduced (according to the R package hmmm).
b) create the latent model using hmmm.model.X() and the observed model
using hmmm.model.T().
c) estimate the parameters of the model using Uhmmm.MLfit(), whose output 
contains all the elements needed for further analyses.
}
\value{
An object representing the marginal model for the latent variables, 
to be passed to the Uhmmm.MLfit function.
}
\references{
Colombi, R., Giordano, S., & Tutz, G. (2021). 
A Rating Scale Mixture Model to Account for the Tendency to Middle and Extreme 
Categories. Journal of Educational and Behavioral Statistics, 46(6), 682-716. 
https://doi.org/10.3102/1076998621992554 

Colombi, R., Giordano, S., & Cazzaro, M. (2014).
hmmm: An R Package for Hierarchical Multinomial Marginal Models. 
Journal of Statistical Software, 59(11), 1–25. https://doi.org/10.18637/jss.v059.i11
}
\author{
Colombi R., Giordano S.
}

\examples{
# run the following code to reproduce example of Section 6 
# of the paper Colombi, R., Giordano, S., & Tutz, G. (2021). JEBS
library(matrixcalc)
library(car)
library(MASS)
library(maxLik)
library(hmmm)
library(nleqslv)
library(Matrix)
library(gtools)
library(CompQuadForm)
remove(list=ls())

source("margmod.R") # it works with R package hmmm
source("Uhmmm_mainsource.R") # it contains the main function Uhmmm.MLfit

options(warn=-1)

# Data
data<-read.csv("datietnia.csv",header=T, sep = ";") ## data from GSS 
head(data)
names(data)

############################ Responses -- see Section 6 of JEBS 
# 1 = "strongly disagree", 2 = disagree, 3 = "Neither agree nor disagree"
# 4 = "agree", 5 = "strongly agree"
E2 <- as.factor(data$ethadapt) # R1 
E5 <- as.factor(data$ethignor) # R2                

############### Covariates ############################
edu<-as.factor(data$edu)         # <= 13 is  1, > 13 is 2
race <- as.factor(data$race)    # 1 = other-black, 2 = white
politics <- as.factor(data$politics) # SREP = 1, NREP = 2, IND = 3, NDEM = 4, SDEM = 5 		

## create dataframe
y<-as.data.frame(cbind(E2,E5,politics,race,edu))
ytab<-matrix(table(y),25,20)
ytab
remove(data)

labelrisp<-c("R1","R2")
labelfac<-c("PO","RC","ED")       
str<-c(5,2,2)   

rc<-5  # num. of response categories

#########################   models hmmmlu --- use package hmmm 
# Refer to Colombi, R., Giordano, S., & Cazzaro, M. (2014).
# hmmm: An R Package for Hierarchical Multinomial Marginal Models. 
# Journal of Statistical Software, 59(11), 1–25. https://doi.org/10.18637/jss.v059.i11

marg<-marg.list(c("l-m","m-l","l-l"), mflag="m")  # see package hmmm

# Formula for latent components (under alternative assumptions)

# Flatnoeffind<-list(LA=~1,LB=~1,LA.LB="zero") # no effects of covariates/only 
# intercepts (~1) in the logit models for the two latent components, 
# and independence (LA.LB="zero")

# Flatnoeff<-list(LA=~1,LB=~1,LA.LB=~1) # no effects of covariates/only 
# intercepts (~1) in the logit models for the two latent components, 
# and association between the two latent var

Flat<-list(LA=~RC+ED,LB=~RC+ED,LA.LB=~1) # covariates RC, ED in the logit models
# for the two latent components, and association between the two latent var

# Flatind<-list(LA=~RC+ED,LB=~RC+ED,LA.LB="zero") # covariates RC, ED in the 
# logit models for the two latent components, and independence between the two 
# latent var

# Specify the model for the latent components 
# using an appropriate formulation from those previously introduced, 
# in line with the assumed framework

# modellat1<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flatnoeff, fnames=labelfac)

modellat2<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flat, fnames=labelfac)

# modellat3<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flatind, fnames=labelfac)

# modellat4<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flatnoeffind, fnames=labelfac)


## The function is currently defined as
function (marg, lev, names, Formula = NULL, strata = 1, fnames = NULL, 
    cocacontr = NULL, ncocacontr = NULL, replace = TRUE) 
{
    str <- prod(strata)
    model <- hmmm.model(marg = marg, lev = lev, names = names, 
        strata = str, cocacontr = cocacontr)
    model <- create.XMAT(model, Formula = Formula, strata = strata, 
        fnames = fnames, cocacontr = cocacontr, ncocacontr = ncocacontr, 
        replace = replace)
  }
}
%
