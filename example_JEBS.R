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
source("Uhmmm_mainsource.R")

options(warn=-1)

# Data
dati<-read.csv("datietnia.csv",header=T, sep = ";") ## dati GSS 
head(dati)
names(dati)

############################ Responses-- see Section 6 of JEBS 

# 1 = "strongly disagree", 2 = disagree, 3 = "Neither agree nor disagree"
# 4 = "agree", 5 = "strongly agree"

E2 <- as.factor(dati$ethadapt) # R1 
E5 <- as.factor(dati$ethignor) # R2                

############### Covariates ############################

edu<-as.factor(dati$edu)         # <= 13 is  1, > 13 is 2
race <- as.factor(dati$race)    # 1 = other-black, 2 = white
politics <- as.factor(dati$politics) # SREP = 1, NREP = 2, IND = 3, NDEM = 4, SDEM = 5 		

## create dataframe
y<-as.data.frame(cbind(E2,E5,politics,race,edu))
ytab<-matrix(table(y),25,20)
# ytab is a tabular data structure, one column for each stratum given by combinations of covariate levels
remove(dati)

labelrisp<-c("R1","R2") 
labelfac<-c("PO","RC","ED")       
str<-c(5,2,2)   

rc<-5  # num. categories for each response

#########################   models hmmmlu --- use package hmmm 
# Refer to Colombi, R., Giordano, S., & Cazzaro, M. (2014).
# hmmm: An R Package for Hierarchical Multinomial Marginal Models. 
# Journal of Statistical Software, 59(11), 1–25. https://doi.org/10.18637/jss.v059.i11

marg<-marg.list(c("l-m","m-l","l-l"), mflag="m")

# the next structure is valid for two responses, for three and more responses it is different (see help pages)

######## Formula model for observed responses

Formula<-list(R1=~LA*PO+LA*R1,R2=~LB*PO+LB*R2) 
# LA corresponds to the latent variable for R1, and LB for R2.
# In this specification, PO is included as a covariate in the logit models for both responses
# (see Eqs. 16 and 17 in JEBS).

######## Formula for latent components (under alternative assumptions)

Flatnoeffind<-list(LA=~1,LB=~1,LA.LB="zero") # no effects of covariates/only intercepts (~1) in the logit models for the two latent components, and independence (LA.LB="zero")

Flatnoeff<-list(LA=~1,LB=~1,LA.LB=~1) # no effects of covariates/only intercepts (~1) in the logit models for the two latent components, and association between the two latent var.

Flat<-list(LA=~RC+ED,LB=~RC+ED,LA.LB=~1) # covariates RC, ED in the logit models for the two latent components, and association between the two latent var.

Flatind<-list(LA=~RC+ED,LB=~RC+ED,LA.LB="zero") # covariates RC, ED in the logit models for the two latent components, and independence between the two latent var.

# Specify the model for the latent components 
# using an appropriate formulation from those previously introduced, 
# in line with the assumed framework

modellat1<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flatnoeff, fnames=labelfac)

modellat2<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flat, fnames=labelfac)

modellat3<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flatind, fnames=labelfac)

modellat4<-hmmm.model.X(marg=marg,lev=c(2,2),names=c("LA","LB"), strata=str, Formula=Flatnoeffind, fnames=labelfac)


control<-list(allowSingular=FALSE)

################################# Fit model 

modelobs <- hmmm.model.T(
                          marg,
                          lev = c(5, 5), # Number of levels for each response
                          names = labelrisp, # Labels for the responses
                          strata = str,
                          Formula = Formula, # Formula specifying the responses 
                          fnames = labelfac, # Labels for covariates
                          replace = TRUE, # If TRUE, allows estimation of parameters
                          INDIP = "CST4", # Specifies 4 association parameters for answering behavior (pairwise association), other options are CST2 and IND
                          UNC = "TUTZ1",  # Scores: sr = 1.5, 0.5, -0.5, -1.5 for r = 1,2,3,4. Other options possible.
                          MOD = "A"  # Model type as in the JEBS paper; other model types allowed by the source file
                        )


fitmodel <- Uhmmm.MLfit(
                        ytab,
                        modelobs = modelobs,
                        modellat = modellat2, 
                        method = "BFGS",  # Optimization method: BFGS (quasi-Newton), alternatives: NR (Fisher Scoring), BFGS (with first derivative). Note: BFGS can be unstable; NM is stable but slower.
                        method.inv = "Newton",  # Inversion method from eta to theta. See the 'nleqslv' package for implementation
                        numeric = FALSE,  # If TRUE, uses numerical derivatives; if FALSE, uses analytic derivatives
                        print.level = 1,  # Controls iteration output: 1 = print details per iteration, 0 = no details
                        finalHessian = TRUE,  # If TRUE, compute Fisher information matrix and standard errors
                        iterlim = 1500,  # Maximum number of iterations for the main optimization method
                        NMiter = 1000,  # Maximum number of iterations for the initial 'startmth' algorithm
                        tol = 1e-7,  # Convergence tolerance for the optimizer
                        startmth = "NM",  # Algorithm used for starting values ('NM' = Nelder-Mead)
                        recalc = TRUE  # If TRUE, recalculate starting parameters (can optionally provide 'start = stpar')
                      )


round(fitmodel$vecpar, 4) # this output contains estimates, standard errors, residuals 

# The estimated parameters reported as output are presented in Table 1 of the paper. 
# Their use in computing the coefficients of the content-driven and EMRS logits models 
# is illustrated in the Excel file Results.xls.

# The output of this script has been saved in the file "output", located in the folder.
# You can therefore work directly with the results by loading it via load("output").
# The main results (estimates, standard errors, and residuals) are stored in
# fitmodel$vecpar and fitmodel$vecres.
# for the interpretation of the parameters see the file Results.xls,
# for the boxplots of residuals see the script residuals.R.


