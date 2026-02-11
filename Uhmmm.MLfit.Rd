\name{Uhmmm.MLfit}
\alias{Uhmmm.MLfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit the model proposed by Colombi C., Giordano S.and Tutz G. (JEBS, 2021)
}
\description{
The sequence of the estimation procedure is as follows:
a) create a data matrix where each column contains the joint response for 
a given combination of covariate levels. Continuous covariates 
are not considered. The order in which the categories change depends on the 
order in which the variables are introduced (according to the R package hmmm).
b) create the latent model using hmmm.model.X() and the observed model
using hmmm.model.T().
c) estimate the parameters of the model using Uhmmm.MLfit(), whose output 
contains all the elements needed for further analyses.
}
\usage{
Uhmmm.MLfit(table, modelobs = NULL, filemodelobs = NULL, modellat = NULL, filemodellat = NULL, modelcomp2 = NULL, start = NULL, startvec = NULL, startmth = "SANN", method = "NM", method.inv = "Broyden", print.level = 0, finalHessian = TRUE, tol = 1e-09, iterlim = 1000, NMiter = 0, numeric = FALSE, fixed = NULL, compare = FALSE, recalc = FALSE, control = list(allowSingular = FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
Data matrix where each column contains the joint response for a combination of 
covariate levels. Continuous covariates are not considered. 
The order in which the categories change depends on the order in which 
the variables were introduced in hmmm.model.X() and hmmm.model.T().
}
  \item{modelobs}{
An object created by hmmm.model.T()}
  \item{filemodelobs}{
A file containing modelobs, if it has been created previously
}
  \item{modellat}{
An object created by hmmm.model.X()
}
  \item{filemodellat}{
A file containing modellat, if it has been created previously
}
  \item{modelcomp2}{
Do not use
}
  \item{start}{
An optional vector of starting values for the parameters 
}
  \item{startvec}{
not used
}
  \item{startmth}{
Typically 'NM' = Nelder-Mead, it specifies the optimization method used 
during the first NMiter iterations.
}
  \item{method}{
Optimization Method (BFGS recomended) to be used after the first 
NMiter iterations 
}
  \item{method.inv}{
An option from the nleqslv package: "Broyden" or "Newton".
}
  \item{print.level}{
Controls iteration output: 1 = print details per iteration, 0 = no details
}
  \item{finalHessian}{
If TRUE, compute Fisher information matrix and standard errors
}
  \item{tol}{
Convergence criterium
}
  \item{iterlim}{
Max number of iterations of the main optimization algorithm (usually BFGS)
}
  \item{NMiter}{
Maximum number of iterations for the initial algorithm (NM recommended)
}
  \item{numeric}{
Do not use 
}
  \item{fixed}{
Do not use }
  \item{compare}{
Do not use 
}
  \item{recalc}{
If FALSE (default), joint probabilities are recomputed only when updating 
the log-likelihood. If TRUE, they are also recomputed when evaluating the 
score functions
}
  \item{control}{
See the package nleqslv
}
}
\details{
BFGS and BFGSR often produce an error of the form:

"Error in nleqslv(x = start, myfun, method = method.inv, jac = myder, 
control = control, : non-finite value(s) returned by jacobian (row=1,col=1)"

which occurs when the nleqslv function is called to compute joint probabilities from the parameters (Appendix B of the supplementary material).

In this case, the solution of the nonlinear system of equations 
f(x)b using the package nleqslv is replaced by a minimization problem 
min∥f(x)−b∥ solved with optim. This situation is signaled by an error message
from nleqslv and typically occurs frequently during the first iterations of 
the maximization algorithm, after which it disappears.
The problem depends on the initial values of the parameters and is less likely 
to occur when a sufficient number of initial EM iterations is used.
}
\value{
A list containing the results of the fitting procedure:

object$vecpar: parameter estimates and standard errors

object$fit$maximum: log-likelihood.

Other elements of object include:

"VarMatrices": variance-covariance matrices of estimates

"vecPtr": estimated probabilities for latent variables

"vecPobs": estimated probabilities for every combination of covariate levels
(for content-driven and EMRS models)

"vecres": standardized residuals for each stratum.

In particular, "fit" is an object of class "maxLik" and can be used accordingly.
}
\references{
Colombi, R., Giordano, S., & Tutz, G. (2021). 
A Rating Scale Mixture Model to Account for the Tendency to Middle and Extreme 
Categories. Journal of Educational and Behavioral Statistics, 46(6), 682-716. 
https://doi.org/10.3102/1076998621992554 
}
\author{
Colombi Roberto, Giordano Sabrina
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
source("Uhmmm_mainsource.R") # it contains the function Uhmmm.MLfit

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

rc<-5  # num. categories

#########################   models hmmmlu --- use package hmmm 
# Refer to Colombi, R., Giordano, S., & Cazzaro, M. (2014).
# hmmm: An R Package for Hierarchical Multinomial Marginal Models. 
# Journal of Statistical Software, 59(11), 1–25. https://doi.org/10.18637/jss.v059.i11

marg<-marg.list(c("l-m","m-l","l-l"), mflag="m")

# the next structure is valid for two responses, for three and more responses it is different (see help pages)

# Formula model for observed responses

Formula<-list(R1=~LA*PO+LA*R1,R2=~LB*PO+LB*R2) 
# LA corresponds to the latent variable for R1, and LB for R2.
# In this specification, PO is included as a covariate in the logit models for both responses
# (see Eqs. 16 and 17 in JEBS).

# Formula for latent components (under alternative assumptions)

# Flatnoeffind<-list(LA=~1,LB=~1,LA.LB="zero") # no effects of covariates/only 
# intercepts (~1) in the logit models for the two latent components, 
# and independence (LA.LB="zero")

# Flatnoeff<-list(LA=~1,LB=~1,LA.LB=~1) # no effects of covariates/only 
# intercepts (~1) in the logit models for the two latent components, 
# and association between the two latent var

Flat<-list(LA=~RC+ED,LB=~RC+ED,LA.LB=~1) # covariates RC, ED in the logit models
# for the two latent components, and association between the two latent var.

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
                          UNC = "TUTZ1",  # Weights: sr = 1.5, 0.5, -0.5, -1.5 for r = 1,2,3,4. Other options possible.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
