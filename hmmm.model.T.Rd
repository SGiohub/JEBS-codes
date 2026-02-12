\name{hmmm.model.T}
\alias{hmmm.model.T}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model for the responses conditional on the response-style latent variables, 
as proposed by Colombi R., Giordano S., and Tutz G. (JEBS, 2021)
}
\description{
A function to construct the marginal model for the observed variables, 
treating the response-style latent variables as covariates.
}
\usage{
hmmm.model.T(marg, lev, names, Formula = NULL, strata = 1, fnames = NULL, cocacontr = NULL, ncocacontr = NULL, replace = TRUE, INDIP = "CST", UNC = "TUTZ2", MOD = "A")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marg}{
A string specifying the marginals where interactions among 
the observed variables are defined. 
See the hmmm package for reference. 
Univariate and bivariate marginals must be included, 
and interactions based on local logits must be used.
}
  \item{lev}{
Numeric vector specifying the number of categories for the response variables, 
in the same order as they are listed in the names argument.
}
  \item{names}{
Labels for the response variables. 
}
  \item{Formula}{
A list of model formulas specifying the interactions defined by `marg`. It must include, at a minimum, formulas for first-order interactions.
Example:
 `Formula <- list(  
   R1 = ~ LA*G + LA*C + LA*R1,  
   R2 = ~ LB*G + LB*C + LB*R2,  
   R3 = ~ LC*G + LC*C + LC*R3  
 )`

 Note that if the covariates are `G` and `C`, specifying `R1 = ~ G + C` is incorrect, because interactions with the latent variable `LA` associated with `R1` must be included. In particular, the term `LA*R1` is required.
}

  \item{strata}{
Numeric vector specifying the levels of the covariates, 
in the same order as they are listed in fnames.
}
  \item{fnames}{
Labels for the covariates.
}
  \item{cocacontr}{
This option is currently not used.
}
  \item{ncocacontr}{
This option is currently not used.
}
  \item{replace}{
If TRUE, the object describing the model is created; if FALSE, only the design matrices are returned.
}
  \item{INDIP}{
This argument must be specified when Formula includes only first-order interaction terms. It is a list describing the type of uniform association for each pair of response variables. See the paper for a detailed discussion of uniform association conditional on the latent variables.

Recall that each pair of responses depends only on the two corresponding response-style latent variables.

For three response variables, for example:
INDIP = list("CST2", "CST2", "CST2")

For each pair of response variables, the available options are:

INDIP = "IND": the two responses are independent in the bivariate distribution given the latent variables. Only first-order interactions are non-zero.

INDIP = "CST1": second-order interactions in the bivariate distribution are non-zero and equal only in the content-driven–content-driven component.

INDIP = "CST2": second-order interactions are non-zero and equal only in the response-style–response-style and content-driven–content-driven components.

INDIP = "CST4": second-order interactions are non-zero and equal in all four components:
response-style–response-style, content-driven–content-driven, content-driven–response-style, and response-style–content-driven.

INDIP = "CST3": as in CST4, but the content-driven–response-style and response-style–content-driven interactions are constrained to be equal.
}
  \item{UNC}{
Specifies the vector of scores s_r used to model response styles 
(see the paper for details). 
The available options are:

"TUTZ0": scores equal to -1 below the midpoint and 1 above the midpoint.

"TUTZ1": scores based on the distance from the midpoint.

"TUTZ2": scores based on the local logits of a symmetric parabolic distribution.

}
  \item{MOD}{
Leave the default unchanged.
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
An object defining the model of observed variables given the latent variables, 
for use with Uhmmm.MLfit.
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
data<-read.csv("datietnia.csv",header=T, sep = ";") ## dati GSS 
head(data)
names(data)

############################ Responses-- see Section 6 of JEBS 
# Ordered scale for the considered responses 
# 1 = "strongly disagree", 2 = disagree, 3 = "Neither agree nor disagree"
# 4 = "agree", 5 = "strongly agree"

# Questionnaire and data are at https://gssdataexplorer.norc.org/

E1 <- as.factor(data$whtsdiff) #
E2 <- as.factor(data$ethadapt) # 
E5 <- as.factor(data$ethignor) #  

# To simplify
E1<-recode(E1,"c(1,2)=1;c(3)=2;c(4,5)=3")
E2<-recode(E2,"c(1,2)=1;c(3)=2;c(4,5)=3")
E5<-recode(E5,"c(1,2)=1;c(3)=2;c(4,5)=3")

############### Covariates ############################
edu<-as.factor(data$edu)         # <= 13 is  1, > 13 is 2
race <- as.factor(data$race)    # 1 = other-black, 2 = white

## create dataframe
y<-as.data.frame(cbind(E1,E2,E5,race,edu))
ytab<-matrix(table(y),27,4)
# ytab is a tabular data structure, one column for each stratum given by combinations of covariate levels
remove(dati)
rc<-2  # num. categories for each response

labelrisp<-c("R1","R2","R3")
labelfac<-c("G","C")       
str<-c(2,2)   

#########################    --- use package hmmm ---
# Refer to Colombi, R., Giordano, S., & Cazzaro, M. (2014).
# hmmm: An R Package for Hierarchical Multinomial Marginal Models. 
# Journal of Statistical Software, 59(11), 1–25. https://doi.org/10.18637/jss.v059.i11

marg <- marg.list(c("l-m-m", "m-l-m", "m-m-l", "l-l-m", "l-m-l", "m-l-l", "l-l-l"), mflag="m")

# model for observed responses

# Formula<-list(R1=~LA*G+LA*C+LA*R1,R2=~LB*G+LB*C+LB*R2,R3=~LC*G+LC*C+LC*R3,
                     # R1.R2=~LA*LB*R1.R2,
                     # R1.R3=~LA*LC*R1.R3,
                     # R2.R3=~LB*LC*R2.R3) # parallel covariate effect, homogeneous association

# parallel covariate effect, association to be specified through INDIP in hmmm.model.T
Formula<-list(R1=~LA*G+LA*C+LA*R1,R2=~LB*G+LB*C+LB*R2,R3=~LC*G+LC*C+LC*R3) 

modelobs<-hmmm.model.T(marg,lev=c(3,3,3), names=labelrisp, strata=str,
                       Formula=Formula,
                       fnames=labelfac, replace=TRUE,
                       INDIP=list("CST4","CST4","CST4"), 
                       UNC="TUTZ1", MOD="A"  
                        )

# model for latent components

Flatnoeff<-list(LA=~1,LB=~1,LC=~1,LA.LB="zero",LA.LC="zero",LB.LC="zero",LA.LB.LC="zero") 

modellat<-hmmm.model.X(marg=marg,lev=c(2,2,2),names=c("LA","LB","LC"),strata=str,Formula=Flatnoeff,
                       fnames=labelfac)

# fit model 
fit<-Uhmmm.MLfit(ytab, 
                 modelobs=modelobs, 
                 modellat=modellat,
                 method="BFGS", 
                 method.inv="Newton",
                 numeric=FALSE, 
                 print.level=1, 
                 finalHessian=TRUE, 
                 iterlim=30, 
                 NMiter=30,
                 recalc=TRUE,
                 tol=1e-12,
                 startmth="NM", 
                 control=control   
                 
)

# Estimates and Standard Errors
fit$vecpar

## The function is currently defined as
function (marg, lev, names, Formula = NULL, strata = 1, fnames = NULL, 
    cocacontr = NULL, ncocacontr = NULL, replace = TRUE, INDIP = "CST", 
    UNC = "TUTZ2", MOD = "A") 
{
    if (length(UNC) == 1) {
        UNC <- rep(UNC, length(lev))
    }
    labelfac <- fnames
    labellat <- paste("L", LETTERS[1:length(names)], sep = "")
    model <- hmmm.model(marg = marg, lev = lev, names = names)
    descr <- summary(model, printflag = FALSE)
    label <- names
    intcode <- descr[, 1]
    intnames <- descr[, 2]
    f <- Formula
    num.lat <- length(labellat)
    label1 <- paste(labellat, "1", sep = "")
    llcst <- paste(rep(label1, each = num.lat), rep(label1, num.lat), 
        sep = ":")
    llcst <- matrix(llcst, num.lat, num.lat)
    llcst <- vech(as.matrix(llcst[-1, -num.lat]))
    if (is.null(Formula)) {
        f <- list()
        for (i in 1:length(label)) {
            f[[i]] <- as.formula((paste("~", paste(paste(labellat[i], 
                labelfac, sep = "*"), collapse = "+"), "+", labellat[i], 
                "*", label[i], sep = "")))
        }
    }
    num.lat <- length(labellat)
    ll <- paste(rep(labellat, each = num.lat), rep(labellat, 
        num.lat), sep = "*")
    LL <- matrix(ll, num.lat, num.lat)
    names.lint <- vech(as.matrix(LL[-1, -num.lat]))
    length.lint <- length(label) + choose(length(label), 2)
    if ((length(Formula) == length(label)) & (length(label) < 
        length.lint)) {
        for (i in (1 + length(label)):length.lint) {
            ii <- i - length(label)
            if (INDIP[[ii]] == "IND") {
                f[[i]] = "zero"
            }
            else {
                f[[i]] <- as.formula((paste("~", intnames[i], 
                  ":", "(", names.lint[ii], ")", sep = "")))
            }
        }
    }
    if (length.lint < length(intnames)) {
        for (i in (1 + length.lint):length(intnames)) {
            f[[i]] = "zero"
        }
    }
    names(f) <- intnames
    fnames <- c(fnames, labellat)
    os <- sum(strata - 1)
    strataold <- strata
    strata <- c(strata, rep(2, length(names)))
    str <- prod(strata)
    if (strata[1] == 1) {
        strata = strata[-1]
    }
    X <- create.XMAT(model, Formula = f, strata = strata, fnames = fnames, 
        cocacontr = cocacontr, ncocacontr = ncocacontr, replace = FALSE)
    XFULL <- X
    nvar <- length(lev)
    l <- length(X)
    c <- dim(X[[l]])[2]
    px <- 0
    npar <- 0
    Xnames = "zero"
    npar <- as.numeric(descr[, "npar"])
    for (i in 1:nvar) {
        if (UNC[[i]] == "TUTZ0") {
            k <- lev[i]
            llg <- rep(-1, k - 1)
            if (floor(k/2) == k/2) {
                llg[1:(k/2 - 1)] <- 1
                llg[k/2] = 0
            }
            if (floor(k/2) != k/2) {
                M <- (k + 1)/2
                llg[1:(M - 1)] <- 1
            }
            idn <- grep(labellat[i], colnames(X[[i]]))
            ncol <- dim(X[[i]])[2] - length(idn) + 1
            X[[i]] <- X[[i]][, -((ncol - lev[i] + 3):ncol)]
            X[[i]][, 1] <- (1 - X[[i]][, 2]) * llg
            X[[i]][, 3:(ncol - lev[i] + 2)] <- X[[i]][, 3:(ncol - 
                lev[i] + 2)] * X[[i]][, 1]
            LA <- X[[i]][, 2]
            X[[i]][LA == 0, ] <- X[[i]][LA == 0, ] + X[[i]][LA == 
                1, ]
            if (MOD == "B") {
                X[[i]][, (dim(X[[i]])[2] - k + 2):dim(X[[i]])[2]] <- X[[i]][, 
                  (dim(X[[i]])[2] - k + 2):dim(X[[i]])[2]] * 
                  LA
            }
            colnames(X[[i]]) <- c("phi1", colnames(X[[i]])[-1])
        }
        if (UNC[[i]] == "TUTZ1") {
            k <- lev[i]
            if (floor(k/2) == k/2) {
                m <- k/2
                llg <- -((1:(k - 1)) - m)
            }
            if (floor(k/2) != k/2) {
                m <- floor(k/2) + 1
                llg <- -((1:(k - 1)) - m + 0.5)
            }
            idn <- grep(labellat[i], colnames(X[[i]]))
            ncol <- dim(X[[i]])[2] - length(idn) + 1
            X[[i]] <- X[[i]][, -((ncol - lev[i] + 3):ncol)]
            X[[i]][, 1] <- (1 - X[[i]][, 2]) * llg
            X[[i]][, 3:(ncol - lev[i] + 2)] <- X[[i]][, 3:(ncol - 
                lev[i] + 2)] * X[[i]][, 1]
            LA <- X[[i]][, 2]
            X[[i]][LA == 0, ] <- X[[i]][LA == 0, ] + X[[i]][LA == 
                1, ]
            if (MOD == "B") {
                X[[i]][, (dim(X[[i]])[2] - k + 2):dim(X[[i]])[2]] <- X[[i]][, 
                  (dim(X[[i]])[2] - k + 2):dim(X[[i]])[2]] * 
                  LA
            }
            colnames(X[[i]]) <- c("phi1", colnames(X[[i]])[-1])
        }
        if (UNC[[i]] == "TUTZ2") {
            m <- lev[i]
            x <- 1:m
            p <- 6 * (x) * (m + 1 - x)/(m + 2)/(m + 1)/(m)
            llg <- log(p[-1]/p[-length(p)])
            idn <- grep(labellat[i], colnames(X[[i]]))
            ncol <- dim(X[[i]])[2] - length(idn) + 1
            X[[i]] <- X[[i]][, -((ncol - lev[i] + 3):ncol)]
            X[[i]][, 1] <- (1 - X[[i]][, 2]) * llg
            X[[i]][, 3:(ncol - lev[i] + 2)] <- X[[i]][, 3:(ncol - 
                lev[i] + 2)] * X[[i]][, 1]
            LA <- X[[i]][, 2]
            X[[i]][LA == 0, ] <- X[[i]][LA == 0, ] + X[[i]][LA == 
                1, ]
            if (MOD == "B") {
                X[[i]][, (dim(X[[i]])[2] - k + 2):dim(X[[i]])[2]] <- X[[i]][, 
                  (dim(X[[i]])[2] - k + 2):dim(X[[i]])[2]] * 
                  LA
            }
            colnames(X[[i]]) <- c("phi1", colnames(X[[i]])[-1])
        }
        px <- c(px, dim(X[[i]])[2])
        Xnames <- c(Xnames, colnames(X[[i]]))
    }
    if (length(names) > 1) {
        if (length(Formula) == length(label)) {
            for (i in (nvar + 1):length.lint) {
                ii <- i - length(label)
                if (f[[i]] != "zero") {
                  c <- dim(X[[i]])[1]
                  c1 <- c/str
                  if (INDIP[[ii]] == "CST4") {
                    idnLA1LB1 <- grep(llcst[ii], colnames(X[[i]]))
                    splitlab <- strsplit(llcst[ii], ":")
                    idnLA1 <- grep(splitlab[[1]][1], colnames(X[[i]]))
                    idnLB1 <- grep(splitlab[[1]][2], colnames(X[[i]]))
                    idnLA2 <- grep(gsub("1", "2", splitlab[[1]][1]), 
                      colnames(X[[i]]))
                    X[[i]] <- cbind(rowSums(X[[i]][, idnLA1]) - 
                      rowSums(X[[i]][, idnLB1]) + rowSums(X[[i]][, 
                      idnLA1LB1]), rowSums(X[[i]][, idnLA2]) - 
                      rowSums(X[[i]][, idnLA1LB1]), rowSums(X[[i]][, 
                      idnLB1]) - 2 * rowSums(X[[i]][, idnLA1LB1]), 
                      rowSums(X[[i]][, idnLA1LB1]))
                    colnames(X[[i]]) <- c("ass00", "ass10", "ass01", 
                      "ass11")
                  }
                  if (INDIP[[ii]] == "CST3") {
                    idnLA1LB1 <- grep(llcst[ii], colnames(X[[i]]))
                    splitlab <- strsplit(llcst[ii], ":")
                    idnLA1 <- grep(splitlab[[1]][1], colnames(X[[i]]))
                    idnLB1 <- grep(splitlab[[1]][2], colnames(X[[i]]))
                    idnLA2 <- grep(gsub("1", "2", splitlab[[1]][1]), 
                      colnames(X[[i]]))
                    X[[i]] <- cbind(rowSums(X[[i]][, idnLA1]) - 
                      rowSums(X[[i]][, idnLB1]) + rowSums(X[[i]][, 
                      idnLA1LB1]), (rowSums(X[[i]][, idnLA2]) - 
                      rowSums(X[[i]][, idnLA1LB1]) + rowSums(X[[i]][, 
                      idnLB1]) - 2 * rowSums(X[[i]][, idnLA1LB1])), 
                      rowSums(X[[i]][, idnLA1LB1]))
                    colnames(X[[i]]) <- c("ass00", "ass10-01", 
                      "ass11")
                  }
                  if (INDIP[[ii]] == "CST2") {
                    idnLA1LB1 <- grep(llcst[ii], colnames(X[[i]]))
                    X[[i]] <- cbind(X[[i]][, 1] - rowSums(X[[i]][, 
                      idnLA1LB1]), rowSums(X[[i]][, idnLA1LB1]))
                    colnames(X[[i]]) <- c("assU", "assA")
                  }
                  if (INDIP[[ii]] == "CST1") {
                    idnLA1LB1 <- grep("LA1:LB1", colnames(X[[i]]))
                    X[[i]] <- matrix(rowSums(X[[i]][, idnLA1LB1]))
                    colnames(X[[i]]) <- "assA"
                  }
                  if (class(INDIP[[ii]]) == "formula") {
                    print("associazion formula used")
                    va <- intcode[i]
                    va <- as.numeric(unlist(strsplit(va, split = character(0))))
                    lev1 <- lev[va[1]] - 1
                    lev2 <- lev[va[2]] - 1
                    x <- factor(rep(1:lev1, lev2))
                    xx <- factor(rep(1:lev2, each = lev1))
                    da <- data.frame(x, xx)
                    colnames(da) <- c("R", "C")
                    disegno <- model.matrix(INDIP[[ii]], da)
                    X[[i]] <- matrix(rep(c(0, 1), c(str - prod(strataold), 
                      prod(strataold)))) \%x\% disegno
                    colnames(X[[i]]) <- colnames(disegno)
                  }
                  px <- c(px, dim(X[[i]])[2])
                  Xnames <- c(Xnames, colnames(X[[i]]))
                }
            }
        }
        else {
            lbl <- paste(labellat, "1", sep = "")
            ll <- paste(rep(lbl, each = num.lat), rep(lbl, num.lat), 
                sep = ":")
            LL <- matrix(ll, num.lat, num.lat)
            names.lint2 <- vech(as.matrix(LL[-1, -num.lat]))
            for (i in (nvar + 1):length.lint) {
                ii <- i - length(label)
                if (any(X[[i]] != "zero")) {
                  idn <- grep(names.lint2[ii], colnames(X[[i]]))
                  npro <- colnames(X[[i]])[idn]
                  X[[i]] <- X[[i]][, idn]
                  if (!is.matrix(X[[i]])) {
                    X[[i]] <- matrix(X[[i]])
                    colnames(X[[i]]) <- npro
                  }
                  px <- c(px, dim(X[[i]])[2])
                  Xnames <- c(Xnames, colnames(X[[i]]))
                }
            }
        }
    }
    px <- px[-1]
    Xnames <- Xnames[-1]
    XX <- matrix(0, sum(npar) * prod(strata), sum(px))
    or <- rep(rep((1:dim(descr)[1]), times = npar), prod(strata))
    pstart <- 0
    for (i in 1:dim(descr)[1]) {
        if (f[[i]] != "zero") {
            XX[or == i, (pstart + 1):(pstart + px[i])] <- X[[i]]
            pstart <- pstart + px[i]
        }
    }
    colnames(XX) <- Xnames
    if (replace) {
        model <- hmmm.model(marg = marg, lev = lev, names = names, 
            strata = str, X = XX)
        model$levcov <- strataold
        model$fnames <- fnames
        model$marglist <- marg
        model
    }
    else {
        list(X = X, XFULL = XFULL, XX = XX, Formula = Formula)
    }
  }
}
% 
