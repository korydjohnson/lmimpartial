context("Impartial Linear Models")
library(impartialReg)

# set up ------------------------------------------------------------------
createCrimeData = function() {
  theData = read.csv("/media/kory/DATA/Data/crime/crime.csv",
                     header=F, sep=",", na.strings = "?")
  colNames = read.csv("/media/kory/DATA/Data/crime/colNames.csv",
                      header=T,sep=",")
  colnames(theData) = colNames[,1]
  hasMissing = apply(theData,2, function(x) sum(is.na(x)))>0
  theData = theData[,!hasMissing]
  theData = theData[,-4]
  yCols = 102
  sCols = c(5:8, 24:28)
  wCols = c(15,18,22,23)
  xCols = setdiff(1:102, c(sCols, wCols, yCols))
  lmCols = 103
  # corMat = cor(theData[, sCols], theData)
  # which(apply(corMat, 2, function(x) sum(abs(x))) > 4)  # determine W
  # W covariates: medIncome pctWInvInc  medFamInc  perCapInc
  colList = list(Y = yCols, S=sCols, X=xCols, W=wCols, fullLM=lmCols)

  fullLM = lm(ViolentCrimesPerPop~., data=theData)
  theData[["fullLM"]] = predict(fullLM)
  theData <<- theData
  colList <<- colList
}
createCrimeData()

source("~/Dropbox/Research/FATML/Code/impartial-fnc.R")
# tests -------------------------------------------------------------------

test_that("makeModelData", {
  fm = makeImMod(theData, colList)
  newNames = lapply(makeModelData(theData, fm$colList), colnames)
  expect_identical(fm$dataNames, newNames)
})

test_that("feo", {
  fm = makeImMod(theData, colList)
  est1 = predict(fm, theData)
  est2 = makeImpartialOld(theData[,colList$Y], theData[,colList$S],
                          X=theData[,colList$X], W=theData[,colList$W])
  expect_equal(est1, est2)
})

test_that("seo", {
  colList$X=NULL
  fm = makeImMod(theData, colList)
  newNames = lapply(makeModelData(theData, fm$colList), colnames)
  expect_identical(fm$dataNames, newNames)
  est1 = predict(fm, theData)
  est2 = makeImpartialOld(Y=theData[,colList$Y], S=theData[,colList$S],
                          W=theData[,colList$W])
  expect_equal(est1,est2)
})
