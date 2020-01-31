context("Impartial Linear Models")
library(impartialReg)

# set up ------------------------------------------------------------------
n = 100
p = 1
S1 = matrix(rnorm(n*p), nrow=n)
X1 = matrix(rnorm(n*p), nrow=n)
W1 = matrix(rnorm(n*p), nrow=n)
Y1 = rowSums(cbind(S1, X1, W1)) + rnorm(n)
data1 = data.frame(Y1, S1, X1, W1)  # vector arguments/results
colnames(data1) = paste0(c("y", "s", "x", "w"), 1)
colList1 = list("Y" = 1, "S" = 1+1:p, "X" = 1+p+1:p, "W" = 1+2*p+1:p)
p = 2
S2 = matrix(rnorm(n*p), nrow=n)
X2 = matrix(rnorm(n*p), nrow=n)
W2 = matrix(rnorm(n*p), nrow=n)
Y2 = rowSums(cbind(S2, X2, W2)) + rnorm(n)
data2 = data.frame(Y2, S2, X2, W2)  # matrix arguments/results
colnames(data2) = c("y", t(outer(c("s", "x", "w"), 1:2, paste0)))
colList2 = list("Y" = 1, "S" = 1+1:p, "X" = 1+p+1:p, "W" = 1+2*p+1:p)
S2 = data.frame(s1=rnorm(n), s2=as.factor(sample.int(2, n, replace=TRUE)))
X2 = matrix(rnorm(n*p), nrow=n)
W2 = matrix(rnorm(n*p), nrow=n)
data2.2 = data.frame(Y2, S2, X2, W2)  # S has factors
colnames(data2.2) = c("y", t(outer(c("s", "x", "w"), 1:2, paste0)))

# tests -------------------------------------------------------------------

test_that("makeModelData", {
  # checks S mean zero and appropriate column names
  theData1 = makeModelData(data1, colList1)
  expect_equal(mean(theData1$S), 0)
  expect_equal(lapply(colList1, function(cols) colnames(data1[, cols, drop=FALSE])),
               lapply(theData1, colnames))

  theData2 = makeModelData(data2, colList2)
  expect_equal(as.vector(colSums(theData2$S)), rep(0, 2))
  expect_equal(lapply(colList2, function(cols) colnames(data2[, cols, drop=FALSE])),
               lapply(theData2, colnames))

  theData2.2 = makeModelData(data2.2, colList2)
  expect_equal(as.vector(colSums(theData2.2$S)), rep(0, 3))
  expect_equal(cor(theData2.2$S), cor(model.matrix(~.-1, data=S2)))
})

test_that("feo_YSXW", {   # lm w ~ s x, est w; then y ~ hat w s x, est y
  # vector
  im = lm_impartial(data1, colList1)
  est1 = predict(im)
  est2 = predict(im, data1)
  theData1 = makeModelData(data1, colList1)
  modW = lm(W~S+X, theData1)  # full model
  theData1$estW = theData1$W - theData1$S %*% modW$coefficients[2]
  modY = lm(Y~S+X+estW, theData1)
  est3 = cbind(1, theData1$X, theData1$estW)%*%modY$coefficients[c(1, 3, 4)]
  expect_equal(est1, est2, est3)

  # matrix
  im = lm_impartial(data2, colList2)
  est1 = predict(im)
  est2 = predict(im, data2)
  theData2 = makeModelData(data2, colList2)
  modW = lm(W~S+X, theData2)  # full model
  theData2$estW = theData2$W - theData2$S %*% modW$coefficients[c(2,3),]
  modY = lm(Y~S+X+estW, theData2)
  est3 = cbind(1, theData2$X, theData2$estW)%*%modY$coefficients[c(1, 4:7)]
  expect_equal(est1, est2, est3)

  # factor
  im = lm_impartial(data2.2, colList2)
  est1 = predict(im)
  est2 = predict(im, data2.2)
  theData2.2 = makeModelData(data2.2, colList2)
  modW = lm(W~S+X, theData2.2)  # full model
  theData2.2$estW = theData2.2$W - theData2.2$S %*% rbind(modW$coefficients[c(2,3),],0)
  modY = lm(Y~S+X+estW, theData2.2)
  est3 = cbind(1, theData2.2$X, theData2.2$estW)%*%modY$coefficients[c(1, 5:8)]
  expect_equal(est1, est2, est3)
})

test_that("seo_YSW", {
  # vector
  im = lm_impartial(data1, colList1[c(1,2,4)])
  est1 = predict(im)
  est2 = predict(im, data1)
  theData1 = makeModelData(data1, colList1[c(1,2,4)])
  modW = lm(W~S, theData1)
  theData1$estW = theData1$W - theData1$S %*% modW$coefficients[2]
  modY = lm(Y~S+estW, theData1)
  est3 = cbind(1, theData1$estW)%*%modY$coefficients[c(1, 3)]
  expect_equal(est1, est2, est3)

  # matrix
  im = lm_impartial(data2, colList2[c(1,2,4)])
  est1 = predict(im)
  est2 = predict(im, data2)
  theData2 = makeModelData(data2, colList2[c(1,2,4)])
  modW = lm(W~S, theData2)  # full model
  theData2$estW = theData2$W - theData2$S %*% modW$coefficients[c(2,3),]
  modY = lm(Y~S+estW, theData2)
  est3 = cbind(1, theData2$estW)%*%modY$coefficients[c(1, 4:5)]
  expect_equal(est1, est2, est3)

  # factor
  im = lm_impartial(data2.2, colList2[c(1,2,4)])
  est1 = predict(im)
  est2 = predict(im, data2.2)
  theData2.2 = makeModelData(data2.2, colList2[c(1,2,4)])
  modW = lm(W~S, theData2.2)  # full model
  theData2.2$estW = theData2.2$W - theData2.2$S %*% rbind(modW$coefficients[c(2,3),],0)
  modY = lm(Y~S+estW, theData2.2)
  est3 = cbind(1, theData2.2$estW)%*%modY$coefficients[c(1, 5:6)]
  expect_equal(est1, est2, est3)
})

