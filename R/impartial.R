# impartialReg ---------------------------------------------------------------

#' @name impartialReg
#' @title Main function for Creating Impartial Regression Estimates
#'
#' @description The function impartialReg takes arguments for the response (y),
#' sensitive covariates (S), legitimate covariates (X), and suspect covariates
#' (W). Given these groupings, it creates impartial estimates according to
#' group fairness. Note that W can include estimates from other models
#' which are not constrained to be fair or impartial.
#'
#' @details The order of the variables in the model is very important:
#' lm(Y~S+X+W) vs lm(Y~S+W+X). R will keep the first variables and drop the
#' latter if they are collinear. Therefore, if W = S + X, X is kept in the first
#' model while W is kept in the second. However, this has a huge impact on the
#' coefficient for S. In fairness cases, it is better to keep X, since this is
#' treated as explaining variability, while W wont be. Therefore use Y~ S+X+W.
#' Put S first because don't want to explain away S (and then drop it). But if
#' explained by an X variable --> probably just means that it should be
#' considered a legitimate variable. This can also cause problems because only
#' need to subtract a component of W if it was actually used in the fitted
#' value. Therefore, it is easier to write the fair estimate constructively or
#' use the original procedure that used hat W in the full model and only removes
#' the S term.
#'
#' @param theData matrix or data.frame of covariates.
#' @return A list which includes the following components: \item{y}{response.}
#'   \item{X}{model matrix from final model.} \item{formula}{final model
#'   formula.}  \item{features}{list of interactions included in formula.}
#'   \item{summary}{if save=TRUE, contains information on each test made by the
#'   algorithm.} \item{time}{run time.} \item{options}{options given to RAI:
#'   alg, searchType, poly, r, startDeg, alpha, omega, m.} \item{subData}{subset
#'   of columns from theData that are used in the final model.}
#'   \item{model}{linear model object using selected model} Summary and predict
#'   methods are provided in order to generate further output and graphics.
#' @examples
#'   data("CO2")
#'   theResponse = CO2$uptake
#'   theData = CO2[ ,-5]
#'   rai_out = rai(theData, theResponse)
#'   summary(rai_out)  # summary information including graphs
#' @importFrom stats lm model.matrix pt qt resid var sd .lm.fit

#' @export
# Create the impartial model. The resulting object includes the appropriate
# model coefficients, in-sample predictions, the colList specifications for
# group membership, and the names of the expanded model matrix. colList is
# included as this obviously needs to stay the same when calling the predict
# function. It includes everything except raw data, which I didn't want to make
# a copy of. The names are included to help model checking. colList includes
# some subset of c("Y","S","X","W"). Note that an element named Y is always
# required.
makeImMod = function(data, colList) {
  data = makeModelData(data, colList)
  dataNames = lapply(data, colnames)
  if (is.null(colList$W)) {  # only Y,S,X
    fm = feoMod(data, theResponse="Y")
    fm$coeffWS = NULL
  } else {
    # fm gives final output, fmw is intermediate model W~S+X
    fmw = feoMod(data, theResponse="W")
    tildeW = with(data, W - S%*%fmw$coeffWS)
    data$X = cbind(data$X,tildeW)
    fm = feoMod(data, theResponse="Y")
    fm$coeffWS = fmw$coeffWS
  }
  fm$colList = colList
  fm$dataNames = dataNames
  class(fm) = "impartialLM"
  fm
}

################# Things to update
# I need a better way of asking for predictions out of sample.
# I think the best way to do this is probably to actually inherit methods from
# standard s3 lm objects. My guess is that if I just store the appropriate
# coefficients, that everything else will work out just fine. This would have the
# benefit of allowing modeling dianostics and such.

########## I actually need to think hard about how to do out of sample prediction.
# The main problem is that my covariates have changed. Given a new observation,
# I need to project the W component off of the S component... I don't know how to do that
# with a single observation
### really doing this is expectation. So the projection is x*cov(X)^-1*cor(X,W).
# I use the training data to compute the cov and cor. Then plug in x from the observation?
###### whoa whoa whoa, thinking about it way to hard. Know what the prediction of W is
# using X in the new model. from that can get the new estimate hat w.

####################### new function
# The new function will produce an S3 object. The important
# method it needs is predict. Either called without new data (giving fitted values)
# or called with new data (giving new predictions). Note I then don't give a test
# option
#
# In order to do these things it needs to save fitted values and coefficients.
# coefficients are needed in both the model w ~ s + x and y ~ s + x + tilde w.
# Instead of giving an option for YisW or whatever, just set
# tildeW = \hat W + resid(W), where resid is from the lm fit for w ~ .
# It saves more information, but not a ridi# fm gives final output, fmw is intermediate model W~S+culous amount more, just to save both
# linear model objects. This only gets inefficient when X is huge. set model=FALSE
# in order to not save the raw data used.

# Inside the function there will be s3 objects. These objects are actually just
# coefficient vectors with an associated predict method. When predict is called
# on the coefficient vector, is gives the appropriate estimate.

# First, a helper function for transforming the data to model matrices. This
# also separates the single data frame into a list with named components. Note
# that empty components are just turned into duplicates of the intercept and
# will be dropped in model fitting.
makeModelData = function(data, colList) {
  if (is.null(colList$Y)) stop("A column specification for Y is required.")
  n = nrow(data)
  dataGroups = c("Y","S","X","W")
  data = lapply(dataGroups,
                function(group)
                  if (is.null(colList[[group]])) {
                    matrix(1, nrow=n, dimnames=list(NULL, "IntSXW"))
                  } else {
                    model.matrix(~.-1, data.frame(data[, colList[[group]], drop=FALSE]))
                  }
  )
  names(data) = dataGroups
  data$S = scale(data$S, scale=FALSE)  # so removing S does not change mean
  data
}

# Core lm function: theResponse with either be y or W depending on how the
# function is called
# coeffWS correspond to W~S+X; included in final output only if W is included
feoMod = function(data, theResponse) {
  fm = list()
  mod = lm(as.formula(paste(theResponse,"~S+X")), data)  # full model
  indS = 1+1:dim(data$S)[2]
  indX = setdiff(1:(1+dim(data$S)[2]+dim(data$X)[2]), indS)  # intercept & X

  if (dim(data[[theResponse]])[2] > 1) {  # coefficients are in a matrix
    fm$coeffWS = mod$coefficients[indS,]
    fm$coeffX = mod$coefficients[indX,]
  } else {  # coefficients are in a vector
    fm$coeffWS = mod$coefficients[indS]
    fm$coeffX = mod$coefficients[indX]
  }
  fm$coeffX[is.na(fm$coeffX)] = 0 # NA for collinear/dropped variables
  fm$coeffWS[is.na(fm$coeffWS)] = 0 # NA for collinear/dropped variables
  fm$fEst = cbind(1,data$X)%*%fm$coeffX
  fm
}

# When called without new data, __predict_im__ just returns the estimates from
# the training model. Note that the data from which to predict does *not*
# include an intercept column. Before computing predictions on new data, a check
# for equivalent model matrices is made.

#' @name impartialReg
#' @export
predict.impartialLM = function(object, newdata=NULL) {
  fm = object
  if (is.null(newdata)) {
    fEst = fm$fEst
  } else {
    if (max(unlist(fm$colList)) > ncol(newdata)) stop("Data & colList mismatch.")
    newdata = makeModelData(newdata, fm$colList)
    newNames = lapply(newdata, colnames)
    if (!identical(newNames, fm$dataNames)) {
      stop("Names for new data do not match those of training data.")
    }
    if (is.null(fm$coeffWS)) {  # model doesn't include W
      fEst = cbind(1, newdata$X)%*%fm$coeffX
    }
    else {
      tildeW = with(newdata, W - S%*%fm$coeffWS)
      fEst = cbind(1, newdata$X, tildeW)%*%fm$coeffX
    }
  }
  fEst
}
