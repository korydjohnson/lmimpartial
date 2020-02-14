# lm_impartial ---------------------------------------------------------------

#' @name lm_impartial
#' @title Creating Impartial Regression Estimates
#'
#' @description The function lm_impartial takes a data frame or matrix as well
#'   as a list specifying which columns correspond to various covariate groups:
#'   response, sensitive, legitimate, and suspect. It then produces estimates of
#'   the response which are impartial with respect to the sensitive covariates.
#'   Note that W can include estimates from other models which are not
#'   constrained to be fair or impartial.
#' @details The order of the variables in the model is very important:
#'   lm(Y~S+X+W) vs lm(Y~S+W+X). R will keep the first variables and drop the
#'   latter if they are collinear. Therefore, if W = S + X, X is kept in the
#'   first model while W is kept in the second. This can have a large impact on
#'   the coefficient for S. In fairness cases, it is better to keep X, since
#'   this is treated as legitimate variability, while W wont be. For the predict
#'   function, the theData from which to predict does not include an intercept
#'   column. Before computing predictions on newdata, a check for equivalent
#'   model matrices is made. This requires column names in newdata to match
#'   those of the original theData. Factor columns also need to include levels for
#'   all categories else equivalent model matrices will not be created.
#' @param theData matrix or data.frame.
#' @param colList A list of named column specifications. List elements can
#'   either be column numbers or column names. The theData frame merely needs to be
#'   able to be subset according to an element of this list. colList includes
#'   some subset of c("Y","S","X","W"). An element named Y is always required,
#'   and S must be included if W is.
#' @param object an object of class impartialLM; expected to be the list output
#'   from the lm_impartial function.
#' @param newdata an optional matrix or data frame in which to look for
#'   variables with which to predict. If omitted, the fitted values are used.
#' @param ... further arguments passed to or from other methods.
#' @return A list which includes the following components:
#'   \item{coeffWS}{coefficients from W~S+X; included in final output only if W
#'   is included.} \item{coeffX}{coefficients from Y~S+X+W.}
#'   \item{imEst}{impartial estimates given colList.} \item{colList}{covariate
#'   groupings.} \item{dataNames}{column names after transforming theData; used for
#'   checking in predict.impartialLM.}
#'
#' @examples
#'   data(swiss)
#'   colList1 = list("Y" = 1, "S" = 5, "X" = setdiff(1:6, c(1,5)))
#'   lmOut1 = lm_impartial(swiss, colList1)
#'   colList2 = list("Y" = 1, "S" = 5, "X" = setdiff(1:6, c(1,4,5)), "W" = 4)
#'   lmOut2 = lm_impartial(swiss, colList2)
#' @importFrom stats lm model.matrix pt qt resid var sd .lm.fit as.formula

#' @export
lm_impartial = function(theData, colList) {
  theData = makeModelData(theData, colList)
  dataNames = lapply(theData, colnames)
  if (is.null(colList$W)) {  # only Y,S,X
    im = feoMod(theData, theResponse="Y")
    im$coeffWS = NULL
  } else {
    imw = feoMod(theData, theResponse="W")  # imw is intermediate model W~S+X
    tildeW = with(theData, W - S%*%imw$coeffWS)
    theData$X = cbind(theData$X,tildeW)
    im = feoMod(theData, theResponse="Y")  # im gives final output
    im$coeffWS = imw$coeffWS
  }
  im$colList = colList
  im$dataNames = dataNames
  class(im) = "lm_impartial"
  im
}

# First, a helper function for transforming the theData to model matrices. This
# also separates the single theData frame into a list with named components. Note
# that empty components are just turned into duplicates of the intercept and
# will be dropped in model fitting.
makeModelData = function(theData, colList) {
  if (is.null(colList$Y)) stop("A column specification for Y is required.")
  if (is.null(colList$S) && !is.null(colList$W)) stop("Cannot have W without S.")
  n = nrow(theData)
  dataGroups = c("Y","S","X","W")
  theData = lapply(dataGroups,
                function(group)
                  if (is.null(colList[[group]])) {
                    matrix(1, nrow=n, dimnames=list(NULL, "IntSXW"))
                  } else {
                    model.matrix(~.-1, data.frame(theData[, colList[[group]], drop=FALSE]))
                  }
  )
  names(theData) = dataGroups
  theData$S = scale(theData$S, scale=FALSE)  # so removing S does not change mean
  theData
}

# Core function: theResponse is either Y or W
feoMod = function(theData, theResponse) {
  im = list()
  mod = lm(as.formula(paste(theResponse,"~S+X")), theData)  # full model
  indS = 1+1:ncol(theData$S)
  indX = setdiff(1:(1+ncol(theData$S)+ncol(theData$X)), indS)  # intercept & X
  if (ncol(theData[[theResponse]]) > 1) {  # coefficients are in a matrix
    im$coeffWS = mod$coefficients[indS,]
    im$coeffX = mod$coefficients[indX,]
  } else {  # coefficients are in a vector
    im$coeffWS = mod$coefficients[indS]
    im$coeffX = mod$coefficients[indX]
  }
  im$coeffX[is.na(im$coeffX)] = 0  # NA for collinear/dropped variables
  im$coeffWS[is.na(im$coeffWS)] = 0  # NA for collinear/dropped variables
  im$imEst = cbind(1,theData$X)%*%im$coeffX
  im
}

#' @name lm_impartial
#' @export
predict.lm_impartial = function(object, newdata=NULL, ...) {
  im = object
  if (is.null(newdata)) {
    imEst = im$imEst
  } else {
    # check that all elements in colList are actually columns
    c1 = max(unlist(im$colList)) > ncol(newdata)  # if string, returns T
    c2 = all(unlist(im$colList) %in% colnames(newdata))
    if (max(unlist(im$colList)) > ncol(newdata)) stop("newdata & colList mismatch.")
    newdata = makeModelData(newdata, im$colList)
    newNames = lapply(newdata, colnames)
    if (!identical(newNames, im$dataNames)) {
      stop("Names for new theData do not match those of training theData.")
    }
    if (is.null(im$coeffWS)) {  # model doesn't include W
      imEst = cbind(1, newdata$X)%*%im$coeffX
    }
    else {
      tildeW = with(newdata, W - S%*%im$coeffWS)
      imEst = cbind(1, newdata$X, tildeW)%*%im$coeffX
    }
  }
  imEst
}
