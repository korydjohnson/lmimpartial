# impartialReg ---------------------------------------------------------------

#' @name impartialReg
#' @title Creating Impartial Regression Estimates
#'
#' @description The function lm_impartial takes a data frame or matrix as well
#'   as a list specifying which columns correspond to various covariate groups:
#'   response, sensitive, legitimate, and suspect. It then produces estimates of
#'   the response which are impartial with respect to the sensitive covariates.
#'   Note that W can include estimates from other models which are not
#'   constrained to be fair or impartial.
#'
#' @details The order of the variables in the model is very important:
#'   lm(Y~S+X+W) vs lm(Y~S+W+X). R will keep the first variables and drop the
#'   latter if they are collinear. Therefore, if W = S + X, X is kept in the
#'   first model while W is kept in the second. This can have a large impact on
#'   the coefficient for S. In fairness cases, it is better to keep X, since
#'   this is treated as legitimate variability, while W wont be. Therefore use
#'   Y~S+X+W. For the predict function the data from which to predict does not
#'   include an intercept column. Before computing predictions on new data, a
#'   check for equivalent model matrices is made.
#'
#' @param data matrix or data.frame of covariates.
#' @param colList A list of named column specifications. List elements can
#'   either be column numbers or column names. The data frame merely needs to be
#'   able to be subsetted according to an element. colList includes some subset
#'   of c("Y","S","X","W"). Note that an element named Y is always required.
#' @param object an object of class impartialLM; expected to be the list output
#'   from the lm_impartial function.
#' @param newdata an optional data frame in which to look for variables with
#'   which to predict. If omitted, the fitted values are used.
#' @return A list which includes the following components:
#'   \item{coeffWS}{coefficients from W~S+X; included in final output only if W
#'   is included.} \item{coeffX}{coefficients from Y~S+X+W.} \item{fEst}{fair
#'   estimates given colList.} \item{colList}{covariate groupings.}
#'   \item{dataNames}{column names after transforming data; used for checking in
#'   predict.impartialLM.}
#'
#' @examples
#'   data("CO2")
#'   colList = list("Y" = 5, "S" = 3, "X" = c(1,2,4))
#'   lm_impartial(CO2, colList)
#' @importFrom stats lm model.matrix pt qt resid var sd .lm.fit

#' @export
lm_impartial = function(data, colList) {
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

# Core function: theResponse is either Y or W
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
