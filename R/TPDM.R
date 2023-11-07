
#' Wrapper function
#'
#' Handles all steps for estimation of EPI from raw-data:
#' 1) Preprocessing into Frechet-Margins
#' 2) Estimation of TPDM
#' 3) Calculation of Principal Components
#' 4) Estimation of EPI
#' @param X A t x n dimensional Data-matrix with t: Number of time steps and n: Number of grid points/stations
#' @param Y Optional: Sames as X but for second variable: If Y!=NULL, cross-TPDM instead of TPDM and SVD instead of PCA is computed
#' @param q Threshold for computation of TPDM. Only data above the 'q'-quantile will be used for estimation. Choose such that 0 < q < 1.
#' @param anz_cores Number of cores for parallel computing (default: 5)
#' @param clust Optional_ Uf clust = NULL, no declustering is performed. Else, declustering according to cluster-length 'clust'
#' @param m Numeric vector: Containing the principal components/expansion coefficients (in case of Y!=NULL) from which the EPI shall be computed (default: modes = c(1:10), calculates the EPI on first ten principle Components)
#' @param thr_EPI Only if Y!=NULL: Threshold for computation of TEPI. Expansion-coefficients that exceed the 'q'-quantile will be used for estimation. Choose such that 0 < q < 1.
#' @return In case of Y =NULL: A list containing:
#' \itemize{
#'     \item basis: The Eigenvectors of TPDM
#'     \item pc: The principal components of TPDM
#'     \item extremal.basis: The Eigenvectors of TPDM but transformed in positive reals with \link{trans}
#'     \item EPI: Extremal pattern index
#' }
#'
#' In case of Y !=NULL: A list containing:
#' \itemize{
#'     \item U, V: The left- and right singular Vectors of cross-TPDM
#'     \item extr.U, extr.V: The left- and right singular vectors of cross-TPDM, but transformed in positive reals with \link{trans}
#'     \item pcU, pcV: The left- and right expansion coefficients of cross-TPDM
#'     \item EPI: Extremal pattern index
#'     \item TEPI: Threshold-based extremal pattern index
#'}
#' @references Szemkus & Friederichs 2023
#' @examples
#' data    <- precipGER
#' \donttest{
#' result  <- wrapper.EPI(data$pr, m = 1:50)
#'
#' rbPal <- colorRampPalette(c('blue', 'white','red'))
#' Col <- rbPal(10)[as.numeric(cut(result$basis[,2],breaks = 10))]
#' plot(data$lat, data$lon,col=Col)
#' plot(data$date, result$EPI, type='l')}
#' @export
wrapper.EPI <- function(X,
                        Y         = NULL,
                        q         = .98,
                        anz_cores = 1,
                        clust     = NULL,
                        m         = 1:10,
                        thr_EPI   = NULL
                        ){
  # 1) #Preprocessing into Frechet
  X       <- to.alpha.2(X)
  # 2) # Estimation of TPDM
  Sigma   <- est.tpdm(X,Y=NULL,anz_cores,clust ,q=q)
  if(is.null(Y)){   # Estimate TPDM and compute PCA
    # 3) # PCA
    res.pca <- pca.tpdm(Sigma, X) #Performing Principal Component Analysis
    # 4) # Estimate EPI
    res.pca$EPI <- compute.EPI(res.pca, m = m)
    return(res.pca)
  }else{        # Estimate cross-TPDM and compunte SVD
    # 3) # SVD
    res.svd <- svd.tpdm(Sigma, X) #Performing Principal Component Analysis
    # 4) # Estimate EPI
    res.svd$EPI  <- compute.EPI(res.svd, m = m)
    if(!is.null(thr_EPI)){
       res.svd$TEPI <- compute.EPI(res.svd, m = m, q=thr_EPI)
    }
    return(res.svd)
  }
}


#' Estimation of EPI
#'
#' Estimates the extremal pattern index (EPI) from either the 'm' principle components after a PCA
#' or left- and right expansion coefficients after an SVD. In case of a SVD,
#' the threshold-based EPI (TEPI) can optionally be calculated.
#' @param coeff A list, containing the t x n dimensional principle components/expansion coefficients of TPDM.
#' Can also be output of function 'est.tpdm'.
#' @param q Optional: A threshold for computation of TEPI
#' @param m numeric vector: Containing the Principle Components from which EPI shall be computed (e.g. with modes = c(1:10), the EPI is calculated on first ten principle components)
#' @return An array of length t, containing EPI. TEPI is computed if if q > 0.
#' @references Szemkus & Friederichs (2023)
#' @details Given the first 'm' modes of principle components u and eigenvalues after a PCA, the EPI is given as:
#' \deqn{EPI_t^{u} = \sqrt{\sum_{k=1}^m (u_{t,k}^2)/\sum_{j=1}^m e_j}.}
#' Given the first 'm' modes of expansion coefficients u and v and singular values e after a SVD, the EPI and TEPI are given as:
#' \deqn{EPI_t^{u, v} = \sqrt{\sum_{k=1}^m (u_{t,k}^2 + v_{t,k}^2)/\sum_{j=1}^m e_j}.}
#' \deqn{TEPI_t^{u, v} = \sqrt{(\sum_{k=1}^m (u_{t,k}^2 + v_{t,k}^2)/\sum_{j=1}^m e_j) |_{(|u_{t,k}| > q_u , |v{t,k}| > q_v)}}.}
#' @importFrom stats quantile
#' @examples
#' data    <- precipGER
#' \donttest{
#' data.alpha2  <- to.alpha.2(data$pr)
#' Sigma   <- est.tpdm(data.alpha2,anz_cores =1)
#' res.pca <- pca.tpdm(Sigma, data.alpha2)
#' EPI <- compute.EPI(res.pca, m = 1:10)
#'
#' plot(data$date, EPI, type='l')}
#' @export
compute.EPI <- function(coeff,
                        m = 1:10,
                        q = 0.98){

  if(!is.list(coeff)){
    stop("ERROR: 'coeff' must be a list")
  }
  names <- attributes(coeff)$names
  if("pc" %in% names & "val" %in% names){
    pc     <- coeff$pc
    val    <- coeff$val
    pc.sec <- NULL
  }else if("pcU" %in% names & "pcV" %in% names & "val" %in% names){
    pc     <- coeff$pcU
    pc.sec <- coeff$pcV
    val    <- coeff$val
  }else{
    message("ERROR: 'coeff' must be the output of pca.tpdm or svd.tpdm")
  }
  norm_vec <- function(x) sqrt(sum(x))

  ## Estimation of EPI from 'm' of Principal components
  if(is.null(pc.sec)){
    id.thr <- which(abs(pc) < quantile(abs(pc),q, na.rm=TRUE))
    pc[id.thr] <- 0
    EPI   <- apply(pc[,m]^2, c(1), norm_vec)/sqrt(sum(val[m]))
  }else{
    id.thr <- which(abs(pc) < quantile(abs(pc),q, na.rm=TRUE)& abs(pc.sec) < quantile(abs(pc.sec),q, na.rm=TRUE))
    pc[id.thr]     <- 0
    pc.sec[id.thr] <- 0
    EPI  <- apply((pc[,m]^2 + pc.sec[,m]^2), c(1), norm_vec)/sqrt(sum(val[m]))
  }
  return(EPI)
}


#' Principal Component Analysis for TPDM
#'
#' Calculates principal component analysis (PCA) of given TPDM
#' @param Sigma A n x n data array, containing the TPDM, can be output of \link{est.tpdm}.
#' @param data A t x n dimensional, numeric Data-matrix with t: Number of time steps and n: Number of grid points/stations.
#' @return list containing
#' \itemize{
#'     \item pc: The Principal Components of TPDM
#'     \item basis: The Eigenvectors of TPDM
#'     \item extremal.basis: The Eigenvectors of TPDM but transformed in positive reals with \link{trans}
#' }
#' @author Yuing Jiang, Dan Cooley
#' @importFrom Matrix nearPD
#' @references Jiang & Cooley (2020) <doi:10.1175/JCLI-D-19-0413.1>
#' @export
pca.tpdm <- function(Sigma, data) {
  eigen.Sigma   <- eigen(Sigma)
  if (min(eigen.Sigma$values) < 0) {
    Sigma       <- nearPD(Sigma)[[1]]
    eigen.Sigma <- eigen(Sigma)
  }

  U   <- eigen.Sigma$vectors
  val <- eigen.Sigma$values

  if (U[1, 1] < 0) {
    U <- -U
  }
  U.ext <- trans(U)                       # Eigenvectors in the positive orthant

  # estimate PCs component-wise to avoid problems with missing values
  data.inv <- invTrans(data)
  V        <- array(NA, dim= dim(data.inv))
  for ( i in 1 : nrow(data.inv)) {
    id.na <- is.na(data.inv[i, ])
    V[i, ] <-  t(U[!id.na, ]) %*% data.inv[i, !id.na]
  }
  return(list(basis = U, pc = V, extremal.basis = U.ext, val = val))
}



#' Singular Value decomposition for cross-TPDM
#'
#' Calculates singular value decomposition (SVD) of given cross-TPDM
#' @param Sigma A n x n data array, containing the cross-TPDM, can be output of \link{est.tpdm}.
#' @param X A t x n dimensional, numeric Data-matrix with t: Number of time steps and n: Number of grid points/stations.
#' @param Y Same as X but for second variable.
#' @return List containing
#' \itemize{
#'     \item pcU, pcV: The left- and right expansion coefficients of cross-TPDM
#'     \item U, V: The left- and right singular Vectors of cross-TPDM
#'     \item extr.U, extr.V: The left- and right singular vectors of cross-TPDM, but transformed in positive reals with \link{trans}
#' }
#' @export
svd.tpdm <- function(Sigma, X, Y) {
  eigen.Sigma <- svd(Sigma)

  U <- eigen.Sigma$u
  V <- eigen.Sigma$v
  if (U[1, 1] < 0) {
    U <- -U
    V <- -V
  }
  val   <- eigen.Sigma$d
  U.ext <- trans(U)
  V.ext <- trans(V)         # Eigenvectors in the positive orthant
  data.invX <- invTrans(X)
  data.invY <- invTrans(Y)

  # estimate PCs component-wise (to avoid problems with missing values)
  Ex.V     <- array(NA, dim= dim(data.invX))
  for ( i in 1 : nrow(data.invX)) {
    id.na <- is.na(data.invX[i, ])
    Ex.V[i, ] <-  t(V[!id.na, ]) %*% data.invX[i, !id.na]
  }

  Ex.U     <- array(NA, dim= dim(data.invY))
  for ( i in 1 : nrow(data.invY)) {
    id.na <- is.na(data.invY[i, ])
    Ex.U[i, ] <-  t(U[!id.na, ]) %*% data.invY[i, !id.na]
  }

  list(U = U, V=V, extr.U = U.ext, extr.V = V.ext, pcU = Ex.U, pcV = Ex.V, val = val )
}



#' Probability integral transformation
#'
#'Performs transformation to make all of the margins follow a Frechet distribution with tail-index alpha = 2.
#' @param data A t x n dimensional, numeric Data-matrix with t: Number of time steps and n: Number of grid points/stations
#' @param orig If known: original distribution of data (currently implemented: 'normal' or 'gamma'), else: NULL
#' @return Data-matrix of same dimension as 'data', but in Frechet-margins with tail-index 2
#' @importFrom  MASS fitdistr
#' @importFrom stats pgamma pnorm
#' @export
to.alpha.2 = function(data, orig=NULL){
  ## Transform marginal distribution into standart-Frechet (alpha = 2)
  P     <- ncol(data)
  X_new <- array(NA, dim=dim(data))

  for(i in 1:P){
    x           <- data[,i]
    x[is.na(x)] <- 0

    ## handles gamma-distribution, normal-distribution and rank-transform so far:
    if(is.null(orig)){
      cdf <- rank(x, na.last = "keep") / (sum(!is.na(x)) + 1)      # +1 to avoid cdf = 1 for highest value
    }else if(orig =='gamma'){
      # No 0-values allowed in gamma-distribution
      x[x <= 0] <- 1e-06
      parameter <- fitdistr(x, "gamma")
      cdf = pgamma(x, parameter$es[1], parameter$es[2])
    }else if(orig =='normal'){
      parameter <- fitdistr(x, "normal")
      cdf = pnorm(x, parameter$es[1], parameter$es[2])
    }
    X_new[,i] <- sqrt(1/(-log(cdf)))
  }
  return(X_new)
}


#' Estimation of TPDM
#'
#' Estimation of tail pairwise dependence matrix (TPDM)
#' @param X A t x n dimensional, numeric data-matrix with t: Number of time steps and n: Number of grid points/stations
#' @param Y Optional: Same as X but for second variable. If Y!=NULL, cross-TPDM is computed
#' @param anz_cores Number of cores for parallel computing (default:1); Be careful not to overload your computer!
#' @param q Threshold for computation of TPDM. Only data above the 'q'-quantile will be used for estimation. Choose such that 0<q<1.
#' @param clust Optional: If clust = NULL, no declustering is performed. Else, declustering according to cluster-length 'clust'.
#' @return An n x n matrix, containing the estimate of the TPDM
#' @import parallel
#' @import doParallel
#' @import foreach
#' @references Jiang & Cooley (2020) <doi:10.1175/JCLI-D-19-0413.1>; Szemkus & Friederichs (2023)
#' @details Given a random vector X with components \eqn{x_{t,i}, x_{t,j}} with \eqn{i,j = 1, \ldots, n} and it's radial component \eqn{r_{t,ij} = \sqrt{x_{t,i}^2 + x_{t,j}^2}} and angular components \eqn{w_{t,i} = x_{t,i}/r_{t,ij}} and \eqn{w_{t,j} = x_{t,j}/r_{t,ij}}, the i'th,j'th element of the TPDM is estimated as:
#' \deqn{\hat{\sigma}_{ij} = 2 n_{ij,exc}^{-1} \sum_{t=1}^{n} w_{t,i} w_{t,j} |_{(r_{t,ij} > r_{0,ij})} }.
#' Given two random vectors X and Y with components \eqn{x_{t,i}, y_{t,j}} with \eqn{i,j = 1, \ldots, n}, and it's radial component \eqn{ r_{t,ij} = \sqrt{x_{t,i}^2 + y_{t,j}^2}} and angular components \eqn{ w_{t,i}^x = \frac{x_{t,i}}{r_{t,ij}} ; w_{t,j}^y = \frac{y_{t,j}}{r_{t,ij}}}, the i'th,j'th element of the cross-TPDM is estimated as:
#' \deqn{\hat{\sigma}_{ij} = 2 n^{-1}_{exc} \sum_{t=1}^{n} w^x_{t,i} w^y_{t,j} |_{(r_{t,ij} > r_{0,ij})} }.
#' @examples
#' data    <- precipGER
#' \donttest{
#' data.alpha2       <- to.alpha.2(data$pr)
#' Sigma   <- est.tpdm(data.alpha2,anz_cores =1)}
#' @export
est.tpdm = function(X,Y=NULL,anz_cores=1,clust =NULL ,q=0.98){
  anz_gridp_x = ncol(X)

  cl <- makeCluster(anz_cores) #not to overload your computer
  clusterExport(cl, c("est.row.tpdm", "est.element.tpdm"))
  registerDoParallel(cl)
  i <- NULL
  final_sigma <- foreach(i= 1:anz_gridp_x, .combine=cbind) %dopar% {
    x_ti    = X[,i]
    if(length(Y)>0){
      sigma_i <- est.row.tpdm(x_ti, Y,clust=clust ,q=q)
    }else{
      sigma_i <- est.row.tpdm(x_ti, X,clust=clust ,q=q)
    }
    sigma_i
  }

  stopCluster(cl)

  return(final_sigma)
}



#' Sub-Routine of function \link{est.tpdm}. Calculates one row of the TPDM
#' @rdname est.tpdm
#' @param x Array of length t, where t is the number of time steps
#' @param Y A t x n dimensional, numeric Data-matrix with t: Number of time steps and n: Number of grid points/stations
#' @param q Threshold for computation of TPDM. Only data above the 'q'-quantile will be used for estimation. Choose such that 0<q<1.
#' @param clust Optional: If clust = NULL, no declustering is performed. Else, declustering according to cluster-length 'clust'.
#' @return Array containing the estimate of one row of the TPDM.
#' @export
est.row.tpdm = function(x,Y,clust = NULL ,q =0.98){
  anz_gridp_y <- ncol(Y)
  sigma_i     <- vector()
  for (j in 1:anz_gridp_y){
    y_tj     = Y[,j]
    sigma_i <- append(sigma_i,est.element.tpdm(x ,y_tj, clust,q))
  }
  return(sigma_i)
}



#' Element-wise Estimation of TPDM
#'
#' Sub-Routine of \link{est.row.tpdm}. Calculates one element of the TPDM
#' @rdname est.tpdm
#' @param x Array of length t, where t is the number of time steps
#' @param y Same as x
#' @param q Threshold for computation of TPDM. Only data above the 'q'-quantile will be used for estimation. Choose such that 0<q<1.
#' @param clust Optional: If clust = NULL, no declustering is performed. Else, declustering according to cluster-length 'clust'.
#' @return Value containing the estimate of one element of the TPDM.
#' @importFrom stats quantile
#' @export
est.element.tpdm = function(x,y ,clust=NULL ,q = 0.98){

  if(is.null(clust)){declust =FALSE}else{declust = TRUE}
  r     <- sqrt(x^2 + y^2)
  quant <- quantile(r, q, na.rm=TRUE)
  if(declust ==TRUE){
    ind   <- decls(r, quant, clust)
  }else{
    ind   <- which(r > quant & is.na(r) ==FALSE)
  }
  n        <- length(ind)
  wt_i     <- x[ind]/r[ind]
  wt_j     <- y[ind]/r[ind]
  sigma_ij <- (2/n * sum(wt_i*wt_j, na.rm=TRUE))

  return(sigma_ij)
}





#' transformation function
#'
#' Applies the transformation \eqn{t(x) = \log(1+\exp{(x)})}
#' @param x Real vector
#' @author Yuing Jiang, Dan Cooley
#' @return Real, positive vector, containing the result of transformation function.
#' @details Transformation from real vector in real, positive vector under preservation of Frechet-distribution.
#' @references Cooley & Thibaud (2019) <doi:10.1093/biomet/asz028>
#' @seealso \link{svd.tpdm}, \link{pca.tpdm}
trans <- function(x)
{
  ##because it takes an exponential, this function flakes out if x is too big
  ##hence for big values of x, we return x
  v                <- log(1 + exp(x))
  id               <- which(x < -20)
  v[!is.finite(v)] <- x[!is.finite(v)]
  v[id]            <- exp(x[id])
  return(v)
}



#' Transformation function
#'
#' Applies the inverse transformation \eqn{t^{-1}(v) = \log(\exp{(v)} -1)}
#' @param v Real, positive vector
#' @author Yuing Jiang, Dan Cooley
#' @return Real vector, containing the result of inverse transformation function.
#' @details Transformation from real, positive vector in real vector under preservation of frechet-distribution.
#' @references Cooley & Thibaud (2019) <doi:10.1093/biomet/asz028>
#' @seealso \link{svd.tpdm}, \link{pca.tpdm}
invTrans <- function(v)
{
  ##same trickeration for big values of v
  ##still returns -Inf if v is machine zero
  x <- log(exp(v) - 1)
  x[!is.finite(x) & v > 1 & !is.na(x)] <- v[!is.finite(x) & v > 1 &
                                              !is.na(x)]
  return(x)
}

#' Declustering
#'
#' Declustering routine, which will can be applied on radial component r in estimation of the TPDM. Subroutine of \link{est.tpdm}.
#' @param x Real vector
#' @param th Threshold
#' @param k Cluster length
#' @references Jiang & Cooley (2020) <doi:10.1175/JCLI-D-19-0413.1>
#' @return numeric vector of declustered threshold exceedances
#' @seealso \link{est.tpdm}
#' @importFrom utils tail
#' @author Yuing Jiang, Dan Cooley
decls <- function(x, th, k) {
  ## Ordinary decluster.
  id.big <- which(x > th)
  id.dif <- diff(id.big)
  tick   <- which(id.dif >= k)
  start  <- id.big[c(1, tick + 1)]              # Where a new cluster begins
  end    <- c(id.big[tick], tail(id.big, 1))
  n      <- length(start)
  id.res <- rep(0, n)
  for ( i in 1 : n) {
    temp <- x[start[i] : end[i]]
    id.res[i] <- which(temp == max(temp, na.rm = TRUE))[1] + start[i] - 1
  }
  id.res

}




