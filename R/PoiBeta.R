#' @export
PoiBeta <- function(k, a, b, n = 100, d = 0) {
  if (a<=0 || b<=0 || k<=0 || n<=0) {
    stop("All parameters must be positive!")
  }
  #Check if we have vectors and that they are the same length
  kLen <- length(k)
  aLen <- length(a)
  bLen <- length(b)
  if ((kLen==1 && aLen==1 && bLen==1) || (kLen>1 && aLen==1 && bLen==1) || (kLen==1 && aLen>1 && bLen==1) || (kLen==1 && aLen==1 && bLen>1) || (kLen==bLen && aLen==bLen) || (kLen==bLen && aLen==1) || (kLen==aLen && bLen==1) || (kLen==1 && aLen==bLen)) {
    m <- max(c(kLen, aLen, bLen))
    x <- mat.or.vec(m, n)
    if (kLen==1) { k <- rep(k, m) }
    if (aLen==1) { a <- rep(a, m) }
    if (bLen==1) { b <- rep(b, m) }
    for (i in 1:m) {
      #First generate Beta random variables
      y <- rbeta(n, a[i], b[i])
      #Then use for the Poisson intensities
      x[i,] <- rpois(n, k[i]*y)
      if (d>0) { #simulate drop-outs
        mu <- k*a/(a+b)
        x[i,as.logical(rbinom(n, 1, 1-mu/(d+mu)))] <- 0
      }
    }
    return(x)
  }
  else {
    stop("Array lengths of input parameters inconsistent!")
  }
}

#' @export
GeneratePoiBetaSamples <- function(ks, as, bs, mult, nGenes = 100, nCells = 50, d = 0, meanFixed = F) {
  genes <- PoiBeta(ks, as, bs, nCells, d)
  ks2 <- ks
  as2 <- as
  bs2 <- bs
  for (i in 1:nGenes) {
    if (meanFixed) {
      as2[i] <- as[i]*mult[i]
      bs2[i] <- bs[i]*mult[i]
    }
    else {
      if (i<nGenes/3) { ks2[i] <- ks[i]*mult[i] }
      else if (i<2*nGenes/3) { as2[i] <- as[i]*mult[i] }
      else { bs2[i] <- bs[i]*mult[i] }
    }
  }
  genes2 <- PoiBeta(ks2, as2, bs2, nCells, d)
  return(list("sample1" = genes, "sample2" = genes2, "ks"=ks, "as"=as, "bs"=bs, "ks2"=ks2, "as2"=as2, "bs2"=bs2, "mult"=mult, "d"=d))
}

#' @export
#' @importFrom moments moment
PoiBetaMMFit <- function(x) {
  #Calculate the central moments
  mu <- mean(x)
  mu2 <- var(x)
  mu3 <- moments::moment(x, order=3, central=T)
  #Estimate the parameters
  k <- 2*mu - (mu3*(mu^2 + mu2))/(- 2*mu2^2 + mu*mu3)
  a <- (2*mu*(mu^2*mu2 + mu3*mu - mu2^2))/(- mu3*mu^2 + 4*mu*mu2^2 + mu3*mu2)
  b <- -(2*mu2*(mu3 + 2*mu*mu2)*(mu^2*mu2 + mu3*mu - mu2^2))/((- 2*mu2^2 + mu*mu3)*(- mu3*mu^2 + 4*mu*mu2^2 + mu3*mu2))
  return(list("k" = k, "a" = a, "b" = b))
}

#' @export
#' @importFrom hypergeo genhypergeo
#' @importFrom orthopolynom lpochhammer
PoiBetaLogLikelihood <- function(x, k = NaN, a = NaN, b = NaN, nMC = 1e3, mcLim = 1e2) {
  if (is.nan(k) | is.nan(a) | is.nan(b)) {
    res <- PoiBetaMMFit(x)
    if (is.nan(k)) { k <- res$k }
    if (is.nan(a)) { a <- res$a }
    if (is.nan(b)) { b <- res$b }
  }
  n <- max(c(ncol(x), nrow(x)))
  logLike <- rep(n, 0)
  for (i in 1:n) {
    if (x[i]<mcLim & x[i]>0) {
      logLike[i] <- x[i]*log(k) - k + orthopolynom::lpochhammer(a, x[i]) -
          lgamma(x[i]) - orthopolynom::lpochhammer(a + b, x[i]) +
          log(hypergeo::genhypergeo(a, a + b + x[i], k))
    }
    else {
      q <- rbeta(nMC, a, b)
      p <- dpois(x[i], k*q)
      logLike[i] <- log(mean(p))
    }
  }
  return(sum(logLike))
}

