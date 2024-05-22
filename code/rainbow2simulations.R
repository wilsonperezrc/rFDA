
foutliers2 <-function(data, method = c("robMah", "lrt", "depth.trim", "depth.pond", "HUoutliers"),
           dfunc = depth.mode, nb = 200, suav = 0.05, trim = 0.1, order = 2, lambda = 3.29, ...)
  {
    #browser()
    method = match.arg(method)
    if (method == "lrt"){
      output = outliers.lrt(data = data, nb = nb, suav = suav, trim = trim, ...)
    }
    if (method == "depth.trim"){
      output = outliers_depth_trim(data = data, dfunc = dfunc, nb = nb, suav = suav, trim = trim, ...)
    }
    if (method == "depth.pond"){
      output = outliers.depth.pond(data = data, dfunc = dfunc, nb = nb, suav = suav, ...)
    }
    if (method == "HUoutliers"){
      k = rapca2(data$y, order = order)
      result = matrix(NA, ncol(data$y), 1)
      for(i in 1:ncol(data$y)){
        result[i,] = sum((data$y[,i] - k$basis %*% as.matrix(k$coef[i,]))^2)
      }
      s = median(result)
      crit = s+sqrt(s) * lambda
      if(is.null(colnames(data$y)))
      {
        stop("Please assign column name for the data matrix.")
      }
      out = colnames(data$y)[which(ifelse(result <= crit,1,0) == 0)]
      output = list(outliers = out)
    }
    if (method == "robMah"){
      sco = PCAproj(t(data$y))$scores
      rownames(sco) = as.numeric(colnames(data$y))
      s = cbind(sco, rep(1, ncol(data$y)))
      output = robout(s, 1, "mcd")
    }
    return(output)
  }

####################
rapca2<- function (x, FUN = Qn, order = 4, mean = TRUE)
{
  #browser()
  if (order < 1)
    stop("Order must be positive.")
  X <- t(x)
  n <- nrow(X)
  p <- ncol(X)
  if (mean) {
    med <- colMeans(X)
    xx <- sweep(X, 2, med)
  }
  else xx <- X
  tmp <- La.svd(xx)
  eigen_value = tmp$d^2
  print("Suma de la variabilidad por componentes: ")
  variab<-eigen_value/sum(eigen_value)
  print(variab)
  
  eigen_value = tmp$d^2  # -- Aqui se calcula los eigenValues
  r = sum(tmp$d > (max(n, p) * max(tmp$d) * 1e-12))
  P <- t(tmp$vt)[, 1:r]
  tmp2 <- rstep2(t(xx %*% P), order = order, r = tmp$r, mean = mean)
  tmp <- P %*% tmp2$basis
  if (mean) {
    med <- c(med + tmp[, 1])
    xx <- sweep(X, 2, med)
    basis <- cbind(med, tmp[, (1:order) + 1])
    coef <- cbind(rep(1, n), xx %*% basis[, -1])
  }
  else {
    basis <- tmp
    coef <- xx %*% basis
  }
  return(list(basis = basis, coeff = coef, X = xx))
}



rstep2 <- function (x, FUN = Qn, order = 4, r = matrix.rank(x), mean = TRUE)
{
  #browser()
  if (order < 1)
    stop("Order must be positive.")
  X <- t(x)
  p <- ncol(X)
  n <- nrow(X)
  p1 <- min(order, r, floor(n/2))
  S <- numeric(p1)
  Bnorm <- numeric(n)
  V <- eig <- matrix(0, p, p1)
  Transfo <- diag(p)
  if (mean) {
    med <- L1median2(X, method = "hoss")
    xxx <- xx <- sweep(X, 2, med)
  }
  else xxx <- xx <- X
  for (l in 1:p1) {
    B <- xxx
    for (i in 1:n) Bnorm[i] <- norm(B[i, ], 2)
    Bnormr <- Bnorm[Bnorm > 1e-12]
    B <- B[Bnorm > 1e-12, ]
    A <- diag(1/Bnormr) %*% B
    Y <- xxx %*% t(A)
    s <- colQn(Y)
    j <- order(s, decreasing = TRUE)[1]
    S[l] <- s[j]
    V[l:p, l] <- A[j, ]
    Base <- diag(p - l + 1)
    ndiff <- norm(Base[, 1] - V[l:p, l], Inf)
    if (ndiff > 1e-12) {
      if (sum(V[l:p, l] * Base[, 1]) < 0)
        V[l:p, l] <- -V[l:p, l]
      u <- matrix(Base[, 1] - V[l:p, l], ncol = 1) / c(norm(Base[,
                                                                 1] - V[l:p, l]))
      U <- Base - 2 * repmat(t(u) %*% Base, p - l + 1,
                             1) * repmat(u, 1, p - l + 1)
    }
    else U <- Base
    eig[, l] <- Transfo %*% V[, l]
    if (l < p1) {
      Edge <- diag(p)
      Edge[l:p, l:p] <- U
      Transfo <- Transfo %*% Edge
      xxx <- xxx %*% U
      xxx <- as.matrix(xxx[, -1])
    }
  }
  coef <- xx %*% eig
  if (mean) {
    basis <- cbind(med, eig)
    coef <- cbind(rep(1, n), coef)
  }
  else basis <- eig
  return(list(basis = basis, coeff = coef, X = xx))
}





#####################################

L1median2 <- function (X, tol = 1e-06, maxstep = 200, na.rm = TRUE, method = c("hossjercroux",
                                                                               "coordinate"))
{
  method <- match.arg(method)
  if (method == "coordinate")
    return(apply(X, 2, median.default, na.rm = na.rm))
  else return(hossjercroux(X, tol = tol, maxstep = maxstep,
                           na.rm = na.rm))
}

colQn <- function (Z)
{
  return(apply(Z, 2, Qn))
}

Qn <- function (x)
{
  n <- length(x)
  diffs <- outer(x, x, "-")
  diffs <- diffs[!lower.tri(diffs, diag = TRUE)]
  qn <- 2.2219 * quantile(abs(diffs), 0.25)
  if (n == 2)
    dn <- 0.399
  else if (n == 3)
    dn <- 0.994
  else if (n == 4)
    dn <- 0.512
  else if (n == 5)
    dn <- 0.844
  else if (n == 6)
    dn <- 0.611
  else if (n == 7)
    dn <- 0.857
  else if (n == 8)
    dn <- 0.669
  else if (n == 9)
    dn <- 0.872
  else if (n %% 2 == 1)
    dn <- n / (n + 1.4)
  else dn <- n / (n + 3.8)
  return(dn * qn)
}

repmat = function (A, m, n = m)
{
  A <- as.matrix(A)
  tmp <- matrix(rep(t(A), m), nrow = m * nrow(A), byrow = TRUE)
  return(matrix(rep(tmp, n), ncol = n * ncol(tmp)))
}

norm <- function (A, p = 2)
{
  A <- as.matrix(A)
  if (min(dim(A)) == 1)
    A <- t(A)
  if (p == 1)
    return(as.matrix(max(colSums(abs(A)))))
  else if (p == 2) {
    A.sv <- La.svd(A)$d
    return(as.matrix(max(A.sv)))
  }
  else if (p > 1e+09)
    return(as.matrix(max(rowSums(abs(A)))))
  else stop("Unknown norm")
}

mrobj <- function (X, m)
{
  return(sum(norme(sweep(X, 2, m))))
}

norme <- function (X)
{
  return(sqrt(rowSums(X^2, na.rm = TRUE)))
}

hossjercroux <- function (X, tol = 1e-06, maxstep = 100, na.rm = TRUE)
{
  n <- nrow(X)
  p <- ncol(X)
  m = apply(X, 2, median.default, na.rm = na.rm)
  hctol <- max(1, min(abs(m), na.rm = na.rm)) * tol
  for (k in 1:maxstep) {
    mold <- m
    XX <- sweep(X, 2, m)
    dx <- norme(XX)
    if (min(abs(dx)) > tol)
      w <- 1/dx
    else {
      w <- rep(0, n)
      w[dx > tol] <- 1/dx[dx > tol]
    }
    delta <- colSums(XX * repmat(w/sum(w), 1, p), na.rm = na.rm)
    nd <- sqrt(sum(delta^2))
    maxhalf <- ifelse(nd < hctol, 0, log2(nd/hctol))
    m <- mold + delta
    nstep <- 0
    oldmobj <- mrobj(X, mold)
    while ((mrobj(X, m) > oldmobj) & (nstep <= maxhalf)) {
      nstep <- nstep + 1
      m <- mold + delta / (2^nstep)
    }
    if (nstep > maxhalf)
      return(mold)
  }
  return(mold)
}

repmat = function (A, m, n = m)
{
  A <- as.matrix(A)
  tmp <- matrix(rep(t(A), m), nrow = m * nrow(A), byrow = TRUE)
  return(matrix(rep(tmp, n), ncol = n * ncol(tmp)))
}