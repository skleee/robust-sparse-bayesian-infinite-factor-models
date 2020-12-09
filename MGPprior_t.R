# ARGUMENTS: X: Data matrix (n x p); 
#            nrun: number of iterations;
#            burn: burn-in period;
#            thin: thinning interval;
#            prop: proportion of elements in each column less than epsilon in magnitude cutoff;
#            epsilon: tolerance;
#            kinit: initial value for the number of factors;
#            adapt: logical. Whether or not to adapt number of factors across sampling
#            output: output type, a vector including some of:
#            c("covMean", "covSamples", "factSamples", "sigSamples", "numFactors");
#            verbose: logical. Show progress bar?
#            dump: logical. Save output object during sampling?
#            filename: if dump, filename for output
#            buffer: if dump, frequency of saving
#            augment: additional sampling steps as an expression

generateData = function( num.vars=3, num.factors = num.vars, n=100, seed=NULL, student_t = FALSE, nu = 10, sparse = TRUE) {
  if(is.null(seed)){
    set.seed(runif(1,0,100))
  } else {
    set.seed(seed)
  }
  L = matrix(runif(num.vars*num.factors,-1,1), num.vars,num.factors)
  if(sparse){
    if(num.vars >= 1000){
      numeff = min(6*num.factors + sample(1:num.factors, 1), num.vars)
    }else if(num.vars >= 500){
      numeff = min(5*num.factors + sample(1:num.factors, 1), num.vars)
    }else{
      numeff = min(3*num.factors + sample(1:num.factors, 1), num.vars)
    }
    for(h in 1:num.factors){
      idx = sample(1:num.vars,numeff, replace = FALSE)
      L[-idx,h] = 0.0
    }
  }
  cov = L%*%t(L) + diag(rep(0.05, num.vars))
  L = t(chol(cov))
  
  z = rnorm(n*num.factors,0,1)
  z = matrix(z, n, num.vars)
  
  if(student_t){
    sqrtgam = diag(1/sqrt(rgamma(shape = nu/2, rate = nu/2, n)))
    data = data.frame(sqrtgam%*%z%*%t(L))
  }else{
    data = data.frame(z%*%t(L))
  }
  
  var.names <- paste("VAR",seq(1,ncol(data)), sep="")
  names(data) <- var.names
  true_cov = cov
  if(student_t){
    true_cov = cov*(nu/(nu-2))
  }
  return( list(data, true_cov) )
}

generateMultipleData = function(num.vars=3, num.factors = num.vars, n=100, rep = 50, seed=NULL, student_t = FALSE, nu = 10 , sparse = TRUE ) {
  if(is.null(seed)){
    set.seed(runif(1,0,100))
  } else {
    set.seed(seed)
  }
  
  N = n*rep
  L = matrix(runif(num.vars*num.factors,-1,1), num.vars,num.factors)
  if(sparse){
    if(num.vars >= 1000){
      numeff = min(6*num.factors + sample(1:num.factors, 1), num.vars)
    }else if(num.vars >= 500){
      numeff = min(5*num.factors + sample(1:num.factors, 1), num.vars)
    }else{
      numeff = min(3*num.factors + sample(1:num.factors, 1), num.vars)
    }
    for(h in 1:num.factors){
      idx = sample(1:num.vars,numeff, replace = FALSE)
      L[-idx,h] = 0.0
    }
  }
  cov = L%*%t(L) + diag(rep(0.05, num.vars))
  L = t(chol(cov))
  
  z = rnorm(N*num.factors, 0,1)
  z = matrix(z, N, num.vars)
  
  if(!student_t){
    data = data.frame(z%*%t(L))
  }else{
    sqrtgam = diag(sqrt(rgamma(shape = nu/2, rate = nu/2, N)))
    data = data.frame(sqrtgam%*%z%*%t(L))
  }
  
  var.names <- paste("VAR",seq(1,ncol(data)), sep="")
  names(data) <- var.names
  
  true_cov = cov
  if(student_t){
    true_cov = cov*(nu/(nu-2))
  }
  return( list(data, L%*%t(L)) )
}


linearMGSP_t = function(X, nrun, burn, thin = 1, prop = 0.7, epsilon = 1e-3, nu = 5, 
                        kinit = NULL, adapt = TRUE, output = c("covMean", "covSamples", 
                                                               "factSamples", "sigSamples", 
                                                               "numFactors"), 
                        verbose = TRUE, dump = FALSE, filename = "samps.Rds",
                        buffer = 10000, augment = NULL, 
                        epseta = 0.1, a1sigma = 5, a2sigma = 0.1){
  
  if(nrun <= burn) stop("nrun must be larger than burn")
  if(!is.matrix(X)) stop("X must be a matrix")
  if(any(is.na(X))) stop("X cannot contain missing data")
  if(!is.null(augment)) if(!is.expression(augment)) stop("augment must be an expression (see expression())")
  
  cm = any(output %in% "covMean")
  cs = any(output %in% "covSamples")
  fs = any(output %in% "factSamples")
  ss = any(output %in% "sigSamples")
  nf = any(output %in% "numFactors")
  
  p = ncol(X)
  n = nrow(X)
  
  as = 1                          # gamma hyperparameters for residual precision
  bs = 0.3                        
  df = 3                          # gamma hyperparameters for t_{ij}
  ad1 = 2.1
  bd1 = 1                         # gamma hyperparameters for delta_1
  ad2 = 3.1
  bd2 = 1                         # gamma hyperparameters delta_h, h >= 2
  adf = 1
  bdf = 1 
  b0 = 1
  b1 = 0.0005
  
  if(is.null(kinit)) kinit = floor(log(p)*3)
  
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  
  VY= apply(X, 2, var)                  # explicitly preserve scale
  scaleMat = sqrt((VY) %*% t(VY))
  X = scale(X)
  num = 0
  k=kinit                               # no. of factors to start with
  
  # --- Initial values --- #
  ps = rgamma(p, shape = as, rate = bs)                           # Sigma = diagonal residual covariance
  lambda = matrix(1, nrow = p, ncol = k)           # latent factor distribution = standard normal
  psijh = matrix(rgamma(p*k, shape = df/2, rate = df/2), nrow = p, ncol = k)     # local shrinkage coefficients
  delta = c(rgamma(1,shape = ad1, rate = bd1), rgamma(k-1, shape = ad2,rate = bd2))       # global shrinkage coefficients multilpliers
  tauh = cumprod(delta)                                       # global shrinkage coefficients
  Plam = t(t(psijh) * (tauh))                                     # precision of loadings rows
  #gam = rgamma(n, shape = nu/2, rate = nu/2)
  #eta = matrix(0, n, k)
  eta = matrix(rnorm(n*k), n, k)
  gam = rep(1,n)
  
  rejectcount1 = 0
  rejectcount2 = 0
  
  if(cm) COVMEAN = matrix(0, nrow = p, ncol = p)
  if(cs) OMEGA = array(dim = c(p, p, sp))
  if(fs) {LAMBDA = list()
  ETA = list()
  GAMMA = matrix(0, nrow = n, ncol = sp)
  A1 = c()
  A2 = c()
  NU = c()
  TREECOUNT = c()
  }
  if(ss) SIGMA = array(dim = c(p, sp))
  if(nf) K = rep(NA, sp)
  ind = 1
  
  at = ceiling(nrun/100)
  if(verbose) {
    pb = txtProgressBar(style = 3)
  }
  
  for(i in 1:nrun){
    treecount = 0
    
    lambda = lam_lin_t(eta, Plam, ps, k, p, X, gam)
    ps = sig_lin_t(lambda, eta, k, p, n, X, as, bs, gam)
    
    eta = eta_NUTS_2(eta, lambda, ps, k, n, p, X, nu, epseta, treecount)
    gam = gam_gibbs_t(lambda, eta, ps, n, k, p, X, nu)
    
    psijh = psi_mg_t(lambda, tauh, ps, k, p, df)
    delta = del_mg_t(lambda, psijh, tauh, delta, k, p, ad1, bd1, ad2, bd2)
    tauh = cumprod(delta)
    Plam = plm_mg_t(psijh, tauh)
    
    ad1 = a1_mg_t(ad1, delta, a1sigma, rejectcount1)
    ad2 = a2_mg_t(ad2, delta, k, a2sigma, rejectcount2)
    
    if(!is.null(augment)) eval(augment)
    
    if(adapt){
      # ----- make adaptations ----#
      prob = 1/exp(b0 + b1*i)                    # probability of adapting
      uu = runif(1)
      lind = colSums(abs(lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
      vec = lind >= prop
      num = sum(vec)                             # number of redundant columns
      
      if(uu < prob) {
        if((i > 20) & (num == 0) & all(lind < 0.995)) {
          k = k + 1
          lambda = cbind(lambda, rep(0,p))
          eta = cbind(eta,rnorm(n))
          psijh = cbind(psijh, rgamma(p,df/2,df/2))
          delta[k] = rgamma(1, ad2,bd2)
          tauh = cumprod(delta)
          Plam = t(t(psijh) * tauh)
        } else {
          if (num > 0) {
            k = max(k - num,1)
            lambda = lambda[,!vec, drop = F]
            psijh = psijh[,!vec, drop = F]
            eta = eta[,!vec, drop = F]
            delta = delta[!vec]
            tauh = cumprod(delta)
            Plam = t(t(psijh) * tauh)
          }
        }
      }
    }
    
    if((i %% thin == 0) & (i > burn)) {
      if(cm | cs) Omega = (nu/(nu-2))*((tcrossprod(lambda) + diag(1/c(ps))) * scaleMat)
      if(cm) COVMEAN = COVMEAN + (Omega / sp)
      if(cs) OMEGA[,,ind] = Omega
      if(fs) {LAMBDA[[ind]] = lambda
      ETA[[ind]] = eta
      GAMMA[,ind] = gam
      A1[ind] = ad1
      A2[ind] = ad2
      NU[ind] = nu
      TREECOUNT[ind] = treecount
      }
      if(ss) SIGMA[,ind] = 1/ps
      if(nf) K[ind] = k
      ind = ind + 1
    }
    
    if(verbose & (i %% at == 0)) setTxtProgressBar(pb, i / nrun)
    
    if(dump & (i %% buffer == 0)){
      out = list()
      if(cm) out = c(out, list(covMean = COVMEAN))
      if(cs) out = c(out, list(omegaSamps = OMEGA))
      if(fs) out = c(out, list(lambdaSamps = LAMBDA, 
        etaSamps = ETA, 
        gammaSamps = GAMMA, a1Samps = A1, a2Samps = A2, 
        nuSamps = NU, 
        treecount = TREECOUNT))
      if(ss) out = c(out, list(sigmaSamps = SIGMA))
      if(nf) out = c(out, list(numFacts = K))
      saveRDS(out, filename, compress = FALSE)
    }
  }
  
  if(verbose) close(pb)
  
  out = list()
  if(cm) out = c(out, list(covMean = COVMEAN))
  if(cs) out = c(out, list(omegaSamps = OMEGA))
  if(fs) out = c(out, list(lambdaSamps = LAMBDA, etaSamps = ETA, gammaSamps = GAMMA, a1Samps = A1, a2Samps = A2, nuSamps = NU, treecount = TREECOUNT))
  if(ss) out = c(out, list(sigmaSamps = SIGMA))
  if(nf) out = c(out, list(numFacts = K))
  
  cat("a1 rejection ratio : ", (100*rejectcount1)/nrun, "%",'\n')
  cat("a2 rejection ratio : ", (100*rejectcount2)/nrun, "%",'\n')
  
  cat("a1 rejection count : ", rejectcount1, "/",nrun,'\n')
  cat("a2 rejection count : ", rejectcount2, "/",nrun,'\n')
  return(out)
}

