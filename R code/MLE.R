"
Functions included in the file:
* mle.copula(): main function, outputs maximum likelihood estimates for the parameter
* For estimating the mle parameters for the marginals:
    * mle.marginals.wrpcauchy()
    * mle.marginals.katojones()
    * mle.copula.params()
    
The functions included in the file 'functions for trivariate wrapped Cauchy copula.R' are also required for this to run.
"

mle.copula = function(x,ninipar=NULL,max=NULL,min=NULL, marginals = rep('uniform', 3),  
                               verbose = FALSE, return_data = FALSE){ 
  
  ## Checks for the inputs
  
  if(!is.matrix(x)) { stop("Argument 'x' must be a numeric matrix")}
  if(dim(x)[2]!=3) { stop("Argument 'x' must be a matrix with three columns")}
  
  if (!is.numeric(ninipar)) {
    warning("Argument 'ninipar' must be a positive integer number. Default value of 'ninipar' was used")
    ninipar=100
  }
  if ((length(ninipar) != 1) | (ninipar%%1 != 0) | (ninipar <= 0)){
    warning("Argument 'ninipar' must be a positive integer number Default value of 'ninipar' was used")
    ninipar=100
  }
  
  if(is.null(max)) max = 10
  if (!is.numeric(max)) {stop("Argument 'max' must be a real number")}
  if (length(max) != 1) {stop("Argument 'max' must be a real number")} 
  
  if(is.null(min)) min = -10
  if (!is.numeric(min)) {stop("Argument 'min' must be a real number")}
  if (length(min) != 1) {stop("Argument 'min' must be a real number")}
  
  marginals <- tolower(marginals)
  if(length(marginals) != 3){stop("'marginals' should be a vector of length 3 with elements 'uniform', 'wrapped cauchy', 'cardioid', 
         'vonMises' or 'katojones' for the circular part and 'weibull' for the linear part")}
  if(any(!(marginals %in% c('uniform', 'wrapped cauchy', 'cardioid', 'vonmises', 'katojones', 'weibull')))){
    stop("'marginals' should be a vector of length 3 with elements 'uniform', 'wrapped cauchy', 'cardioid', 
         'vonMises' or 'katojones' for the circular part and 'weibull' for the linear part")
  }
  
  require(Rsolnp)
  # require(R.utils)
  require(CircStats)
  require(Directional)
  require(VGAM)
  require(circular)
  
  ## First calculate parameters for the marginals
  
  pars1 <- list(); pars2 <- list(); pars3 <- list()
  if (any(marginals == 'uniform')){
    index <- which(marginals == 'uniform')
    for(i in 1:length(index)){
      eval(parse(text = paste('uniform_x', index[i], ' = x[,', index[i], ']', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = 1', sep ='')))
    }
  } 
  if (any(marginals == 'wrapped cauchy')){
    index <- which(marginals == 'wrapped cauchy')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = mle.marginals.wrpcauchy(x[,', index[i], '], ninipar = ninipar, verbose = verbose)', sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'wrapped cauchy', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dwrpcauchy(x[,', index[i],'],pars', index[i],'$mu,pars', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'cardioid')){
    index <- which(marginals == 'cardioid')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = cardio.mle(x[,', index[i], '], rads = TRUE)', sep ='')))
      eval(parse(text = paste('pars', index[i], " = list('mu' = as.numeric(pars", index[i], "$param[1]), 'rho' = as.numeric(pars", index[i], "$param[2]))", sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'cardioid', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dcard(x[,', index[i],'],pars', index[i],'$mu,pars', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'vonmises')){
    index <- which(marginals == 'vonmises')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = mle.vonmises(circular(x[,', index[i], ']))', sep ='')))
      eval(parse(text = paste('pars', index[i], " = list('mu' = as.numeric(pars", index[i], "$mu), 'kappa' = pars", index[i], "$kappa)", sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'vonmises', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dvonmises(circular(x[,', index[i],']),circular(pars', index[i],'$mu),pars', index[i], '$kappa)', sep ='')))
    }
  } 
  if(any(marginals == 'katojones')){
    index <- which(marginals == 'katojones')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = mle.marginals.katojones(x[,', index[i], '], ninipar = ninipar, verbose = verbose)', sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'katojones', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dkatojones(x[,', index[i],'],pars', index[i],'$mu,pars', index[i], '$gamma,pars', index[i], '$rho,pars', index[i], '$lambda)', sep ='')))
    }
  }
  if(any(marginals == 'weibull')){
    require(EnvStats)
    index <- which(marginals == 'weibull')
    for(i in 1:length(index)){
      eval(parse(text = paste('pars', index[i], ' = eweibull(x[,', index[i], '], method = "mle")', sep ='')))
      eval(parse(text = paste('pars', index[i], " = list('shape' = as.numeric(pars", index[i], "$parameters[1]), 'scale' = as.numeric(pars", index[i], "$parameters[2]))", sep ='')))
      eval(parse(text = paste('uniform_x', index[i], ' = vector_to_uniform(x[,', index[i], "], marginals = 'weibull', params = pars", index[i], ')', sep ='')))
      eval(parse(text = paste('dens', index[i], ' = dweibull(x[,', index[i],'],pars', index[i],'$shape,pars', index[i], '$scale)', sep ='')))
    }
  }
  if(marginals[1] != 'uniform') names(pars1) <- paste(names(pars1), 1, sep = '')
  if(marginals[2] != 'uniform') names(pars2) <- paste(names(pars2), 2, sep = '')
  if(marginals[3] != 'uniform') names(pars3) <- paste(names(pars3), 3, sep = '')
  
  dens_marginals <- dens1*dens2*dens3
  
  uniform_x <- cbind(uniform_x1, uniform_x2, uniform_x3)
  
  copula_params <- mle.copula.params(uniform_x = uniform_x, ninipar = ninipar, max = max, 
                                     min = min, verbose = verbose)
  
  resultmle=list()
  resultmle$rho12 = copula_params$final_rho12
  resultmle$rho13 = copula_params$final_rho13
  resultmle$rho23 = copula_params$final_rho23
  resultmle$params1 = pars1
  resultmle$params2 = pars2
  resultmle$params3 = pars3
  # resultmle <- append(resultmle, append(pars1, append(pars2, pars3, length(pars2)), length(pars1)), length(resultmle))
  
  resultmle$LL = copula_params$final_LL + sum(log(dens_marginals))
  
  num_params <- 2 + length(resultmle$params1) + length(resultmle$params2) + length(resultmle$params3)
  resultmle$AIC <- -2*resultmle$LL + 2*num_params
  resultmle$BIC <- -2*resultmle$LL + num_params*log(dim(x)[1])
  resultmle$no_params <- num_params
  
  if(return_data){
    resultmle$u = uniform_x
    resultmle$density_marginals = sum(log(dens_marginals))
  }
  
  return(resultmle)
  
}

mle.marginals.wrpcauchy <- function(x, ninipar=NULL, verbose = FALSE){
  if(is.null(ninipar)) ninipar = 100
  
  llpar <- function(par){
    -sum(log(dwrpcauchy(x, mu = par[1], rho = par[2])))
  }
  
  LB = c(0, 0); UB = c(2*pi, 1)
  
  valllfin=Inf
  paramfin=numeric()
  
  for(iterip in 1:ninipar){
    if(verbose){print(iterip)}
    
    inipar=numeric()
    
    inipar[1] = runif(1,0,2*pi)
    inipar[2] = runif(1,0,1)
    
    require(Rsolnp)
    
    paramtot <- try(solnp(pars=inipar, fun=llpar, LB = LB, UB = UB,
                          control = list(trace=0)), silent = T)
    
    if(class(paramtot)=="try-error"){
      valll=Inf
    }else{
      valll=paramtot$values[length(paramtot$values)]
      
      if(valll<valllfin){
        paramfin=paramtot$pars
        valllfin=valll
      }
    }
  }
  
  resultmle=list()
  resultmle$mu = paramfin[1]
  resultmle$rho = paramfin[2]
  # resultmle$LL = valllfin
  
  return(resultmle)
}

mle.marginals.katojones <- function(theta, ninipar=NULL, verbose = FALSE){
  if(is.null(ninipar)) ninipar = 10
  
  likeli=function(parame,theta){
    n=length(theta)
    alpha2m=parame[2]+(1-parame[2])*parame[4]*cos(parame[3])
    beta2m=(1-parame[2])*parame[4]*sin(parame[3])
    ll=sum(log(1+2*parame[2]*(cos(theta-parame[1])-alpha2m)/(1+alpha2m^2+beta2m^2-2*(alpha2m*cos(theta-parame[1])+beta2m*sin(theta-parame[1])))))
    return(-ll+n*log(2*pi))
  }
  
  mini = Inf
  for(k in 1:ninipar){
    for(l in 1:ninipar){
      for(m in 1:ninipar){
        for(n in 1:ninipar){
          para=c(2*pi*k/ninipar,0.9*l/ninipar,2*pi*m/ninipar,0.9*n/ninipar)
          inte=nlminb(para,function(x) likeli(x,theta),lower=c(-1e+7,0,-1e+7,0),upper=c(1e+7,1-1e-3,1e+7,1-1e-3))
          if (mini>inte[[2]]){
            mini=inte[[2]]
            parameter=inte[[1]]
          }
        }
      }
    }
  }
  
  mu_ml=Arg(exp((1i)*parameter[1]))
  gamma_ml=parameter[2]
  phi_ml=Arg(exp((1i)*parameter[3]))
  omega_ml=parameter[4]
  alpha2_ml=gamma_ml^2+gamma_ml*(1-gamma_ml)*omega_ml*cos(phi_ml)
  beta2_ml=gamma_ml*(1-gamma_ml)*omega_ml*sin(phi_ml)
  
  rho_ml = 1/gamma_ml * sqrt(alpha2_ml^2 + beta2_ml^2) 
  lambda_ml = Arg(complex(real = alpha2_ml, imaginary = beta2_ml))
  
  resultmle=list()
  resultmle$mu = mu_ml
  resultmle$gamma = gamma_ml
  resultmle$rho = rho_ml
  resultmle$lambda = lambda_ml
  # resultmle$LL = valllfin
  
  return(resultmle)
}

mle.copula.params <- function(uniform_x,ninipar=NULL,max=NULL,min=NULL, verbose = FALSE){
  LB = c(-1,-Inf); UB = c(1,Inf)
  
  ## (k,l) = (1,2)
  
  ## - Log-likelihood because the function used finds the minimum of the input function
  llpar12=function(par){
    rho23 = par[2]
    rho13 = (1 + sqrt(1 + 4*abs(rho23)^3)) / (2*abs(rho23)^2) * 1 / par[1]
    rho12 = 1/(rho23*rho13)
    
    -sum(log(dtri_copula(uniform_x, rho12 = rho12, rho13 = rho13, rho23 = rho23)))
  }
  
  valllfin12=Inf
  paramfin12=numeric()
  
  for(iterip in 1:ninipar){
    if(verbose){print(iterip)}
    
    inipar=numeric()
    
    inipar[1] = runif(1,-1,1)
    inipar[2] = runif(1,min,max)
    
    paramtot <- try(solnp(pars=inipar, fun=llpar12, LB = LB, UB = UB,
                          control = list(trace=0)), silent = T)
    
    if(class(paramtot)=="try-error"){
      valll12=Inf
    }else{
      valll12=paramtot$values[length(paramtot$values)]
      
      if(valll12<valllfin12){
        paramfin12=paramtot$pars
        valllfin12=valll12
      }
    }
  }
  
  
  
  ## (k,l) = (2,3)
  
  ## - Log-likelihood because the function used finds the minimum of the input function
  llpar23=function(par){
    # cat(par[1],par[2],'\n')
    rho13 = par[2]
    rho12 = (1 + sqrt(1 + 4*abs(rho13)^3)) / (2*abs(rho13)^2) * 1 / par[1]
    rho23 = 1/(rho12*rho13)
    
    -sum(log(dtri_copula(uniform_x, rho12 = rho12, rho13 = rho13, rho23 = rho23)))
  }
  
  valllfin23=Inf
  paramfin23=numeric()
  
  for(iterip in 1:ninipar){
    if(verbose){print(iterip)}
    
    inipar=numeric()
    
    inipar[1] = runif(1,-1,1)
    inipar[2] = runif(1,min,max)
    
    paramtot <- try(solnp(pars=inipar,fun=llpar23,LB = LB, UB = UB,
                          control = list(trace=0)), silent=T)
    
    if(class(paramtot)=="try-error"){
      valll23=Inf
    }else{
      valll23=paramtot$values[length(paramtot$values)]
      
      if(valll23<valllfin23){
        paramfin23=paramtot$pars
        valllfin23=valll23
      }
    }
  }
  
  
  ## (k,l) = (3,1)
  
  ## - Log-likelihood because the function used finds the minimum of the input function
  llpar13=function(par){
    rho12 = par[2]
    rho23 = (1 + sqrt(1 + 4*abs(rho12)^3)) / (2*abs(rho12)^2) * 1 / par[1]
    rho13 = 1/(rho23*rho12)
    
    -sum(log(dtri_copula(uniform_x, rho12 = rho12, rho13 = rho13, rho23 = rho23)))
  }
  
  valllfin13=Inf
  paramfin13=numeric()
  
  for(iterip in 1:ninipar){
    if(verbose) {print(iterip)}
    
    inipar=numeric()
    
    inipar[1] = runif(1,-1,1)
    inipar[2] = runif(1,min,max)
    
    paramtot <- try(solnp(pars=inipar,fun=llpar13,LB = LB, UB = UB,
                          control = list(trace=0)), silent=T)
    
    if(class(paramtot)=="try-error"){
      valll13=Inf
    }else{
      valll13=paramtot$values[length(paramtot$values)]
      
      if(valll13<valllfin13){
        paramfin13=paramtot$pars
        valllfin13=valll13
      }
    }
  }
  
  if((valllfin12 < valllfin13) && (valllfin12 < valllfin23)){
    paramfin = paramfin12
    final_rho23 = paramfin[2]
    final_rho13 = (1 + sqrt(1 + 4*abs(paramfin[2])^3)) / (2*abs(paramfin[2])^2) * 1/paramfin[1]
    final_rho12 = 1/(final_rho23*final_rho13)
    
    final_LL=-valllfin12 ## Estimated LL value
    
  } else if((valllfin23 < valllfin12) && (valllfin23 < valllfin13)){
    paramfin = paramfin23
    final_rho13 = paramfin[2]
    final_rho12 = (1 + sqrt(1 + 4*abs(paramfin[2])^3)) / (2*abs(paramfin[2])^2) * 1/paramfin[1]
    final_rho23 = 1/(final_rho12*final_rho13)
    
    final_LL=-valllfin23 ## Estimated LL value
    
  } else {
    paramfin = paramfin13
    final_rho12 = paramfin[2]
    final_rho23 = (1 + sqrt(1 + 4*abs(paramfin[2])^3)) / (2*abs(paramfin[2])^2) * 1/paramfin[1]
    final_rho13 = 1/(final_rho23*final_rho12)
    
    final_LL=-valllfin13 ## Estimated LL value
  }
  
  return(list('final_rho12' = final_rho12, 'final_rho23' = final_rho23, 'final_rho13' = final_rho13, 'final_LL' = final_LL))
}


## Examples
# The functions from 'functions for trivariate wrapped Cauchy copula.R' are required
source("~/functions for trivariate wrapped Cauchy copula.R")

set.seed(2)
x <- rtri(1500, rho12 = 1, rho13 = 0.25, rho23 = 4, 
          marginals = c('wrapped Cauchy', 'uniform', 'weibull'), 
          params1 = list(mu = 0, rho = 0.2), params3 = list(shape = 1, scale = 1))

# Estimating the maximum likelihood estimates of the parameters
parameters <- mle.copula(x, ninipar = 100, max = 50, min = -50,
                         marginals = c('wrapped Cauchy', 'uniform', 'weibull'))
# Copula parameters
parameters[c('rho12', 'rho13', 'rho23')]
parameters$params1
parameters$params2
parameters$params3

