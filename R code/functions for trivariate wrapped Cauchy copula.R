"
Functions included in the file:
* dtri(): function for calculating the density of given points given the parameters
* multiply_copula_denisty(): function for calculating the product of marginal densities. Needs to be multiplied with dtri_copula() to get the denisty of the copula.
* dtri_copula(): function for calculating the value of the density c(). Needs to be multiplied with multiply_copula_denisty() to get the denisty of the copula.
* rtri(): function for random variate generation
* check_parameters(): Function for checking that all inputs are allowed
* cdf_wrpcauchy(): calculating the cdf of wrapped Cauchy distribution given theta, mu, rho
* inv_cdf_wrpcauchy(): calculating the inverse cdf of wrapped Cauchy distribution given p, mu, rho
* dkatojones(): calculating the pdf of Kato-Jones distribution given theta, mu, gamma, rho, lambda
* cdf_katojones(): calculating the cdf of Kato-Jones distribution given theta, mu, gamma, rho, lambda
* vector_to_uniform(): finding the value of F^{-1}(theta) for all three components according to the input of marginals
"

dtri <- function(x, rho12, rho13, rho23, marginals, params1 = NULL, params2 = NULL, params3 = NULL){
  dtri_copula(x, rho12, rho13, rho23, marginals = marginals, params1 = params1, params2 = params2, params3 = params3) * multiply_copula_denisty(x, marginals = marginals, params1 = params1, params2 = params2, params3 = params3)
}

multiply_copula_denisty <- function(x, marginals, params1 = NULL, params2 = NULL, params3 = NULL){
  dens1 <- 1; dens2 <- 1; dens3 <- 1
  if (any(marginals == 'wrapped cauchy')){
    index <- which(marginals == 'wrapped cauchy')
    for(i in 1:length(index)){
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dwrpcauchy(x[,', index[i],'],params', index[i],'$mu,params', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'cardioid')){
    require(VGAM)
    index <- which(marginals == 'cardioid')
    for(i in 1:length(index)){
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dcard(x[,', index[i],'],params', index[i],'$mu,params', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'vonmises')){
    require(circular)
    index <- which(marginals == 'vonmises')
    for(i in 1:length(index)){
      suppressWarnings(
        eval(parse(text = paste('dens', index[i], ' = (2*pi) * dvonmises(circular(x[,', index[i],']),circular(params', index[i],'$mu),params', index[i], '$kappa)', sep ='')))
      )
    }
  } 
  if (any(marginals == 'katojones')){
    index <- which(marginals == 'katojones')
    for(i in 1:length(index)){
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dkatojones(x[,', index[i],'],params', index[i],'$mu,params', index[i], '$gamma,params', index[i], '$rho,params', index[i], '$lambda)', sep ='')))
    }
  }
  if (any(marginals == 'weibull')){
    index <- which(marginals == 'weibull')
    for(i in 1:length(index)){
      if(any(is.nan(dweibull(x[,3],params3$shape,params3$scale)))){
        temp <- which(is.nan(dweibull(x[,3],params3$shape,params3$scale)))
        cat(x[temp,3],params3$shape,params3$scale, '\n')
      }
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dweibull(x[,', index[i],'],params', index[i],'$shape,params', index[i], '$scale)', sep ='')))
    }
  }
  
  density <- dens1 * dens2 * dens3
  
  return(density)
}

## Note that this function returns the value of the copula part, without multiplying by the marginals
dtri_copula <- function(x, rho12, rho13, rho23, marginals = rep('uniform', 3), params1 = NULL, params2 = NULL, params3 = NULL){
  if(!is.matrix(x)) { stop("Argument 'x' must be a numeric matrix")}
  if(dim(x)[2]!=3) { stop("Argument 'x' must be a matrix with three columns")}
  
  cond1 <- abs(rho23) < abs(rho12 * rho13) / (abs(rho12) + abs(rho13))
  cond2 <- abs(rho13) < abs(rho12 * rho23) / (abs(rho12) + abs(rho23))
  cond3 <- abs(rho12) < abs(rho23 * rho13) / (abs(rho23) + abs(rho13))
  if(!(cond1 | cond2 | cond3)){stop("abs(rhojk) < abs(rhojl * rhokl) / (abs(rhojl) + abs(rhokl)) must
                                  be satisfied for some combination of j,k,l=1,2,3")}
  if(rho12*rho13*rho23 < 0){stop("rho12*rho13*rho23 must be greater than 0")}
  if(rho12*rho13*rho23 == 0){stop("rho12*rho13*rho23 must be greater than 0")}
  
  marginals <- tolower(marginals)
  check_parameters(marginals = marginals, params1 = params1, params2 = params2, params3 = params3)
  
  if(any(marginals != 'uniform')){
    not_unif <- which(marginals != 'uniform')
    for(i in 1:length(not_unif)){
      x[,not_unif[i]] <- vector_to_uniform(x[,not_unif[i]], marginals = marginals[not_unif[i]], params = eval(parse(text = paste('params', not_unif[i], sep =''))))
    }
  }
  c1 = rho12*rho13/rho23 + rho12*rho23/rho13 + rho13*rho23/rho12
  c2 = (1/(2*pi)^3) * sqrt((rho12*rho13/rho23)^2 + (rho12*rho23/rho13)^2 + (rho13*rho23/rho12)^2 - 
                             2*(rho12^2 + rho13^2 + rho23^2))
  density = c2 / (c1 + 2*(rho12*cos(x[,1]-x[,2]) + rho13*cos(x[,1]-x[,3]) + 
                            rho23*cos(x[,2]-x[,3])))
  
  return(density)
}

## Draw sample from trivariate copula
# n is the number of observations
# rho and mu are parameters for the marginals
rtri <- function(n, rho12, rho13, rho23, marginals = rep('uniform', 3), 
                          params1 = NULL, params2 = NULL, params3 = NULL){
  
  marginals <- tolower(marginals)
  check_parameters(marginals = marginals, params1 = params1, params2 = params2, params3 = params3)
  
  cond1 <- abs(rho23) < abs(rho12 * rho13) / (abs(rho12) + abs(rho13))
  cond2 <- abs(rho13) < abs(rho12 * rho23) / (abs(rho12) + abs(rho23))
  cond3 <- abs(rho12) < abs(rho23 * rho13) / (abs(rho23) + abs(rho13))
  if(!(cond1 | cond2 | cond3)){stop("abs(rhojk) < abs(rhojl * rhokl) / (abs(rhojl) + abs(rhokl)) must 
                                    be satisfied for some combination of j,k,l=1,2,3")}
  if(rho12*rho13*rho23 < 0){stop("rho12*rho13*rho23 must be greater than 0")}
  if(rho12*rho13*rho23 == 0){stop("rho12*rho13*rho23 must be greater than 0")}
  
  c2 = (1/(2*pi)^3) * sqrt((rho12*rho13/rho23)^2 + (rho12*rho23/rho13)^2 + (rho13*rho23/rho12)^2 - 2*(rho12^2 + rho13^2 + rho23^2))
  
  phi12 <- (rho13*rho23/rho12 - rho12*rho13/rho23 - rho12*rho23/rho13 - (2*pi)^3*c2) / (2*rho12)
  eta12 <- Arg(phi12)
  delta12 <- Mod(phi12)
  
  omega1 <- runif(n, 0, 1)
  omega2 <- runif(n, 0, 1)
  omega3 <- runif(n, 0, 1)
  
  u1 <- 2*pi*omega1
  
  u2 <- u1 + eta12 + 2*atan((1-delta12)/(1+delta12) * tan(pi*(omega2 - 0.5)))
  
  phi312 <- -rho12 * (rho23^(-1)*exp(complex(imaginary = u1)) + rho13^(-1)*exp(complex(imaginary = u2)))
  delta312 <- Mod(phi312)
  eta312 <- Arg(phi312)
  
  u3 <- eta312 + 2*atan((1-delta312)/(1+delta312) * tan(pi*(omega3 - 0.5)))
  
  u1 <- u1 %% (2*pi); u2 <- u2 %% (2*pi); u3 <- u3 %% (2*pi)
  
  suppressWarnings({
    if (any(marginals == 'wrapped cauchy')){
      index <- which(marginals == 'wrapped cauchy')
      for(i in 1:length(index)){
        eval(parse(text = paste('u', index[i], ' = inv_cdf_wrpcauchy(u', index[i],',params', index[i],'$mu,params', index[i], '$rho)', sep ='')))
        eval(parse(text = paste('u', index[i], ' = u', index[i], '%% (2*pi)', sep ='')))
      }
    } 
    if (any(marginals == 'cardioid')){
      require(VGAM)
      index <- which(marginals == 'cardioid')
      for(i in 1:length(index)){
        eval(parse(text = paste('u', index[i], ' = qcard(u', index[i],'/(2*pi),params', index[i],'$mu,params', index[i], '$rho)', sep ='')))
        eval(parse(text = paste('u', index[i], ' = u', index[i], '%% (2*pi)', sep ='')))
      }
    } 
    if (any(marginals == 'vonmises')){
      require(circular)
      index <- which(marginals == 'vonmises')
      for(i in 1:length(index)){
        eval(parse(text = paste('u', index[i], ' = qvonmises(u', index[i],'/(2*pi),circular(params', index[i],'$mu),params', index[i], '$kappa)', sep ='')))
        eval(parse(text = paste('u', index[i], ' = u', index[i], '%% (2*pi)', sep ='')))
      }
    } 
    if (any(marginals == 'katojones')){
      stop('Data generation is not supported for Kato-Jones marginals')
    }
    if (any(marginals == 'weibull')){
      index <- which(marginals == 'weibull')
      for(i in 1:length(index)){
        eval(parse(text = paste('u', index[i], ' = qweibull(u', index[i],'/(2*pi),params', index[i],'$shape,params', index[i], '$scale)', sep ='')))
      }
    }
  
  u <- matrix(c(u1, u2, u3), ncol = 3)
  return(u)
  })
}

## Function for checking that all inputs are allowed
check_parameters <- function(marginals = rep('uniform', 3), params1 = NULL, params2 = NULL, params3 = NULL){
  if(length(marginals) != 3){stop("'marginals' should be a vector of length 3 with elements 'uniform', 'wrapped cauchy', 'cardioid',
         'vonMises' or 'katojones' for the circular part and 'weibull' for the linear part")}
  if(any(!(marginals %in% c('uniform', 'wrapped cauchy', 'cardioid', 'vonmises', 'katojones', 'weibull')))){
    stop("'marginals' should be a vector of length 3 with elements 'uniform', 'wrapped cauchy', 'cardioid', 
         'vonMises' or 'katojones' for the circular part and 'weibull' for the linear part")
  }
  
  if(any(marginals != 'uniform')){
    
    ind_temp = ifelse(is.null(params3), 2, 3)
    
    for(i in 1:ind_temp){
      if(marginals[i] == 'wrapped cauchy'){
        if(!is.list(eval(parse(text = paste('params', i, sep =''))))){stop("'params1', 'params2', 'params3' should be a list containing the parameters of the marginal distribution. Wrapped Cauchy and cardioid require rho, mu. von Mises requires mu, kappa. Kato-Jones requires mu, gamma, rho, lambda. Weibull requires shape and scale.")}
        
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$rho)){stop("The list of parameters for wrapped cauchy marginals should include 'rho'")}
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$mu)){stop("The list of parameters for wrapped cauchy marginals should include 'mu'")}
      } else if (marginals[i] == 'cardioid'){
        if(!is.list(eval(parse(text = paste('params', i, sep =''))))){stop("'params1', 'params2', 'params3' should be a list containing the parameters of the marginal distribution. Wrapped Cauchy and cardioid require rho, mu. von Mises requires mu, kappa. Kato-Jones requires mu, gamma, rho, lambda. Weibull requires shape and scale.")}
        
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$rho)){stop("The list of parameters for cardioid marginals should include 'rho'")}
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$mu)){stop("The list of parameters for cardioid marginals should include 'mu'")}
      } else if (marginals[i] == 'vonmises'){
        if(!is.list(eval(parse(text = paste('params', i, sep =''))))){stop("'params1', 'params2', 'params3' should be a list containing the parameters of the marginal distribution. Wrapped Cauchy and cardioid require rho, mu. von Mises requires mu, kappa. Kato-Jones requires mu, gamma, rho, lambda. Weibull requires shape and scale.")}
        
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$mu)){stop("The list of parameters for von Mises marginals should include 'mu'")}
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$kappa)){stop("The list of parameters for von Mises marginals should include 'kappa'")}
      } else if (marginals[i] == 'katojones'){
        if(!is.list(eval(parse(text = paste('params', i, sep =''))))){stop("'params1', 'params2', 'params3' should be a list containing the parameters of the marginal distribution. Wrapped Cauchy and cardioid require rho, mu. von Mises requires mu, kappa. Kato-Jones requires mu, gamma, rho, lambda. Weibull requires shape and scale.")}
        
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$mu)){stop("The list of parameters for Kato-Jones marginals should include 'mu'")}
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$gamma)){stop("The list of parameters for Kato-Jones marginals should include 'gamma'")}
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$rho)){stop("The list of parameters for Kato-Jones marginals should include 'rho'")}
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$lambda)){stop("The list of parameters for Kato-Jones marginals should include 'lambda'")}
      } else if (marginals[i] == 'weibull'){
        if(!is.list(eval(parse(text = paste('params', i, sep =''))))){stop("'params1', 'params2', 'params3' should be a list containing the parameters of the marginal distribution. Wrapped Cauchy and cardioid require rho, mu. von Mises requires mu, kappa. Kato-Jones requires mu, gamma, rho, lambda. Weibull requires shape and scale.")}
        
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$shape)){stop("The list of parameters for Weibull marginals should include 'shape'")}
        if(!is.numeric(eval(parse(text = paste('params', i, sep ='')))$scale)){stop("The list of parameters for Weibull marginals should include 'scale'")}
      }
    }
  }
}


## calculating the cdf of wrapped Cauchy distribution given theta, mu, rho
cdf_wrpcauchy <- function(theta, mu, rho){
  require(CircStats)
  
  if((rho >= 1) | (rho < 0)){stop('rho should be in the interval [0,1)')}
  
  theta <- theta %% (2*pi)
  mu <- mu %% (2*pi)
  
  pdf <- function(x){dwrpcauchy(x, mu, rho)}
  cdf <- numeric(length(theta))
  
  for(i in 1:length(theta)){
    # cat('theta', theta[i], 'rho', rho, 'mu', mu, '\n')
    tryCatch(
      {
        cdf[i] <- integrate(pdf, 0, theta[i])$value
      }, 
      error=function(e) {
        message('An Error Occurred')
        cat('theta', theta[i], 'rho', rho, 'mu', mu, '\n')
        print(e)
      }
    )
  }
  
  return(cdf)
}

## calculating the inverse cdf of wrapped Cauchy distribution given p, mu, rho
inv_cdf_wrpcauchy <- function(p, mu, rho){
  if((rho >= 1) | (rho < 0)){stop('rho should be in the interval [0,1)')}
  
  mu <- mu %% (2*pi)
  
  inv_cdf <- mu + 2*atan((1-rho)/(1+rho)*tan(p/2))
  
  return(inv_cdf)
}

## calculating the pdf of Kato-Jones distribution given theta, mu, gamma, rho, lambda
dkatojones <- function(theta, mu, gamma, rho, lambda){
  alpha2 = rho*gamma*cos(lambda)
  beta2 = rho*gamma*sin(lambda)
  value = (1+2/rho*(alpha2*cos(theta-mu-lambda)-beta2*sin(theta-mu-lambda)-rho*alpha2)/(1+rho^2-2*rho*cos(theta-mu-lambda)))/(2*pi)
  return(value)
}

## calculating the cdf of Kato-Jones distribution given theta, mu, gamma, rho, lambda
cdf_katojones <- function(theta, mu, gamma, rho, lambda){
  a=((1-gamma*cos(lambda)/rho)*(theta-mu+pi)+gamma*sin(lambda)/rho*log((1+rho^2+2*rho*cos(lambda))/(1+rho^2-2*rho*cos(theta-mu-lambda)))+2*gamma*cos(lambda)/rho*(atan((1+rho)/(1-rho)*tan((theta-mu-lambda)/2))-atan((1+rho)/(1-rho)/tan(lambda/2))))/(2*pi)
  check <- tan((theta-mu-lambda)/2)< tan((-pi-lambda)/2)
  if(any(check)){
    a[check] = a[check] + gamma*cos(lambda)/rho
  }
  return(a)
}

## Function to find the value of F^{-1}(theta) for all three components according to the input of marginals
vector_to_uniform <- function(x, marginals = rep('wrapped cauchy', 3), params = NULL){
  x <- x %% (2*pi)
  if(marginals == 'wrapped cauchy'){
    uniform_x <- cdf_wrpcauchy(x, params$mu, params$rho)
  } else if (marginals == 'cardioid'){
    require(VGAM)
    uniform_x <- pcard(x, params$mu, params$rho)
    
  } else if (marginals == 'vonmises'){
    require(circular)
    suppressWarnings(
      uniform_x <- pvonmises(circular(x), circular(params$mu), params$kappa, from = circular(0))
    )
    
  } else if(marginals == 'katojones'){
    uniform_x <- cdf_katojones(x, params$mu, params$gamma, params$rho, params$lambda)
    
  } else if(marginals == 'weibull'){
    uniform_x <- pweibull(x, params$shape, params$scale)
  } 
  uniform_x <- uniform_x * 2*pi
  
  return(uniform_x)
}



## Examples
# Simulate from the trivariate wrapped Cauchy copula with marginals wrapped Cauchy, uniform, weibull
set.seed(2)
x <- rtri(1000, rho12 = 1, rho13 = 0.25, rho23 = 3, 
          marginals = c('wrapped Cauchy', 'uniform', 'weibull'), 
          params1 = list(mu = 1, rho = 0.2), params3 = list(shape = 1, scale = 2))

# Calculating the density for each of the 1000 datapoints generated
density <- dtri(x, rho12 = 1, rho13 = 0.25, rho23 = 3, 
               marginals = c('wrapped Cauchy', 'uniform', 'weibull'), 
               params1 = list(mu = 1, rho = 0.2), params3 = list(shape = 1, scale = 2))
density


