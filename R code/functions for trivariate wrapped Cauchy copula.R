copula_density <- function(x, rho12, rho13, rho23, marginals, params1 = NULL, params2 = NULL, params3 = NULL){
  dtri_cylinder(x, rho12, rho13, rho23, marginals = marginals, params1 = params1, params2 = params2, params3 = params3) * multiply_copula_denisty(x, marginals = marginals, params1 = params1, params2 = params2, params3 = params3)
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
      eval(parse(text = paste('dens', index[i], ' = (2*pi) * dvonmises(circular(x[,', index[i],']),circular(params', index[i],'$mu),params', index[i], '$kappa)', sep ='')))
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
dtri <- function(x, rho12, rho13, rho23, marginals = rep('uniform', 3), params1 = NULL, params2 = NULL, params3 = NULL){
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
  
  if (any(marginals == 'wrapped cauchy')){
    index <- which(marginals == 'wrapped cauchy')
    for(i in 1:length(index)){
      eval(parse(text = paste('u', index[i], ' = inv_cdf_wrpcauchy(u', index[i],',params', index[i],'$mu,params', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'cardioid')){
    require(VGAM)
    index <- which(marginals == 'cardioid')
    for(i in 1:length(index)){
      eval(parse(text = paste('u', index[i], ' = qcard(u', index[i],'/(2*pi),params', index[i],'$mu,params', index[i], '$rho)', sep ='')))
    }
  } 
  if (any(marginals == 'vonmises')){
    require(circular)
    index <- which(marginals == 'vonmises')
    for(i in 1:length(index)){
      eval(parse(text = paste('u', index[i], ' = qvonmises(u', index[i],'/(2*pi),circular(params', index[i],'$mu),params', index[i], '$kappa)', sep ='')))
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
  
  u1 <- u1 %% (2*pi); u2 <- u2 %% (2*pi); u3 <- u3 %% (2*pi)
  
  u <- matrix(c(u1, u2, u3), ncol = 3)
  return(u)
}