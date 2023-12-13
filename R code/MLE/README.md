This folder contains codes that can be used to fit the trivariate wrapped Cauchy copula on data using MLE, with user-specified marginals. The choices for the marginals are 
* Uniform distribution
* wrapped Cauchy distribution
* von Mises distribution
* cardioid distribution
* Kato-Jones distribution
* weibull distribution (for linear components)

Functions included in the file:
* mle.copula(): main function, outputs parameter estimates
* For estimating the mle parameters for the marginals:
    * mle.marginals.wrpcauchy()
    * mle.marginals.katojones()
    * mle.copula.params()
* cdf_wrpcauchy(): calculating the cdf of wrapped Cauchy distribution given theta, mu, rho
* inv_cdf_wrpcauchy(): calculating the inverse cdf of wrapped Cauchy distribution given p, mu, rho
* dkatojones(): calculating the pdf of Kato-Jones distribution given theta, mu, gamma, rho, lambda
* cdf_katojones(): calculating the cdf of Kato-Jones distribution given theta, mu, gamma, rho, lambda
* vector_to_uniform(): finding the value of F^{-1}(theta) for all three components according to the input of marginals
* check_parameters(): Function for checking that all inputs are allowed

The functions included in the file "functions for trivariate wrapped Cauchy copula" are also required for this to run.