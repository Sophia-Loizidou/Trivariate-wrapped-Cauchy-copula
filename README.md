# Trivariate wrapped Cauchy copula

This GitHub repo contains codes that can be used to fit the trivariate wrapped Cauchy copula on data using MLE, with user-specified marginals. The choices for the marginals are 
* Uniform distribution
* wrapped Cauchy distribution
* von Mises distribution
* cardioid distribution
* Kato-Jones distribution
* weibull distribution (for linear components)

Any combination of the above distributions is possible.
Data generation is currently not supported for Kato-Jones marginals.


The trivariate wrapped Cauchy copula (TWCC) is proposed in S. Kato, C. Ley, S. Loizidou and K.V.Mardia (2024): `A versatile trivariate wrapped Cauchy copula with applications to toroidal and cylindrical data' (https://arxiv.org/abs/2401.10824) and is a copula on the three dimensional torus that is given by:

$$ c(u_1,u_2,u_3) =  c_2 \Bigl[ c_1 + 2 \left( \rho_{12} \cos (u_1 - u_2) + \rho_{13} \cos (u_1 - u_3) + \rho_{23} \cos (u_2 - u_3) \right) \Bigr]^{-1}$$

for $0 \leq u_1,u_2, u_3 < 2\pi$
where $\rho_{12},\rho_{13},\rho_{23} \in \mathbb{R} \setminus \{0\} $, $\rho_{12}\rho_{13} \rho_{23} >0$.

The constants $c_1$ and $c_2$ are given by
$$c_1 = \frac{\rho_{12} \rho_{13}}{ \rho_{23}} + \frac{\rho_{12} \rho_{23}}{\rho_{13}} + \frac{ \rho_{13} \rho_{23} }{\rho_{12}}$$
and
$$c_2 = \frac{1}{(2\pi)^3} \left[ \left( \frac{\rho_{12} \rho_{13}}{ \rho_{23}} \right)^2 + \left( \frac{\rho_{12} \rho_{23}}{\rho_{13}} \right)^2 + \left(\frac{\rho_{13} \rho_{23}}{\rho_{12}} \right)^2 - 2 \rho_{12}^2 - 2 \rho_{13}^2 - 2 \rho_{23}^2 \right]^{1/2}.$$

We are assuming that there exists a permutation of $(1,2,3)$, $(j,k,\ell)$, such that $|\rho_{k \ell}| < |\rho_{jk} \rho_{j \ell}| / ( |\rho_{jk}| + |\rho_{j \ell}|)$, where $\rho_{kj} = \rho_{jk}$ for $1 \leq j < k \leq 3$.

This repository includes two R files. The functions included in the 'MLE.R' file require the functions included in the 'functions for trivariate wrapped Cauchy copula.R' to run
