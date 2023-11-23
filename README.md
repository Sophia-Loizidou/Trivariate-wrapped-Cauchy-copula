# Trivariate-wrapped-Cauchy-copula

The trivariate wrapped Cauchy copula is proposed in LINK TO ARXIV PAPER and is a copula on the three dimensional torus that is given by:
$$ c(u_1,u_2,u_3) =  c_2 \Bigl\[ c_1 + 2 \left( \rho_{12} \cos (u_1 - u_2) + \rho_{13} \cos (u_1 - u_3) + \rho_{23} \cos (u_2 - u_3) \right) \Bigr\]^{-1}$$
for $0 \leq u_1,u_2, u_3 < 2\pi$
where $\rho_{12},\rho_{13},\rho_{23} \in \mathbb{R} \setminus \{0\} $, $\rho_{12}\rho_{13} \rho_{23} >0$.

The constants $c_1$ and $c_2$ are given by
$$c_1 = \frac{\rho_{12} \rho_{13}}{ \rho_{23}} + \frac{\rho_{12} \rho_{23}}{\rho_{13}} + \frac{ \rho_{13} \rho_{23} }{\rho_{12}}$$
and
$$c_2 = \frac{1}{(2\pi)^3} \left[ \left( \frac{\rho_{12} \rho_{13}}{ \rho_{23}} \right)^2 + \left( \frac{\rho_{12} \rho_{23}}{\rho_{13}} \right)^2 + \left(\frac{\rho_{13} \rho_{23}}{\rho_{12}} \right)^2 - 2 \rho_{12}^2 - 2 \rho_{13}^2 - 2 \rho_{23}^2 \right]^{1/2}.$$

We are assuming that there exists a permutation of $(1,2,3)$, $(j,k,\ell)$, such that $|\rho_{k \ell}| < |\rho_{jk} \rho_{j \ell}| / ( |\rho_{jk}| + |\rho_{j \ell}|)$, where $\rho_{kj} = \rho_{jk}$ for $1 \leq j < k \leq 3$.

## Properties of the copula

