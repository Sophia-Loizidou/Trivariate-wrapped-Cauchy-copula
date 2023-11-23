# Trivariate-wrapped-Cauchy-copula

Let
\begin{align}
	\lefteqn{ c(u_1,u_2,u_3) =  c_2 \Bigl[ c_1 + 2 \left\{ \rho_{12} \cos (u_1 - u_2) + \rho_{13} \cos (u_1 - u_3) + \rho_{23} \cos (u_2 - u_3) \right\} \Bigr]^{-1}, } \label{eq:tri_density} \hspace{11cm} \\
	&  0 \leq u_1,u_2, u_3 < 2\pi, \nonumber
\end{align}
where $\rho_{12},\rho_{13},\rho_{23} \in \mathbb{R} \setminus \{0\} $, $\rho_{12}\rho_{13} \rho_{23} >0$,
\begin{equation}
c_1 = \frac{\rho_{12} \rho_{13}}{ \rho_{23}} + \frac{\rho_{12} \rho_{23}}{\rho_{13}} + \frac{ \rho_{13} \rho_{23} }{\rho_{12}}, \label{eq:c1}
\end{equation}
and
\begin{equation} \label{eq:c2}
c_2 = \frac{1}{(2\pi)^3} \left\{ \left( \frac{\rho_{12} \rho_{13}}{ \rho_{23}} \right)^2 + \left( \frac{\rho_{12} \rho_{23}}{\rho_{13}} \right)^2 + \left(\frac{\rho_{13} \rho_{23}}{\rho_{12}} \right)^2 - 2 \rho_{12}^2 - 2 \rho_{13}^2 - 2 \rho_{23}^2 \right\}^{1/2}.
\end{equation}
Suppose that there exists a permutation of $(1,2,3)$, $(j,k,\ell)$, such that $|\rho_{k \ell}| < |\rho_{jk} \rho_{j \ell}| / ( |\rho_{jk}| + |\rho_{j \ell}|)$, where $\rho_{kj} = \rho_{jk}$ for $1 \leq j < k \leq 3$.
Then the function (\ref{eq:tri_density}) is a probability density function on the three-dimensional torus $[0,2 \pi)^3$.