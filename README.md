# Spline-MHE
[![](https://img.shields.io/badge/MATLAB-R2024a-blue)](https://www.mathworks.com/campaigns/products/trials.html?prodcode=ML&s_tid=prod_getml)
[![](https://img.shields.io/badge/DOI-Paperlink-orange)](https://www.mathworks.com/campaigns/products/trials.html?prodcode=ML&s_tid=prod_getml)
![](Spline-MHE/docs/Lorenz_animation.gif)

`Spline-MHE` is a **Smoothing Spline Nonlinear Moving Horizon Estimation** framework for fast and consistently convergent parameter identification and initial-state estimation of nonlinear dynamic systems of the form,
```math
\begin{aligned}\dot{x}( t) =f( x( t) ,u(t),\theta ),\\ y( t) =g( x( t) ,u(t),\theta ), \end{aligned}
```
where:
* $x(t)\in \mathbb{R}^{n_{x}}$ is the $n_{x}$ dimensional state vector,
* $y(t)\in \mathbb{R}^{n_{y}}$ is the $n_{y}$ dimensional output vector, and
* $\theta \in \mathbb{R}^{n_{\theta }}$ is the $n_{\theta}$ dimensional output vector
