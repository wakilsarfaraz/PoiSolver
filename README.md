# PoiSolver
PoiSolver is a Finite Element Solver written in R for Poisson equation. PoiSolver employes piecewise linear basis functions on uniform triangulation. It solves a Poisson equation on a two dimensional rectangular domain with homogeneous boundary conditions of Dirichlet type. We illustrate PoiSolver through an example in which we present both the relevant theory of the finite element method a long with the implementation for a given example.
## Example problem
Consider a two dimensional bounded convex domain in the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega=\[0\1\]\times\[0\1\]\subset\mathbb{R}^2"/> on the positive real <img src="https://latex.codecogs.com/svg.latex?\Large&space;x-y"/> plane. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y)"/> satisfy the probelm 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\,u(x,y)=f(x,y)"/> 
in <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/> 
and 
<img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y)=0"/> on <img src="https://latex.codecogs.com/svg.latex?\Large&space;\partial\Omega"/>. 
