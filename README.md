# PoiSolver
PoiSolver is a Finite Element Solver written in R for Poisson equation. PoiSolver employes piecewise linear basis functions on uniform triangulation. It solves a Poisson equation on a two dimensional rectangular domain with homogeneous boundary conditions of Dirichlet type. We illustrate PoiSolver through an example in which we present both the relevant theory of the finite element method a long with the implementation for a given example.
## Finite element method overview
Consider a two dimensional open bounded [convex](https://en.wikipedia.org/wiki/Convex_set) domain denoted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega\subset\mathbb{R}^2"/> on the positive real <img src="https://latex.codecogs.com/svg.latex?\Large&space;x-y"/> plane. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y)"/> satisfy the probelm <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\,u(x,y)=f(x,y)"/> in <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y)=0"/> on <img src="https://latex.codecogs.com/svg.latex?\Large&space;\partial\Omega"/>. Let the solution and test spaces be denoted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}=\{w\in\,H^1:\|w(x,y)\|^2<\infty,\,\|\nabla\,w(x,y)\|^2<\infty\}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}=\{w\in\,H^1:\,w(x,y)=0,\;\,x,y\in\partial\Omega\}"/> respectively, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are a special class of [Hilbert spaces](https://en.wikipedia.org/wiki/Hilbert_space) known as [Sobelov spaces](https://en.wikipedia.org/wiki/Hilbert_space). Assuming that the problem satisfies all the conditions stated by [Lax-Milgram Theorem](http://mathworld.wolfram.com/Lax-MilgramTheorem.html), we multiply both sides of the equation by a test function <img src="https://latex.codecogs.com/svg.latex?\Large&space;w(x,y)"/> and integrate by parts with application of [Green's Theorem](https://en.wikipedia.org/wiki/Green%27s_theorem), then the weak formulation of the problem is to find <img src="https://latex.codecogs.com/svg.latex?\Large&space;u\in\,H^1_{\Omega}"/> such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega}\nabla\!u\cdot\nabla\!w\,d\Omega=\int_{\Omega}fw\,d\Omega"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;w\in\,H^1_{\Omega_0}"/>. Given that the problem is [well-posed](https://en.wikipedia.org/wiki/Well-posed_problem) in Hilbert space of functions of which both <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are subsets, then we can define finite dimensional solution and test spaces namely <img src="https://latex.codecogs.com/svg.latex?\Large&space;V^h\subset\,H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;V_0^h\subset\,H^1_{\Omega_0}"/> respectively. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/> be a triangulated approximation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, containing <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> discrete number of nodes and <img src="https://latex.codecogs.com/svg.latex?\Large&space;K"/> triangles, then using the finite dimentional solution and test spaces, the finite element formulation of the problem is to find <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h\in\,V^h"/> such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega^h}\nabla\!u^h\cdot\nabla\!w^h\,d\Omega=\int_{\Omega^h}fw^h\,d\Omega"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;\,w^h\in\,V_0^h"/>.  We expand <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h"/> in terms of its finite i.e. <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> dimensional [basis functions](https://en.wikipedia.org/wiki/Basis_function) in the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h=\sum_{i=1}^N\,U_i\phi_i"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h=\sum_{j=1}^N\,\phi_j"/>. Substituting these in the finite element formulation we obtain a discrete set of linear equations of the from <img src="https://latex.codecogs.com/svg.latex?\Large&space;S\,U=L"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> is the stiffness matrix, <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> is the vector of unknowns containing the approximate solution values at each node and <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> is the right hand-side load vector. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> are given by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[S\]_{i,j}=\int_{\Omega^h}\,\nabla\phi_i\cdot\nabla\phi_j\,d\Omega^h"/>, <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[F\]_{j}=\int_{\Omega^h}\,f\phi_j\,d\Omega^h"/>. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> are the finite element approximate solution values to be found from the linear system. 
## Code Implementation
Using a uniform scheme we proceed to descritise  <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>. This is achieved by creating an <img src="https://latex.codecogs.com/svg.latex?\Large&space;N\times\,N"/> grid points within <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, which provides a uniform quadrilateral mesh. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathcal{T}"/> denote the uniform triangulation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, consisting of triangles. The algorithm is programmed to compute locally on each triangle <img src="https://latex.codecogs.com/svg.latex?\Large&space;K\in\mathcal{T}"/>, the entries of all the matrices <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and vector <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/>. The entries of the local <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times3"/> matrix <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and  <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times1"/> vector <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> are constructed by the formulea <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[S\]_{i,j}=\sum_{K}\int_{\Omega^h}\nabla\phi_i\cdot\nabla\phi_j\,d\Omega"/>  and <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[L\]_j=\sum_K\int_{\Omega^h}f\phi_j\,d\Omega"/>.
Consider a two dimensional open bounded [convex](https://en.wikipedia.org/wiki/Convex_set) domain denoted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega\subset\mathbb{R}^2"/> on the positive real <img src="https://latex.codecogs.com/svg.latex?\Large&space;x-y"/> plane. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y)"/> satisfy the probelm <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\,u(x,y)=f(x,y)"/> in <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y)=0"/> on <img src="https://latex.codecogs.com/svg.latex?\Large&space;\partial\Omega"/>. Let the solution and test spaces be denoted by <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}=\{w\in\,H^1:\|w(x,y)\|^2<\infty,\,\|\nabla\,w(x,y)\|^2<\infty\}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}=\{w\in\,H^1:\,w(x,y)=0,\;\,x,y\in\partial\Omega\}"/> respectively, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are a special class of [Hilbert spaces](https://en.wikipedia.org/wiki/Hilbert_space) known as [Sobelov spaces](https://en.wikipedia.org/wiki/Hilbert_space). Assuming that the problem satisfies all the conditions stated by [Lax-Milgram Theorem](http://mathworld.wolfram.com/Lax-MilgramTheorem.html), we multiply both sides of the equation by a test function <img src="https://latex.codecogs.com/svg.latex?\Large&space;w(x,y)\in\H^1_{\Omega_0}"/> and integrate by parts with application of [Green's formula](https://math.stackexchange.com/questions/1936916/greens-formula-integration-in-parts), then the weak formulation of the problem is to fine <img src="https://latex.codecogs.com/svg.latex?\Large&space;u\in\,H^1"/> such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega}\nabla\!u\cdot\nabla\!w\,d\Omega=\int_{\Omega}fw\,d\Omega"/> for all  <img src="https://latex.codecogs.com/svg.latex?\Large&space;w\in\,H_{\Omega_0}^1"/>. Given that the problem is [well-posed](https://en.wikipedia.org/wiki/Well-posed_problem) in Hilbert space of functions of which both <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;H^1_{\Omega_0}"/> are subsets, then we can define finite dimensional solution and test spaces namely <img src="https://latex.codecogs.com/svg.latex?\Large&space;V^h\subset\,H^1_{\Omega}"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;V_0^h\subset\,H^1_{\Omega_0}"/> respectively. Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/> be a triangulated approximation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, containing <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> discrete number of nodes and <img src="https://latex.codecogs.com/svg.latex?\Large&space;K"/> triangles, then using the finite dimentional solution and test spaces, the finite element formulation of the problem is to find <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h\in\,V^h"/> such that <img src="https://latex.codecogs.com/svg.latex?\Large&space;\int_{\Omega^h}\nabla\!u^h\cdot\nabla\!w^h\,d\Omega=\int_{\Omega^h}fw^h\,d\Omega"/> for all <img src="https://latex.codecogs.com/svg.latex?\Large&space;\,w^h\in\,V_0^h"/>.  We expand <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h"/> in terms of its finite i.e. <img src="https://latex.codecogs.com/svg.latex?\Large&space;N"/> dimensional [basis functions](https://en.wikipedia.org/wiki/Basis_function) in the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;u^h=\sum_{i=1}^N\,U_i\phi_i"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;w^h=\sum_{j=1}^N\,\phi_j"/>. Substituting these in the finite element formulation we obtain a discrete set of linear equations of the from <img src="https://latex.codecogs.com/svg.latex?\Large&space;S\,U=L"/>, where <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> is the stiffness matrix, <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> is the vector of unknowns containing the approximate solution values at each node and <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> is the right hand-side load vector. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> are given by <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[S\]_{i,j}=\int_{\Omega^h}\,\nabla\phi_i\cdot\nabla\phi_j\,d\Omega^h"/>, <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[F\]_{j}=\int_{\Omega^h}\,f\phi_j\,d\Omega^h"/>. The entries of <img src="https://latex.codecogs.com/svg.latex?\Large&space;U"/> are the finite element approximate solution values to be found from the linear system. 
For code implementation we start with a scheme to descritise  <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>. This is achieved by creating an <img src="https://latex.codecogs.com/svg.latex?\Large&space;N\times\,N"/> grid points within <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, which provides a uniform quadrilateral mesh. Then we triangulate the quadrilateral mesh by drawing a line of slope -1 diangonally through each square in order to divide every quadrilateral element into two adjacent triangles, which leads to the construction of a uniform triangulated domain <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/>.  Let <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mathcal{T}"/> denote the uniform triangulation of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>, consisting of triangles. The algorithm is programmed to compute locally on each triangle <img src="https://latex.codecogs.com/svg.latex?\Large&space;K\in\mathcal{T}"/>, the entries of the matrix <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and the load vector <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/>. The entries of the local <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times3"/> matrix <img src="https://latex.codecogs.com/svg.latex?\Large&space;S"/> and  <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times1"/> vector <img src="https://latex.codecogs.com/svg.latex?\Large&space;L"/> are constructed by the formulea <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[S\]_{i,j}=\sum_{K}\int_{\Omega^h}\nabla\phi_i\cdot\nabla\phi_j\,d\Omega"/>  and <img src="https://latex.codecogs.com/svg.latex?\Large&space;\[L\]_j=\sum_K\int_{\Omega^h}f\phi_j\,d\Omega"/>.
# Documentation
We go through the implementation process for the problem <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Delta\,u(x,y)=1"/> that is posed on <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega=\[0\;1\]\times\[0\;1\]"/>. It satisfies homogeneous Dirichlet type boundary conditions of the form <img src="https://latex.codecogs.com/svg.latex?\Large&space;u(x,y)=0\quad\,x,y\in\partial\Omega"/>.
The first 11 lines of the code mainly assign numerical values for variables to be used. The variables that user may require to interfere with in PoiSolver are 
``` r 
L = 1
```
, which stores the value for the side length of <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega"/>. The current code is implemented so that it only works for rectangular domains that have the same side lengths i.e. square only. The variable 
``` r
N = 50
``` 
is the mesh refinement controlling parameter i.e. it is the number of points that discretises `L`. The line 
``` r
X = seq(0,L,len=N+1)
```
is equivalent to `linspace` function in MATLAB. It outputs the actual values of the `N+1` grid points. We can creat a quadrilateral grid that consists of <img src="https://latex.codecogs.com/svg.latex?\Large&space;N^2"/> points each of which is equipped with the <img src="https://latex.codecogs.com/svg.latex?\Large&space;(x,y)"/> numerical values for the coordinates in two matrices namely `x` and `y`, which were obtained using the following built-in functions of R.   
``` r
m = length(X); n=length(X);
x = matrix(rep(X,each=n),nrow=n);
y = matrix(rep(X,m),nrow=n)
```
Now we reshape the matrices <img src="https://latex.codecogs.com/svg.latex?\Large&space;x,y"/> into column vectors, which is accomplished by  
``` r
x = c(x)
y = c(y)
```  
The variables `GNodes` is meant to stand for Global Nodes and it is instantiated in terms of `N`. The variable `NumTRI` stores the number of triangles that a discretisation of this type outputs. If we have <img src="https://latex.codecogs.com/svg.latex?\Large&space;N^2"/> quadrilateral elements and each square is divided by two triangles we obtain <img src="https://latex.codecogs.com/svg.latex?\Large&space;2\times\,N^2"/> triangles. The array `LocNodes` stores the connectivity array of the vertices of all the triangles, hence it is a matrix with `NumTRI` rows and 3 columns.
``` r
GNodes = (N+1)^2
NumTRI = 2*N^2
LocNodes = matrix(0,NumTRI,3)
```
Triangulation of the domain is obtained such that the vertices of each triangle is locally counted in anti-clockwise orientation and such that the local counting of vertices of each triangle fills the global connectivity array `LocNodes` in the correct order. To achieve this we need a nested for loop which is given by 
``` r
for (i in 1:N){
	for (j in 1:N){
		    LocNodes[i+2*(j-1)*N,1] = i+(j-1)*(N+1)
        LocNodes[i+2*(j-1)*N,2] = i+j*(N+1)
        LocNodes[i+2*(j-1)*N,3] = (i+1)+(j-1)*(N+1)
        LocNodes[i+N+2*(j-1)*N,1] = i+1+j*(N+1)
        LocNodes[i+N+2*(j-1)*N,2] = (i+1)+(j-1)*(N+1)
        LocNodes[i+N+2*(j-1)*N,3] = i+j*(N+1)
    }
}
```
We need to introduce two types of global arrays namely a system matrix `Sparsity` and a right hand-side vector `LoadVect`, which must be of size <img src="https://latex.codecogs.com/svg.latex?\Large&space;(N+1)^2\times\,(N+1)^2"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;(N+1)^2\times\,1"/> respectively. These are introduced by 
``` r
Sparsity <- matrix(0, GNodes, GNodes)
LoadVect = matrix(0, GNodes, 1)
```
The next segment of code is a for loop that executes a series of computations on each one of the triangles. It starts with 
``` r
for (n in 1:NumTRI)
	{
```
The first part within the loop defines the local position vectors `r1`, `r2`, and `r3` in terms of `x` and `y` coordinate values through which we refer to different vertices of each triangle. 
``` r
r1 = matrix(c(x[LocNodes[n,1]],y[LocNodes[n,1]]),nrow=2, byrow=FALSE)
r2 = matrix(c(x[LocNodes[n,2]],y[LocNodes[n,2]]),nrow=2, byrow=FALSE);
r3 = matrix(c(x[LocNodes[n,3]],y[LocNodes[n,3]]),nrow=2, byrow=FALSE);
```
The second part within the loop defines a <img src="https://latex.codecogs.com/svg.latex?\Large&space;2\times2"/> jacobian matrix `J` for the mapping from an arbitrary triangle in <img src="https://latex.codecogs.com/svg.latex?\Large&space;\Omega^h"/> to a reference triangle that has vertices in order <img src="https://latex.codecogs.com/svg.latex?\Large&space;(0,0),\!(1,0)"/> and <img src="https://latex.codecogs.com/svg.latex?\Large&space;(0,1)"/>.
``` r
J = matrix(c(r2[1]-r1[1],r2[2]-r1[2]
    ,r3[1]-r1[1],r3[2]-r1[2]), nrow=2, byrow=TRUE) 
```
The Jacobian of the mapping namely `J` serves to reduce the computational cost by a significant amount, particularly due to a property of integration for computing the integral of a function on a reference domain with a given mapping between  the arbitrary domain and the reference domain. Further details on this topic can be found on [Integral domain transformation](http://www.iue.tuwien.ac.at/phd/nentchev/node58.html). 
The variable name `Astiff` is the local <img src="https://latex.codecogs.com/svg.latex?\Large&space;3\times3"/> stiffness matrix, which is computed and looped over all the triangles.
``` r
 Astiff = (1/(2*det(J)))*matrix(c(sum((r2-r3)*(r2-r3))
          ,sum((r2-r3)*(r3-r1)),sum((r2-r3)*(r1-r2))
          ,sum((r2-r3)*(r3-r1)),sum((r3-r1)*(r3-r1))
          ,sum((r3-r1)*(r1-r2)),sum((r2-r3)*(r1-r2))
          ,sum((r3-r1)*(r1-r2)),sum((r1-r2)*(r1-r2))), nrow=3, byrow=TRUE);
```



