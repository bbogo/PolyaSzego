# Computing upper bounds for interpolation constants using Morley finite elements

Based on the paper: *On the interpolation constants over triangular elements* by K. Kobayashi.

Requires a working installation of Intlab, since interval arithmetics is used in the certification process. 

The interpolation constant of interest verifies 
$$ \|\nabla u - \nabla \Pi_h(u)\|_{L^2(\Omega)} \leq C |u|_{H^2(\Omega)}$$
where $\Pi_h$ is the $P_1$ Lagrange interpolation operator. The optimal constant $C$ is the solution of an eigenvalue problem involving the bi-Laplacian operator. The paper cited above shows that a discrete eigenvalue related to Morley finite elements can give a certified upper bound for the interpolation constant, accurate enough even when using relatively small meshes. The constant $C$ scales linearly with the size of the triangle.

`ComputeConstants` computes the constants for the triangles needed regarding the Polya-Szego conjecture. 

`P1InterpolationConstant` computes a certified upper bound for the interpolation constant on the triangle with vertices $(0,0),(1,0),(a,b)$.

`AssemblyMorleyTriInterval` assembles the rigidity matrices for the first and second derivatives of Morley elements, following the ideas written in the paper B. Bogosel, D. Bucur, *Polygonal Faber-Krahn inequality: local minimality via validated computing*. 