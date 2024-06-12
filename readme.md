# Polygonal Faber-Krahn inequality: local minimality via validated computing

This repository contains the Matlab codes which provide a validated computing proof of a conjecture by Polya and Szego.

> **Polya Szego Conjecture** The regular $n$-gon minimizes the first Dirichlet Laplace eigenvalue among $n$ gons having fixed area.

The code accompanies the paper *The polygonal Faber-Krahn inequality: local minimality via validated computing* by Beniamin Bogosel and Dorin Bucur. 

The local minimality depends on the positivity of the eigenvalues of a particular Hessian matrix. The code computes these eigenvalues using validated computing (interval arithmetic) and outputs certified interval enclosures for these eigenvalues. 

The code is writte in Matlab and requires a working installation of the interval arithmetic library [INTLAB](https://www.tuhh.de/ti3/intlab/). 

Before running the functions listed below make sure that the subfolders are in the Matlab path. Run the command `addpath(genpath(pwd))` inside this folder or run the startup file provided. 

Main functions:

`PolyaHessInterval.m`: computes interval enclosures for the finite element approximations for the first two eigenvalues and the first eigenfunction for the discrete Laplace problem with Dirichlet boundary conditions on the regular $n$-gon. Usage:

`PolyaHessInterval(n,m)`, where `n` is the number of vertices and `m` is the number of division points of the mesh contained on a ray of the polygon. The mesh has $nm(m+1)/2$ points. Choose $n \leq 10$ and $m \in [100,400]$.  

`PolyaHessIntervalU(n,m)`: computes Hessian eigenvalues and corresponding interval enclosures. If $2n-4$ of the resulting intervals are contained in $(0,\infty)$ then the validation of local minimality succeeds. This should be run after the previous function since it uses information about the first two eigenvalues and the first eigenfunction. 

The certification of local minimality is successful for $n=5$, $m=250$ and $n=6$, $m=380$. The routines may take a long time, or even fail to execute for $m \geq 400$.

The `Tools` subfolder contains all necessary sub-routines. Some of them are of potential interest for other applications:

`verifyeig_BB.m`: a modification of `verifyeig` from INTLAB allowing to handle larger matrices (a matrix inverse is replaced with a few linear systems). `verifyeig_FEM.m` goes even further, replacing a linear system resolution which fails for large matrices with an iterative one. 