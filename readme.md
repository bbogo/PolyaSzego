# Polygonal Faber-Krahn inequality: local minimality via validated computing

This repository contains the Matlab codes which provide a validated computing proof of a conjecture by Polya and Szego.

> **Polya Szego Conjecture** The regular $n$-gon minimizes the first Dirichlet Laplace eigenvalue among $n$ gons having fixed area.

The code accompanies the paper *The polygonal Faber-Krahn inequality: local minimality via validated computing* by Beniamin Bogosel and Dorin Bucur. 

The local minimality depends on the positivity of the eigenvalues of a particular Hessian matrix. The code computes these eigenvalues using validated computing (interval arithmetic) and outputs certified interval enclosures for these eigenvalues. 

The code is writte in Matlab and requires a working installation of the interval arithmetic library [INTLAB](https://www.tuhh.de/ti3/intlab/). 