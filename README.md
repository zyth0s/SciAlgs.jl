# SciAlgs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://zyth0s.github.io/SciAlgs.jl/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://zyth0s.github.io/SciAlgs.jl/latest)
[![Build Status](https://travis-ci.com/zyth0s/SciAlgs.jl.svg?branch=master)](https://travis-ci.com/zyth0s/SciAlgs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/zyth0s/SciAlgs.jl?svg=true)](https://ci.appveyor.com/project/zyth0s/SciAlgs-jl)
[![Coveralls](https://coveralls.io/repos/github/zyth0s/SciAlgs.jl/badge.svg?branch=master)](https://coveralls.io/github/zyth0s/SciAlgs.jl?branch=master)

SciAlgs contains a set of common scientific algorithms written in a 21st century high-level programming language.
Ideally it would be like a Numerical Recipes of our age. However, acomplishing that task would be daunting and is not necesarily the goal.

Some of the aims would be:
* To rewrite some old programs to avoid their obsolescence.
* Provide simple functions to perform basic operations, such that an average researcher can understand.
* Modern development practices: tests, continuous integration, collaborative development, ...
* Make algorithms accessible to researchers without expertise in computer programming.
* Enrich the ecosystem of scientific packages in Julia planting a seed.


Some of the topics that are going to be covered are:
* Crystallography
* Quantum Chemistry
* Condensed Matter Physics
* Quantum Mechanics models
* Numerical methods


Algorithms
===========

- [x] Linear least squares minimization http://mathworld.wolfram.com/LeastSquaresFitting.html
- [x] Theil's method to fit data to straight line. https://pubs.acs.org/doi/pdf/10.1021/ed082p1472.2
- [x] Heron's method to find roots of function. https://press.princeton.edu/books/hardcover/9780691146867/numerical-analysis
- [x] Cartesian to polar conversion.
- [x] Find closest atom
* Numerical quadrature (NumQuad)
  - [x] 1D quadratures:
    + [x] Trapezoidal quadrature
    + [x] Euler-McLaurin quadrature
    + [x] Clenshaw-Curtis quadrature
    + [x] Gauss-Chebyshev 1st kind
    + [x] Gauss-Chebyshev 2nd kind
    + [x] Gauss-Legendre
    + [x] Pérez-Jordá
    + [x] Multiexp quadrature
    + [x] Maxwell quadrature
- [x] Finite difference solution of 1D single particle Schrödinger equation with 
  + [x] Two point central difference formula
  + [x] Matrix Numerov method https://aapt.scitation.org/doi/10.1119/1.474
- [x] Kronig-Penney model
- [x] Brusselator
* Crystallography (Xtal)
  - [x] Calculate metric tensor G from cell parameters
  - [x] Calculate unit cell volume from cell parameters or metric tensor
  - [x] Calculate lattice vectors from cell parameters
  - [x] Calculate reciprocal cell parameters
  - [x] Calculate reciprocal vectors
  - [x] Calculate interplanar spacing
  - [x] Convert between cartesian and fractionary coordinates
  - [x] Convert mapping to transformation matrix+vector (Seitz symbol)
  - [x] Charazterize a crystallographic symmetry operation
  - [x] Apply crystallographic symmetry operation to a point 
  - [x] Listing of planes that diffract X-rays in non triclinic systems

