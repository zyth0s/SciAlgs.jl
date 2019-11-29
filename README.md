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
* Reproduce figures/data of published papers, mainly educative ones.


Some of the topics that are going to be covered are:
* Crystallography
* Quantum Chemistry
* Condensed Matter Physics
* Quantum Mechanics models
* Numerical methods


Algorithms
===========

* Statistics
  - [x] Linear least squares minimization http://mathworld.wolfram.com/LeastSquaresFitting.html
  - [x] Theil's method to fit data to straight line. 
        [J. Chem. Educ. 2005, 82, 10, 1472](https://pubs.acs.org/doi/10.1021/ed082p1472.2)
* Optimization
  - [x] Heron's method to find roots of function. 
        [Numerical Analysis by L. Ridgway Scott](https://press.princeton.edu/books/hardcover/9780691146867/numerical-analysis)
* Geometry
  - [x] Cartesian to polar conversion.
  - [x] Find closest atom
* Numerical Quadrature (NumQuad)
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
* Quantum Mechanics
  - [x] Finite difference solution of 1D single particle Schrödinger equation with 
    + [x] Two point central difference formula 
          [J. Chem. Educ. 2017, 94, 6, 813-815](https://pubs.acs.org/doi/10.1021/acs.jchemed.7b00003)
    + [x] Matrix Numerov method [American Journal of Physics 80, 1017 (2012)](https://aapt.scitation.org/doi/10.1119/1.474)
  - [x] Kronig-Penney model
  - [x] Linear variational calculation 
        [Eur. J. Phys. 31 (2010) 101–114](https://iopscience.iop.org/article/10.1088/0143-0807/31/1/010/meta)
* Electrochemistry
  - [x] Linear Sweep Voltammogram 
        [J. Chem. Educ. 2019, 96, 10, 2217-2224](https://pubs.acs.org/doi/abs/10.1021/acs.jchemed.9b00542)
  - [x] "Lifelike" Linear Sweep Voltammogram
        [J. Chem. Educ. 2000, 77, 1, 100](https://pubs.acs.org/doi/10.1021/ed077p100)
* Chemical Kinetics
  - [x] Brusselator
* Crystallography (Xtal) 
  [Symmetry Relationships between Crystal Structures by Ulrich Müller]( https://global.oup.com/academic/product/symmetry-relationships-between-crystal-structures-9780198807209?cc=de&lang=en&)
  - [x] Calculate metric tensor G from cell parameters and lattice vectors
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

