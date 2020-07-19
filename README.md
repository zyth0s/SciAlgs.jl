# SciAlgs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://zyth0s.github.io/SciAlgs.jl/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://zyth0s.github.io/SciAlgs.jl/latest)
[![Build Status](https://travis-ci.com/zyth0s/SciAlgs.jl.svg?branch=master)](https://travis-ci.com/zyth0s/SciAlgs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/zyth0s/SciAlgs.jl?svg=true)](https://ci.appveyor.com/project/zyth0s/SciAlgs-jl)
[![Coveralls](https://coveralls.io/repos/github/zyth0s/SciAlgs.jl/badge.svg?branch=master)](https://coveralls.io/github/zyth0s/SciAlgs.jl?branch=master)

SciAlgs is a compilation of fundamental scientific algorithms.
The collection is not intended to be exhaustive but to offer clear and concise implementations.

A disgression of the reasons for choosing Julia follows.
Translating formulas to code is a non trivial task. This process was first addressed
by Fortran (that is, Formula translator) with great success, and thus Fortran has found widespread 
application until nowadays. However, some of its characteristics (not interactive; low level string manipulation;
lack of standard library; attached to bad programming practices) make it less attractive to young scientists
and professional programmers. In particular, both the compilation/execution
process and the verbosity of the language slow down the implementation and distract the programmer
from scientific to technical questions.
Interactive and modern languages allow a more agile development and a more concise
implementation of the algorithms.
Matlab has such functionality but it is not open source. On the other hand,
Octave offers similar capabilties but it is too focused on array handling like Matlab, with less
attention to other general purpose programming tasks and the developement of new packages has stagnated.
Python has been the alternative but numerical arrays are not first class objects in the language and
optimization usually lends to rewriting code in other language.
Julia offers an easy and natural array language like Fortran (no need of extra libraries to handle numerical
arrays; same 1 index convention) whereas it also includes modern programming language capabilities.

Ideally, the implementation would be in only one language (here Julia). In case that a routine
is not fast enough, the kernel could be written in Fortran (0-indexation in C could be a source of trouble).

Some of the aims would be:
* Concise and pedagogical implementations
* Make algorithms accessible to researchers without expertise in computer programming.
* Rewrite legacy programs more clearly (e.g. no spaguetti code)
* Provide useful functions to perform basic operations
* Integrate modern development practices: tests, continuous integration, collaborative development, ...
* Enrich Julia's package ecosystem with new functionalities
* Reproduce figures/data of published papers (mainly educative ones).


Some of the topics that are going to be covered are:
* Crystallography
* Quantum Chemistry
* Condensed Matter Physics
* Quantum Mechanics models
* Numerical methods


Algorithms
===========

* [Algorithms from THE BOOK](https://bookstore.siam.org/ot168/bonus) by Kenneth Lange
* [Computational Fluid Dynamics Course](https://github.com/surajp92/CFD_Julia)
* Statistics [Stat]
  - [x] Linear least squares minimization http://mathworld.wolfram.com/LeastSquaresFitting.html
  - [x] Theil's method to fit data to straight line. 
        [J. Chem. Educ. 2005, 82, 10, 1472](https://pubs.acs.org/doi/10.1021/ed082p1472.2)
  - [x] Sample size needed to have certain maximum error.
* Optimization, Analysis
  - [x] Heron's method to find roots of function. 
        [Numerical Analysis by L. Ridgway Scott](https://press.princeton.edu/books/hardcover/9780691146867/numerical-analysis)
  - [x] Illustrated Fourier Transform [J. Appl. Cryst. (2007). 40, 1153–1165](https://doi.org/10.1107/S0021889807043622)
* Geometry
  - [x] Cartesian to polar conversion.
  - [x] Find closest atom
  - [x] Cartesian to internal coordinates (Z-matrix)
  - [x] Read a XYZ formatted file
* Numerical Quadrature [NumQuad]
  - [x] 1D quadratures:
    + [x] Trapezoidal quadrature
    + [x] Euler-McLaurin quadrature
    + [x] Clenshaw-Curtis quadrature
    + [x] Gauss-Chebyshev 1st kind
    + [x] Gauss-Chebyshev 2nd kind
    + [x] Gauss-Legendre
    + [x] Pérez-Jordá
          [Comput. Phys. Commun., 77, 1, 1993, 46-56](https://www.sciencedirect.com/science/article/pii/001046559390035B?via%3Dihub)
    + [x] Multiexp quadrature
          [J. Comput. Chem. 24 (2003) 732-740](https://www.onlinelibrary.wiley.com/doi/abs/10.1002/jcc.10211)
    + [x] Maxwell quadrature
* Quantum Mechanics [QM]
  - [x] Finite difference solution of 1D single particle Schrödinger equation with 
          (i) infinite potential well
          (ii) finite potential well
          (iii) double finite well (unequal depth)
          (iv) harmonic well
          (v) Morse well, and
          (vi) Kronig-Penney finite wells, using
    + [x] Two point central difference formula 
          [J. Chem. Educ. 2017, 94, 6, 813-815](https://pubs.acs.org/doi/10.1021/acs.jchemed.7b00003)
    + [x] Matrix Numerov method [American Journal of Physics 80, 1017 (2012)](https://aapt.scitation.org/doi/10.1119/1.474)
  - [x] Kronig-Penney model
  - [x] Linear variational calculation 
        [Eur. J. Phys. 31 (2010) 101–114](https://iopscience.iop.org/article/10.1088/0143-0807/31/1/010/meta)
  - [x] Tight-binding 2/3 centers
  - [x] Tight-binding 1D homoatomic chain/ring + impurity. Surface state.
  - [x] Tight-binding 1D heteroatomic chain/ring with s, and s&p orbitals
  - [x] Tight-binding in 2D/3D homoatomic
  - [x] Green's function for a free particle
  - [x] Green's function for a two level system
  - [x] Green's function method for H2 [Revista Brasileira de Ensino de Fisica, vol. 39, no 1, e1303, 2017](http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1806-11172017000100402&lng=en&tlng=en)
  - [x] Green's function method for an infinite chain [Revista Brasileira de Ensino de Fisica, vol. 39, no 1, e1303, 2017](http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1806-11172017000100402&lng=en&tlng=en)
  - [x] Green's function method (surface-bulk recursive) for an semi-infinite chain [Revista Brasileira de Ensino de Fisica, vol. 39, no 1, e1303, 2017](http://www.scielo.br/scielo.php?script=sci_arttext&pid=S1806-11172017000100402&lng=en&tlng=en)
  - [x] Ising chain (1D) model with direct diagonalization/exact diagonalization/fullCI
* Electronic Structure [ElStruct]
  - [x] Slater-Koster Tight-binding of Sr2RuO4 [PRL 116, 197003 (2016)](https://link.aps.org/doi/10.1103/PhysRevLett.116.197003)
  - [x] McMurchie-Davidson molecular integrals evaluation scheme
  - [x] Simple plane-wave Density Functional Theory (DFT)
  - [x] Simple Hartree-Fock (HF) with GTOs (H2,He) [Szabo-Ostlund. Modern Quantum Chemistry.](https://store.doverpublications.com/0486691861.html)
  - [x] Simple Hartree-Fock (HF) with GTOs (H2,He) [J. Thijssen. Computational Physics](https://www.cambridge.org/core/books/computational-physics/BEE73B0139D4A9993193B57CDC62096E)
  - [x] Hartree-Fock (HF) with PySCF AO integrals
  - [x] Moller-Plesset order 2 (MP2) with PySCF AO integrals
  - [x] Davidson diagonalization [J. Comput. Phys. 17, 87–94](https://doi.org/10.1016/0021-9991(75)90065-0)
* Spectroscopy
  - [x] X-ray Photoelectron Spectrocopy (XPS) [J. Chem. Educ. 2019, 96, 7, 1502-1505](https://pubs.acs.org/doi/10.1021/acs.jchemed.9b00236)
* Electrochemistry
  - [x] Linear Sweep Voltammogram 
        [J. Chem. Educ. 2019, 96, 10, 2217-2224](https://pubs.acs.org/doi/abs/10.1021/acs.jchemed.9b00542)
  - [x] "Lifelike" Linear Sweep Voltammogram
        [J. Chem. Educ. 2000, 77, 1, 100](https://pubs.acs.org/doi/10.1021/ed077p100)
* Chemical Kinetics
  - [x] Brusselator
* Crystallography [Xtal]
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
  - [x] Probe wavelengths
* Epidemiology [Epidemics]
  - [x] SIS model
  - [x] SIS Discrete Time Markov Chain model
  - [x] SIS Continuous Time Markov Chain model
  - [x] SIR model (final size also)
  - [x] SIR Discrete State Discrete Time Markov Chain (Chain Binomial) model
  - [x] SIR Discrete State Continuous Time Markov Chain model
  - [x] SEIR model
  - [x] SEIR model (Erlang)
  - [x] Microparasite scaling model [Nature volume 379, 720–722(1996)](http://www.nature.com/articles/379720a0)
  - [x] SEIRC model
  - [x] SEIRV model [A mathematical model for the novel coronavirus epidemic in Wuhan, China](http://www.aimspress.com/article/10.3934/mbe.2020148)
  - [x] SIR with age segregation
  - [x] SQLIHUHURF model for COVID-19 in Spain
* Health
  - [x] NUTRI-SCORE see [Public health panorama, 03 (04), 712 - 725](https://apps.who.int/iris/handle/10665/325207)
  - [x] Dog-to-human age [bioRxiv doi: https://doi.org/10.1101/829192](https://www.biorxiv.org/content/10.1101/829192v2)
* Astro (-nomy, -physics) [Astro] 
  - [x] Gauss' Easter algorithm [Mathematics Magazine, 92:2, 91-98](https://doi.org/10.1080/0025570X.2019.1549889)
  - [x] Meesus-Jones-Butcher Easter algorithm. [New Scientist, 9 (228): 828. (30 March 1961).](https://books.google.co.uk/books?id=zfzhCoOHurwC)
  - [ ] Doomsday algorithm
  - [x] Julian date and its inverse
  - [x] Julian century
  - [x] Astronomical units and conversions
  - [x] Earth's atmosphere refraction with altitude
  - [x] Conic sections visualization
  - [x] Kepler's equation for the eccentric anomaly
  - [x] Solar system simulation from orbital elements. [Murray-Dermott. Solar System Dynamics](http://ssdbook.maths.qmul.ac.uk/)
  - [x] Arenstorf orbits to reach the Moon
* Chaos
  - [x] Bifurcation diagram of the logistic map
* Plots [plots]
  - [x] Simple PyPlot figure
  - [x] Subplots with PyPlot
  - [x] Some mplstyles
* Mixed language
  - [x] Examples of Fortran interfaces
  - [x] Examples of C interfaces
  - [x] Import Python modules

