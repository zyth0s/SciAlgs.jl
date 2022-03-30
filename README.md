# SciAlgs

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://zyth0s.github.io/SciAlgs.jl/)
[![Build Status](https://travis-ci.com/zyth0s/SciAlgs.jl.svg?branch=master)](https://travis-ci.com/zyth0s/SciAlgs.jl)
[![codecov](https://codecov.io/gh/zyth0s/SciAlgs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/zyth0s/SciAlgs.jl)


SciAlgs is a compilation of fundamental scientific algorithms.  This collection
does not intended to be exhaustive but to offer clear and concise
implementations.

## [Why Julia?](https://github.com/zyth0s/SciAlgs.jl/wiki/Why-Julia%3F)

## Purpose

Some of the aims would be:
* Concise and pedagogical implementations
* Make algorithms accessible to researchers without expertise in computer programming.
* Rewrite legacy programs more clearly (e.g. no spaghetti code)
* Provide useful functions to perform basic operations
* Integrate modern development practices: tests, continuous integration, collaborative development, ...
* Enrich Julia's package ecosystem with new functionalities
* Reproduce figures/data of published papers (mainly educative ones).


Some of the topics will be covered are:
* Crystallography
* Quantum Chemistry
* Condensed Matter Physics
* Quantum Mechanics models
* Numerical methods

NOTE: Tutorials, `tutorial_<name>.jl`, can be translated with `jupytext --sync tutorial_<name>.jl`
 to Jupyter notebooks `tutorial_<name>.ipynb`, and be executed in the browser.

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
    + [x] Lebedev-Laikov sphere quadrature.
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
  - [x] Simulation of a quantum teleportation circuit [Andy Matuschak and Michael A. Nielsen, “How Quantum Teleportation Works”,San Francisco (2019)](https://quantum.country/teleportation)
* Electronic Structure [ElStruct]
  - [x] [Psi4Julia](https://github.com/zyth0s/psi4julia) shows how to interact with [Psi4](https://psicode.org).
    + [x] Hartree-Fock (HF)
    + [x] SCF with DIIS [Chem. Phys. Lett. 73, 393-398 (1980)](https://doi.org/10.1016/0009-2614(80)80396-4)
    + [x] Density Fitting/Resolution of Identity approximation
    + [x] Density Functional Theory (DFT) wit LDA, GGA, and VV10 functionals
    + [x] Moller-Plesset order 2 (MP2)
    + [x] Coupled-Electron Pair Approximation (CEPA0) and Coupled-Cluster Doubles (CCD)
    + [x] Configuration Interaction Singles (CIS)
    + [x] Orbital-Optimized Moller-Plesset 2 (OMP2)
    + [x] One-electron integrals over GTOs with Obara-Saika scheme.
  - [x] [PySCF.jl](https://github.com/zyth0s/PySCF.jl) shows how to interact with [PySCF](https://pyscf.org).
    + [x] Hartree-Fock (HF)
    + [x] Moller-Plesset order 2 (MP2)
    + [x] Coupled Cluster Singles and Doubles (CCSD) [J. Chem. Phys. 94, 4334](https://doi:10.1063/1.460620)
    + [x] SCF with DIIS [Chem. Phys. Lett. 73, 393-398 (1980)](https://doi.org/10.1016/0009-2614(80)80396-4)
    + [x] AO-based general population analysis (special cases: Mulliken and Löwdin analysis)
  - [x] Slater-Koster Tight-binding of Sr2RuO4 [PRL 116, 197003 (2016)](https://link.aps.org/doi/10.1103/PhysRevLett.116.197003)
  - [x] Chadi-Cohen Tight-binding of Si [Phys. Stat. Sol. (b) 68, 405 (1975)](https://onlinelibrary.wiley.com/doi/abs/10.1002/pssb.2220680140)
  - [x] McMurchie-Davidson molecular integrals evaluation scheme
  - [x] Simple plane-wave Density Functional Theory (DFT)
  - [x] Simple Hartree-Fock (HF) with GTOs (H2,He) [Szabo-Ostlund. Modern Quantum Chemistry.](https://store.doverpublications.com/0486691861.html)
  - [x] Simple Hartree-Fock (HF) with GTOs (H2,He) [J. Thijssen. Computational Physics](https://www.cambridge.org/core/books/computational-physics/BEE73B0139D4A9993193B57CDC62096E)
  - [x] Simple Hartree-Fock (HF) with two functions (He) [Am. J. Phys. 89 (4)](https://aapt.scitation.org/doi/10.1119/10.0002644)
  - [x] Davidson diagonalization [J. Comput. Phys. 17, 87–94](https://doi.org/10.1016/0021-9991(75)90065-0)
  - [x] Madelung sum with Pickard's algorithm [Phys. Rev. Materials 2, 013806](https://link.aps.org/doi/10.1103/PhysRevMaterials.2.013806)
  - [x] Madelung sum with Tavernier's algorithm [J. Phys. Chem. Lett.](https://pubs.acs.org/doi/10.1021/acs.jpclett.0c01684)
  - [x] Madelung sum with Ewald summation
  - [x] Integrate DOS in an energy range
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
  - [x] Magnitude of radiated energy when two black holes merge
* Chaos
  - [x] Bifurcation diagram of the logistic map
* Plots [plots]
  - [x] Simple PyPlot figure
  - [x] Subplots with PyPlot
  - [x] Some mplstyles
* Mixed language
  - [x] Examples of Fortran interfaces
  - [x] Examples of C interfaces
  - [x] Examples of Rust interfaces
  - [x] Examples of C++ interfaces with CxxWrap.jl
  - [x] Examples of calling LLVM intrinsics
  - [x] Import Python modules
* Util
  - [x] Validate Spanish DNI digits.
  - [x] Validate ISBN.

