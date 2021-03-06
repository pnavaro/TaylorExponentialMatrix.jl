# TaylorExponentialMatrix.jl

This is a julia translation of the matlab code available [here](http://www.gicas.uji.es/Research/MatrixExp.html).

**Computing the Matrix Exponential with an Optimized Taylor Polynomial Approximation**

**Philipp Bader** (Departament de Matemàtiques, Universitat Jaume I, Castellón, Spain), 
**Sergio Blanes** (Instituto de Matemática Multidisciplinar, Universitat Politècnica de València, Spain) 
and **Fernando Casas** (IMAC and Departament de Matemàtiques, Universitat Jaume I, Castellón, Spain)


```julia
julia> using Pkg

julia> pkg" add https://github.com/pnavaro/TaylorExponentialMatrix.jl"

julia> A = rand(5,5)
5×5 Array{Float64,2}:
 0.0224285  0.160116   0.504822  0.370332   0.203693
 0.861772   0.156394   0.178399  0.645844   0.229411
 0.0630692  0.584537   0.358806  0.763173   0.410573
 0.320181   0.391341   0.78607   0.619399   0.055634
 0.450914   0.0945151  0.277274  0.0576302  0.560325

julia> exp(A) # version from LinearAlgebra
5×5 Array{Float64,2}:
 1.45688   0.636229  1.1295    1.13607   0.591914
 1.41956   1.77159   1.2015    1.70089   0.734244
 0.918259  1.29852   2.42545   2.00871   1.02066
 1.04361   1.20838   1.88584   2.99255   0.655514
 0.848529  0.454553  0.838875  0.638862  2.02498

julia> using TaylorExponentialMatrix

julia> expm2(A) # Version using Taylor polynomial aproximation (simple algorithm)
5×5 Array{Float64,2}:
 1.45688   0.636229  1.1295    1.13607   0.591914
 1.41956   1.77159   1.2015    1.70089   0.734244
 0.918259  1.29852   2.42545   2.00871   1.02066
 1.04361   1.20838   1.88584   2.99255   0.655514
 0.848529  0.454553  0.838875  0.638862  2.02498

julia> expm3(A) # Version using Taylor polynomial aproximation (sophisticated algorithm)
5×5 Array{Float64,2}:
 1.45688   0.636229  1.1295    1.13607   0.591914
 1.41956   1.77159   1.2015    1.70089   0.734244
 0.918259  1.29852   2.42545   2.00871   1.02066
 1.04361   1.20838   1.88584   2.99255   0.655514
 0.848529  0.454553  0.838875  0.638862  2.02498

```

# See also

- [ExponentialUtilities.jl](https://github.com/JuliaDiffEq/ExponentialUtilities.jl): Utility functions used by the exponential integrators in [OrdinaryDiffEq.jl](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl)
- [Expokit.jl](https://github.com/acroy/Expokit.jl): Julia implementation of EXPOKIT routines
- [ExpMV.jl](https://github.com/matteoacrossi/ExpmV.jl): Julia package to compute the result of `expm(t*A)*v` when A is a sparse matrix, without computing `expm(t*A)`.
