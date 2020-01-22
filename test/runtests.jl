using TaylorExponentialMatrix
using Test

@testset "TaylorExponentialMatrix.jl" begin

A = rand( 5, 5)

@show expm2(A) .- exp(A)

@test all( expm2(A) .≈ exp(A))

end
