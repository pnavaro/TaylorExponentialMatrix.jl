using ExponentialTaylor
using Test

@testset "ExponentialTaylor.jl" begin

A = rand( 5, 5)

@show expm2(A) .- exp(A)

@test all( expm2(A) .≈ exp(A))

end
