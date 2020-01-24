using TaylorExponentialMatrix
using Test

@testset " simplest implementation " begin

A = rand( 5, 5)

@test all( expm2(A) .≈ exp(A))

end

@testset " sophisticated implementation " begin

A = rand( 5, 5)

@test all( expm3(A) .≈ exp(A))

end
