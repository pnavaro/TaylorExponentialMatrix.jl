using TaylorExponentialMatrix
using Test

@testset " simplest implementation " begin

A = rand( 100, 100)

@test expm2(A) ≈ exp(A)

end

@testset " sophisticated implementation " begin

A = rand( 100, 100)

@test expm3(A) ≈ exp(A)

end
