using DrWatson
@quickactivate

using Graphs
using LinearAlgebra, KrylovKit
using Plots

# function analyze(A)
#     _, R, _ = eigsolve(transpose(A), 2, :LR); 
#     return scatter(abs.(R[1]), R[2])
# end

# function P(A)
#     Λ, R, _ = eigsolve(A, 1, :LR, tol = 1e-20)
#     Λ, R = real(Λ[1]), real.(R[1])
#     P = spdiagm(1 ./R)*A*spdiagm(R)./Λ
#     return P
# end

function embedding_2D(P; weights = nothing)
    Λ, R, _ = eigsolve(P, 3, :LR)
    X, Y = R[2].*R[1]./sqrt.(Λ[1] - Λ[2]), R[3].*R[1]./sqrt.(Λ[1] - Λ[3])
    scatter(X./X[1], Y./Y[1], marker_z = weights == nothing ? log.(R[1]) : weights, label = false, alpha = .9)
end

function embedding_3D(P; weights = nothing)
    Λ, R, _ = eigsolve(P, 4, :LR)
    X, Y, Z = R[2]./sqrt.(Λ[1] - Λ[2]), R[3]./sqrt.(Λ[1] - Λ[3]), R[4]./sqrt.(Λ[1] - Λ[4])
    scatter(X, Y, Z, marker_z = weights == nothing ? log.(R[1]) : weights, label = false)
end