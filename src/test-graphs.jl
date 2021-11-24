include("../src/retrospective-process.jl")

using Graphs, SparseArrays 
using DifferentialEquations

lattice = Graphs.grid([50, 50], periodic = true)

ϕ = rand(length(vertices(lattice)))
m = spdiagm(ϕ) - laplacian_matrix(lattice)

function f!(du, u, p, t)
    mul!(du, p[:m], u)
end



p = @dict m 

u0 = zeros(length(vertices(lattice)))
u0[rand(vertices(lattice))] = 1
pb = ODEProblem(f!, u0, (0., 500.), p)
sol = DifferentialEquations.solve(pb)

normalize(v) = v/sum(v)
@gif for i  = 1:length(sol.t)
    embedding_2D(m; weights = normalize(sol.u[i]))
end

embedding_2D(m; weights = normalize(ϕ))