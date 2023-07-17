using Plots 
include("/home/pol/cell_crawling/funcions.jl")
using .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions, .PhaseField

#  ----- RUN SIMULATION ----- # 

N = 4
nx = ny = 85
stoptime = 50.0
repulsion = 0.4

# ξ, ϕ, ∇ϕ, ∇²ϕ  = initialization(N, nx,ny);
# @time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, N, params, repulsion, stoptime)

params = Params(dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2)
k = 1.0
σ_c = 0.5
α = 0.1
D = 1
diffusion = Diffusion(k,σ_c,α,D)

ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ  = initialization_wdiffusion(N, nx,ny);
@time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, c, ∇c, ∇²c,N, params, diffusion, repulsion, stoptime)