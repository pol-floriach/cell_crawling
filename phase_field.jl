include("/home/pol/cell_crawling/funcions.jl")
using Plots, .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions, .PhaseField

#  ----- RUN SIMULATION ----- # 

N = 1
nx = ny = 83
stoptime = 200.0
repulsion = 0.5

# ξ, ϕ, ∇ϕ, ∇²ϕ  = initialization(N, nx,ny);
# @time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, N, params, repulsion, stoptime)

params = Params(dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2)
k = 0.0 # 1 
σ_c = 0.0 # 0.3
α = 5 # 0.3
D = 0.0005
diffusion = Diffusion(k,σ_c,α,D)

ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ  = initialization_wdiffusion(N, nx,ny);
@time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, c, ∇c, ∇²c,N, params, diffusion, repulsion, stoptime) 