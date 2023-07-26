# Code to run the simulation. Functions in 
include("/home/pol/cell_crawling/funcions.jl")
using Plots, .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions, .PhaseField, ProgressBars

# Mutable simulation parameters
N = 18;
nx = 110;
ny = 230
stoptime = 500.0;
repulsion = 2.3;

params = Params(dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2);

# External diffusion parameters
k = 1.4;# 1 
σ_c = 0.4; # 0.3
α = 1.4; # 0.3
D = 1; # 0.0005
diffusion = Diffusion(k,σ_c,α,D);

ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ  = initwdiff(N, nx,ny);

@time phasefield!(ϕ, ∇ϕ, ∇²ϕ, ξ, c, ∇c, ∇²c, N, params, diffusion, repulsion, stoptime) 

