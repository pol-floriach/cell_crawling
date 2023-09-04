# Code to run the simulation. Functions in 
include("/home/pol/cell_crawling/funcions.jl")
using Plots, .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions, .PhaseField, ProgressBars

# Mutable simulation parameters
N = 4;
nx = ny = 55
stoptime = 400.0;
repulsion = 2.3;

params = Params(dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2);

# External diffusion parameters
k = 1.4;# 1 
σ_c = 0.4; # 0.3
α = 1.4; # 0.3
D = 1; # 0.0005
diffusion = Diffusion(k,σ_c,α,D);

ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ  = initwdiff(N, nx,ny);

A = 1.5;
B = 2.3
attr_rep = AttRep(A,B)
@time phasefield2!(ϕ, ∇ϕ, ∇²ϕ, ξ, N, params, attr_rep, stoptime) 

