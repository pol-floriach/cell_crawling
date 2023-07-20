include("/home/pol/cell_crawling/funcions.jl")
using Plots, .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions, .PhaseField, ProgressBars
#  ----- RUN SIMULATION ----- # 

# Mutable simulation parameters
N = 4;
nx = ny = 83;
stoptime = 300.0;
repulsion = 1.7;

params = Params(dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2);

# External diffusion parameters
k = 1.9;# 1 
σ_c = 0.2; # 0.3
α = 1.0; # 0.3
D = 1; # 0.0005
diffusion = Diffusion(k,σ_c,α,D);

ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ  = initwdiff(N, nx,ny);
@time phasefield!(ϕ, ∇ϕ, ∇²ϕ, ξ, c, ∇c, ∇²c, N, params, diffusion, repulsion, stoptime) 

