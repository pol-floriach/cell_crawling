include("/home/pol/cell_crawling/funcions.jl")
using Plots, .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions, .PhaseField, ProgressBars
#  ----- RUN SIMULATION ----- # 

N = 4;
nx = ny = 83;
stoptime = 100.0;
repulsion = 1.3;

# ξ, ϕ, ∇ϕ, ∇²ϕ  = initialization(N, nx,ny);
# @time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, N, params, repulsion, stoptime)

params = Params(dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2);
k = 1.9;# 1 
σ_c = 0.2; # 0.3
α = 0.5; # 0.3
D = 1; # 0.0005

diffusion = Diffusion(k,σ_c,α,D);

ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ  = initialization_wdiffusion(N, nx,ny);
@time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, c, ∇c, ∇²c, N, params, diffusion, repulsion, stoptime) 

