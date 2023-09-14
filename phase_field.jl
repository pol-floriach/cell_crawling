#!/home/pol/julia-1.9.2/bin/julia
ENV["GKSwstype"] = "nul"
# Code to run the simulation. Functions in 
include("/home/pol/cell_crawling/funcions.jl")
using Plots, .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions, .PhaseField, ProgressBars

# Mutable simulation parameters
N = 4;
nx = ny = 55;
stoptime = 400.0;
# repulsion = 2.3;

params = Params(dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2);

# External diffusion parameters
k = 1.4;# 1 
σ_c = 0.4; # 0.3
α = 1.4; # 0.3
D = 1; # 0.0005
diffusion = Diffusion(k,σ_c,α,D);

# ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ  = initwdiff(N, nx,ny);
ξ, ϕ, ∇ϕ, ∇²ϕ  = init(N, nx,ny);

rang_integracio = parse(Int,ARGS[1])

if rang_integracio == 1 
    A_range = 0.5:0.5:2
elseif rang_integracio == 2
    A_range = 2.5:0.5:4
elseif rang_integracio == 3
    A_range = 4.5:0.5:6
elseif rang_integracio == 4
    A_range = 6.5:0.5:8
elseif rang_integracio == 5
    A_range = 8.5:0.5:10
end

println(A_range)
# # for A in 0.5:0.5:10, B in 0.5:0.5:5 
for A in A_range, B in 0.5:0.5:5
    rep = RepField(A, B)
    @time phasefield2!(ϕ, ∇ϕ, ∇²ϕ, ξ, N, params, rep, stoptime);
end


