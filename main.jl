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
for i in 8:16
    if rang_integracio == i 
        global A_range = i+0.25:0.25:Int(i+1)
    end
end

println(A_range)

for A in A_range, B in 0.25:0.25:1.25
    rep = RepField(A, B)
    # @time phasefield3!(ϕ, ∇ϕ, ∇²ϕ, ξ, N, params, rep, stoptime)
    @time phasefield_rigged!(ϕ, ∇ϕ, ∇²ϕ, ξ, N, A, B, stoptime)
end
