# Used constants and struct definitions for clarity (currently not used)
module PhaseFieldConstants
    export dt, dx, rodx, vol, ϵ, γ, τ, β, σ2, τ_ξ, k, σ_c, α, repulsion, D
    export Params, Diffusion, RepField

    # Integration constants 
    const dt = 1e-3;
    const dx = 0.15;

    # Phase Field
    rodx = 9.0
    const ro = rodx*dx;
    const vol = π*rodx^2;

    const ϵ = 0.750; # 0.750
    const γ = 2.0;
    const τ = 2.0;
    const β = 0.50;
    repulsion = 0.70; 
    # Noise 
    const σ2 = 0.0150;
    const τ_ξ = 10.0;

    # External diffusion
    k = 1.0
    σ_c = 1.0
    α = 0.01
    D = 0.1

    struct Params
        dt::Float64
        dx::Float64
        rodx::Float64
        vol::Float64
        ϵ::Float64
        γ::Float64
        τ::Float64
        β::Float64
        τ_ξ::Float64
        σ2::Float64
    end

    struct Diffusion
        k::Float64
        σ_c::Float64
        α::Float64
        D::Float64
    end

    struct RepField
        A::Float64
        B::Float64
    end
end