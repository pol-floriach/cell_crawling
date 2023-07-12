using Plots 
# include("/home/pol/cell_crawling/codis/funcions.jl");
include("/home/pol/phase_field/funcions.jl")
using .PhaseFieldConstants, .Numerical, .Initialize, .OtherFunctions
# ENV["GKSwstype"]="nul"

begin
    Mat = Matrix{Float64}
    ArrMat = Vector{Mat}
    # Function for integration, vectorized
    function PhaseField!(ϕ::Mat, ∇ϕ::Mat, ∇²ϕ::Mat, ξ::Mat, params)
        dt, dx, stoptime, rodx, ro,  vol, α, ϵ, γ, τ, β, repulsion, τ_ξ, σ2 = params
        # Initialization of phase field
        generation_phi!(ϕ,30,50)
        ϕ_tot = vol
        amplitude = sqrt(2*σ2*dt)/dx

        @gif for timestep in 1:Int(stoptime/dt)
            # Update values for laplacian and gradient
            laplacian!(ϕ,∇²ϕ,dx); 
            gradient!(ϕ,∇ϕ,dx) 

            # Next step
            dϕ =  γ/τ * (∇²ϕ + Gd(ϕ)/(ϵ^2)) - ma/τ * (ϕ_tot-vol) * ∇ϕ + ξ.*∇ϕ  
            ϕ = ϕ + dt*dϕ
            ξ = ξ + dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            ϕ_tot = sum(ϕ)        
            # Plot
            heatmap(ϕ, title = "time = $(round((timestep*dt),digits = 1))", colormap = :Accent_4, colorbar = false, size = (800,800));
        end every 500
    end
    # Function for integration for multiple cells
    function PhaseField!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, params, repulsion, stoptime)
        dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2 = params
        nx, ny = size(ϕ[1])
        # Initialization of phase field
        # generation_phi!(ϕ[1],20,20) ; generation_phi!(ϕ[2],80,80)
        generation_phi_equis!(ϕ, N, rodx, dx, ϵ)

        ϕ_tot = vol*ones(N)
        ∇ϕ_tot = sum(∇ϕ)
        amplitude = sqrt(2*σ2*dt)/dx
        dϕ = [zeros(nx,ny) for i in 1:N]

        anim = @animate for timestep in 1:Int(stoptime/dt)
            # Update values for laplacian and gradient
            @. laplacian!(ϕ,∇²ϕ,dx, nx, ny); 
            @. gradient!(ϕ,∇ϕ,dx, nx, ny);

            # Next step
            ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            
            for k in 1:N
                # ∇ϕ_iter = @view ∇ϕ[k]
                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕ[k] + ξ*∇ϕ[k]   - repulsion*∇ϕ[k]*(∇ϕ_tot - ∇ϕ[k])  
                # dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - ma/τ *(ϕ_tot[k]-vol)*∇ϕ_iter + ξ*∇ϕ_iter   - repulsion*∇ϕ_iter*(∇ϕ_tot - ∇ϕ_iter)  

                ϕ[k] = ϕ[k] + dt*dϕ[k]
                ϕ_tot[k] = sum(ϕ[k])
            end
            ϕ_all = sum(ϕ)
            ∇ϕ_tot = sum(∇ϕ)
            # Plot
            heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 0))", colormap = :Accent_3, colorbar = false, size = (800,800))
        end every 500
        gif(anim, "pf_$(N)_cells_w_repulsion_$(repulsion).gif", fps = 15);
    end
    # Function for integration of multiple cells, with external concentration diffusion
    function PhaseField!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, c::Mat,params, repulsion, stoptime)
        dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2 = params
        nx, ny = size(ϕ[1])
        # Initialization of phase field
        # generation_phi!(ϕ[1],20,20) ; generation_phi!(ϕ[2],80,80)
        generation_phi_equis!(ϕ, N, rodx, dx, ϵ)
        ϕ_tot = vol*ones(N)
        ∇ϕ_tot = sum(∇ϕ)
        amplitude = sqrt(2*σ2*dt)/dx
        dϕ = [zeros(nx,ny) for i in 1:N]

        anim = @animate for timestep in 1:Int(stoptime/dt)
            # Update values for laplacian and gradient

            @. gradient!(ϕ,∇ϕ,dx);
            @. laplacian!(ϕ,∇²ϕ,dx); 
            gradient!(c,∇c,dx);
            laplacian!(c,∇²c,dx);
            
            # Next step
            ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            dċ = k*c*∇ϕ_tot - σ*c + ∇²c
            c += dt*()
            for k in 1:N
                # ∇ϕ_iter = @view ∇ϕ[k]
                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕ[k] + ξ*∇ϕ[k]   - repulsion*∇ϕ[k]*(∇ϕ_tot - ∇ϕ[k])  + α*c*∇c[k]
                # dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - ma/τ *(ϕ_tot[k]-vol)*∇ϕ_iter + ξ*∇ϕ_iter   - repulsion*∇ϕ_iter*(∇ϕ_tot - ∇ϕ_iter)  

                ϕ[k] = ϕ[k] + dt*dϕ[k]
                ϕ_tot[k] = sum(ϕ[k])
            end
            ϕ_all = sum(ϕ)
            ∇ϕ_tot = sum(∇ϕ)

            # Plot
            heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 0))", colormap = :Accent_3, colorbar = false, size = (800,800))
        end every 500
        gif(anim, "pf_$(N)_cells_w_repulsion_$(repulsion).gif", fps = 15);
    end
end

#  ----- RUN SIMULATION ----- # 

N = 4
nx = ny = 100
ξ, ϕ, ∇ϕ, ∇²ϕ  = initialization(N, nx,ny);
stoptime = 10.0
repulsion = 1.0

@time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, params, repulsion, stoptime)


