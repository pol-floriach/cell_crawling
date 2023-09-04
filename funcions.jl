# Constants and struct definitions for organizing
module PhaseFieldConstants
    export dt, dx, rodx, vol, ϵ, γ, τ, β, σ2, τ_ξ, k, σ_c, α, repulsion, D
    export Params, Diffusion, AttRep

    # Integration constants 
    const dt = 0.001;
    const dx = 0.15;

    # Phase Field
    rodx = 10.0
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

    struct AttRep
        A::Float64
        B::Float64
    end
end

# Gradient and laplacian functions 
module Numerical
    export gradient!, laplacian!
    # Gradient
    function gradient!(c,∇c,dx, nx, ny)
        @fastmath @inbounds for j in 2:ny-1, i in 2:nx-1
            ∇c[i,j] = sqrt( (c[i+1,j] - c[i-1,j])^2 +  (c[i,j+1] - c[i,j-1])^2)/(2*dx)
        end
        @fastmath @inbounds for j in 2:ny-1
            ∇c[1, j]  = sqrt((c[2, j] - c[nx, j])^2 + (c[1, j+1] - c[1, j-1])^2) / (2*dx)
            ∇c[nx, j] = sqrt((c[1, j] - c[nx-1, j])^2 + (c[nx, j+1] - c[nx, j-1])^2) / (2*dx)
        end
        @fastmath @inbounds for i in 2:nx-1
            ∇c[i, 1] = sqrt((c[i+1, 1] - c[i-1, 1])^2 + (c[i, 2] - c[i, ny])^2) / (2*dx)
            ∇c[i, ny] = sqrt((c[i+1, ny] - c[i-1, ny])^2 + (c[i, 1] - c[i, ny-1])^2) / (2*dx)
        end
        # Compute ∇c for the four corner cells using periodic conditions
        @fastmath @inbounds ∇c[1, 1] = sqrt((c[2, 1] - c[nx, 1])^2 + (c[1, 2] - c[1, ny])^2) / (2*dx)
        @fastmath @inbounds ∇c[1, ny] = sqrt((c[2, ny] - c[nx, ny])^2 + (c[1, 1] - c[1, ny-1])^2) / (2*dx)
        @fastmath @inbounds ∇c[nx, 1] = sqrt((c[1, 1] - c[nx-1, 1])^2 + (c[nx, 2] - c[nx, ny])^2) / (2*dx)
        @fastmath @inbounds ∇c[nx, ny] = sqrt((c[1, ny] - c[nx-1, ny])^2 + (c[nx, 1] - c[nx, ny-1])^2) / (2*dx)
        nothing
    end
    # Laplacian with finite differences
    function laplacian!(c,∇²c,dx, nx, ny)
        @fastmath @inbounds for j in 2:ny-1, i in 2:nx-1
            ∇²c[i,j] = (c[i+1,j] + c[i-1,j] + c[i,j+1] + c[i,j-1] - 4*c[i,j]) / dx^2
        end
        @fastmath @inbounds for j in 2:ny-1
            ∇²c[1, j] = (c[2, j] + c[nx, j] + c[1, j+1] + c[1, j-1] - 4*c[1, j]) / dx^2
            ∇²c[nx, j] = (c[1, j] + c[nx-1, j] + c[nx, j+1] + c[nx, j-1] - 4*c[nx, j]) / dx^2
        end
        @fastmath @inbounds for i in 2:nx-1
            ∇²c[i, 1] = (c[i+1, 1] + c[i-1, 1] + c[i, 2] + c[i, ny] - 4*c[i, 1]) / dx^2
            ∇²c[i, ny] = (c[i+1, ny] + c[i-1, ny] + c[i, 1] + c[i, ny-1] - 4*c[i, ny]) / dx^2
        end

        @fastmath @inbounds ∇²c[1, 1] = (c[2, 1] + c[nx, 1] + c[1, 2] + c[1, ny] - 4*c[1, 1]) / dx^2
        @fastmath @inbounds ∇²c[1, ny] = (c[2, ny] + c[nx, ny] + c[1, 1] + c[1, ny-1] - 4*c[1, ny]) / dx^2
        @fastmath @inbounds ∇²c[nx, 1] = (c[1, 1] + c[nx-1, 1] + c[nx, 2] + c[nx, ny] - 4*c[nx, 1]) / dx^2
        @fastmath @inbounds ∇²c[nx, ny] = (c[1, ny] + c[nx-1, ny] + c[nx, 1] + c[nx, ny-1] - 4*c[nx, ny]) / dx^2
        nothing
    end
end

# Generation of phase field and array initialization functions
module Initialize
    export gen_phi!, gen_phi_equis!, gen_phi_equis_2groups!, init, initwdiff
    # Initial conditions - gives circle 
    function gen_phi!(phi,io,jo)
        for j = 1:ny, i = 1:nx
            r = dx*sqrt((i-io)^2+(j-jo)^2) ;
            phi[i,j] = 0.5+0.5*tanh((ro-r)/(dx*ϵ));
        end
    end
    # Initial conditions - gives N equispaced circles 
    function gen_phi_equis!(ϕ,N,rodx,dx,ϵ)
        ro = rodx * dx
        nx, ny = size(ϕ[1])
        io = jo = rodx + 10
        ro2 = 2*rodx
        k = 0
        for ky in 1:Int(sqrt(N))
                io = rodx+1
                for kx in 1:Int(sqrt(N))
                    k += 1; 
                    for j = 1:ny,  i = 1:nx
                        r = dx*sqrt((i-io)^2+(j-jo)^2) ;
                        ϕ[k][i,j] = 0.5+0.5*tanh((ro-r)/(dx*ϵ));
                    end
                io += ro2 + 1
                end
            jo += ro2 + 1
        end
    end
    #Generation of phase field, 2 groups of N cells
    function gen_phi_equis_2groups!(ϕ,N,rodx,dx,ϵ)
        ro = rodx * dx
        nx, ny = size(ϕ[1])
        io1 = jo1 = rodx + 10
        io2 = copy(io1)
        jo2 = jo1 + Int(round(ny/2 + 5))
    
        ro2 = 2*rodx
        # Group 1 
        k = 0
        for ky in 1:Int(sqrt(N/2))
                io1 = rodx+10
                io2 = rodx+10
                for kx in 1:Int(sqrt(N/2))
                    k += 1; 
                    for j = 1:Int(ny/2),  i = 1:nx
                        r = dx*sqrt((i-io1)^2+(j-jo1)^2) ;
                        ϕ[k][i,j] = 0.5+0.5*tanh((ro-r)/(dx*ϵ));
                    end
                    k+=1
                    for j = Int(ny/2+1):ny,  i = 1:nx
                        r = dx*sqrt((i-io2)^2+(j-jo2)^2) ;
                        ϕ[k][i,j] = 0.5+0.5*tanh((ro-r)/(dx*ϵ));
                    end
                io1 += ro2 + 1
                io2 += ro2 + 1
                end
            jo1 += ro2 + 1
            jo2 += ro2 + 1
        end
    end

    # Vector initialization
    function init(N, nx, ny)
        ξ = randn(nx,ny);
        if N == 1 
            # One cell
            ϕ   = zeros(nx,ny);     # Phase field
            ∇ϕ  = zeros(nx,ny);     # Gradient
            ∇²ϕ = zeros(nx,ny);     # Laplacian
        else
            # Multiple cells
            ϕ   = [zeros(nx,ny) for i in 1:N];
            ∇ϕ  = [zeros(nx,ny) for i in 1:N];
            ∇²ϕ = [zeros(nx,ny) for i in 1:N];
        end
        return ξ, ϕ, ∇ϕ, ∇²ϕ
    end
    # Vector initialization with external concentration diffusion
    function initwdiff(N, nx, ny)
        ξ = randn(nx,ny);
        c = rand(nx,ny);
        ∇c = zeros(nx,ny;)
        ∇²c = zeros(nx,ny);
        if N == 1
            # One cell
            ϕ   = zeros(nx,ny);     # Phase field
            ∇ϕ  = zeros(nx,ny);     # Gradient
            ∇²ϕ = zeros(nx,ny);     # Laplacian
        else
            # Multiple cells
            ϕ   = [zeros(nx,ny) for i in 1:N];
            ∇ϕ  = [zeros(nx,ny) for i in 1:N];
            ∇²ϕ = [zeros(nx,ny) for i in 1:N];
        end
        return ξ, c, ∇c, ∇²c, ϕ, ∇ϕ, ∇²ϕ
    end
end

# Other used functions
module OtherFunctions
    export Gd
    # G'(ϕ)
    Gd(ϕ::Matrix{Float64}) = @. 72.0*(ϕ * (1-ϕ)*(ϕ-0.5))
    Gd(ϕ::Float64) = @. 72.0*(ϕ * (1-ϕ)*(ϕ-0.5))
end


# Main functions, used for the simulation
module PhaseField 
    export phasefield!, phasefield2!
    using Plots,ProgressBars, ..PhaseFieldConstants, ..Initialize, ..Numerical, ..OtherFunctions 

    Mat = Matrix{Float64}
    ArrMat = Vector{Mat}

    # Function for integration for a single cell
    function phasefield!(ϕ::Mat, ∇ϕ::Mat, ∇²ϕ::Mat, ξ::Mat, params::Params, stoptime)
        # dt, dx, stoptime, rodx, ro,  vol, α, ϵ, γ, τ, β, repulsion, τ_ξ, σ2 = params
        (; dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2) = params

        # Initialization of phase field
        gen_phi!(ϕ,30,50)
        ϕ_tot = vol
        amplitude = sqrt(2*σ2*dt)/dx

        @gif for timestep in 1:Int(stoptime/dt)
            # Update values for laplacian and gradient
            laplacian!(ϕ,∇²ϕ,dx); 
            gradient!(ϕ,∇ϕ,dx); 

            # Next step
            dϕ =  γ/τ * (∇²ϕ + Gd(ϕ)/(ϵ^2)) - β/τ * (ϕ_tot-vol) * ∇ϕ + ξ.*∇ϕ  
            ϕ = ϕ + dt*dϕ
            ξ = ξ + dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            ϕ_tot = sum(ϕ)        
            # Plot
            heatmap(ϕ, title = "time = $(round((timestep*dt),digits = 1))", colormap = :Accent_4, colorbar = false, size = (800,800));
        end every 500
    end
    # Function for integration for multiple cells
    function phasefield!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, N, params::Params, repulsion, stoptime)
        (; dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2) = params
        nx, ny = size(ϕ[1])
        # Initialization of phase field
        # generation_phi!(ϕ[1],20,20) ; generation_phi!(ϕ[2],80,80)
        gen_phi_equis!(ϕ, N, rodx, dx, ϵ)
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
                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕ[k] + ξ*∇ϕ[k]   - repulsion*∇ϕ[k]*(∇ϕ_tot - ∇ϕ[k])  
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
    function phasefield!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, c::Mat, ∇²c::Mat, N, params::Params, diffusion::Diffusion, repulsion, stoptime)
        (; dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2) = params
        (; k, σ_c, α, D) = diffusion
        nx, ny = size(ϕ[1])
        # Initialization of phase field
        gen_phi_equis!(ϕ, N, rodx, dx, ϵ)
        ϕ_tot = vol*ones(N)
        ∇ϕ_tot = sum(∇ϕ)
        amplitude = sqrt(2*σ2*dt)/dx
        dϕ = [zeros(nx,ny) for i in 1:N]
        ∇ϕ_k = zeros(nx,ny)
        # ϕ_tot_plot = sum(ϕ)

        anim = @animate for timestep in ProgressBar(1:Int(round(stoptime/dt)))
            # Update values for laplacian and gradient

            @. gradient!(ϕ,∇ϕ,dx,nx,ny);
            @. laplacian!(ϕ,∇²ϕ,dx,nx,ny); 
            laplacian!(c,∇²c,dx,nx,ny);
            
            # Next step
            ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            c += dt*(k*∇ϕ_tot - σ_c*c + D*∇²c)
            @fastmath @inbounds for k in 1:N
                ∇ϕ_k = ∇ϕ[k]
                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕ_k + ξ*∇ϕ_k   - repulsion*∇ϕ_k*(∇ϕ_tot - ∇ϕ_k) + α/τ*c*∇ϕ_k #c*∇c  

                ϕ[k] = ϕ[k] + dt*dϕ[k]
                ϕ_tot[k] = sum(ϕ[k])
            end
            ϕ_all = sum(ϕ)
            ∇ϕ_tot = sum(∇ϕ)

            # Plot
                # Passar a 1 si >1 pel plot
                # for i in eachindex(ϕ_tot_plot)
                #     ϕ_all[i] > 1 ? ϕ_tot_plot[i] = 1 : ϕ_tot_plot[i] = ϕ_all[i]
                # end
            heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 0))", colormap = :Accent_4, colorbar = false, size = (800,800))
            # heatmap(c, title = "time = $(round((timestep*dt), digits = 0))", colorbar = false, size = (800,800))
        end every 2500
        gif(anim, "pf_2groups_$(N)_rep_$(repulsion)_k_$(k)_attr_$(α)_deg_$(σ_c)_gen_$(α)_diff_$(D).gif", fps = 15);
    end
# Function for integration of multiple cells with quadratic attraction/repulsion potential
    function phasefield2!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, N, params::Params, attr_rep::AttRep, stoptime)
        (; dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2) = params
        nx, ny = size(ϕ[1])
        (; A, B) = attr_rep
        # Initialization of phase field
        # generation_phi!(ϕ[1],20,20) ; generation_phi!(ϕ[2],80,80)
        gen_phi_equis!(ϕ, N, rodx, dx, ϵ)
        ϕ_tot = vol*ones(N)
        ∇ϕ_tot = sum(∇ϕ)
        amplitude = sqrt(2*σ2*dt)/dx
        dϕ = [zeros(nx,ny) for i in 1:N]

        anim = @animate for timestep in ProgressBar(1:Int(round(stoptime/dt)))
            # Update values for laplacian and gradient
            @. laplacian!(ϕ,∇²ϕ,dx, nx, ny); 
            @. gradient!(ϕ,∇ϕ,dx, nx, ny);

            # Next step
            ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            
            for k in 1:N
                ∇ϕ_k = ∇ϕ[k]
                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕ_k + ξ*∇ϕ_k  + A*∇ϕ_k*(∇ϕ_tot-∇ϕ_k) - B*(∇ϕ_k*(∇ϕ_tot-∇ϕ_k))^2
                ϕ[k] = ϕ[k] + dt*dϕ[k]
                ϕ_tot[k] = sum(ϕ[k])
            end
            ϕ_all = sum(ϕ)
            ∇ϕ_tot = sum(∇ϕ)
            # Plot
            heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 0))", colormap = :Accent_4, colorbar = false, size = (800,800))
        end every 500
        gif(anim, "pf_$(N)_cells_w_repulsion_$(repulsion).gif", fps = 15);
    end


end



 # Other functions not used now 
#= 
    # Function for integration, with a loop
    function PhaseFieldLoop!(ϕ, ∇ϕ, ∇²ϕ, params)
        nx, ny, dt, dx, stoptime, rodx, ro, threshold, vol, α, ϵ, γ, τ, ma, mb = params
        generation_index!()
        generation_phi!(ϕ,30,50)
        ϕ_tot = vol
        sqrtΔt = sqrt(dt)
        # Initialization of phase field
        actualtime = 0.0
        plt = plot()
        for timestep in 1:Int(stoptime/dt)
            actualtime = round(actualtime+dt,digits=1)
            # Update values for laplacian and gradient
            laplacian!(ϕ,∇²ϕ,dx); 
            gradient!(ϕ,∇ϕ,dx) 

            for i in 1:nx, j in 1:ny
                dϕ = γ/τ * (∇²ϕ[i,j] + Gd(ϕ[i,j])/(ϵ^2)) - ma/τ * (ϕ_tot-vol) * ∇ϕ[i,j] + randn()*sqrtΔt*∇ϕ[i,j]      
                ϕ[i,j] = ϕ[i,j] + dt*dϕ
                # if ϕ[i,j] < 0.0 
                #     ϕ[i,j] = 0.0
                # elseif ϕ[i,j] > 1.0 
                #     ϕ[i,j] = 1.0
                # end
            end
            ϕ_tot = sum(ϕ)
            #if timestep % 500 == 0
            #   heatmap(ϕ, title = "time = $actualtime", colormap = :Accent_4);
            #   savefig(p1,"fig_$time.png")
            #end
        end #every 500
    end
    # For boundary conditions
    function generation_index!()
        for i in 1:nx-1
            ixyp1[i] = i+1
        end
        ixyp1[nx] = 1
        for i in 2:nx
            ixym1[i] = i-1
        end
        ixym1[1] = nx-1
        nothing
    end
=#