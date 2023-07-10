using Plots 
include("/home/pol/cell_crawling/codis/funcions.jl");
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
# Function for integration, vectorized
function PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, params)
    nx, ny, dt, dx, stoptime, rodx, ro, threshold, vol, α, ϵ, γ, τ, ma, mb, τ_ξ = params
    # Initialization of phase field
    generation_index!()
    generation_phi!(ϕ,30,50)
    ϕ_tot = vol
    # sqrtΔt = sqrt(dt)
    amplitude = sqrt(2*σ2*dt)/dx
    actualtime = 0.0

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
        heatmap(ϕ, title = "time = $(round((timestep*dt),digits = 1))", colormap = :Accent_4, colorbar = false);
    end every 500
end

# Function for integration for two cells
function PhaseFieldMultiple!(ϕ, ∇ϕ, ∇²ϕ, ξ, params)
    nx, ny, dt, dx, stoptime, rodx, ro, threshold, vol, α, ϵ, γ, τ, ma, mb = params
    # Initialization of phase field
    generation_index!()
    generation_phi!(ϕ[1],25,25)
    generation_phi!(ϕ[2],75,75)
    ϕ_tot = vol*ones(2)

    sqrtΔt = 10*sqrt(dt)
    actualtime = 0.0
    dϕ = [zeros(nx,ny),zeros(nx,ny)]
    rand_term = zeros(nx,ny)
    @gif for timestep in 1:Int(stoptime/dt)
        # Update values for laplacian and gradient
        @. laplacian!(ϕ,∇²ϕ,dx); 
        @. gradient!(ϕ,∇ϕ,dx);

        # Next step
        ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude

        for k in 1:2
            dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕ[k])/(ϵ^2)) - ma/τ *(ϕ_tot[k]-vol)*∇ϕ[k] + ξ*∇ϕ[k]     
            ϕ[k] = ϕ[k] + dt*dϕ[k]
            ϕ_tot[k] = sum(ϕ[k])
        end
        ϕ_all = sum(ϕ)
        # Plot
        heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 1))", colormap = :Accent_3, colorbar = false);
    end every 500
end

# ----- CONSTANTS ----- #
# Integration constants 
const nx = 110;
const ny = 110;
const dt = 0.001;
const dx = 0.15;
stoptime = 10.0

# Phase Field
const rodx = 20
const ro = rodx*dx
const threshold = 0.001
const vol = π*rodx^2

const α = 3.0
const ϵ = 0.750
const γ = 2.0
const τ = 2.0
const ma = 0.50
const mb = 1.0
# Noise 
σ2 = 0.0150;
const τ_ξ = 10.0;
const amplitude = sqrt(2*σ2*dt)/dx

params = (nx, ny, dt, dx, stoptime, rodx, ro, threshold, vol, α, ϵ, γ, τ, ma, mb, τ_ξ)

# ----- VECTOR INICIALIZATION ----- #

# One cell
ϕ = zeros(nx,ny);    # Phase field
∇ϕ = zeros(nx,ny);   # Gradient
∇²ϕ = zeros(nx,ny);  # Laplacian

ξ = randn(nx,ny);

# Two cells

ϕ   = [zeros(nx,ny),zeros(nx,ny)];
∇ϕ  = [zeros(nx,ny),zeros(nx,ny)];
∇²ϕ = [zeros(nx,ny),zeros(nx,ny)];

# Multiple cells

ϕ   = [zeros(nx,ny), for i in 1:N];
∇ϕ  = [zeros(nx,ny), for i in 1:N];
∇²ϕ = [zeros(nx,ny), for i in 1:N];


ixyp1 = zeros(Int,nx);
ixym1 = zeros(Int,nx);

#  ----- RUN SIMULATION ----- # 

# @time PhaseField!(ϕ, ∇ϕ, ∇²ϕ, ξ, params)

@time PhaseFieldTwo!(ϕ, ∇ϕ, ∇²ϕ, ξ,params)

@time PhaseFieldMultiple!(ϕ, ∇ϕ, ∇²ϕ, ξ,params)
