using Plots

# ----- FUNCTIONS ----- #
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
# Gradient
function gradient!(c,∇c)
    for i in 1:nx, j in 1:ny
        ∇c[i,j] = sqrt( (c[ixyp1[i],j] - c[ixym1[i],j])^2 +  (c[i,ixyp1[j]] - c[i,ixym1[j]])^2)/(2*dx)
    end
end
# Laplacian with finite differences
function laplacian!(c,dx2c)
    for i in 1:nx, j in 1:ny
        dx2c[i,j] = (c[ixyp1[i],j] + c[ixym1[i],j] + c[i,ixyp1[j]] + c[i,ixym1[j]] - 4*c[i,j]) / dx^2
    end
end
# Initial conditions - gives circle 
function generation_phi!(phi,io,jo)
    for i = 1:nx, j = 1:ny
        r = dx*sqrt((i-io)^2+(j-jo)^2) ;
        phi[i,j] = 0.5+0.5*tanh((ro-r)/(dx*ϵ));
    end
end
# G'(ϕ)
Gd(ϕ) = 72.0*(ϕ * (1-ϕ)*(ϕ-0.5))
# Function for integration
function phase_field(params)
    nx, ny, dt, dx, stoptime, rodx, ro, threshold, vol, α, ϵ, γ, τ, ma, mb = params
    generation_index!()
    generation_phi!(ϕ,30,50)
    ϕ_tot = vol
    sqrtΔt = sqrt(dt)
    # Initialization of phase field
    actualtime = 0.0
    plt = plot()
    @gif for timestep in 1:Int(stoptime/dt)
        actualtime = round(actualtime+dt,digits=1)
        # Update values for laplacian and gradient
        laplacian!(ϕ,∇²ϕ); 
        gradient!(ϕ,∇ϕ) 

        for i in 1:nx, j in 1:ny
            # dϕ = γ/τ * (∇²ϕ + Gd(ϕ[i,j])/(ϵ^2))
            #     - mb/τ*(ϕ[i,j]*ϕ_t[i,j]*∇ϕ[i,j]) 
            #     - ma/τ * (ϕ_tot-vol) * ∇ϕ[i,j] 
            #     + α/τ * bct[i,j]*ϕ[i,j] * ∇ϕ[i,j]

            dϕ = γ/τ * (∇²ϕ[i,j] + Gd(ϕ[i,j])/(ϵ^2)) - ma/τ * (ϕ_tot-vol) * ∇ϕ[i,j] + randn()*sqrtΔt*∇ϕ[i,j]
                  
            ϕ[i,j] = ϕ[i,j] + dt*dϕ

            if ϕ[i,j] < 0.0 
                ϕ[i,j] = 0.0
            elseif ϕ[i,j] > 1.0 
                ϕ[i,j] = 1.0
            end
        end

        ϕ_tot = sum(ϕ)

        #if timestep % 500 == 0
        heatmap(ϕ, title = "time = $actualtime", colormap = :Accent_4);
        #    savefig(p1,"fig_$time.png")
        #end
    end every 500
end
# ----- CONSTANTS ----- #

# Integration constants 
const nx = 100;
const ny = 100;
const dt = 0.001;
const dx = 0.15;
stoptime = 50.0

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

params = (nx, ny, dt, dx, stoptime, rodx, ro, threshold, vol, α, ϵ, γ, τ, ma, mb)

# ----- VECTOR INICIALIZATION ----- #

ϕ = zeros(nx,ny);    # Phase field
∇ϕ = zeros(nx,ny);   # Gradient
∇²ϕ = zeros(nx,ny);  # Laplacian



ixyp1 = zeros(Int,nx);
ixym1 = zeros(Int,nx);

#  ----- RUN SIMULATION ----- # 

phase_field(params)