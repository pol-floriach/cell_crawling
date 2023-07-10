# ----- FUNCTIONS ----- #

# Gradient
function gradient!(c,∇c,dx)    
    for i in 2:nx-1, j in 2:ny-1
        ∇c[i,j] = sqrt( (c[i+1,j] - c[i-1,j])^2 +  (c[i,j+1] - c[i,j-1])^2)/(2*dx)
    end
    for j in 2:ny-1
        ∇c[1, j]  = sqrt((c[2, j] - c[nx, j])^2 + (c[1, j+1] - c[1, j-1])^2) / (2*dx)
        ∇c[nx, j] = sqrt((c[1, j] - c[nx-1, j])^2 + (c[nx, j+1] - c[nx, j-1])^2) / (2*dx)
    end
    for i in 2:nx-1
        ∇c[i, 1] = sqrt((c[i+1, 1] - c[i-1, 1])^2 + (c[i, 2] - c[i, ny])^2) / (2*dx)
        ∇c[i, ny] = sqrt((c[i+1, ny] - c[i-1, ny])^2 + (c[i, 1] - c[i, ny-1])^2) / (2*dx)
    end
    # Compute ∇c for the four corner cells using periodic conditions
    ∇c[1, 1] = sqrt((c[2, 1] - c[nx, 1])^2 + (c[1, 2] - c[1, ny])^2) / (2*dx)
    ∇c[1, ny] = sqrt((c[2, ny] - c[nx, ny])^2 + (c[1, 1] - c[1, ny-1])^2) / (2*dx)
    ∇c[nx, 1] = sqrt((c[1, 1] - c[nx-1, 1])^2 + (c[nx, 2] - c[nx, ny])^2) / (2*dx)
    ∇c[nx, ny] = sqrt((c[1, ny] - c[nx-1, ny])^2 + (c[nx, 1] - c[nx, ny-1])^2) / (2*dx)
    nothing
end
# Laplacian with finite differences
function laplacian!(c,∇²c,dx)
    for i in 2:nx-1, j in 2:ny-1
        ∇²c[i,j] = (c[i+1,j] + c[i-1,j] + c[i,j+1] + c[i,j-1] - 4*c[i,j]) / dx^2
    end
    for j in 2:ny-1
        ∇²c[1, j] = (c[2, j] + c[nx, j] + c[1, j+1] + c[1, j-1] - 4*c[1, j]) / dx^2
        ∇²c[nx, j] = (c[1, j] + c[nx-1, j] + c[nx, j+1] + c[nx, j-1] - 4*c[nx, j]) / dx^2
    end
    for i in 2:nx-1
        ∇²c[i, 1] = (c[i+1, 1] + c[i-1, 1] + c[i, 2] + c[i, ny] - 4*c[i, 1]) / dx^2
        ∇²c[i, ny] = (c[i+1, ny] + c[i-1, ny] + c[i, 1] + c[i, ny-1] - 4*c[i, ny]) / dx^2
    end

    ∇²c[1, 1] = (c[2, 1] + c[nx, 1] + c[1, 2] + c[1, ny] - 4*c[1, 1]) / dx^2
    ∇²c[1, ny] = (c[2, ny] + c[nx, ny] + c[1, 1] + c[1, ny-1] - 4*c[1, ny]) / dx^2
    ∇²c[nx, 1] = (c[1, 1] + c[nx-1, 1] + c[nx, 2] + c[nx, ny] - 4*c[nx, 1]) / dx^2
    ∇²c[nx, ny] = (c[1, ny] + c[nx-1, ny] + c[nx, 1] + c[nx, ny-1] - 4*c[nx, ny]) / dx^2
    nothing
end

# Initial conditions - gives circle 
function generation_phi!(phi,io,jo)
    for i = 1:nx, j = 1:ny
        r = dx*sqrt((i-io)^2+(j-jo)^2) ;
        phi[i,j] = 0.5+0.5*tanh((ro-r)/(dx*ϵ));
    end
end
# G'(ϕ)
Gd(ϕ::Matrix{Float64}) = @. 72.0*(ϕ * (1-ϕ)*(ϕ-0.5))
Gd(ϕ::Float64) = @. 72.0*(ϕ * (1-ϕ)*(ϕ-0.5))


 

# Other functions not used now 
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