using BenchmarkTools, Plots, LoopVectorization

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

function generation_phi!(phi,io,jo)
    for i = 1:nx, j = 1:ny
        r = dx*sqrt((i-io)^2+(j-jo)^2) ;
        phi[i,j] = 0.5+0.5*tanh((ro-r)/(dx*ϵ));
    end
end
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
function gradient_1loop!(c,∇c,dx)
    for i in 1:nx, j in 1:ny
        ∇c[i,j] = sqrt( (c[ixyp1[i],j] - c[ixym1[i],j])^2 +  (c[i,ixyp1[j]] - c[i,ixym1[j]])^2)/(2*dx)
    end
    nothing
end
function gradient_1loop_turbo!(c,∇c,dx)
    @turbo for i in 1:nx, j in 1:ny
        ∇c[i,j] = sqrt( (c[ixyp1[i],j] - c[ixym1[i],j])^2 +  (c[i,ixyp1[j]] - c[i,ixym1[j]])^2)/(2*dx)
    end
    nothing
end
function gradient_standard!(c,∇c,dx)    
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
function gradient_turbo!(c,∇c,dx)    
    @turbo for i in 2:nx-1, j in 2:ny-1
        ∇c[i,j] = sqrt( (c[i+1,j] - c[i-1,j])^2 +  (c[i,j+1] - c[i,j-1])^2)/(2*dx)
    end
    @turbo for j in 2:ny-1
        ∇c[1, j]  = sqrt((c[2, j] - c[nx, j])^2 + (c[1, j+1] - c[1, j-1])^2) / (2*dx)
        ∇c[nx, j] = sqrt((c[1, j] - c[nx-1, j])^2 + (c[nx, j+1] - c[nx, j-1])^2) / (2*dx)
    end
    @turbo for i in 2:nx-1
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
function laplacian1!(c,dx2c,dx)
    for i in 1:nx, j in 1:ny
        dx2c[i,j] = (c[ixyp1[i],j] + c[ixym1[i],j] + c[i,ixyp1[j]] + c[i,ixym1[j]] - 4*c[i,j]) / dx^2
    end
    nothing
end


# ----- VECTOR INICIALIZATION ----- #

# One cell
ϕ = zeros(nx,ny);    # Phase field
∇ϕ = zeros(nx,ny);   # Gradient
∇²ϕ = zeros(nx,ny);  # Laplacian

ixyp1 = zeros(Int,nx);
ixym1 = zeros(Int,ny)

generation_phi!(ϕ,30,50)
generation_index!()



# @benchmark laplacian1!(ϕ,∇²ϕ,dx)

@benchmark gradient_1loop!(ϕ,∇ϕ,dx)
∇ϕ = zeros(nx,ny);   
@benchmark gradient_1loop_turbo!(ϕ,∇ϕ,dx)
∇ϕ = zeros(nx,ny);   
@benchmark gradient_standard!(ϕ,∇ϕ,dx)
∇ϕ = zeros(nx,ny);   
@benchmark gradient_turbo!(ϕ,∇ϕ,dx)
