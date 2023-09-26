# Functions to compute the gradient and laplacian with pbc
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

# Function to compute derivative of G function from phase field equation
module OtherFunctions
    export Gd
    # G'(ϕ)
    Gd(ϕ::Matrix{Float64}) = @. 72.0*(ϕ * (1-ϕ)*(ϕ-0.5))
    Gd(ϕ::Float64) = @. 72.0*(ϕ * (1-ϕ)*(ϕ-0.5))
end