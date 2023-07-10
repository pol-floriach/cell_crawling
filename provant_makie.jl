using GLMakie

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

# One cell
ϕ = zeros(nx,ny);    # Phase field
∇ϕ = zeros(nx,ny);   # Gradient
∇²ϕ = zeros(nx,ny);  # Laplacian

ξ = randn(nx,ny);

Observable(ϕ)




fig = lines(xs, ys_1, color = :blue, linewidth = 4,
    axis = (title = @lift("t = $(round($time, digits = 1))"),))
scatter!(xs, ys_2, color = :red, markersize = 15)

framerate = 30
timestamps = range(0, 2, step=1/framerate)

record(fig, "time_animation.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end