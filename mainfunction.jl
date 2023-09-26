# Main function definitions for different modifications of the phase field model
module PhaseField 
    export phasefield!, phasefield2!, phasefield3!, phasefield_rigged!
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
# Function for integration of multiple cells with quadratic attraction/repulsion force
    function phasefield!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, N, params::Params, rep::RepField, stoptime)
        (; dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2) = params
        (; A, B) = rep
        nx, ny = size(ϕ[1])
        # Initialization of phase field
        gen_phi_equis!(ϕ, N, rodx, dx, ϵ)
        ϕ_tot = vol*ones(N)
        ∇ϕ_tot = sum(∇ϕ)
        ϕ_all = sum(ϕ)
        amplitude = sqrt(2*σ2*dt)/dx* 0.4
        dϕ = [zeros(nx,ny) for i in 1:N]

        anim = @animate for timestep in 1:Int(round(stoptime/dt))
        # anim = @animate for timestep in ProgressBar(1:Int(round(stoptime/dt)))

            # Update values for laplacian and gradient
            @. laplacian!(ϕ,∇²ϕ,dx, nx, ny); 
            @. gradient!(ϕ,∇ϕ,dx, nx, ny);

            # Next step
            ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            
            for k in 1:N
                ∇ϕₖ = ∇ϕ[k]
                ϕₖ = ϕ[k]
                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕₖ)/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕₖ + ξ*∇ϕₖ + A*(-ϕₖ*(ϕ_all-ϕₖ) + B*ϕₖ^2*(ϕ_all-ϕₖ)^2)
                ϕ[k] = ϕ[k] + dt*dϕ[k]
                ϕ_tot[k] = sum(ϕ[k])
            end
            ϕ_all = sum(ϕ)
            ∇ϕ_tot = sum(∇ϕ)
            # Plot
            heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 0))", colormap = :Accent_4, colorbar = false, size = (800,800))
        end every 1000
        gif(anim, "/home/pol/figs2/pf_$(N)_A_$(A)_B_$(B).gif", fps = 15);

    end

    # Rigged force, quadratic only if too close (repulsion), constant > 0 if not to far (attraction), 0 if too far
    function phasefield_rigged!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, N, params::Params, rep::RepField, stoptime)
        (; dt, dx, rodx, vol, ϵ, γ, τ, β, τ_ξ, σ2) = params
        (; A, B) = rep
        nx, ny = size(ϕ[1])
        # Initialization of phase field
        gen_phi_equis!(ϕ, N, rodx, dx, ϵ)
        ϕ_tot = vol*ones(N)
        ∇ϕ_tot = sum(∇ϕ)
        ϕ_all = sum(ϕ)
        amplitude = sqrt(2*σ2*dt)/dx
        dϕ = [zeros(nx,ny) for i in 1:N]

        F_l = zeros(nx,ny)
        cond = 0.0
        anim = @animate for timestep in 1:Int(round(stoptime/dt))
        # anim = @animate for timestep in ProgressBar(1:Int(round(stoptime/dt)))

            # Update values for laplacian and gradient
            @. laplacian!(ϕ,∇²ϕ,dx, nx, ny); 
            @. gradient!(ϕ,∇ϕ,dx, nx, ny);

            # Next step
            ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            
            F_tot = zeros(nx,ny)
            for k in 1:N
                ∇ϕₖ = ∇ϕ[k]
                ϕₖ = ϕ[k]
                # Attraction / repulsion force bc of other cells
                    for l in 1:N
                	    if l != k
                            for m in 1:nx, n in 1:ny
                                ϕₗ = ϕ[l]
                                cond = ϕₖ[m,n]*ϕₗ[m,n]
                                if cond < 0.05
                                    F_l[m,n] = 0.0
                                elseif cond >= 0.25
                                    F_l[m,n] += -A*cond
                                else 
                                    F_l[m,n] = B
                                end
                            end
                        end
                        F_tot += F_l
                    end

                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕₖ)/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕₖ + ξ*∇ϕₖ + F_tot
                ϕ[k] += dt*dϕ[k]
                ϕ_tot[k] = sum(ϕ[k])
            end
            ϕ_all = sum(ϕ)
            ∇ϕ_tot = sum(∇ϕ)

            # Stop if NaN
            any(isnan,ϕ_all) ? (println("Algun valor és NaN"); break) : nothing

            # Plot
            heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 0))", colormap = :Accent_4, colorbar = false, size = (800,800))
        end every 1000
        gif(anim, "/home/pol/figs_rigged/pf_$(N)_A_$(A)_B_$(B).gif", fps = 15);
    end
    # mateixa pero sense els structs per comparar performance . + @inbounds i @fastmath
    function phasefield_rigged!(ϕ::ArrMat, ∇ϕ::ArrMat, ∇²ϕ::ArrMat, ξ::Mat, N, A, B, stoptime)
        nx, ny = size(ϕ[1])
        # Initialization of phase field
        gen_phi_equis!(ϕ, N, rodx, dx, ϵ)
        ϕ_tot = vol*ones(N)
        ∇ϕ_tot = sum(∇ϕ)
        ϕ_all = sum(ϕ)
        amplitude = sqrt(2*σ2*dt)/dx
        dϕ = [zeros(nx,ny) for i in 1:N]

        F_l = zeros(nx,ny)
        cond = 0.0
        anim = @animate for timestep in 1:Int(round(stoptime/dt))
        # anim = @animate for timestep in ProgressBar(1:Int(round(stoptime/dt)))

            # Update values for laplacian and gradient
            @. laplacian!(ϕ,∇²ϕ,dx, nx, ny); 
            @. gradient!(ϕ,∇ϕ,dx, nx, ny);

            # Next step
            ξ += dt*(-ξ/τ_ξ) + randn(nx,ny)*amplitude
            
            F_tot = zeros(nx,ny)
            @inbounds for k in 1:N
                ∇ϕₖ = ∇ϕ[k]
                ϕₖ = ϕ[k]
                # Attraction / repulsion force bc of other cells
                    for l in 1:N
                	    if l != k
                            for m in 1:nx, n in 1:ny
                                ϕₗ = ϕ[l]
                                cond = ϕₖ[m,n]*ϕₗ[m,n]
                                if cond < 0.05
                                    F_l[m,n] = 0.0
                                elseif cond >= 0.25
                                    F_l[m,n] += -A*cond
                                else 
                                    F_l[m,n] = B
                                end
                            end
                        end
                        F_tot += F_l
                    end

                dϕ[k] = @.  γ/τ * (∇²ϕ[k] + Gd(ϕₖ)/(ϵ^2)) - β/τ *(ϕ_tot[k]-vol)*∇ϕₖ + ξ*∇ϕₖ + F_tot
                ϕ[k] += dt*dϕ[k]
                ϕ_tot[k] = sum(ϕ[k])
            end
            ϕ_all = sum(ϕ)
            ∇ϕ_tot = sum(∇ϕ)

            # Stop if NaN
            any(isnan,ϕ_all) ? (println("Algun valor és NaN"); break) : nothing

            # Plot
            heatmap(ϕ_all, title = "time = $(round((timestep*dt),digits = 0))", colormap = :Accent_4, colorbar = false, size = (800,800))
        end every 1000
        gif(anim, "/home/pol/figs_rigged_long/pf_$(N)_A_$(A)_B_$(B).gif", fps = 15);
    end
end