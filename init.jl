# Functions to initialize arrays and generate initial conditions
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