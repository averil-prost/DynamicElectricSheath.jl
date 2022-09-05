using DynamicElectricSheath
using ProgressMeter
using Plots 

function run_simulation(Nx, Nv)

    # ioparams = IOparameters()

    physics = Physics()

    T = physics.T
    xmin, xmax = physics.xmin, physics.xmax
    vmin, vmax = physics.vmin, physics.vmax

    numerics = Discretization(physics, Nx=Nx, Nv=Nv)

    Nt = numerics.Nt
    dt = numerics.dt
    CFL_x = numerics.CFL_x
    CFL_v = numerics.CFL_v
    Nx = numerics.Nx
    dx = numerics.dx
    Nv = numerics.Nv
    dv = numerics.dv

    println("Parameters : ")
    println("T = $T, Nt = $Nt, dt = $dt (CFL_x = $CFL_x, CFL_v = $CFL_v)")
    println("[xmin,xmax] = [$xmin,$xmax], Nx = $Nx, dx = $dx")
    println("[vmin,vmax] = [$vmin,$vmax], Nv = $Nv, dv = $dv")

    xx = collect(LinRange(xmin, xmax, Nx + 1)) # vector (vertical) of spatial grid points 
    vv = collect(LinRange(vmin, vmax, Nv + 1)) # vector (vertical) of speeds 
    vv_plus = vv .* (vv .> 0.0)
    vv_minus = vv .* (vv .< 0.0) # signed positive and negative parts

    # Initialisation
    EE = E0.(xx)        # electric field
    EExact = copy(EE)   # should be stationary
    f = f_0.(xx, vv')   # particle density
    f_old = 1.0 * f     # to optimize parallel version
    source = 0.0 * f    # no source for the advection
    ρ = zeros(Nx + 1)   # charge density
    ρexact = -1.0 * (-1.0 .< xx .< 1.0)

    # Boundary conditions (all 0, never updated afterwards)
    f[begin,:] .= 0.0
    f[end,  :] .= 0.0 
    f[:, begin] .= 0.0
    f[:, end] .= 0.0 

    @showprogress 1 for n = 1:Nt # loop over time
        compute_charge!(ρ, -f, dv)
        EE_minus, EE_plus = compute_e!(EE, ρ, 1.0, Nx, dx) # λ=1
        # advection!(f, vv_plus, vv_minus, EE_plus, EE_minus, source, dx, dv, dt)
        advection!(f, f_old, vv_plus, vv_minus, EE_plus, EE_minus, source, dx, dv, dt)
    end 

    
    # Error on ρ
    err = abs.(ρ - ρexact)
    errL∞ = maximum(err) / maximum(abs.(ρexact))
    errL¹ = sum(err) / sum(abs.(ρexact))
    println("Erreur ρ - Nx : $Nx, Erreur : L∞ $errL∞, L¹ $errL¹")
    
    # # Error on E
    # err = abs.(EE - EExact)
    # errL∞ = maximum(err) / maximum(abs.(EExact))
    # errL¹ = sum(err) / sum(abs.(EExact))
    # println("Erreur E - Nx : $Nx, Erreur : L∞ $errL∞, L¹ $errL¹")

    return errL∞, errL¹ 
end

#########################
# Main 
#########################

# Nxs = [100, 200, 400, 800]
Nxs = [2400]
Nvs = [2048]
errors = zeros(length(Nxs),2)

@time for (iNx, Nx) in enumerate(Nxs)
    errors[iNx,:] .= run_simulation(Nx, 2*Nx)
end 

println("Errors : \n", errors)