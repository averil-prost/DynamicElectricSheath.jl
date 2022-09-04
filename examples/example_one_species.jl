using DynamicElectricSheath
using ProgressMeter
using Plots 

function run_simulation(Nx)

    # ioparams = IOparameters()

    physics = Physics()

    T = physics.T
    xmin, xmax = physics.xmin, physics.xmax
    vmin, vmax = physics.vmin, physics.vmax

    numerics = Discretization(physics, Nx=Nx, Nv=2*Nx+1)

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
    source = 0.0 * f    # no source for the advection
    ρ = zeros(Nx + 1)   # charge density

    # Boundary conditions (all 0, never updated afterwards)
    f[begin, :] .= 0.0
    f[end, :] .= 0.0 # speed distribution is almost 0 
    f[:, begin] .= 0.0
    f[:, end] .= 0.0 # non-emmiting wall

    # @showprogress 1 for n = 1:Nt # loop over time
    # @gif for n = 1:Nt # loop over time
    for n = 1:Nt # loop over time

        compute_charge!(ρ, -f, dv)

        EE_minus, EE_plus = compute_e!(EE, ρ, 1.0, Nx, dx) # λ=1

        advection!(f, vv_plus, vv_minus, EE_minus, EE_plus, source, dx, dv, dt)

        # if (mod(n,ceil(Int,Nt/20))==0)
        #     errL∞ = maximum(abs.(EE - EExact))
        #     errL¹ = sum(abs.(EE-EExact)) * dx
        #     println("Nx : $Nx, Erreur à $n / $Nt : L∞ $errL∞, L¹ $errL¹")
        # end

        # plot(xx, EE)
    end # every 20


    errL∞ = maximum(abs.(EE - EExact))
    errL¹ = sum(abs.(EE-EExact)) * dx
    println("Nx : $Nx, Erreur : L∞ $errL∞, L¹ $errL¹")
    return errL∞, errL¹ 
end

#########################
# Main 
#########################

# Nxs = [50, 100, 200, 400]
Nxs = 100:10:300
# Nxs = [250]
errors = zeros(length(Nxs),2)

@time for (iNx, Nx) in enumerate(Nxs)
    errors[iNx,:] .= run_simulation(Nx)
end

println("Errors : \n", errors)