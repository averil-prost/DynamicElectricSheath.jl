using DynamicElectricSheath
using ProgressMeter
using Plots 
using DelimitedFiles

function run_simulation(Nx, Nv)

    IOparam = IOparameters()

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

    # f = exp.(-100 .* ((xx .- 0.5).^2 .+ (vv').^2))
    f = f_0.(xx, vv')   # particle density
    # f_old = 1.0 * f     # to optimize parallel version
    source = 0.0 * f    # no source for the advection
    ρ = zeros(Nx + 1)   # charge density
    ρexact = -1.0 * (-1.0 .< xx .< 1.0)

    # Boundary conditions (all 0, never updated afterwards)
    f[begin,:] .= 0.0
    f[end,  :] .= 0.0 
    f[:, begin] .= 0.0
    f[:, end] .= 0.0 

    # errors : max over n of |ρ - ρₑₓ| (resp. |E - Eₑₓ|)
    errL∞ρ = errL¹ρ = errL∞E = errL¹E = 0.0

    @showprogress 1 for n = 1:Nt # loop over time
        compute_charge!(ρ, -f, dv)
        EE_minus, EE_plus = compute_e!(EE, ρ, 1.0, Nx, dx) # λ=1
        # advection!(f, vv_plus, vv_minus, EE_plus, EE_minus, source, dx, dv, dt)
        # EE_plus = max.(EE, 0.0); EE_minus = min.(EE, 0.00)
        advection!(f, vv_plus, vv_minus, EE_plus, EE_minus, source, dx, dv, dt)
        # advection!(f, f_old, vv_plus, vv_minus, EE_plus, EE_minus, source, dx, dv, dt)

        # if (mod(n,ceil(Int,Nt/20))==0)
        #     display(heatmap(f',title="n=$n")); readline()
        # end

        # Error on ρ
        errρ = abs.(ρ - ρexact)
        errL∞ρ = max(errL∞ρ, maximum(errρ))
        errL¹ρ = max(errL¹ρ, sum(errρ) * dx)
        # Error on E
        errE = abs.(EE - EExact)
        errL∞E = max(errL∞E, maximum(errE))
        errL¹E = max(errL¹E, sum(errE) * dx)
    end 

    # err = abs.(f - exp.(-100 .* ((xx).^2 .+ (vv' .+ 0.5).^2)))
    # errL∞ = maximum(err)
    # errL¹ = sum(err) * dx * dv
    # println("Nx : $Nx, Erreur : L∞ $errL∞, L¹ $errL¹")

    # Saving datas
    try 
        println("Saving f, E and ρ...")
        for (tag, data) in zip(["f","EE","rho"],[f,EE,ρ])
            writedlm("$(IOparam.Output_folder)/$(tag)_Nx$(Nx)_Nv$(Nv)_Nt$(Nt).dat", data)
        end
        println("... Done.")
    catch exception 
        println("Error while saving f, E and ρ.")
        throw(exception)
    end
    
    println("Erreur ρ ---- Nx : $Nx, Erreur : L∞ $errL∞ρ, L¹ $errL¹ρ")
    println("Erreur E ---- Nx : $Nx, Erreur : L∞ $errL∞E, L¹ $errL¹E")

    return Nx, Nv, Nt, T, errL∞E, errL¹E, errL∞ρ, errL¹ρ
end

#########################
# Main 
#########################

# Nxs = [100, 200, 400, 800]
# Nxs = 200:100:700
# Nxs = 700:100:1200
# Nvs = [200]
Nxs = [100, 200, 400, 800]
Nvs = [4097 for _ ∈ Nxs]
errors = zeros(length(Nxs),4)

errfile = open("data/one_species/errs.txt", "a")

@time for (irun, (Nx, Nv)) in enumerate(zip(Nxs, Nvs))
    oNx, oNv, oNt, dt, T, errL∞E, errL¹E, errL∞ρ, errL¹ρ = run_simulation(Nx, Nv)
    println(errfile, "$oNx\t$oNv\t$oNt\t$T\t$errL∞E\t$errL¹E\t$errL∞ρ\t$errL¹ρ")
    errors[:,irun] .= [errL∞E, errL¹E, errL∞ρ, errL¹ρ]
end 

close(errfile)

println("Errors : \n", errors)