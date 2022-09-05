using DynamicElectricSheath
using ProgressMeter
using Plots 
using DelimitedFiles

function run_simulation(Nx, Nv)

    IOparam = IOparameters()

    physics = Physics()

    ν = physics.ν
    λ = physics.λ
    T = physics.T
    μ = physics.μ
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
    println("ν = $ν, λ = $λ, μ = $μ")
    println("T = $T, Nt = $Nt, dt = $dt (CFL_x = $CFL_x, CFL_v = $CFL_v)")
    println("[xmin,xmax] = [$xmin,$xmax], Nx = $Nx, dx = $dx")
    println("[vmin,vmax] = [$vmin,$vmax], Nv = $Nv, dv = $dv")

    xx = LinRange(xmin, xmax, Nx + 1) # vector (vertical) of spatial grid points (shared by e and i)
    vv = LinRange(vmin, vmax, Nv + 1) # vector (vertical) of speeds (shared by e and i)
    vv_plus = vv .* (vv .> 0.0)
    vv_minus = vv .* (vv .< 0.0) # signed positive and negative parts

    # Initialisation
    EE = E0.(xx)         # electric field

    if (IOparam.Init_source=="analytical")
        fi = fi_0.(xx, vv')  # electron density 
        fe = fe_0.(xx, vv')  # ion density
    elseif (IOparam.Init_source=="analytical")
        finame = "$(IOparam.Init_folder)/fi_Nx$Nx_Nv$Nv.dat"
        fename = "$(IOparam.Init_folder)/fe_Nx$Nx_Nv$Nv.dat"
        try
            fi = readdlm(finame)
            fe = readdlm(fename)
        catch exception
            println("Error while loading fi and fe. Please ensure that $finame and $finame exist.")
            throw(exception)
        end
    else
        throw(ErrorException("unknown Ioparam.Init_source $(IOparam.Init_source)"))
    end
    ρ = zeros(Nx + 1)    # charge density
    # ρi = zeros(Nx + 1)   # ion charge density
    # ρe = zeros(Nx + 1)   # electron charge density

    # Boundary conditions (all 0, never updated afterwards)
    fi[begin, :] .= 0.0
    fi[end, :] .= 0.0 # speed distribution is almost 0 
    fi[:, begin] .= 0.0
    fi[:, end] .= 0.0 # non-emmiting wall
    fe[begin, :] .= 0.0
    fe[end, :] .= 0.0 # speed distribution is almost 0 
    fe[:, begin] .= 0.0
    fe[:, end] .= 0.0 # non-emmiting wall

    source_fe = 0.0 * fe # no source term for the electrons
    fi_old = 1.0 * fi # container for parallel
    fe_old = 1.0 * fe # container for parallel

    @showprogress 1 for n = 1:Nt # loop over time
    # @gif for n = 1:Nt # loop over time

        # compute_charge!(ρi, fi, dv)
        # compute_charge!(ρe, fe, dv)
        # ρ = ρi - ρe
        compute_charge!(ρ, fi .- fe, dv)

        # J_l, J_r = compute_current(fi, fe, vv, dv)

        EE_minus, EE_plus = compute_e!(EE, ρ, λ, Nx, dx) # λ=1
        # EE_minus, EE_plus = compute_e!(EE, ρ, λ, J_l, J_r, dx, dt)

        # advection!(fi, fe, vv_plus, vv_minus, EE_plus, EE_minus, ν, μ, dx, dv, dt)

        # (a / μ)₊ = (- a)₊ * (- μ) = - a₋ * (- μ) = a₋ * μ

        advection!(fe, fe_old, vv_plus, vv_minus, -EE_minus/μ, -EE_plus/μ, source_fe, dx, dv, dt)
        advection!(fi, fi_old, vv_plus, vv_minus, EE_plus, EE_minus, ν * fe, dx, dv, dt)

        # if (mod(n,ceil(Nt/10))==0)
        #     # display(plot(xx, EE,legend=false,title="t = $(n*dt)")); readline()
        #     # display(heatmap(min.(1.0,fi'),title="fᵢ, t = $(n*dt)")); readline()
        #     display(heatmap(fe',title="fₑ, t = $(n*dt)")); readline()
        # end

    end #every 50

    # Saving datas
    try 
        println("Saving fᵢ, fₑ, E and ρ...")
        for (tag, data) in zip(["fi","fe","EE","rho"],[fi,fe,EE,ρ])
            writedlm("$(IOparam.Output_folder)/$(tag)_Nx$(Nx)_Nv$(Nv).dat", data)
        end
        println("... Done.")
    catch exception 
        println("Error while saving fᵢ, fₑ, E and ρ.")
        throw(exception)
    end

end

#########################
# Main 
#########################

Nx = 200
Nv = 400
@time run_simulation(Nx, Nv)
