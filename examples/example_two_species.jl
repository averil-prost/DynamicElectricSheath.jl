using DynamicElectricSheath
using ProgressMeter
using Plots 
using DelimitedFiles

function run_simulation(Nx, Nv; savedata=false)

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
    i0 = ceil(Int, Nv/2) + 1

    println("Parameters : ")
    println("ν = $ν, λ = $λ, μ = $μ")
    println("T = $T, Nt = $Nt, dt = $dt (CFL_x = $CFL_x, CFL_v = $CFL_v)")
    println("[xmin,xmax] = [$xmin,$xmax], Nx = $Nx, dx = $dx, i0=$i0")
    println("[vmin,vmax] = [$vmin,$vmax], Nv = $Nv, dv = $dv")

    xx = LinRange(xmin, xmax, Nx + 1) # vector (vertical) of spatial grid points (shared by e and i)
    vv = LinRange(vmin, vmax, Nv + 1) # vector (vertical) of speeds (shared by e and i)
    vv_plus = vv .* (vv .> 0.0)
    vv_minus = vv .* (vv .< 0.0) # signed positive and negative parts
    

    # Initialisation
    EE = E0.(xx)         # electric field

    if (IOparam.Init_source=="analytical")
        fi = fi_0.(xx, vv', μ=μ)  # electron density 
        fe = fe_0.(xx, vv', μ=μ)  # ion density
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
    fi[begin,i0:end] .= 0.0 # non-emmiting wall
    fi[end,begin:i0] .= 0.0 
    fi[:, begin] .= 0.0 # speed distribution is almost 0 
    fi[:, end] .= 0.0 
    fe[begin,i0:end] .= 0.0 # non-emmiting wall
    fe[end,begin:i0] .= 0.0 
    fe[:, begin] .= 0.0 # speed distribution is almost 0 
    fe[:, end] .= 0.0 

    source_fe = 0.0 * fe # no source term for the electrons
    # fi_old = 1.0 * fi # container for parallel
    # fe_old = 1.0 * fe # container for parallel

    @showprogress 1 for n = 1:Nt # loop over time
    # @gif for n = 1:Nt # loop over time

        # compute_charge!(ρi, fi, dv)
        # compute_charge!(ρe, fe, dv)
        # ρ = ρi - ρe
        compute_charge!(ρ, fi .- fe, dv)

        # J_l, J_r = compute_current(fi, fe, vv, dv)

        EE_minus, EE_plus = compute_e!(EE, ρ, λ, Nx, dx)
        # EE_minus, EE_plus = compute_e!(EE, ρ, λ, J_l, J_r, dx, dt)

        # advection!(fi, fe, vv_plus, vv_minus, EE_plus, EE_minus, ν, μ, dx, dv, dt)

        # (a / μ)₊ = (- a)₊ * (- μ) = - a₋ * (- μ) = a₋ * μ

        advection!(fe, vv_plus, vv_minus, -EE_minus/μ, -EE_plus/μ, source_fe, dx, dv, dt, i0)
        advection!(fi, vv_plus, vv_minus, EE_plus, EE_minus, ν * fe, dx, dv, dt, i0)
        # advection!(fe, fe_old, vv_plus, vv_minus, -EE_minus/μ, -EE_plus/μ, source_fe, dx, dv, dt)
        # advection!(fi, fi_old, vv_plus, vv_minus, EE_plus, EE_minus, ν * fe, dx, dv, dt)

        # if (mod(n,ceil(Nt/10))==0)
        #     # display(plot(xx, EE,legend=false,title="t = $(n*dt)")); readline()
        #     # display(heatmap(min.(1.0,fi'),title="fᵢ, t = $(n*dt)")); readline()
        #     display(heatmap(fe',title="fₑ, t = $(n*dt)")); readline()
        # end

        # display(heatmap(fe',title="fₑ, t = $(n*dt)")) #; readline()
        # heatmap(fe',title="fₑ, t = $(n*dt)") #; readline()
        # heatmap(fi',title="fᵢ, t = $(n*dt)") #; readline()
        # plot(xx, EE,legend=false,title="t = $(n*dt)")

    end # every 50

    # Saving datas
    if (savedata)
        try 
            # thefolder = "$(IOparam.Output_folder)/run_comp_short_time_2sp_Nx$(Nx)_Nv$(Nv)_Nt$(Nt)"
            thefolder = "$(IOparam.Output_folder)/run_comp_long_time_2sp_Nx$(Nx)_Nv$(Nv)_Nt$(Nt)"
            mkdir("$thefolder")
            # raw outputs
            println("Saving datas in $thefolder...")
            for (tag, data) in zip(["fi","fe","EE","rho"],[fi,fe,EE,ρ])
                writedlm("$thefolder/$(tag).dat", data)
            end
            # diags 
            imin = floor(Int, Nv * (-8.0 - vmin) / (vmax - vmin)) + 1
            imax = floor(Int, Nv * ( 8.0 - vmin) / (vmax - vmin)) + 1
            # png(plot(xx, EE, title="Electric field at T=$T", legend=false), "$thefolder/E.png")
            # png(plot(xx, ρ, title="Density ρ at T=$T", legend=false), "$thefolder/rho.png")
            # png(plot(heatmap(xx, vv, fe',clim=(0.0,0.04), title="Electron density at T=$T")), "$thefolder/fe.png")
            # png(plot(heatmap(xx, vv[imin:imax], fi'[imin:imax,:],clim=(0.0,0.5), title="Ion density at T=$T")), "$thefolder/fi.png")
            png(plot(xx, EE, legend=false), "$thefolder/E.png")
            png(plot(xx, ρ, legend=false), "$thefolder/rho.png")
            png(plot(heatmap(xx, vv, fe', clim=(-0.002,0.04),color=:gnuplot)), "$thefolder/fe.png")
            png(plot(heatmap(xx, vv[imin:imax], fi'[imin:imax,:],clim=(0.0,0.5),color=:gnuplot)), "$thefolder/fi.png")
            # parameters
            params = [Nx, Nv, Nt, T, λ, μ, ν, xmin, xmax, vmin, vmax, imin, imax]
            writedlm("$thefolder/params.dat", params)
            println("... Done.")
        catch exception 
            println("Error while saving datas.")
            throw(exception)
        end
    end

end

function generate_diags(Nx, Nv, Nt)
    # thefolder = "data/two_species/run_comp_short_time_2sp_Nx$(Nx)_Nv$(Nv)_Nt$(Nt)"
    thefolder = "data/two_species/run_comp_long_time_2sp_Nx$(Nx)_Nv$(Nv)_Nt$(Nt)"

    Nx, Nv, Nt, T, λ, μ, ν, xmin, xmax, vmin, vmax, imin, imax = readdlm("$thefolder/params.dat")
    Nx = round(Int,Nx); Nv = round(Int,Nv); Nt = round(Int,Nt); imin=round(Int,imin); imax=round(Int,imax)
    xx = LinRange(xmin,xmax,Nx+1)
    vv = LinRange(vmin,vmax,Nv+1)

    # EE = readdlm("$thefolder/EE.dat")
    # # display(plot(xx, EE, title="Electric field at T=$T", legend=false))
    # png(plot(xx, EE, title="Electric field at T=$T", zlims=(-5.0,5.0), legend=false), "$thefolder/E.png")

    # rho = readdlm("$thefolder/rho.dat")
    # # display(plot(xx, rho, title="Density ρ at T=$T", legend=false))
    # png(plot(xx, rho, title="Density ρ at T=$T", zlims=(0.0,1.6), legend=false), "$thefolder/rho.png")

    fi = readdlm("$thefolder/fi.dat")
    # display(plot(heatmap(xx,vv,fi')))
    # display(plot(heatmap(xx, vv[imin:imax], fi'[imin:imax,:])))
    # png(plot(heatmap(xx, vv[imin:imax], fi'[imin:imax,:], clim=(0.0,0.5), title="Ion density at T=$T")), "$thefolder/fi.png")
    png(plot(heatmap(xx, vv[imin:imax], fi'[imin:imax,:], clim=(0.0,2.2), color=:gnuplot)), "$thefolder/fi.png")

    # fe = readdlm("$thefolder/fe.dat")
    # # display(plot(heatmap(xx, vv, fe')))
    # # png(plot(heatmap(xx, vv, fe', clim=(0.0,0.04), title="Electron density at T=$T")), "$thefolder/fe.png")
    # png(plot(heatmap(xx, vv, fe', clim=(0.0,0.04), color=:gnuplot)), "$thefolder/fe.png")

end

#########################
# Main 
#########################

# # # Nxs = [100, 200, 400, 800]
# # # Nvs = [100, 200, 400, 800]
# Nxs = [100]
# Nvs = [800]

# for Nv in Nvs
#     for Nx in Nxs
#         # println("\n\n\nNx = $Nx, Nv = $Nv\n\n")
#         @time run_simulation(Nx, Nv, savedata=true)
#     end
# end

# Nxs = [100,100,100,100,200,200,200,200,400,400,400,400,800,800,800,800,1000]
# Nvs = [100,200,400,800,100,200,400,800,100,200,400,800,100,200,400,800,2000]
# Nts = [625,625,1000,2000,1250,1250,1250,2000,2500,2500,2500,2500,5000,5000,5000,5000,6250]

Nxs = [200]
Nvs = [400]
Nts = [2500000]

for (Nx, Nv, Nt) in zip(Nxs, Nvs, Nts)
    generate_diags(Nx, Nv, Nt)
end