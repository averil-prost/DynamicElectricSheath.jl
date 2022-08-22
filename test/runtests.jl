using Test
using DynamicElectricSheath

@testset "Initialization" begin

    println("Welcome in main.")

    physics = Physics()

    ν = physics.ν
    λ = physics.λ
    T = physics.T
    μ = physics.μ
    xmin, xmax = physics.xmin, physics.xmax
    vmin, vmax = physics.vmin, physics.vmax

    numerics = Discretization(physics)

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
    fi = fi_0.(xx, vv')    # electron density 
    fe = fe_0.(xx, vv')    # ion density
    ρ = zeros(Nx + 1)       # charge density
    ρi = zeros(Nx + 1)        # ion charge density
    ρe = zeros(Nx + 1)        # electron charge density

    t = 0.0
    iplot = 0

    # Boundary conditions (all 0, never updated afterwards)
    fi[begin, :] .= 0.0
    fi[end, :] .= 0.0 # speed distribution is almost 0 
    fi[:, begin] .= 0.0
    fi[:, end] .= 0.0 # non-emmiting wall
    fe[begin, :] .= 0.0
    fe[end, :] .= 0.0 # speed distribution is almost 0 
    fe[:, begin] .= 0.0
    fe[:, end] .= 0.0 # non-emmiting wall


    compute_charge!(ρi, fi, dv)
    compute_charge!(ρe, fe, dv)
    compute_charge!(ρ, fi .- fe, dv)

    J_l, J_r = compute_current(fi, fe, vv, dv)

    EE_minus, EE_plus = compute_e!(EE, ρ, λ, J_l, J_r, dx, dt)

    advection!(fi, fe, vv_plus, vv_minus, EE_plus, EE_minus, ν, μ, dx, dv, dt)

    @test true

end
