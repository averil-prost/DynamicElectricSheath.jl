using Parameters

export Physics

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
@with_kw struct Physics
    "Collision frequency"
    ν::Float64 = 20.0
    "Debye length"
    λ::Float64 = 0.5
    "Time horizon"
    T::Float64 = 10.0
    "Thermal speed"
    μ::Float64 = 0.5
    "Lower bound space "
    xmin::Float64 = -1.0
    "Upper bound space "
    xmax::Float64 = 1.0
    "Lower bound speed "
    vmin::Float64 = -10.0 #-5.0/sqrt(μ); 
    "Uupper bound speed "
    vmax::Float64 = 10.0 #5.0/sqrt(μ); 
end

export Discretization

"""
$(TYPEDEF)

Grid in velocity corresponds to the electronic one which contains the 
support of boths ions and electron density function

$(TYPEDFIELDS)
"""
struct Discretization
    "Number of points of the space mesh "
    Nx::Int
    "Number of points of the speed mesh "
    Nv::Int
    "Space step"
    dx::Float64
    "Speed step"
    dv::Float64
    "CFL in x"
    CFL_x::Float64
    "CFL in v"
    CFL_v::Float64
    "Number of time steps"
    Nt::Int
    "Time step"
    dt::Float64

    function Discretization(physics::Physics; Nx = 600, Nv = 1001)

        dx = (physics.xmax - physics.xmin) / Nx
        dv = (physics.vmax - physics.vmin) / Nv
        CFL_x = 0.5 * dx / physics.vmax
        CFL_v = physics.μ * 0.1dv
        Nt = floor(Int, physics.T / min(CFL_x, CFL_v)) + 1
        dt = physics.T / Nt

        new(Nx, Nv, dx, dv, CFL_x, CFL_v, Nt, dt)

    end
end
