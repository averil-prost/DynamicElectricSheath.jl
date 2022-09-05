using Parameters

export E0 

"""
$(SIGNATURES)

Electric field 
```math
    E(x) = - \\partial \\phi (x) = - x
```
"""
E0(x) = min.(1.0, max.(-1.0, -x))
# E0(x) = -x

export f_0

"""
$(SIGNATURES)

Particles distribution
"""
function f_0(x, v)
    coeff = 1.0 - x^2 - v^2
    1 / (pi * sqrt(abs(coeff))) * (coeff > 1e-8)
end

export Physics

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
@with_kw struct Physics
    "Time horizon"
    T::Float64 = 0.1
    "Lower bound space"
    xmin::Float64 = -1.5
    "Upper bound space"
    xmax::Float64 = 1.5
    "Lower bound speed"
    vmin::Float64 = -2.0 #-5.0/sqrt(μ); 
    "Upper bound speed"
    vmax::Float64 = 2.0 #5.0/sqrt(μ); 
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

    function Discretization(physics::Physics; Nx = 500, Nv = 1001)

        dx = (physics.xmax - physics.xmin) / Nx
        dv = (physics.vmax - physics.vmin) / Nv
        CFL_x = 0.5 * dx / physics.vmax
        CFL_v = 0.5 * dv / 10
        dt = min(0.8 * min(CFL_x, CFL_v), physics.T / 10)
        Nt = ceil(Int, physics.T / dt)

        # Nt = floor(Int, physics.T / min(CFL_x, CFL_v)) + 1
        # dt = physics.T / Nt

        new(Nx, Nv, dx, dv, CFL_x, CFL_v, Nt, dt)

    end
end

export IOparameters

"""
$(TYPEDEF)

Folders for inputs and outputs.

$(TYPEDFIELDS)
"""
@with_kw struct IOparameters 
    "Output folder"
    Output_folder::String = "data/one_species"
end