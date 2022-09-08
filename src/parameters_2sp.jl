using Parameters

export E0 

"""
$(SIGNATURES)

Electric field
"""
E0(x) = 0.0

"""
$(SIGNATURES)

Initial mask
"""
mask(x) = 0.5 * (tanh((x + 0.8) ./ 0.1) - tanh((x - 0.8) ./ 0.1))

export fi_0

"""
$(SIGNATURES)

Ions distribution
"""
fi_0(x, v; μ = 0.5) = mask(x) * exp(-0.5 .* v^2) / sqrt(2π)

export fe_0

"""
$(SIGNATURES)

Electrons distribution
""" 
fe_0(x, v; μ = 0.5) = sqrt(μ) * mask(x) * exp(-0.5 * μ * v^2) / sqrt(2π)

export Physics

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
@with_kw struct Physics
    "Collision frequency"
    ν::Float64 = 20.0
    "Debye length"
    λ::Float64 = 0.5 # physical : 0.1
    "Time horizon"
    # T::Float64 = 0.025 # really short time
    # T::Float64 = 0.025 # short time
    T::Float64 = 200.0 # long time
    "Thermal speed"
    μ::Float64 = 0.01 # physical : 1/3000
    "Lower bound space "
    xmin::Float64 = -1.0
    "Upper bound space "
    xmax::Float64 = 1.0
    "Lower bound speed "
    vmin::Float64 = -50.0 #-5.0/sqrt(μ); 
    "Uupper bound speed "
    vmax::Float64 = 50.0 #5.0/sqrt(μ); 
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

    function Discretization(physics::Physics; Nx = 50, Nv = 100)

        dx = (physics.xmax - physics.xmin) / Nx
        dv = (physics.vmax - physics.vmin) / Nv
        CFL_x = 0.5 * dx / physics.vmax
        CFL_v = physics.μ * 0.5 * dv / 10 # |E| supposed bounded by 10
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
    "Initial data source"
    Init_source::String = "analytical" # "loaded"
    "Initial data folder"
    Init_folder::String = "data/two_species"
    "Output folder"
    Output_folder::String = "data/two_species"
end