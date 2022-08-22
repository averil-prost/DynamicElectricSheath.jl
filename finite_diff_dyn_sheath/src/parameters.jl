using Parameters

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
@with_kw struct PhysicalParameters
    "Collision frequency"
    ν :: Float64  =  20.0 
    "Debye length"
    λ :: Float64  =  0.5 
    "Time horizon"
    T :: Float64  =  10. 
end

"""
$(TYPEDEF)
    grid in velocity corresponds to the electronic one which contains the support of boths ions and electron density function
$(TYPEDFIELDS)
"""
@with_kw struct DiscretizationParameters
    "Number of points of the space mesh "
    Nx::Int      =  600; 
    "Number of points of the speed mesh "
    Nv::Int      =  1001; 
    "Lower bound space mesh "
    xmin::Float64    = -1.0 
    "Upper bound space mesh "
    xmax::Float64    =  1.0; 
    "Lower bound speed mesh "
    vmin::Float64    = -10. #-5.0/sqrt(mu); 
    "Uupper bound speed mesh "
    vmax::Float64    =  10. #5.0/sqrt(mu); 
    "Space step"
    dx::Float64      =  (xmax - xmin)/Nx; # space step
    "Speed step"
    dv::Float64      =  (vmax - vmin)/Nv; # speed step
    "CFL in x"
    CFL_x::Float64   =  0.5*dx/vmax;
    "CFL in v"
    CFL_v::Float64   =  mu * dv/10;
    "Number of time steps"
    Nt::Float64      =  floor(Int,T/min(CFL_x,CFL_v))+1
    "Time step"
    dt::Float64      =  T/Nt;
end
