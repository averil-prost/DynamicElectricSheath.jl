export compute_charge!

"""
$(SIGNATURES)

compute charge density
"""
function compute_charge!(rho, f, dv)
    # rho .= vec(sum(f, dims = 2)) * dv # Rectangles
    rho .= vec(sum(f, dims = 2) .- 0.5 .* (f[:,begin] + f[:,end])) * dv # Trapezoids
end

export compute_current

"""
$(SIGNATURES)

compute current at boundaries : rectangle integration of 

```math
\\int_{v} v * (fi(t,\\pm 1,v) - fe(t,\\pm 1,v)) dv
```
"""
function compute_current(fi, fe, vv, dv)

    @views J_l = dv * sum(vv .* (fi[begin, :] - fe[begin, :]))
    @views J_r = dv * sum(vv .* (fi[end, :] - fe[end, :]))
    return J_l, J_r

end

export compute_e!

"""
$(SIGNATURES)

Update electric field under assumption E(0)=0 : use formula 

```math
E(x) = \\frac{1}{\\lambda^2} \\int_{0}^{x} \\rho(y)dy 
```
"""
function compute_e!(EE, rho, lambda, Nx, dx)
    i0 = floor(Int, Nx/2)+1
    factor = 0.5 * dx / (lambda*lambda)
    if (2*i0 == Nx+2)
        # Boundary condition
        EE[i0]=0.
        # Poisson solved on two subdomains
        for i in i0+1:Nx+1
            EE[i] = EE[i-1] + factor * (rho[i] + rho[i-1])
        end
        for i in i0-1:-1:1
            EE[i] = EE[i+1] - factor * (rho[i+1] + rho[i])
        end
    else
        throw(DomainError("0 is not part of the space mesh"))
    end
    EE_plus = max.(EE, 0.0)
    EE_minus = min.(EE, 0.0)
    return EE_minus, EE_plus
end

export advection!

# """
# $(SIGNATURES)

# ```math
# \\frac{d f_i}{dt} + v \\frac{d f_i}{dx} + E \\frac{d f_i}{dv} = \\nu * f_e
# ```

# ```math
# \\frac{d f_e}{dt} + v \\frac{d f_e}{dx} - \\frac{1}{\\mu} E \\frac{d f_e}{dv} = 0
# ```
# """
# function advection!(fi, fe, vv_plus, vv_minus, EE_plus, EE_minus, nu, mu, dx, dv, dt)

#     Nx = size(fi, 1) - 1
#     Nv = size(fi, 2) - 1

#     @views fi[2:Nx, 2:Nv] .+= (
#         -(dt / dx) .* (
#             vv_plus[2:Nv]' .* (fi[2:Nx, 2:Nv] .- fi[1:Nx-1, 2:Nv]) .+
#             vv_minus[2:Nv]' .* (fi[3:Nx+1, 2:Nv] .- fi[2:Nx, 2:Nv])
#         ) .-
#         (dt / dv) .* (
#             EE_plus[2:Nx] .* (fi[2:Nx, 2:Nv] .- fi[2:Nx, 1:Nv-1]) .+
#             EE_minus[2:Nx] .* (fi[2:Nx, 3:Nv+1] .- fi[2:Nx, 2:Nv])
#         ) .+ dt * nu .* fe[2:Nx, 2:Nv]
#     )
#     @views fe[2:Nx, 2:Nv] .+= (
#         -(dt / dx) .* (
#             vv_plus[2:Nv]' .* (fe[2:Nx, 2:Nv] .- fe[1:Nx-1, 2:Nv]) .+
#             vv_minus[2:Nv]' .* (fe[3:Nx+1, 2:Nv] .- fe[2:Nx, 2:Nv])
#         ) .+
#         (dt / (dv * mu)) .* (
#             EE_minus[2:Nx] .* (fe[2:Nx, 2:Nv] .- fe[2:Nx, 1:Nv-1]) .+
#             EE_plus[2:Nx] .* (fe[2:Nx, 3:Nv+1] .- fe[2:Nx, 2:Nv])
#         )
#     )

# end

"""
$(SIGNATURES)

Upwind scheme for the advection equation. 

```math
\\frac{\\partial f_s}{\\partial t} + a_{x} \\frac{\\partial f_s}{\\partial x} + a_{v} \\frac{\\partial f_s}{\\partial v} = \\text{source}
```
"""
function advection!(fs, axp, axm, avp, avm, source, dx, dv, dt)
    Nx = size(fs, 1) - 1
    Nv = size(fs, 2) - 1

    # center of the domain
    @views fs[2:Nx, 2:Nv] .+= (
        -(dt / dx) .* (
            axp[2:Nv]' .* (fs[2:Nx, 2:Nv] .- fs[1:Nx-1, 2:Nv]) .+
            axm[2:Nv]' .* (fs[3:Nx+1, 2:Nv] .- fs[2:Nx, 2:Nv])
        ) .-
        (dt / dv) .* (
            avp[2:Nx] .* (fs[2:Nx, 2:Nv] .- fs[2:Nx, 1:Nv-1]) .+
            avm[2:Nx] .* (fs[2:Nx, 3:Nv+1] .- fs[2:Nx, 2:Nv])
        ) .+ dt * source[2:Nx, 2:Nv]
    )

    i0 = findlast(axp .< 1e-8) # last index such that v <= 0.0
    # region (x=-1, v<=0)
    @views fs[1,2:i0] .+= (
        -(dt / dx) .* (
            axm[2:i0] .* (fs[2, 2:i0] .- fs[1, 2:i0])
        ) .-
        (dt / dv) .* (
            avp[1] .* (fs[1, 2:i0] .- fs[1, 1:i0-1]) .+
            avm[1] .* (fs[1, 3:i0+1] .- fs[1, 2:i0])
        ) .+ dt * source[1, 2:i0]
    )

    # region (x=1, v>0)
    @views fs[end,i0:Nv] .+= (
        -(dt / dx) .* (
            axp[i0:Nv] .* (fs[end, i0:Nv] .- fs[end-1, i0:Nv]) 
        ) .-
        (dt / dv) .* (
            avp[end] .* (fs[end, i0:Nv] .- fs[end, i0-1:Nv-1]) .+
            avm[end] .* (fs[end, i0+1:Nv+1] .- fs[end, i0:Nv])
        ) .+ dt * source[end, i0:Nv]
    )
end