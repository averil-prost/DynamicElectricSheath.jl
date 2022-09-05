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
end


"""
$(SIGNATURES)

Upwind scheme for the advection equation (parallel version).
The matrix fs_old is allocated only once at the beginning of the program. 

```math
\\frac{\\partial f_s}{\\partial t} + a_{x} \\frac{\\partial f_s}{\\partial x} + a_{v} \\frac{\\partial f_s}{\\partial v} = \\text{source}
```
"""
function advection!(fs, fs_old, axp, axm, avp, avm, source, dx, dv, dt)
    Nx = size(fs, 1) - 1
    Nv = size(fs, 2) - 1

    fs_old .= 1*fs # to avoid concurrent access to fs
    for i=2:Nx
        @Threads.threads for j=2:Nv
            fs[i,j] += - (dt / dx) .* (
                axp[j] .* (fs_old[i,j] .- fs_old[i-1,j]) .+
                axm[j] .* (fs_old[i+1,j] .- fs_old[i,j])
            ) .-
            (dt / dv) .* (
                avp[i] .* (fs_old[i,j] .- fs_old[i,j-1]) .+
                avm[i] .* (fs_old[i,j+1] .- fs_old[i,j])
            ) .+ dt * source[i,j]
        end
    end
end