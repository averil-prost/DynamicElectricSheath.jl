using Test
using FiniteDiffDynSheath

@testset "Initialization" begin

println("Welcome in main.")

physics = Physics()

nu = physics.ν
lambda = physics.λ
T = physics.T
mu = physics.μ
xmin,xmax = physics.xmin, physics.xmax
vmin,vmax = physics.vmin, physics.vmax

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
println("nu = $nu, lambda = $lambda, mu = $mu")
println("T = $T, Nt = $Nt, dt = $dt (CFL_x = $CFL_x, CFL_v = $CFL_v)")
println("[xmin,xmax] = [$xmin,$xmax], Nx = $Nx, dx = $dx")
println("[vmin,vmax] = [$vmin,$vmax], Nv = $Nv, dv = $dv")

@test true

end

#=

# The discretized quantities 
xx = collect(LinRange(xmin,xmax,Nx+1)) # vector (vertical) of spatial grid points (shared by e and i)
vv = collect(LinRange(vmin,vmax,Nv+1)) # vector (vertical) of speeds (shared by e and i)
vv_plus = vv .* (vv .> 0.0); vv_minus = vv .* (vv .< 0.0) # signed positive and negative parts

# Initialisation
EE = E0(xx)         # electric field
fi = fi_0(xx,vv)    # electron density 
fe = fe_0(xx,vv)    # ion density
rho = 0.0*xx        # charge density
rhoi = 0.0*xx        # ion charge density
rhoe = 0.0*xx        # electron charge density

t = 0.0
iplot = 0

# Boundary conditions (all 0, never updated afterwards)
fi[begin,:] .*= 0.0; fi[end,:] .*= 0.0; # speed distribution is almost 0 
fi[:,begin] .*= 0.0; fi[:,end] .*= 0.0; # non-emmiting wall
fe[begin,:] .*= 0.0; fe[end,:] .*= 0.0; # speed distribution is almost 0 
fe[:,begin] .*= 0.0; fe[:,end] .*= 0.0; # non-emmiting wall

rho .= vec(sum(fi.-fe,dims=2)) * dv       
rhoi .= vec(sum(fi,dims=2)) * dv       
rhoe .= vec(sum(fe,dims=2)) * dv       


#contourf(xx[2:end-1],vv[fislice],fi[begin+1:end-1,fislice]',colormap=thecmap,title="Ions at t=$t/$T",xlabel="x",ylabel="v",linewidth=0);
#png("jlimages/fi$iplot.png")
contourf(xx[2:end-1],vv[2:end-1],fi[begin+1:end-1,begin+1:end-1]',colormap=thecmap,title="Ions at t=$t/$T",xlabel="x",ylabel="v",linewidth=0);
png("jlimages/fi$iplot.png")
contourf(xx[2:end-1],vv[2:end-1],fe[begin+1:end-1,begin+1:end-1]',colormap=thecmap,title="Electrons at t=$t/$T",xlabel="x",ylabel="v",linewidth=0);
png("jlimages/fe$iplot.png")
plot(xx,EE,legend=false,title="E at t=$t/$T",xlabel="x",ylabel="E") #,ylim=[-10.0,10.0]);
png("jlimages/E$iplot.png")
plot(xx,rho,legend=false,title="rho at t=$t/$T",xlabel="x",ylabel="rho") #,ylim=[-0.05,0.5]);
png("jlimages/initrho$iplot.png")
plot(xx,rhoi,legend=false,title="rhoi at t=$t/$T",xlabel="x",ylabel="rhoi") #,ylim=[-0.05,0.5]);
png("jlimages/initrhoi$iplot.png")
plot(xx,rhoe,legend=false,title="rhoe at t=$t/$T",xlabel="x",ylabel="rhoe") #,ylim=[-0.05,0.5]);
png("jlimages/initrhoe$iplot.png")

open("jlimages/init_rho_$iplot.txt", "w") do io
    writedlm(io, [xx rho rhoe rhoi EE])
end

try 
    global iplot=0;
    # Main temporal loop
    for n in 1:Nt
        global t = n * dt;
        # plots 
        if ((mod(n+1,floor(Int,Nt/8))==0) || (n==Nt))
            println("Iteration $n / $Nt, time t = ", dt * n)

            #contourf(xx[2:end-1],vv[fislice],fi[begin+1:end-1,fislice]',colormap=thecmap,title="Ions at t=$t/$T",xlabel="x",ylabel="v",linewidth=0);
            #png("jlimages/fi$iplot.png")
            contourf(xx[2:end-1],vv[2:end-1],fi[begin+1:end-1,begin+1:end-1]',colormap=thecmap,title="Ions at t=$t/$T",xlabel="x",ylabel="v",linewidth=0);
            png("jlimages/fi$iplot.png")
            save("jlfiles/fi$iplot.jld", "fi", fi)
            contourf(xx[2:end-1],vv[2:end-1],fe[begin+1:end-1,begin+1:end-1]',colormap=thecmap,title="Electrons at t=$t/$T",xlabel="x",ylabel="v",linewidth=0);
            png("jlimages/fe$iplot.png")
            save("jlfiles/fe$iplot.jld", "fe", fe)
            plot(xx,EE,legend=false,title="E at t=$t/$T",xlabel="x",ylabel="E") #,ylim=[-10.0,10.0]);
            png("jlimages/E$iplot.png")
            save("jlfiles/E$iplot.jld", "E", EE)
            plot(xx,rho,legend=false,title="rho at t=$t/$T",xlabel="x",ylabel="rho") #,ylim=[-0.05,0.5]);
            png("jlimages/rho$iplot.png")
            save("jlfiles/rho$iplot.jld", "rho", rho)
            plot(xx,rhoi,legend=false,title="rhoi at t=$t/$T",xlabel="x",ylabel="rhoi") #,ylim=[-0.05,0.5]);
            png("jlimages/rhoi$iplot.png")
            save("jlfiles/rhoi$iplot.jld", "rhoi", rhoi)
            plot(xx,rhoe,legend=false,title="rhoe at t=$t/$T",xlabel="x",ylabel="rhoe") #,ylim=[-0.05,0.5]);
            png("jlimages/rhoe$iplot.png")
            save("jlfiles/rhoe$iplot.jld", "rhoe", rhoe)

			open("jlimages/rho_$iplot.txt", "w") do io
    			writedlm(io, [xx rho rhoe rhoi EE])
			end

            iplot += 1;
        end

        # electronic charge density and total mass
        rho .= vec(sum(fi.-fe,dims=2)) * dv       
        rhoi .= vec(sum(abs.(fi),dims=2)) * dv       
       	rhoe .= vec(sum(abs.(fe),dims=2)) * dv       

        # compute current at boundaries : rectangle integration of \int_{v} v * (fi(t,+-1,v) - fe(t,+-1,v)) dv
        J_l = dv * sum(vv .* (fi[begin,:] - fe[begin,:]))
        J_r = dv * sum(vv .* (fi[end  ,:] - fe[end  ,:]))

        # update electric field : use formula E(x) = 0.5*(E(1) + E(-1)) + nu/2 * ( \int_{-1}^{x} rho(y)dy - \int_{x}^{1} rho(y)dy )
        # with E(t^{n+1},1) + E(t^{n+1},-1) approximated by E(t^n,1) + E(t^n,-1) + dt/lambda^2 * (J_l + J_r) (Ampère boundary conditions)
        EE .= 0.5 * (EE[begin]+EE[end] + dt/(lambda^2) * (J_l + J_r)); # the +- Mass_e cancels out (?????)
        EE[begin+1:end] .+= 0.5/(lambda^2) * 0.5*dx.*cumsum(rho[1:end-1] .+ rho[2:end]);                        # add forward  trapeze integral
        EE[begin:end-1] .-= 0.5/(lambda^2) * 0.5*dx.*cumsum(rho[end-1:-1:1] .+ rho[end:-1:2])[end:-1:begin];    # add backward trapeze integral
        EE_plus = max.(EE,0.0); EE_minus = min.(EE,0.0)

        # Update ion and electron densities
        # dt fi + v dx fi + E dv fi = nu * fe
        fi[2:Nx,2:Nv] += (
             - (dt/dx) .* (vv_plus[2:Nv]' .* (fi[2:Nx,2:Nv].-fi[1:Nx-1,2:Nv]) .+ vv_minus[2:Nv]' .* (fi[3:Nx+1,2:Nv].-fi[2:Nx,2:Nv]))
            .- (dt/dv) .* (EE_plus[2:Nx]  .* (fi[2:Nx,2:Nv].-fi[2:Nx,1:Nv-1]) .+ EE_minus[2:Nx]  .* (fi[2:Nx,3:Nv+1].-fi[2:Nx,2:Nv]))
            .+ dt * nu .* fe[2:Nx,2:Nv]);
        # dt fe + v dx fe - 1/mu E dv fe = 0
        fe[2:Nx,2:Nv] += (
             - (dt/dx)      .* (vv_plus[2:Nv]' .* (fe[2:Nx,2:Nv].-fe[1:Nx-1,2:Nv]) .+ vv_minus[2:Nv]' .* (fe[3:Nx+1,2:Nv].-fe[2:Nx,2:Nv]))
            .+ (dt/(dv*mu)) .* (EE_minus[2:Nx] .* (fe[2:Nx,2:Nv].-fe[2:Nx,1:Nv-1]) .+ EE_plus[2:Nx]   .* (fe[2:Nx,3:Nv+1].-fe[2:Nx,2:Nv])));
    end # while n
catch e
    println(e.msg) # just so that ctrl+c works
end # try-catch

=#
