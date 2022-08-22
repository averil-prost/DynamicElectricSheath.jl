# DynamicElectricSheath.jl

Two-species Vlasov-Poisson solver on $[-1,1] \times \RR$ using
 
 with :
 - non-emitting boundary conditions on the distribution functions
 - Source term for the ions
 - Dynamic electric field at the boundary
 
 
 Simulation of the symmetric sheath problem with ionization term.
 
 Discretization :
 - Upwind for the transport
 - Integration for the Poisson equation -\\lambda^2 dx E = rho with trapeze formula for rho.  The integration  preserves the symmetry of rho if it is even.
 
 
 Initial data :
 
 - Initial data are in the file initial_data.hpp
 
 Author : Mehdi BASDI. 
 Date : 27/07/2022.
 Translation in Julia : Pierre Navaro & Averil Prost.
 Date : 05/08/2022.

```bash
git clone https://github.com/juliavlasov/DynamicElectricSheath.jl
cd DynamicElectricSheath.jl
julia --project
julia> using Pkg
julia> Pkg.instantiate()
julia> include("example.jl")
```
