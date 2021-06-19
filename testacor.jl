# test acor
using Libdl
dir = pwd()
push!(Libdl.DL_LOAD_PATH,dir)
include("MCMC.jl"); using .MCMC
using LinearAlgebra
# MCMC parameters
dt = 0.5
gamma = 2.
# equilibrate using 50,000 steps
N = 50000
z0 = [0., 1]
using Random
# rng = MersenneTwister(1234)
z = MCMC.xact(N, gamma, dt, z0, rng)
z0 = z[:,end]
Nruns = 1000
N = 2 .^collect(13:22)
imax = size(N,1)
# tauS[i] = sum of taus for N[i] over Nruns independent realisations
tauS = Vector{Float64}(undef, imax)  
imin = 1
for run = 1:Nruns
  z = MCMC.xact(N[imax], gamma, dt, z0, rng)
  q = z[1,:]
  # u(q) = He3(q)/sqrt(3!) − He2(q)/sqrt(2!)
  u = q.*(q.*q .- 3.)/sqrt(6.) - (q.*q .- 1.)/sqrt(2.)
  for i = 1:imax
    Ni = N[i]
    result = MCMC.ondiag!(copy(u[1:Ni]))
    if result === nothing
      imin = min(imin, i + 1)
    else
      C0, D, reducs = result
      tauS[i] += D/C0
    end
  end
end  # for run
tauS = tauS/Nruns
tau2 = coth(-0.5*dt*[0 1 0; -2 -gamma 2 ; 0 -1 -2*gamma])[1,1]
tau3 =
 coth(-0.5*dt*[0 1 0 0; -3 -gamma 2 0; 0 -2 -2*gamma 3; 0 0 -1 -3*gamma])[1,1]
xact_tau = (tau2 + tau3)/2
println("N[1:imin] = ", N[1:imin])
relerror = abs.(tauS[imin:end] .- xact_tau)/xact_tau
N = N[imin:end]
using PyPlot
loglog(N, abs.(relerror))
savefig("loglog.png")
