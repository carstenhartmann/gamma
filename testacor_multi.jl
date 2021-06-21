# generate multiple realisations of testacor.jl
include("MCMC.jl"); using .MCMC
using Libdl 
using Random 
using LinearAlgebra
# init
dir = pwd()
push!(Libdl.DL_LOAD_PATH,dir)
rng = MersenneTwister(1234)
# MCMC parameters
dt = 0.5
gamma = 2.
# equilibrate using 50,000 steps
N = 50000
z0 = [0., 1]
z = MCMC.xact(N, gamma, dt, z0, rng)
z0 = z[:,end]
# start Nruns independent realisation
Nruns = 1000
N = 2 .^collect(13:22)
imax = size(N,1)
tauS = Vector{Float64}(undef, imax) 
imin = 1
for run = 1:Nruns
  tauS = zeros(imax)
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
  end # for i 
  open("iACT.txt","a") do io
    println(io,tauS)
end # write file
end  # for run
tau2 = coth(-0.5*dt*[0 1 0; -2 -gamma 2 ; 0 -1 -2*gamma])[1,1]
tau3 =
 coth(-0.5*dt*[0 1 0 0; -3 -gamma 2 0; 0 -2 -2*gamma 3; 0 0 -1 -3*gamma])[1,1]
xact_tau = (tau2 + tau3)/2
open("iACT.txt","a") do io
    println(io,xact_tau*ones(imax))
end # write file
println("N[imin]:N[imax] = ", N[imin], ":", N[imax])
relerror = abs.(tauS[imin:end] .- xact_tau)/xact_tau
N = N[imin:end]
using PyPlot
loglog(N, abs.(relerror))
savefig("loglog.png")
