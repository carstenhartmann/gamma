# optimal gamma for LeMa13 example
include("MCMC.jl"); using .MCMC
using LinearAlgebra
# MCMC parameters
dt = 0.2
# equilibrate using 50,000 steps
gamma = 1.
N = 50000
using Random
rng = MersenneTwister(1234)
q = Vector{Float64}(undef, N); q[1] = 0.
p = Vector{Float64}(undef, N); p[1] = 1.
MCMC.LeMa13!!(q, p, gamma=gamma, dt=dt, seed=235711)
q1, p1 = q[end], p[end]
# compute gramma* from Markov chain
N = 2000000  # 2 M
q = Vector{Float64}(undef, N); q[1] = q1
p = Vector{Float64}(undef, N); p[1] = p1
MCMC.LeMa13!!(q, p, gamma=gamma, dt=dt, seed=235711)
q1, p1 = q[end], p[end]
mean = sum(q[:])/N
sigma2 = sum((q[:] .- mean).^2)/N
gammaStar = 1/sqrt(sigma2)
println("gammaStar = ", gammaStar)
gamma = collect(1:25)/5
imax = size(gamma,1)
taumax = Vector(undef, imax)
N = 10000000  # 10 M 
for i = 1:imax
  # generate Markov chain
  q = Vector{Float64}(undef, N); q[1] = q1
  p = Vector{Float64}(undef, N); p[1] = p1
  MCMC.LeMa13!!(q, p, gamma=gamma[i], dt=dt, seed=235711)
  m = 7
  ub = zeros(m, N)
  for i = 1:m  ub[i,:] = q[:].^i end
  eigenTau, eachTau, reducs0 = MCMC.taus0!(copy(ub))
  eigenTau, eachTau, reducs2 = MCMC.taus2!(copy(ub))
  if reducs2 != reducs0
    println("reducs2 = ", reducs2, ", reducs0 = ", reducs0) end
  taumax[i] = maximum(eigenTau.values)
  println("taumax[", i,"] = ", taumax[i])
end
using PyPlot
plot(gamma./gammaStar, taumax)
savefig("optLeMa.png")
