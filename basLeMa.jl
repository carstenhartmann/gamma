# basis for LeMa13 example
include("MCMC.jl"); using .MCMC
using LinearAlgebra
# MCMC parameters
dt = 0.2
gamma = 1.
# equilibrate using 50,000 steps
N = 50000
using Random
rng = MersenneTwister(1234)
q = Vector{Float64}(undef, N); q[1] = 0.
p = Vector{Float64}(undef, N); p[1] = 1.
MCMC.LeMa13!!(q, p, gamma=gamma, dt=dt, seed=235711)
q1, p1 = q[end], p[end]
# compute gamma* from Markov chain
N = 2000000  # 2 M
q = Vector{Float64}(undef, N); q[1] = q1
p = Vector{Float64}(undef, N); p[1] = p1
MCMC.LeMa13!!(q, p, gamma=gamma, dt=dt, seed=235711)
q1, p1 = q[end], p[end]
mean = sum(q[:])/N
sigma2 = sum((q[:] .- mean).^2)/N
gammaStar = 1/sqrt(sigma2)
# generate Markov chain
N = 10000000  # 10 M
q = Vector{Float64}(undef, N); q[1] = q1
p = Vector{Float64}(undef, N); p[1] = p1
MCMC.LeMa13!!(q, p, gamma=gammaStar, dt=dt, seed=235711)
using PyPlot
# plot unnormalized pdf for −2 ≤ q ≤ 2.
  qs = collect(-2:1/128:2)
  pdf = exp.(-0.25*qs.^4 - sin.(1 .+ 5*qs))
  plot(qs, pdf, ":")
for m = 3:2:7  # m = degree
  println("\ndegree = ", m)
  ub = zeros(m, N)
  for i = 1:m  ub[i,:] = q[:].^i end
  eigenTau, eachTau, reducs0 = MCMC.taus0!(copy(ub))
  eigenTau, eachTau, reducs2 = MCMC.taus2!(copy(ub))
  if reducs2 != reducs0  println("reducs2 !=	reducs0") end
  println("taumax = ", maximum(eigenTau.values))
  x = eigenTau.vectors[:,end]
  qs = collect(-2:1/128:2)
  npts = size(qs,1)
  ubs = zeros(m, npts) 
  for i = 1:m  ubs[i,:] = qs[:].^i end
  # flip sign for m = 3 and 7
  if (m == 3 || m == 7) x = -x end
  plot(qs, ubs'*x, "-")
  # top to bottom: red=7, green=5, orange=3
  # in the same order, the values of tau are 69.5, 67.7, 62.7; N = 200 000
  # in the same order, the values of tau are 68.9, 67.0, 62.8; N = 2000000
  # in the same order, the values of tau are 66.6, 64.8, 60.5; N = 4000000
  # in the same order, the values of tau are 65.8, 63.9, 59.8; N = 10000000
end  # for m
savefig("eigfuncsLeMa.png")
# function changes sign from one mode to theother
