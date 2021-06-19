# basis for 2 Gaussians
using Libdl
dir = pwd()
push!(Libdl.DL_LOAD_PATH,dir)
# 
include("MCMC.jl"); using .MCMC
using LinearAlgebra
# MCMC parameters
dt = 0.5
gamma = 2.0  
# equilibrate using 50,000 steps
N = 50000
using Random
rng = MersenneTwister(1234)
q = Vector{Float64}(undef, 2*N); q[1:2] = [0., 0.]
p = Vector{Float64}(undef, 2*N); p[1:2] = [1., 1.]
d = 4.4
MCMC.threeGauss!!(q, p, gamma=gamma, d=d, dt=dt, seed=235711)
q1, p1 = q[end], p[end]
# compute gammaStar from Markov chain
N = 20000000  # 20 M
q = Vector{Float64}(undef, 2*N); q[1] = q1
p = Vector{Float64}(undef, 2*N); p[1] = p1
MCMC.threeGauss!!(q, p, gamma=gamma, d=d, dt=dt, seed=235711)
q1, p1 = q[end], p[end]
x = q[1:2:end-1]; y = q[2:2:end]
Cov = (1/N)*[(x .- sum(x)/N)'*(x .- sum(x)/N) (x .- sum(x)/N)'*(y .- sum(y)/N);
       (y .- sum(y)/N)'*(x .- sum(x)/N) (y .- sum(y)/N)'*(y .- sum(y)/N)]
using LinearAlgebra
gammaStar = 1/sqrt(maximum(LinearAlgebra.eigvals(Cov)))
println("gammaStar = ", gammaStar)
# generate Markov chain
N = 50000000  # 10 M
q = Vector{Float64}(undef, 2*N); q[1] = q1
p = Vector{Float64}(undef, 2*N); p[1] = p1
MCMC.threeGauss!!(q, p, gamma=gamma, d=d, dt=dt, seed=235711)
x = q[1:2:end-1]; y = q[2:2:end]
for m = 1:2  # m = degree
  println("\ndegree = ", m)
  ub = zeros(m==1 ? 2 : 5, N)
  ub[1,:] = x
  ub[2,:] = y
  if m == 2
    ub[3,:] = x.*x
    ub[4,:] = x.*y
    ub[5,:] = y.*y
  end
  eigenTau, eachTau, reducs0 = MCMC.taus0!(copy(ub))
  eigenTau, eachTau, reducs2 = MCMC.taus2!(copy(ub))
  if reducs2 != reducs0  println("reducs2 !=	reducs0") end
  println("taumax = ", maximum(eigenTau.values))
  # in the same order, the values of tau are 152.0 for m = 1,
  # 153.8 for m = 2; N = 5000000  d = 4
  # in the same order, the values of tau are 2777.8 for m = 1,
  # 2822.9 for m = 2; N = 10 000000  d = 4
  # in the same order, the values of tau are 2066.5 for m = 1,
  # 2100.3 for m = 2; N = 10 000000  d = 4 gammaStar = 0.307
  # in the same order, the values of tau are 25154 for m = 1,
  # 38465 for m = 2; N = 5 000000  d = 4.8 gammaStar = 0.261
  # in the same order, the values of tau are 18774 for m = 1,
  # 19408 for m = 2; N = 10 000000  d = 4.8 gammaStar = 0.261
end  # for m
