# test taumax, taumax1
include("MCMC.jl"); using .MCMC
using LinearAlgebra
# MCMC parameters
dt = 0.5  # 2
N = 2000000  # 1 M
d = 4.8  # 3.
#:::println("\n dt = ", dt, ", d = ", d, ", N = ", N);
println("\n dt = ", dt, ", N = ", N);
println("gamma  E[q]  gamma*Var[q]^(1/2)  tau")
gamma = 1.  #1.:0.0625:1.625  # gamma = 1. is optimal for d = 0. 
  # U(x, y) = 0.5*(x^2 + y^2)
  #         - log(exp(dx-d^2/8) + 2exp(-dx/2)cosh(sqrt(3)dy/2))
  # U(x, y) = 0.5*(x^2 + y^2) - log(xi^2 + (2/xi)(eta+1/eta))
  # where xi = exp(dx/2) and eta = exp(sqrt(3)dy/2)
  # and because (d/dx)xi = (d/2)xi and (d/dy)eta = (sqrt(3)d/2)eta, one has
  # Fx = -x + d(xi^2 - (1/xi)(eta+1/eta))
  #           /(xi^2 + (2/xi)(eta+1/eta))
  # Fy = -y + sqrt(3)d(1/xi)(eta-1/eta)
  #           /(xi^2 + (2/xi)(eta+1/eta))
  # generate Markov chain
  q = Vector{Float64}(undef, 2*N); q[1] = 0.; q[2] = 0.
  p = Vector{Float64}(undef, 2*N); p[1] = 1.; p[2] = 1.
  MCMC.threeGauss!!(q, p, gamma=gamma, d=d, dt=dt, seed=235711)
  x = q[1:2:2*N-1]; y = q[2:2:2*N]
  using PyPlot;
  plot(x[1:80000], y[1:80000], ".")
  savefig("trajectory.png")
#=
  z = [q'; p']
  # end C version
  mean = sum(z[1,:])/N
  sigma = sqrt(sum((z[1,:] .- mean).^2)/N)
  m = 4
  U(z) = [z[1], z[1]^3, z[1]^5, z[1]^7]
  ub = zeros(4, N)
  for n = 1:N
    ub[:,n] = U(z[:,n])
  end
  eigenTau, eachTau, reducs = MCMC.taus0!(copy(ub[1:m,:]))
  tau = maximum(eigenTau.values)
  println(gamma, " ", mean, " ", gamma*sigma, " ", tau)
=#
