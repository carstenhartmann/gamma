module MCMC
# propagator generates time series of state variables
# omit implement of numerical version of Delta
# include("Europhys/MCMC.jl"); using .MCMC

# DELETE following after confirming agreement with ondiag!
# clang -dynamiclib my_acc.c -o my_acc.so
# function my_acor(X)
#   mean = Ref{Cdouble}(0.)
#   sigma = Ref{Cdouble}(0.)
#   tau = Ref{Cdouble}(0.)
#   ccall((:my_acor, "my_acc.so"), Cint,
#         (Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Cint),
#         mean, sigma, tau, X, size(X, 1))
#   mean[], sigma[], tau[]
# end

# clang -dynamiclib dcoeff.c -o dcoeff.so
function ondiag!(X)
  C0 = Ref{Cdouble}(0.)
  D = Ref{Cdouble}(0.)
  reducs = Ref{Cint}(0)
  error = ccall((:ondiag, "dcoeff.so"), Cint,
        (Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Cint, Ref{Cint}),
        C0, D, X, size(X, 1), reducs)
  if error != 0  return end
  C0[], D[], reducs[]
end

function offdiag!!(X, Y, reducs)
  C0 = Ref{Cdouble}(0.)
  D = Ref{Cdouble}(0.)
  ccall((:offdiag, "dcoeff.so"), Cint,
        (Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Cint, Cint),
        C0, D, X, Y, size(X, 1), reducs)
  C0[], D[]
end

# clang -dynamiclib propagate.c -o propagate.so
function twoGauss!!(q, p; gamma=gamma, d=d, dt=dt, seed=seed)
  ccall((:twoGauss, "propagate.so"), Cint,
        (Cint, Ref{Cdouble}, Ref{Cdouble}, Cdouble, Cdouble, Cdouble, Cint),
        size(q, 1), q, p, gamma, d, dt, seed)
end


# clang -dynamiclib propagate.c -o propagate.so
function threeGauss!!(q, p; gamma=gamma, d=d, dt=dt, seed=seed)
  # q = [x1, y1, x2, y2, ..., xN, yN], same for p
  ccall((:threeGauss, "propagate.so"), Cint,
        (Cint, Ref{Cdouble}, Ref{Cdouble}, Cdouble, Cdouble, Cdouble, Cint),
        size(q, 1), q, p, gamma, d, dt, seed)
end

# clang -dynamiclib propagate.c -o propagate.so
function LeMa13!!(q, p; gamma=gamma, dt=dt, seed=seed)
  #    q[1], p[1] are initial values
  ccall((:LeMa13, "propagate.so"), Cint,
        (Cint, Ref{Cdouble}, Ref{Cdouble}, Cdouble, Cdouble, Cint),
        size(q, 1), q, p, gamma, dt, seed)
end

using LinearAlgebra
function xact(N, gamma, dt, z0, rng) # z0 = initial z
  expAt = exp(dt*[0 1; -1 -gamma])
  delta = sqrt(Complex(gamma^2 - 4))
  exp_gt = exp(-gamma*dt)
  hdt = 0.5*delta*dt 
  sinhchdt = Real(delta == 0 ? 1 : sinh(hdt)/hdt)
  coshhdt = Real(cosh(hdt))
  Sigma = (1 - exp_gt)*[1 0; 0 1]
  Sigma -= 0.5*gamma*dt^2*exp_gt*sinhchdt^2*[gamma -2; -2 gamma]
  Sigma += exp_gt*gamma*dt*coshhdt*sinhchdt*[-1 0; 0 1]
  Sigmah = cholesky(Sigma).L
  m = size(z0, 1)
  zs = zeros(m,N)
  z = copy(z0)
  for n = 1:N
    @inbounds zs[:,n] = z
    z = expAt*z
    z = z + Sigmah*randn(rng, 2)
  end
  zs
end

function BAOAB(F, N, gamma, dt, z0, rng)  # actually ABBAO
  expgt = exp(-gamma*dt)
  sigma2 = 1 - expgt*expgt
  sigma = sqrt(sigma2)
  m = size(z0, 1)
  zs = zeros(m,N)
  z = copy(z0)
  for n = 1:N
    @inbounds zs[:,n] = z
    z[1] = z[1] + 0.5*dt*z[2]
    z[2] = z[2] + dt*F(z[1])
    z[1] = z[1] + 0.5*dt*z[2]
    z[2] = expgt*z[2]
    z[2] = z[2] + sigma*randn(rng)
  end
  zs
end

function tauArrays!(ub; reducs=nothing)
  # ub time series of vectors
  m, N = size(ub)
  C0 = zeros(m, m)
  K = zeros(m, m)
  if reducs === nothing
    reducs = 0
    for i = 1:m
      result = ondiag!(copy(ub[i,:]))
      if result === nothing  return end
      reducs = max(reducs, result[3])
      end
  end
  for i = 1:m
    for j = 1:m
      C0[i,j], K[i,j] = offdiag!!(copy(ub[i,:]), copy(ub[j,:]), reducs)
    end
  end
  C0, K, reducs
end

function taus!(ub; reducs=nothing)  # new algorithm uncorrected
  if reducs===nothing
    result = MCMC.tauArrays!(ub, reducs=nothing)
    if result === nothing  return
    else C0, K, reducs  = result end
    eigvals = LinearAlgebra.eigvals(0.5*(K + K'), C0)
  end
  C0, K, reducs = MCMC.tauArrays!(ub,reducs=reducs)
  eigenTau = LinearAlgebra.eigen(0.5*(K + K'), C0)
  eachTau = diag(0.5*(K + K')) ./ diag(C0)  # one by one
  eigenTau, eachTau, reducs
  # eigenTau.values, eigenTau.vectors are eigenvalues, eigenvectors
end

function taus2!(ub; reducs=nothing)  # new algorithm
  if reducs===nothing
    result = MCMC.tauArrays!(ub, reducs=nothing)
    if result === nothing  return
    else C0, K, reducs  = result end
    eigvals = LinearAlgebra.eigvals(0.5*(K + K'), C0)
    if minimum(eigvals) < 0 && reducs > 0
      reducs -= 1
    end
  end
  C0, K, reducs = MCMC.tauArrays!(ub,reducs=reducs)
  eigenTau = LinearAlgebra.eigen(0.5*(K + K'), C0)
  eachTau = diag(0.5*(K + K')) ./ diag(C0)  # one by one
  eigenTau, eachTau, reducs
  # eigenTau.values, eigenTau.vectors are eigenvalues, eigenvectors
end

function taus0!(ub)  # from FaCS20
  # ub time series of vectors
  m, N = size(ub)
  C0 = zeros(m, m)
  K = zeros(m, m)
  maxtau = 0
  for i = 1:m
    result = ondiag!(copy(ub[i,:]))
    if result === nothing  return end
    tau = result[2]/result[1]
    if tau > maxtau  reducs = result[3]; maxtau = tau end
  end
  converged = false
  while !converged
    for i = 1:m
      for j = 1:m
        C0[i,j], K[i,j] = offdiag!!(copy(ub[i,:]), copy(ub[j,:]), reducs)
      end
    end
    x = LinearAlgebra.eigvecs(0.5*(K + K'), C0)[:,end]
    result = ondiag!(ub'*x)
    if result === nothing  return end
    reducs1 = result[3]
    if reducs1 >= reducs  converged = true
    else reducs = reducs1
    end
  end
  C0, K, reducs = MCMC.tauArrays!(ub,reducs=reducs)
  eigenTau = LinearAlgebra.eigen(0.5*(K + K'), C0)
  eachTau = diag(0.5*(K + K')) ./ diag(C0)  # one by one
  eigenTau, eachTau, reducs
  # eigenTau.values, eigenTau.vectors are eigenvalues, eigenvectors
end

end
