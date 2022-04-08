# using Pkg
# Pkg.activate(temp=true)
# Pkg.add("Trixi")

using Trixi
using Test
using LinearAlgebra

using Trixi: ln_mean, inv_ln_mean
import Trixi: flux_ranocha

function flux_ranocha_A(uA_ll, uA_rr, equations::CompressibleEulerEquations1D)
  gamma = equations.gamma
  u_ll = A2cons(uA_ll, equations)
  u_rr = A2cons(uA_rr, equations)
  A_ll, A_rr = uA_ll[4], uA_rr[4]
  A_avg = .5 * (A_ll + A_rr)

  rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
  rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

  rho_mean = ln_mean(rho_ll, rho_rr)
  inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
  v1_avg = 0.5 * (v1_ll + v1_rr)
  p_avg  = 0.5 * (p_ll + p_rr)
  velocity_square_avg = 0.5 * (v1_ll*v1_rr)

  # rhoA_mean = ln_mean(rho_ll * A_ll, rho_rr * A_rr)
  v1A_ll = v1_ll * A_ll 
  v1A_rr = v1_rr * A_rr
  v1A_avg = 0.5 * (v1A_ll + v1A_rr)

  # f1 = rho_mean * v1_avg 
  f1 = rho_mean * v1A_avg
  f2 = f1 * v1_avg + A_ll * p_avg
    f3 = f1 * ( velocity_square_avg + inv_rho_p_mean / (gamma-1) ) + 0.5 * (p_ll*v1A_rr + p_rr*v1A_ll)
  # f3 = f1 * ( velocity_square_avg + inv_rho_p_mean / (gamma-1) ) + 0.5 * (p_ll*v1_rr + p_rr*v1_ll)

  return SVector(f1, f2, f3)
end

function A2cons(uA, equations::CompressibleEulerEquations1D)
  A = uA[4]
  return SVector{3}(uA[1:3] ./ A)
end

function separable_flux(uA_ll, uA_rr, equations::CompressibleEulerEquations1D)
  A_ll = uA_ll[4]
  A_rr = uA_rr[4]
  u_ll = A2cons(uA_ll, equations)
  u_rr = A2cons(uA_rr, equations)
  _, v_ll, p_ll = cons2prim(u_ll, equations)
  _, v_rr, p_rr = cons2prim(u_rr, equations)
  p_avg = 0.5 * A_ll * (p_ll + p_rr)
  pv_mean = 0.5 * (p_ll * A_rr * v_rr + p_rr * A_ll * v_ll)
  return SVector{3}(0.0, p_avg, pv_mean)
end

equations = CompressibleEulerEquations1D(1.4)

# make some random solution state
A = 1.0 .+ rand(4)

u = zeros(SVector{3}, 4)
uA = zeros(SVector{4}, 4)
for i in 1:length(u)
  u_i = prim2cons(SVector{3}(2.0, 1.0, 3.0) .+ .1 * randn(3), equations)
  u[i] = u_i
  uA[i] = SVector{4}((u_i .* A[i])..., A[i])
end

# make an SBP operator
Q = Float64[-1  1 0 0;
            -1  0  1 0;
            0 -1 0 1;
            0 0 -1  1]/2

# skew symmetric part
S = Q - Q'

# primitive variables
q = cons2prim.(u, equations)
rho = getindex.(q, 1)
v1 = getindex.(q, 2)
p = getindex.(q, 3)

# entropy variables
v = cons2entropy.(u, equations)

# # check consistency b/w flux diff and mat-vec formulas
# r = sum(S .* separable_flux.(uA, uA', equations), dims=2)
# @test sum(-rho .* (Q * (A .* v1))) ≈ sum(map((x, y) -> dot(x, y), v, r))

# test if we get the entropy potential back
r = sum(S .* flux_ranocha_A.(uA, uA', equations), dims=2)
@test -sum(Q * (rho .* A .* v1)) ≈ sum(map((x, y) -> dot(x, y), v, r))