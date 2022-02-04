using Revise
using Plots
using StartUpDG
# using NodesAndModes
using LinearAlgebra

# polynomial degree and mesh size
N = 3
K1D = 8

# init ref element and mesh
rd = RefElemData(Tri(), N)
VXY, EToV = uniform_mesh(Tri(), K1D)
md = MeshData(VXY, EToV, rd)

# Define a function by interpolation
@unpack x, y = md
u = @. 2 + .5*exp(-100*(x^2 + y^2))

# Compute derivatives using geometric mapping + chain rule
@unpack Dr, Ds = rd
@unpack rxJ, sxJ, J = md
dudx = (rxJ .* (Dr*u) + sxJ .* (Ds*u)) ./ J