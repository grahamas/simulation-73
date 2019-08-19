Pkg.activate("docs/")
using Test, Documenter, Simulation73
doctest(Simulation73)
