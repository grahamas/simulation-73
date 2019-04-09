using Pkg; Pkg.activate(pwd())

using Documenter, Simulation73

makedocs(
    sitename="Simulation73 Documentation",
    pages=["index.md"]
)

deploydocs(
    repo = "github.com/grahamas/Simulation73.git"
)
