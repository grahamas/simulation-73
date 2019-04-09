using Documenter, Simulation73

makedocs(
    sitename="Simulation73 Documentation",
    pages=["index.md"],
    modules=[Simulation73]
)

deploydocs(
    repo = "github.com/grahamas/Simulation73.git"
)
