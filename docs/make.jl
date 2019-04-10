using Documenter, Simulation73

makedocs(
    sitename="Simulation73 Documentation",
    modules=[Simulation73],
    pages = [
       "Home" => "index.md",
       "Manual" => Any[
           "man/simulation.md"
           # "man/examples.md",
           # "man/syntax.md",
           # "man/doctests.md",
           # "man/latex.md",
           # hide("man/hosting.md", [
           #     "man/hosting/walkthrough.md"
           # ]),
           # "man/other-formats.md",
       ],
       "Library" => Any[
           "Public" => "lib/public.md",
           hide("Internals" => "lib/internals.md",
           Any[joinpath("lib/internals", f) for f in readdir("docs/src/lib/internals")]
           )
       ]
   ]
)

deploydocs(
    repo = "github.com/grahamas/Simulation73.git"
)
