# see documentation at https://juliadocs.github.io/Documenter.jl/stable/

using Documenter, npformation

makedocs(
    modules = [npformation],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "rebekaszabo95",
    sitename = "npformation.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/rebekaszabo95/npformation.jl.git",
    push_preview = true
)
