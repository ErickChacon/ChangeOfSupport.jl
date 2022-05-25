using Literate, Documenter, ChangeOfSupport

# Filenames
jls = filter(
    x -> occursin(r"/[0-9]+-.*\.jl$|/index.jl$", x),
    readdir("literate", join = true)
)

# Create markdown files
# rm(joinpath("docs", "src"), recursive = true, force = true)
repo_url = "https://github.com/ErickChacon/ChangeOfSupport.jl/blob/main"
Literate.markdown.(jls, joinpath("src"), documenter = true,
    repo_root_url = repo_url, credit = false)

makedocs(sitename = "ChangeOfSupport Documentation")


