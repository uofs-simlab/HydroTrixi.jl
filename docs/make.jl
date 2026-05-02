using Documenter
using DocumenterInterLinks
using CairoMakie
using LaTeXStrings
using HydroTrixi

# Fix for https://github.com/trixi-framework/Trixi.jl/issues/668
if (get(ENV, "CI", nothing) != "true") &&
   (get(ENV, "TRIXI_DOC_DEFAULT_ENVIRONMENT", nothing) != "true")
    push!(LOAD_PATH, dirname(@__DIR__))
end

const REPOSITORY_ROOT = dirname(@__DIR__)
const DOCS_ROOT = @__DIR__
const DOCS_SRC = joinpath(DOCS_ROOT, "src")
const README_ASSETS = joinpath(REPOSITORY_ROOT, "assets", "images")
const DOCS_ASSETS = joinpath(DOCS_SRC, "assets", "images")

links = InterLinks("Trixi" => ("https://trixi-framework.github.io/TrixiDocumentation/stable/",
                               "https://trixi-framework.github.io/TrixiDocumentation/stable/objects.inv"))

fallbacks = ExternalFallbacks(; automatic = true)

DocMeta.setdocmeta!(HydroTrixi,
                    :DocTestSetup,
                    :(using HydroTrixi, CairoMakie, LaTeXStrings);
                    recursive = true)

if isnothing(Base.get_extension(HydroTrixi, :HydroTrixiVisualizationExt))
    error("HydroTrixiVisualizationExt did not load in the documentation environment.")
end

function sync_readme_assets!()
    isdir(README_ASSETS) || return nothing

    mkpath(DOCS_ASSETS)

    for (root, _, files) in walkdir(README_ASSETS)
        relative_root = relpath(root, README_ASSETS)
        destination_root = if relative_root == "."
            DOCS_ASSETS
        else
            joinpath(DOCS_ASSETS, relative_root)
        end
        mkpath(destination_root)

        for file in files
            cp(joinpath(root, file), joinpath(destination_root, file); force = true)
        end
    end

    return nothing
end

function sync_readme_homepage!()
    readme_text = read(joinpath(REPOSITORY_ROOT, "README.md"), String)
    write(joinpath(DOCS_SRC, "index.md"), readme_text)
    return nothing
end

sync_readme_assets!()
sync_readme_homepage!()

makedocs(; modules = [HydroTrixi],
         repo = Remotes.GitHub("uofs-simlab", "HydroTrixi.jl"),
         sitename = "HydroTrixi.jl",
         format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://tjbmontoya.com/HydroTrixi.jl",
                                  edit_link = "main",
                                  assets = String[]),
         pages = ["Home" => "index.md",
                  "Reference" => "reference.md",
                  "License" => "license.md"],
         plugins = [links, fallbacks])

deploydocs(; repo = "github.com/uofs-simlab/HydroTrixi.jl.git",
           deploy_repo = "github.com/tristanmontoya/HydroTrixi.jl.git",
           repo_previews = "github.com/tristanmontoya/HydroTrixi.jl.git",
           devbranch = "main",
           push_preview = all(!isempty,
                              (get(ENV, "GITHUB_TOKEN", ""),
                               get(ENV, "DOCUMENTER_KEY", ""))))
