using Documenter
using DocumenterInterLinks
using CairoMakie
using LaTeXStrings
using HydroTrixi
using Literate

# Fix for https://github.com/trixi-framework/Trixi.jl/issues/668
if (get(ENV, "CI", nothing) != "true") &&
   (get(ENV, "TRIXI_DOC_DEFAULT_ENVIRONMENT", nothing) != "true")
    push!(LOAD_PATH, dirname(@__DIR__))
end

const DOCS_ROOT = @__DIR__
const DOCS_SRC = joinpath(DOCS_ROOT, "src")
const DOCS_LITERATE = joinpath(DOCS_ROOT, "literate")
const DOCS_TUTORIALS = joinpath(DOCS_SRC, "tutorials")
const DOCS_GENERATED_ASSETS = joinpath(DOCS_SRC, "assets", "generated")

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

function prepare_generated_docs!()
    rm(DOCS_TUTORIALS; recursive = true, force = true)
    mkpath(DOCS_TUTORIALS)

    rm(DOCS_GENERATED_ASSETS; recursive = true, force = true)
    mkpath(DOCS_GENERATED_ASSETS)

    return nothing
end

function generate_tutorials!()
    ENV["HYDROTRIXI_DOCS_SRC"] = DOCS_SRC
    ENV["HYDROTRIXI_DOCS_LITERATE"] = DOCS_LITERATE

    for source in ("celia_1990.jl",)
        Literate.markdown(joinpath(DOCS_LITERATE, source),
                          DOCS_TUTORIALS;
                          execute = true,
                          flavor = Literate.DocumenterFlavor(),
                          credit = false)
    end

    return nothing
end

prepare_generated_docs!()
generate_tutorials!()

makedocs(; modules = [HydroTrixi],
         repo = Remotes.GitHub("uofs-simlab", "HydroTrixi.jl"),
         sitename = "HydroTrixi.jl",
         format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://tjbmontoya.com/HydroTrixi.jl",
                                  edit_link = "main",
                                  assets = String[]),
         pages = ["Home" => "index.md",
                  "Tutorials" => ["Overview" => "tutorials.md",
                                  "Celia (1990) infiltration problem" =>
                                      "tutorials/celia_1990.md"],
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
