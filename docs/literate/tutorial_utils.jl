module TutorialUtils

function docs_generated_dir(name)
    docs_source_directory = if haskey(ENV, "HYDROTRIXI_DOCS_SRC")
        ENV["HYDROTRIXI_DOCS_SRC"]
    else
        candidate = normpath(joinpath(@__DIR__, "..", "src"))
        isdir(candidate) || error("Set HYDROTRIXI_DOCS_SRC before executing this tutorial.")
        candidate
    end

    return mkpath(joinpath(docs_source_directory, "assets", "generated", name))
end

end
