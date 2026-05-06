module TutorialUtils

function docs_src_dir()
    if haskey(ENV, "HYDROTRIXI_DOCS_SRC")
        return ENV["HYDROTRIXI_DOCS_SRC"]
    end

    candidate = normpath(joinpath(@__DIR__, "..", "src"))
    isdir(candidate) && return candidate

    error("Set HYDROTRIXI_DOCS_SRC before executing this tutorial.")
end

docs_generated_dir(name) = mkpath(joinpath(docs_src_dir(), "assets", "generated", name))

function quietly(f)
    redirect_stdout(devnull) do
        return f()
    end
end

end
