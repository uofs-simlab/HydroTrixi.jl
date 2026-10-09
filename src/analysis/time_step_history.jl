"""
    accepted_step_history(path::AbstractString; step_column = "timestep",
                          time_column = "time", dt_column = "dt",
                          include_initial = false)

Read accepted-step numbers, end times, and step sizes from an analysis file and return
`(; steps, times, dts)`. The step numbers are integers. By default, omit step zero, whose
step size is the initial guess rather than an accepted step. Set `include_initial = true`
to retain that row when copying a complete analysis table.

An analysis callback with `analysis_interval = 1` records every accepted step.
Larger intervals produce a sampled history;
the recorded step sizes are not differences between successive sample times.
Rejected steps are not recorded.

The default column names match Trixi analysis files. Supply explicit names for other
tables with the same commented-header format:
```julia
history = accepted_step_history("analysis.dat")
history = accepted_step_history("conservation.dat";
                                step_column = "accepted_step", time_column = "time_s",
                                dt_column = "dt_s")
```
"""
function accepted_step_history(path::AbstractString; step_column = "timestep",
                                time_column = "time", dt_column = "dt",
                                include_initial = false)
    steps, times, dts = read_analysis_columns(path, (step_column => Int,
                                                    time_column => Float64,
                                                    dt_column => Float64))
    if include_initial
        return (; steps, times, dts)
    end

    accepted = steps .> 0
    return (; steps = steps[accepted], times = times[accepted], dts = dts[accepted])
end
