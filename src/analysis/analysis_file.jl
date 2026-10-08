# Read selected, typed columns from the commented header of a Trixi analysis table.
function read_analysis_columns(path, column_types)
    columns = map(column_types) do (_, value_type)
        Vector{value_type}()
    end

    open(path, "r") do io
        header = nothing
        for line in eachline(io)
            stripped_line = strip(line)
            if isempty(stripped_line)
                continue
            end
            header = stripped_line
            break
        end

        header_columns = split(strip(header[2:end]))
        column_indices = Dict(name => index for (index, name) in pairs(header_columns))
        indices = map(column -> column_indices[first(column)], column_types)

        for line in eachline(io)
            stripped_line = strip(line)
            if isempty(stripped_line) || startswith(stripped_line, "#")
                continue
            end
            values = split(stripped_line)
            for (column, index, (_, value_type)) in zip(columns, indices, column_types)
                push!(column, parse(value_type, values[index]))
            end
        end
    end

    return columns
end
