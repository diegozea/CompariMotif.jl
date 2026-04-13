const ORACLE_RELATIVE_PATH = joinpath("tools", "comparimotif_V3.py")
const ORACLE_NORMALIZED_COLUMNS = [
    "Motif1", "Motif2", "Sim1", "Sim2", "Match",
    "MatchPos", "MatchIC", "NormIC", "CoreIC", "Score", "Info1", "Info2"
]
const ALTERNATION_NORMALIZED_COLUMNS = [
    "Name1", "Name2", "Motif1", "Motif2", "Sim1", "Sim2", "Match",
    "MatchPos", "MatchIC", "NormIC", "CoreIC", "Score", "Info1", "Info2"
]
const CORNERCASE_NORMALIZED_COLUMNS = ALTERNATION_NORMALIZED_COLUMNS

function read_motifs(path::String)
    motifs = String[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        startswith(line, '#') && continue
        cols = split(line, '\t')
        length(cols) == 1 || error("Invalid motif line: $line")
        push!(motifs, cols[1])
    end
    return motifs
end

function write_motif_file(path::String, motifs::Vector{String})
    open(path, "w") do io
        for (index, motif) in enumerate(motifs)
            name = "M$(lpad(string(index), 4, '0'))"
            println(io, name, '\t', motif)
        end
    end
end

function read_table(path::String)
    lines = readlines(path)
    isempty(lines) && error("Empty oracle output: $path")
    header = split(lines[1], '\t')
    rows = Dict{String, String}[]
    for line in Iterators.drop(lines, 1)
        stripped = strip(line)
        isempty(stripped) && continue
        cols = split(stripped, '\t')
        length(cols) == length(header) || error("Malformed oracle output row: $line")
        row = Dict{String, String}()
        for (key, value) in zip(header, cols)
            row[key] = value
        end
        push!(rows, row)
    end
    return header, rows
end

function normalize_rows(
        rows::Vector{Dict{String, String}},
        columns::Vector{String} = ORACLE_NORMALIZED_COLUMNS
)
    for row in rows
        for col in columns
            haskey(row, col) || error("Oracle output missing expected column: $col")
        end
    end
    sort!(rows; by = row -> Tuple(row[col] for col in columns))
    return rows
end

function resolve_oracle_path(env::AbstractDict = ENV)
    slimsuite_root = strip(get(env, "SLiMSuite_PATH", ""))
    if !isempty(slimsuite_root)
        return joinpath(slimsuite_root, ORACLE_RELATIVE_PATH)
    end

    error(
        "Oracle executable not configured. Set `SLiMSuite_PATH` to a local clone of " *
        "https://github.com/slimsuite/SLiMSuite " *
        "(expected executable at `\$SLiMSuite_PATH/$ORACLE_RELATIVE_PATH`)."
    )
end

function _cleanup_oracle_logs!(dir::String)
    for entry in readdir(dir; join = true)
        isfile(entry) || continue
        filename = basename(entry)
        if occursin(r"^comparimotif.*\.log$", filename)
            rm(entry; force = true)
        end
    end
end
