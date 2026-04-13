#!/usr/bin/env julia

import Dates
import Dates: @dateformat_str

const FIXTURES_DIR = @__DIR__
include(joinpath(FIXTURES_DIR, "oracle_utils.jl"))

# Backward-compatible defaults used by tests that include this script.
const MOTIF_SET = joinpath(FIXTURES_DIR, "regression_motifs.tsv")
const NORMALIZED_OUTPUT = joinpath(FIXTURES_DIR, "oracle_regression_normalized.tsv")

const FIXTURE_SPECS = (
    (
        name = "regression",
        motif_set = joinpath(FIXTURES_DIR, "regression_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_regression_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "defaults",
        motif_set = joinpath(FIXTURES_DIR, "default_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_default_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=2", "normcut=0.5", "matchfix=0",
            "mismatches=0", "overlaps=T", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "boundary_wildcard",
        motif_set = joinpath(FIXTURES_DIR, "boundary_wildcard_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_boundary_wildcard_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "query_fixed",
        motif_set = joinpath(FIXTURES_DIR, "matchfix_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_query_fixed_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "matchfix=1", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "search_fixed",
        motif_set = joinpath(FIXTURES_DIR, "matchfix_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_search_fixed_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "matchfix=2", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "both_fixed",
        motif_set = joinpath(FIXTURES_DIR, "matchfix_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_both_fixed_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "matchfix=3", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "coreic",
        motif_set = joinpath(FIXTURES_DIR, "coreic_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_coreic_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "mismatches=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "coreic_mismatch",
        motif_set = joinpath(FIXTURES_DIR, "coreic_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_coreic_mismatch_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "mismatches=1", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "alternation",
        motif_set = joinpath(FIXTURES_DIR, "alternation_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_alternation_probe_normalized.tsv"),
        columns = ALTERNATION_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "wildcard_alias_divergence",
        motif_set = joinpath(FIXTURES_DIR, "wildcard_alias_divergence_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_wildcard_alias_divergence_normalized.tsv"),
        columns = ALTERNATION_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "exact_prefilter",
        motif_set = joinpath(FIXTURES_DIR, "exact_prefilter_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_exact_prefilter_probe_normalized.tsv"),
        columns = ALTERNATION_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "score_tiebreak",
        motif_set = joinpath(FIXTURES_DIR, "score_tiebreak_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_score_tiebreak_probe_normalized.tsv"),
        columns = ALTERNATION_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "variant_order",
        motif_set = joinpath(FIXTURES_DIR, "variant_order_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_variant_order_probe_normalized.tsv"),
        columns = ALTERNATION_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "relationship_tiebreak",
        motif_set = joinpath(FIXTURES_DIR, "relationship_tiebreak_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_relationship_tiebreak_probe_normalized.tsv"),
        columns = ALTERNATION_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "relationship_overlap_tiebreak",
        motif_set = joinpath(FIXTURES_DIR, "relationship_overlap_tiebreak_probe_motifs.tsv"),
        output = joinpath(
            FIXTURES_DIR,
            "oracle_relationship_overlap_tiebreak_probe_normalized.tsv"
        ),
        columns = ALTERNATION_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    ),
    (
        name = "nonuniform",
        motif_set = joinpath(FIXTURES_DIR, "nonuniform_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_nonuniform_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = joinpath(FIXTURES_DIR, "nonuniform_probe.aafreq.tsv"),
        dna = false
    ),
    (
        name = "dna",
        motif_set = joinpath(FIXTURES_DIR, "dna_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_dna_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = true
    ),
    (
        name = "dna_nonuniform",
        motif_set = joinpath(FIXTURES_DIR, "dna_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_dna_nonuniform_probe_normalized.tsv"),
        columns = ORACLE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = joinpath(FIXTURES_DIR, "dna_nonuniform_probe.aafreq.tsv"),
        dna = true
    ),
    (
        name = "cornercases",
        motif_set = joinpath(FIXTURES_DIR, "cornercase_probe_motifs.tsv"),
        output = joinpath(FIXTURES_DIR, "oracle_cornercase_probe_normalized.tsv"),
        columns = CORNERCASE_NORMALIZED_COLUMNS,
        oracle_args = ["minshare=1", "normcut=0", "xgmml=F"],
        aafreq = nothing,
        dna = false
    )
)

function write_normalized(path::String, rows::Vector{Dict{String, String}};
        columns::Vector{String} = ORACLE_NORMALIZED_COLUMNS)
    open(path, "w") do io
        println(io, "# Generated by data/fixtures/generate_oracle_fixture.jl")
        println(io, "# Date: ", Dates.format(Dates.now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "# Source: CompariMotif oracle (black-box)")
        println(io, join(columns, '\t'))
        for row in rows
            println(io, join((row[col] for col in columns), '\t'))
        end
    end
end

function _fixture_spec_map()
    return Dict(spec.name => spec for spec in FIXTURE_SPECS)
end

function _selected_specs(args::Vector{String})
    specs = _fixture_spec_map()
    if isempty(args) || any(arg -> arg == "all", args)
        return collect(FIXTURE_SPECS)
    end

    selected = NamedTuple[]
    for arg in args
        haskey(specs, arg) || error(
            "Unknown fixture set '$arg'. Allowed: " *
            join(sort!(collect(keys(specs))), ", ") * ", all"
        )
        push!(selected, specs[arg])
    end
    return selected
end

function _run_fixture!(oracle::String, spec)
    motifs = read_motifs(spec.motif_set)
    isempty(motifs) && error("No motifs found in fixture set: $(spec.motif_set)")

    mktempdir() do tmpdir
        motif_file = joinpath(tmpdir, "motifs.tsv")
        result_prefix = joinpath(tmpdir, "oracle")
        write_motif_file(motif_file, motifs)

        cmd_parts = [
            "python3", oracle, "motifs=$motif_file", "searchdb=$motif_file",
            "resfile=$result_prefix", "i=-1", "v=0", spec.oracle_args...
        ]
        if spec.dna
            push!(cmd_parts, "dna=T")
        end
        if spec.aafreq !== nothing
            push!(cmd_parts, "aafreq=$(spec.aafreq)")
        end
        cmd = Cmd(cmd_parts)
        @info "Running oracle fixture" fixture=spec.name cmd
        cd(tmpdir) do
            run(cmd)
            _cleanup_oracle_logs!(tmpdir)
        end

        raw_table = "$(result_prefix).compare.tdt"
        isfile(raw_table) || error("Oracle output not found: $raw_table")

        mkpath(FIXTURES_DIR)
        _, rows = read_table(raw_table)
        normalize_rows(rows, spec.columns)
        write_normalized(spec.output, rows; columns = spec.columns)
    end

    @info "Wrote fixtures" fixture=spec.name normalized=spec.output
end

function main(args::Vector{String} = ARGS)
    oracle = resolve_oracle_path()
    isfile(oracle) || error("Oracle executable not found: $oracle")
    _cleanup_oracle_logs!(pwd())

    for spec in _selected_specs(args)
        _run_fixture!(oracle, spec)
    end

    _cleanup_oracle_logs!(pwd())
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
