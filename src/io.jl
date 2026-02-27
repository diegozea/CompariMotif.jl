"""
    _empty_result_columns(nrows::Int) -> NamedTuple

Allocate typed result columns for `nrows` rows.
"""
function _empty_result_columns(nrows::Int)
    (
        query = Vector{String}(undef, nrows),
        search = Vector{String}(undef, nrows),
        normalized_query = Vector{String}(undef, nrows),
        normalized_search = Vector{String}(undef, nrows),
        matched = Vector{Bool}(undef, nrows),
        query_relationship = Vector{String}(undef, nrows),
        search_relationship = Vector{String}(undef, nrows),
        matched_pattern = Vector{String}(undef, nrows),
        matched_positions = Vector{Int}(undef, nrows),
        match_ic = Vector{Float64}(undef, nrows),
        normalized_ic = Vector{Float64}(undef, nrows),
        core_ic = Vector{Float64}(undef, nrows),
        score = Vector{Float64}(undef, nrows),
        query_information = Vector{Float64}(undef, nrows),
        search_information = Vector{Float64}(undef, nrows)
    )
end

"""
    _set_result_row!(columns, row::Int, result::ComparisonResult)

Write one `ComparisonResult` into preallocated column vectors.
"""
function _set_result_row!(columns::NamedTuple, row::Int, result::ComparisonResult)
    columns.query[row] = result.query
    columns.search[row] = result.search
    columns.normalized_query[row] = result.normalized_query
    columns.normalized_search[row] = result.normalized_search
    columns.matched[row] = result.matched
    columns.query_relationship[row] = result.query_relationship
    columns.search_relationship[row] = result.search_relationship
    columns.matched_pattern[row] = result.matched_pattern
    columns.matched_positions[row] = result.matched_positions
    columns.match_ic[row] = result.match_ic
    columns.normalized_ic[row] = result.normalized_ic
    columns.core_ic[row] = result.core_ic
    columns.score[row] = result.score
    columns.query_information[row] = result.query_information
    columns.search_information[row] = result.search_information
    return nothing
end

"""
    to_column_table(results) -> NamedTuple

Convert comparison results into a column-oriented `NamedTuple` where each key is
a column name and each value is a vector column.

- `to_column_table(::ComparisonResult)` returns a one-row table.
- `to_column_table(::AbstractVector{<:ComparisonResult})` adds `result_index`.
- `to_column_table(::AbstractMatrix{<:ComparisonResult})` adds `query_index` and
  `search_index` with one row per matrix cell in deterministic row-major order.

The returned object can be converted to a `DataFrame` or written using `CSV.write`
without requiring either dependency in the package itself.

# Examples
```jldoctest
julia> using CompariMotif, DataFrames

julia> motifs = ["RKLI", "R[KR]L[IV]"];

julia> options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> table = to_column_table(compare(motifs, options));

julia> df = DataFrame(table);

julia> show(select(df, [:query_index, :search_index, :query, :search, :query_relationship]), allrows = true, allcols = true, truncate = 0)
4×5 DataFrame
 Row │ query_index  search_index  query       search      query_relationship
     │ Int64        Int64         String      String      String
─────┼───────────────────────────────────────────────────────────────────────
   1 │           1             1  RKLI        RKLI        Exact Match
   2 │           1             2  RKLI        R[KR]L[IV]  Variant Match
   3 │           2             1  R[KR]L[IV]  RKLI        Degenerate Match
   4 │           2             2  R[KR]L[IV]  R[KR]L[IV]  Exact Match
```

$(_DOC_COMPARE_REF)
$(_DOC_RESULT_REF)
"""
function to_column_table(result::ComparisonResult)
    columns = _empty_result_columns(1)
    _set_result_row!(columns, 1, result)
    return columns
end

"""
    to_column_table(results::AbstractVector{<:ComparisonResult}) -> NamedTuple

Convert a result vector to a column table with `result_index`.
"""
function to_column_table(results::AbstractVector{<:ComparisonResult})
    nrows = length(results)
    columns = merge((result_index = Vector{Int}(undef, nrows),), _empty_result_columns(nrows))
    for (row, result) in enumerate(results)
        columns.result_index[row] = row
        _set_result_row!(columns, row, result)
    end
    return columns
end

"""
    to_column_table(results::AbstractMatrix{<:ComparisonResult}) -> NamedTuple

Convert a result matrix to a column table with `query_index` and `search_index`.
"""
function to_column_table(results::AbstractMatrix{<:ComparisonResult})
    nrows = length(results)
    columns = merge(
        (
            query_index = Vector{Int}(undef, nrows),
            search_index = Vector{Int}(undef, nrows)
        ),
        _empty_result_columns(nrows)
    )

    row = 0
    for i in axes(results, 1)
        for j in axes(results, 2)
            row += 1
            columns.query_index[row] = i
            columns.search_index[row] = j
            _set_result_row!(columns, row, results[i, j])
        end
    end
    return columns
end
