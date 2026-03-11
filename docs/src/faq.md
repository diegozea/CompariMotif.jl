```@meta
CurrentModule = Main.CompariMotif
DocTestSetup = quote
    using CompariMotif
    using CSV
    using DataFrames
end
```

# FAQ / How-To

## How do I use a Julia `Regex` as input?

The package accepts motif strings, not compiled `Regex` objects. When you already 
have a Julia regex, pass its pattern text with `rx.pattern` (please note that this is not 
part of the public API, so it may change in future versions of Julia). You should also 
ensure that the regex is compatible with the selected alphabet, otherwise an error 
will be thrown.

```jldoctest persist_results
julia> options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> rx = r"A[CG]T";

julia> result = compare(rx.pattern, "ACT", options);

julia> result.matched
true

julia> result.query_relationship
"Degenerate Match"
```

## How do I handle IUPAC ambiguity codes such as `N`, `R`, `B` or `Z`?

CompariMotif accepts motifs written with the selected core alphabet, bracket
classes, and wildcard syntax (`x`, `X`, or `.`). If your input uses IUPAC/IUBMB
ambiguity letters, rewrite them into bracket classes or wildcards before
calling the API.

For nucleotide alphabets, use the NC-IUB ambiguity codes as follows:

| Code | `DNAAlphabet()` | `RNAAlphabet()` |
|:-----|:----------------|:----------------|
| `R` | `[AG]` | `[AG]` |
| `Y` | `[CT]` | `[CU]` |
| `W` | `[AT]` | `[AU]` |
| `S` | `[CG]` | `[CG]` |
| `M` | `[AC]` | `[AC]` |
| `K` | `[GT]` | `[GU]` |
| `H` | `[ACT]` | `[ACU]` |
| `B` | `[CGT]` | `[CGU]` |
| `V` | `[ACG]` | `[ACG]` |
| `D` | `[AGT]` | `[AGU]` |
| `N` | `x` or `.` | `x` or `.` |

For the current 20-residue protein alphabet, the relevant IUPAC-IUB ambiguity
codes are:

| Code | `ProteinAlphabet()` |
|:-----|:--------------------|
| `B` | `[DN]` |
| `Z` | `[EQ]` |

These mappings come from the NC-IUB recommendation
[_Nomenclature for incompletely specified bases in nucleic acid sequences_](https://iubmb.qmul.ac.uk/misc/naseq.html)
and the IUPAC-IUB/JCBN amino-acid nomenclature section
[_AA-21.2 The one-letter system_](https://iupac.qmul.ac.uk/AminoAcid/A2021.html).

You can use Julia's `replace` function to rewrite motifs with IUPAC codes into the 
appropriate format for CompariMotif. For example:

```jldoctest persist_results
julia> clean_dna_motif = replace("ARYN", "R" => "[AG]", "Y" => "[CT]", "N" => "x")
"A[AG][CT]x"

julia> clean_rna_motif = replace("AUHB", "H" => "[ACU]", "B" => "[CGU]")
"AU[ACU][CGU]"

julia> clean_protein_motif = replace("ABZX", "B" => "[DN]", "Z" => "[EQ]")
"A[DN][EQ]X"
```jldoctest persist_results

## How do I switch alphabets?

Choose the alphabet when constructing [`ComparisonOptions`](@ref), or when
calling [`normalize_motif`](@ref) directly. The alphabet options are `DNAAlphabet()`, 
`RNAAlphabet()`, and `ProteinAlphabet()` (the default).

```jldoctest persist_results
julia> options = ComparisonOptions(; alphabet = RNAAlphabet());

julia> compare("AxG", "AUG", options).matched
true
```

## How do I turn results into a table?

Use [`to_column_table`](@ref) when you want to convert a matrix of comparison results into 
a column-oriented `NamedTuple` that can be easily converted to a `DataFrame` or written to 
a CSV file. The `CSV.write` function can consume the `NamedTuple` directly, which is 
convenient for exporting results (you need to do `using CSV` to use `CSV.write`). 
When interactively exploring results in a Julia session, converting the table to a 
`DataFrame` is often easier, as `DataFrame` operations provide convenient tools for 
data manipulation, such as sorting and filtering. In order to use `DataFrame`, you need to 
do `using DataFrames` in your Julia session.

```jldoctest persist_results
julia> motifs = ["RKLI", "R[KR]L[IV]", "[ST]P"];

julia> options = ComparisonOptions(; min_shared_positions = 1, normalized_ic_cutoff = 0.0);

julia> matrix = compare(motifs, options);

julia> table = to_column_table(matrix);

julia> df = DataFrame(table);

julia> show(df, allrows = true, allcols = true)
9×17 DataFrame
 Row │ query_index  search_index  query       search      normalized_query  normalized_search  matched  query_relationship  search_relationship  matched_pattern  matched_positions  match_ic  normalized_ic  core_ic   score    query_information  search_information
     │ Int64        Int64         String      String      String            String             Bool     String              String               String           Int64              Float64   Float64        Float64   Float64  Float64            Float64
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │           1             1  RKLI        RKLI        RKLI              RKLI                  true  Exact Match         Exact Match          RKLI                             4   4.0                1.0  1.0           4.0            4.0                 4.0
   2 │           1             2  RKLI        R[KR]L[IV]  RKLI              R[RK]L[IV]            true  Variant Match       Degenerate Match     R[rk]L[iv]                       4   3.53724            1.0  0.884311      4.0            4.0                 3.53724
   3 │           1             3  RKLI        [ST]P       RKLI              [ST]P                false  No Match            No Match                                              0   0.0                0.0  0.0           0.0            0.0                 0.0
   4 │           2             1  R[KR]L[IV]  RKLI        R[RK]L[IV]        RKLI                  true  Degenerate Match    Variant Match        R[rk]L[iv]                       4   3.53724            1.0  0.884311      4.0            3.53724             4.0
   5 │           2             2  R[KR]L[IV]  R[KR]L[IV]  R[RK]L[IV]        R[RK]L[IV]            true  Exact Match         Exact Match          R[RK]L[IV]                       4   3.53724            1.0  0.884311      4.0            3.53724             3.53724
   6 │           2             3  R[KR]L[IV]  [ST]P       R[RK]L[IV]        [ST]P                false  No Match            No Match                                              0   0.0                0.0  0.0           0.0            0.0                 0.0
   7 │           3             1  [ST]P       RKLI        [ST]P             RKLI                 false  No Match            No Match                                              0   0.0                0.0  0.0           0.0            0.0                 0.0
   8 │           3             2  [ST]P       R[KR]L[IV]  [ST]P             R[RK]L[IV]           false  No Match            No Match                                              0   0.0                0.0  0.0           0.0            0.0                 0.0
   9 │           3             3  [ST]P       [ST]P       [ST]P             [ST]P                 true  Exact Match         Exact Match          [ST]P                            2   1.76862            1.0  0.884311      2.0            1.76862             1.76862
```jldoctest persist_results

## How do I persist comparison results and load them in another Julia session?

For a matrix of results, write `to_column_table(results)` to CSV. The matrix
form includes `query_index` and `search_index`, which lets you reconstruct the
original `Matrix{ComparisonResult}` after loading the file again.

Let's see an example step by step. First, we create some motifs and compare them. Then we 
write the results to a CSV file to persist them. For that, we first need to 
do `using CSV` in our Julia session.

```jldoctest persist_results
julia> motifs = ["RKLI", "R[KR]L[IV]", "[ST]P"];

julia> options = ComparisonOptions();

julia> results = compare(motifs, options);

julia> path = tempname();  # you can choose any filename; here we create a temporary file

julia> CSV.write(path, to_column_table(results));
```

Now we can load the file again in another Julia session and reconstruct the original 
matrix of `ComparisonResult` objects. For that, we need to do `using CSV` to read the 
CSV file and `using DataFrames` to convert the table into a `DataFrame`, which is easier 
to work with:

```jldoctest persist_results
julia> loaded = DataFrame(CSV.File(path));
```

Once we have the `DataFrame`, we can create an empty matrix of `ComparisonResult` objects 
with the correct dimensions. The `query_index` and `search_index` columns indicate the 
required size. We then fill the matrix by iterating over the rows of the `DataFrame` and 
constructing `ComparisonResult` objects from each row:

```jldoctest persist_results
julia> restored = Matrix{ComparisonResult}(undef, maximum(loaded.query_index), maximum(loaded.search_index));

julia> for row in eachrow(loaded)
           restored[row.query_index, row.search_index] = ComparisonResult(;
               query = row.query,
               search = row.search,
               normalized_query = row.normalized_query,
               normalized_search = row.normalized_search,
               matched = row.matched,
               query_relationship = row.query_relationship,
               search_relationship = row.search_relationship,
               matched_pattern = coalesce(row.matched_pattern, ""),
               matched_positions = row.matched_positions,
               match_ic = row.match_ic,
               normalized_ic = row.normalized_ic,
               core_ic = row.core_ic,
               score = row.score,
               query_information = row.query_information,
               search_information = row.search_information,
           )
       end

julia> restored[1, 2]
ComparisonResult(
  query               = RKLI
  search              = R[KR]L[IV]
  normalized_query    = RKLI
  normalized_search   = R[RK]L[IV]
  matched             = true
  query_relationship  = Variant Match
  search_relationship = Degenerate Match
  matched_pattern     = R[rk]L[iv]
  matched_positions   = 4
  match_ic            = 3.537243573680481
  normalized_ic       = 1.0
  core_ic             = 0.8843108934201203
  score               = 4.0
  query_information   = 4.0
  search_information  = 3.537243573680481
)
```
