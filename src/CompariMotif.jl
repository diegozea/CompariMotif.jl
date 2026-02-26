module CompariMotif

using AutoPrettyPrinting

export ComparisonOptions, ComparisonResult, MatchFixMode, MatchFixBothFixed,
       MatchFixNone, MatchFixQueryFixed, MatchFixSearchFixed, compare,
       normalize_motif, write_results_tsv

# Module assembly order matters:
# - foundational constants/types first,
# - then option parsing and motif parsing helpers,
# - then relationship/alignment core logic,
# - public API methods and output writers last.
# This keeps includes acyclic and makes data flow easier to follow.
include("constants.jl")
include("types.jl")
include("options.jl")
include("parsing.jl")
include("relationships.jl")
include("alignment.jl")
include("api.jl")
include("io.jl")

end
