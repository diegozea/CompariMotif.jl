abstract type _Alphabet end

struct ProteinAlphabet <: _Alphabet end

struct DNAAlphabet <: _Alphabet end

struct RNAAlphabet <: _Alphabet end

const _AlphabetValue = Union{ProteinAlphabet, DNAAlphabet, RNAAlphabet}

struct _AlphabetSpec
    chars::String
    index::Dict{Char, Int}
    mask::ResidueMask
    log_base::Float64
end

_alphabet_chars(::Type{ProteinAlphabet}) = _ProteinAlphabetString
_alphabet_chars(::Type{DNAAlphabet}) = _DNAAlphabetString
_alphabet_chars(::Type{RNAAlphabet}) = _RNAAlphabetString

@doc """
    ProteinAlphabet

Select the standard protein alphabet for [`ComparisonOptions`](@ref).
Allowed residues: `$(_ProteinAlphabetString)`. Use as `ProteinAlphabet()`.
""" ProteinAlphabet

@doc """
    DNAAlphabet

Select the DNA alphabet for [`ComparisonOptions`](@ref).
Allowed residues: `$(_DNAAlphabetString)`. Use as `DNAAlphabet()`.
""" DNAAlphabet

@doc """
    RNAAlphabet

Select the RNA alphabet for [`ComparisonOptions`](@ref).
Allowed residues: `$(_RNAAlphabetString)`. Use as `RNAAlphabet()`.
""" RNAAlphabet

function _build_alphabet_spec(::Type{A}) where {A <: _Alphabet}
    chars = _alphabet_chars(A)
    index = Dict{Char, Int}()
    for (i, aa) in enumerate(chars)
        index[aa] = i
    end
    mask = ResidueMask((1 << length(chars)) - 1)
    return _AlphabetSpec(chars, index, mask, log(length(chars)))
end

const _PROTEIN_SPEC = _build_alphabet_spec(ProteinAlphabet)
const _DNA_SPEC = _build_alphabet_spec(DNAAlphabet)
const _RNA_SPEC = _build_alphabet_spec(RNAAlphabet)

_alphabet_spec(::Type{ProteinAlphabet}) = _PROTEIN_SPEC
_alphabet_spec(::Type{DNAAlphabet}) = _DNA_SPEC
_alphabet_spec(::Type{RNAAlphabet}) = _RNA_SPEC

_alphabet_spec(::ProteinAlphabet) = _PROTEIN_SPEC
_alphabet_spec(::DNAAlphabet) = _DNA_SPEC
_alphabet_spec(::RNAAlphabet) = _RNA_SPEC
