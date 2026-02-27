"""
    _relationship_type_from_flags(has_variant, has_degenerate, has_complex) -> _RelationshipType

Reduce position-level relationship evidence into a final relationship type.
"""
# Reduce per-position relationship evidence into the final relationship type.
function _relationship_type_from_flags(has_variant::Bool, has_degenerate::Bool, has_complex::Bool)
    if has_complex || (has_variant && has_degenerate)
        return _REL_COMPLEX
    elseif has_variant
        return _REL_VARIANT
    elseif has_degenerate
        return _REL_DEGENERATE
    else
        return _REL_EXACT
    end
end

"""
    _reverse_type(relationship_type::_RelationshipType) -> _RelationshipType

Reverse relation type for query-vs-search direction inversion.
"""
# Reverse relation type for query-vs-search direction inversion.
# Variant and Degenerate swap perspective; others are symmetric.
function _reverse_type(relationship_type::_RelationshipType)
    if relationship_type == _REL_VARIANT
        return _REL_DEGENERATE
    elseif relationship_type == _REL_DEGENERATE
        return _REL_VARIANT
    end
    return relationship_type
end

"""
    _relationship_length(qlen, slen, overlap) -> _RelationshipLength

Classify relationship length based on overlap coverage.
"""
# Length class is based on how much of each motif is covered by overlap.
function _relationship_length(qlen::Int, slen::Int, overlap::Int)
    if overlap == qlen && qlen == slen
        return _LEN_MATCH
    elseif overlap == slen && qlen > slen
        return _LEN_PARENT
    elseif overlap == qlen && slen > qlen
        return _LEN_SUBSEQUENCE
    else
        return _LEN_OVERLAP
    end
end

"""
    _reverse_length(length_type::_RelationshipLength) -> _RelationshipLength

Reverse length class for the opposite comparison direction.
"""
# Reverse length relation for the opposite perspective.
function _reverse_length(length_type::_RelationshipLength)
    if length_type == _LEN_PARENT
        return _LEN_SUBSEQUENCE
    elseif length_type == _LEN_SUBSEQUENCE
        return _LEN_PARENT
    end
    return length_type
end

"""
    _relationship_word(relationship_type, length_type) -> String

Render full relationship words for user-facing result fields.
"""
# Render full relationship words for user-facing result fields.
function _relationship_word(relationship_type::_RelationshipType, length_type::_RelationshipLength)
    return _RELATIONSHIP_TYPE_WORDS[Int(relationship_type) + 1] * " " *
           _RELATIONSHIP_LENGTH_WORDS[Int(length_type) + 1]
end
