using TestItems

@testitem "FAQ docs examples" begin
    using Test
    using CompariMotif
    using Graphs

    # This test mirrors the ELM clustering example in `docs/src/faq.md`.
    # Keep this file and that section in sync when changing the example
    # motifs, thresholds, graph construction, or documented outputs.

    is_cluster_match(result) = result.matched && result.match_ic >= 1.5

    function connected_component_assignments(results)
        n = size(results, 1)
        size(results, 2) == n ||
            throw(ArgumentError("expected a square all-vs-all comparison matrix"))

        graph = SimpleDiGraph(n)
        for i in 1:n, j in 1:n

            i == j && continue
            is_cluster_match(results[i, j]) || continue
            add_edge!(graph, i, j)
        end

        components = weakly_connected_components(graph)
        sort!(components, by = minimum)

        clusters = zeros(Int, n)
        for (cluster_id, component) in enumerate(components)
            clusters[component] .= cluster_id
        end

        return clusters
    end

    @testset "ELM clustering example" begin
        motifs = ["RKLI", "R[KR]L[IV]", "[ST]P", "[ST]Px[KR]"]
        options = ComparisonOptions()
        results = compare(motifs, options)
        clusters = connected_component_assignments(results)

        @test clusters == [1, 1, 2, 2]
        @test maximum(clusters) == 2
    end

    @testset "directed graph rationale" begin
        options = ComparisonOptions(;
            matchfix = :query_fixed,
            min_shared_positions = 1,
            normalized_ic_cutoff = 0.0
        )

        ab = compare("A", "[AC]", options)
        ba = compare("[AC]", "A", options)

        @test !ab.matched
        @test ba.matched
        @test ab.match_ic != ba.match_ic
    end
end
