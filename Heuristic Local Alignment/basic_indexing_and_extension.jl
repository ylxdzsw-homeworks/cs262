using OhMyJulia

immutable EvalMetric
    match::Int
    mismatch::Int
    gap::Int
end

function align_index_and_extend(query::Bytes, ref::Bytes, eval_metric::EvalMetric; index_length::Int=4)
    qlen, rlen = length(query), length(ref)

    index = let
        I = Dict{Bytes, Vector{Int}}()

        for s in 1:rlen-index_length+1
            word = ref[s:s+index_length-1]

            if word in keys(I)
                push!(I[word], s)
            else
                I[word] = [s]
            end
        end

        I
    end

    seeds = @task for s in 1:qlen-index_length+1
        word = query[s:s+index_length-1]

        for i in get(index, word, Int[])
            produce((s, i))
        end
    end
end

###=== tests ===###
if !isinteractive()

using Base.Test

@assert align_index_and_extend(b"GTAAGGTCCAGTAA",
                               b"GTTAGGTCAGTCA",
                               EvalMetric(1,-1,-1))

end
