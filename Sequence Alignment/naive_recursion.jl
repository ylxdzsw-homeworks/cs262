using OhMyJulia

immutable EvalMetric
    mismatch::Int
    gap::Int
end

const EDIT_DISTANCE = EvalMetric(1, 1)

function align_naive(seq1::Bytes, seq2::Bytes, eval_metric::EvalMetric)
    _ = Byte('_')

    function align(s1::Bytes, s2::Bytes)
        if isempty(s1)
            map(x->_, s2), s2, eval_metric.gap * length(s2)
        elseif isempty(s2)
            s1, map(x->_, s1), eval_metric.gap * length(s1)
        else
            min(
                begin # match or replace
                    a1, a2, cost = align(s1[1:end-1], s2[1:end-1])
                    e1, e2 = s1[end], s2[end]
                    a1++e1, a2++e2, cost+(e1!=e2)*eval_metric.mismatch
                end,
                begin # gap on s1
                    a1, a2, cost = align(s1[1:end], s2[1:end-1])
                    a1++_, a2++s2[end], cost+eval_metric.gap
                end,
                begin # gap on s2
                    a1, a2, cost = align(s1[1:end-1], s2[1:end])
                    a1++s1[end], a2++_, cost+eval_metric.gap
                end,
                key=x->x[3]
            )
        end
    end

    align(seq1, seq2)
end

function align_memo(seq1::Bytes, seq2::Bytes, eval_metric::EvalMetric)
    memo = Dict{Tuple{Bytes, Bytes}, Tuple{Bytes, Bytes, Int}}()
    _ = Byte('_')

    function _align(s1::Bytes, s2::Bytes)
        if isempty(s1)
            map(x->_, s2), s2, eval_metric.gap * length(s2)
        elseif isempty(s2)
            s1, map(x->_, s1), eval_metric.gap * length(s1)
        else
            min(
                begin # match or replace
                    a1, a2, cost = align(s1[1:end-1], s2[1:end-1])
                    e1, e2 = s1[end], s2[end]
                    a1++e1, a2++e2, cost+(e1!=e2)*eval_metric.mismatch
                end,
                begin # gap on s1
                    a1, a2, cost = align(s1[1:end], s2[1:end-1])
                    a1++_, a2++s2[end], cost+eval_metric.gap
                end,
                begin # gap on s2
                    a1, a2, cost = align(s1[1:end-1], s2[1:end])
                    a1++s1[end], a2++_, cost+eval_metric.gap
                end,
                key=x->x[3]
            )
        end
    end

    function align(s1::Bytes, s2::Bytes)
        if (s1, s2) in keys(memo)
            memo[(s1, s2)]
        else
            memo[(s1, s2)] = _align(s1, s2)
        end
    end

    align(seq1, seq2)
end

###=== tests ===###
if !isinteractive()

using Base.Test

@assert align_naive(b"AGTA", b"ATA", EDIT_DISTANCE)[3] == 1
# @assert align_naive(b"AGGCTATCACCTGACCTCCAGGCCGATGCCC",
#                     b"TAGCTATCACGACCGCGGTCGATTTGCCCGAC",
#                     EDIT_DISTANCE)[3] == 13

@assert align_memo(b"AGTA", b"ATA", EDIT_DISTANCE)[3] == 1
@assert align_memo(b"AGGCTATCACCTGACCTCCAGGCCGATGCCC",
                   b"TAGCTATCACGACCGCGGTCGATTTGCCCGAC",
                   EDIT_DISTANCE)[3] == 13

end
