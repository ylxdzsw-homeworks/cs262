using OhMyJulia

immutable EvalMetric
    match::Int
    mismatch::Int
    gap::Int
end

typealias AbstractBytes AbstractArray{Byte, 1}

const EDIT_DISTANCE = EvalMetric(0, -1, -1)

function score_lastline(seq1::AbstractBytes, seq2::AbstractBytes, eval_metric::EvalMetric)
    l1, l2 = length(seq1)+1, length(seq2)+1
    prev, curr = Vector{Int}(l2), Vector{Int}(l2)

    curr[:] = eval_metric.gap * (0:l2-1)
    prev[1] = 0

    s(x,y) = seq1[x-1]==seq2[y-1] ? eval_metric.match : eval_metric.mismatch

    for i in 2:l1
        prev, curr = curr, prev
        curr[1] = prev[1] + eval_metric.gap

        for j in 2:l2
            curr[j] = max(
                prev[j-1] + s(i, j),
                prev[j]   + eval_metric.gap,
                curr[j-1] + eval_metric.gap
            )
        end
    end

    curr
end

function align_rec(seq1::AbstractBytes, seq2::AbstractBytes, eval_metric::EvalMetric)
    l1, l2 = length(seq1), length(seq2)

    if l1 < 2 || l2 < 2
        return align_basic(seq1, seq2, eval_metric)
    end

    xmid = l1 ÷ 2

    rev(x::AbstractVector) = sub(x, endof(x):-1:1)
    rev(x::AbstractVector, y) = sub(x, reverse(y))

    lscore = score_lastline(sub(seq1, 1:xmid)            , seq2     , eval_metric)
    rscore = score_lastline(rev(seq1, xmid+1:endof(seq1)), rev(seq2), eval_metric)

    ymid = findmax(lscore .+ rev(rscore)) |> cadr
    ymid -= 1 # adjust index

    la1, la2, ls = align_rec(sub(seq1, 1:xmid)            , sub(seq2, 1:ymid)            , eval_metric)
    ra1, ra2, rs = align_rec(sub(seq1, xmid+1:endof(seq1)), sub(seq2, ymid+1:endof(seq2)), eval_metric)

    la1++ra1, la2++ra2, ls+rs
end

### copied from ../Sequence Alignment ###
function align_basic(seq1::AbstractBytes, seq2::AbstractBytes, eval_metric::EvalMetric)
    const l1, l2 = length(seq1)+1, length(seq2)+1
    const F      = Array{Int}(l1, l2)
    const Ptr    = Array{Byte}(l1, l2)

    const DIAG, UP, LEFT, STOP = 0x01, 0x02, 0x03, 0x00

    F[:,1] = eval_metric.gap * (0:l1-1)
    F[1,:] = eval_metric.gap * (0:l2-1)

    Ptr[:,1] = UP
    Ptr[1,:] = LEFT
    Ptr[1,1] = STOP

    s(x,y) = seq1[x-1]==seq2[y-1] ? eval_metric.match : eval_metric.mismatch

    for i in 2:l1, j in 2:l2
        score, case = findmax([
            F[i-1, j-1] + s(i, j),
            F[i-1, j] + eval_metric.gap,
            F[i, j-1] + eval_metric.gap
        ])

        F[i,j]   = score
        Ptr[i,j] = case
    end

    a1, a2 = let
        a1, a2 = IOBuffer(), IOBuffer()

        function traceback(i, j)
            ptr = Ptr[i, j]

            if ptr == DIAG
                traceback(i-1, j-1)
                a1 << seq1[i-1]
                a2 << seq2[j-1]
            elseif ptr == UP
                traceback(i-1, j)
                a1 << seq1[i-1]
                a2 << '_'
            elseif ptr == LEFT
                traceback(i, j-1)
                a1 << '_'
                a2 << seq2[j-1]
            end

            nothing
        end

        traceback(l1, l2)

        takebuf_array(a1), takebuf_array(a2)
    end

    a1, a2, F[l1, l2]
end

###=== tests ===###
if !isinteractive()

using Base.Test

@assert score_lastline(b"AGTACGCA", b"TATGC", EvalMetric(2, -1, -2)) == [-16, -12, -8, -7, -3, 1]

@assert align_rec(b"AGGCTATCACCTGACCTCCAGGCCGATGCCC",
                  b"TAGCTATCACGACCGCGGTCGATTTGCCCGAC",
                  EDIT_DISTANCE)[3] == -13
end
