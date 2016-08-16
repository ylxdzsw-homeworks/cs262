using OhMyJulia

immutable EvalMetric
    match::Int
    mismatch::Int
    gap::Int
end

function align_local(seq1::Bytes, seq2::Bytes, eval_metric::EvalMetric)
    const l1, l2 = length(seq1)+1, length(seq2)+1
    const F      = Array{Int}(l1, l2)
    const Ptr    = Array{Byte}(l1, l2)

    const DIAG, UP, LEFT, STOP = 0x01, 0x02, 0x03, 0x00

    F[:,1] = 0
    F[1,:] = 0

    Ptr[:,1] = UP
    Ptr[1,:] = LEFT
    Ptr[1,1] = STOP

    s(x,y) = seq1[x-1]==seq2[y-1] ? eval_metric.match : eval_metric.mismatch

    for i in 2:l1, j in 2:l2
        score, case = findmax([
            0,
            F[i-1, j-1] + s(i, j),
            F[i-1, j] + eval_metric.gap,
            F[i, j-1] + eval_metric.gap,
        ])

        F[i,j]   = score
        Ptr[i,j] = case - 1
    end

    opt = findmax(F)

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

        traceback(ind2sub(size(F), cadr(opt))...)

        takebuf_array(a1), takebuf_array(a2)
    end

    a1, a2, car(opt)
end

###=== tests ===###
if !isinteractive()

using Base.Test

@assert align_local(b"AGTA", b"ATA", EvalMetric(1, -1, -1))[3] == 2
@assert align_local(b"AAAACCCCCGGGGTTA",
                    b"TTCCCGGGAACCAACC",
                    EvalMetric(1, -1, -1))[3] == 6

end
