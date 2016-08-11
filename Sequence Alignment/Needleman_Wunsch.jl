using OhMyJulia

immutable EvalMetric
    match::Int
    mismatch::Int
    gap::Int
    extend::Int
end

const EDIT_DISTANCE = EvalMetric(0, -1, -1, -1)

function align_basic(seq1::Bytes, seq2::Bytes, eval_metric::EvalMetric)
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

function align_with_overlap_detection(seq1::Bytes, seq2::Bytes, eval_metric::EvalMetric)
    const l1, l2 = length(seq1)+1, length(seq2)+1
    const F      = Array{Int}(l1, l2)
    const Ptr    = Array{Byte}(l1, l2)

    const DIAG, UP, LEFT, STOP = 0x01, 0x02, 0x03, 0x00

    F[:,1] = 0
    F[1,:] = 0

    Ptr[:,1] = STOP
    Ptr[1,:] = STOP

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

    opt1, opt2 = let
        opt1 = findmax(F[:, l2])
        opt2 = findmax(F[l1, :])
        car(opt1) > car(opt2) ? (cadr(opt1), l2) : (l1, cadr(opt2))
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

        traceback(opt1, opt2)

        takebuf_array(a1), takebuf_array(a2)
    end

    a1, a2, F[opt1, opt2]
end

###=== tests ===###
if !isinteractive()

using Base.Test

@assert align_basic(b"AGTA", b"ATA", EDIT_DISTANCE)[3] == -1
@assert align_basic(b"AGGCTATCACCTGACCTCCAGGCCGATGCCC",
                    b"TAGCTATCACGACCGCGGTCGATTTGCCCGAC",
                    EDIT_DISTANCE)[3] == -13

@assert align_with_overlap_detection(b"AGTA", b"ATA", EvalMetric(2, -1, -3, -1))[3] == 3
@assert align_with_overlap_detection(b"AAAACCCCCGGGGTTA",
                                     b"TTCCCGGGAACCAACC",
                                     EvalMetric(2, -1, -3, -1))[3] == 6

end
