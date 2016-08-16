using OhMyJulia

immutable EvalMetric
    match::Int
    mismatch::Int
    gap::Int
    extend::Int
end

function align_with_affine_gap(seq1::Bytes, seq2::Bytes, eval_metric::EvalMetric)
    const l1, l2 = length(seq1)+1, length(seq2)+1
    const F      = Array{Int}(l1, l2)
    const G      = Array{Int}(l1, l2)
    const H      = Array{Int}(l1, l2)
    const Ptr    = Array{Byte}(l1, l2)

    const FMASK = 0b0011
    const STOP  = 0b0000
    const DIAG  = 0b0001
    const UP    = 0b0010
    const LEFT  = 0b0011

    const GMASK = 0b0100
    const GF    = 0b0000
    const GG    = 0b0100

    const HMASK = 0b1000
    const HF    = 0b0000
    const HH    = 0b1000

    F[2:end,1]   = eval_metric.gap + eval_metric.extend * (0:l1-2)
    G[2:end,1]   = eval_metric.gap + eval_metric.extend * (0:l1-2)
    H[2:end,1]   = 2eval_metric.gap + eval_metric.extend * (-1:l1-3) - 1
    Ptr[2,1]     = UP | GF | HF
    Ptr[3:end,1] = UP | GG | HF

    F[1,2:end]   = eval_metric.gap + eval_metric.extend * (0:l2-2)
    H[1,2:end]   = eval_metric.gap + eval_metric.extend * (0:l2-2)
    G[1,2:end]   = 2eval_metric.gap + eval_metric.extend * (-1:l2-3) - 1
    Ptr[1,2]     = LEFT | GF | HF
    Ptr[1,3:end] = LEFT | GF | HH

    F[1,1]   = 0
    # G[1,1] and H[1,1] are undefined and shouldn't be accessed
    Ptr[1,1] = STOP

    s(x,y) = seq1[x-1]==seq2[y-1] ? eval_metric.match : eval_metric.mismatch

    for i in 2:l1, j in 2:l2
        Gscore, Gcase = findmax((
            F[i-1, j] + eval_metric.gap,
            G[i-1, j] + eval_metric.extend
        ))

        Hscore, Hcase = findmax((
            F[i, j-1] + eval_metric.gap,
            H[i, j-1] + eval_metric.extend
        ))

        Fscore, Fcase = findmax((
            F[i-1, j-1] + s(i, j),
            Gscore,
            Hscore
        ))

        F[i,j] = Fscore
        G[i,j] = Gscore
        H[i,j] = Hscore

        Ptr[i,j] = Fcase % Byte | (GF, GG)[Gcase] | (HF, HH)[Hcase]
    end

    a1, a2 = let
        a1, a2 = IOBuffer(), IOBuffer()

        function traceback(i, j, mask)
            ptr = Ptr[i, j] & mask

            if mask == FMASK
                if ptr == DIAG
                    traceback(i-1, j-1, FMASK)
                    a1 << seq1[i-1]
                    a2 << seq2[j-1]
                elseif ptr == UP
                    traceback(i, j, GMASK)
                elseif ptr == LEFT
                    traceback(i, j, HMASK)
                end
            elseif mask == GMASK
                ptr == GF ? traceback(i-1, j, FMASK) :
                ptr == GG ? traceback(i-1, j, GMASK) :
                            error("unknown ptr")
                a1 << seq1[i-1]
                a2 << '_'
            elseif mask == HMASK
                ptr == HF ? traceback(i, j-1, FMASK) :
                ptr == HH ? traceback(i, j-1, HMASK) :
                            error("unknown ptr")
                a1 << '_'
                a2 << seq2[j-1]
            end

            nothing
        end

        traceback(l1, l2, FMASK)

        takebuf_array(a1), takebuf_array(a2)
    end

    a1, a2, F[l1, l2]
end


###=== tests ===###
if !isinteractive()

using Base.Test

@assert align_with_affine_gap(b"AGTA", b"ATA", EvalMetric(2, -2, -3, -1))[3] == 3
@assert align_with_affine_gap(b"AAAACCCCCGGGGTTA",
                              b"TTCCCGGGAACCAACC",
                              EvalMetric(2, -2, -3, -1))[3] == -10
@assert align_with_affine_gap(b"AGGCTATCACCTGACCTCCAGGCCGATGCCC",
                              b"TAGCTATCACGACCGCGGTCGATTTGCCCGAC",
                              EvalMetric(2, -2, -3, -1))[3] == 21
end
