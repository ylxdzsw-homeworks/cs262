#!/usr/bin/env julia
using OhMyJulia
using RedisAlchemy

conn = RedisConnection(db=2)

chomp!(s::Bytes) = ccall(:jl_array_del_end, Void, (Any, UInt), s, 1)

let R
    for i in eachline("D:/hg19/hg19.fa")
        chomp!(i.data)

        if startswith(i, '>')
            R = RedisBlob(conn, cdr(i))
        else
            R += i
        end
    end
end
