##
## METIS
## Minimum interface to the Metis graph partitioner
##

const METIS_NOPTIONS = 40
const METIS_OK = Int32(1)
const metis_options = - ones(Cint, METIS_NOPTIONS)
const libmetis = "/Users/lcambier/.julia/packages/Homebrew/l8kUw/deps/usr/lib/libmetis.dylib"

function vertex_sep_fast!(n, colptr::Array{Cint,1}, rowval::Array{Cint,1}, part::Array{Cint,1})
    sepSize = zeros(Cint,1)
    n2 = Cint(n)
    err = ccall((:METIS_ComputeVertexSeparator,libmetis), Cint,
                (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},     Ptr{Cint}, Ptr{Cint}),
                Ref(n2),    colptr,    rowval,    C_NULL,    metis_options, sepSize,   part)
    err == METIS_OK || error("METIS_ComputeVertexSeparator returned error code $err")
    return
end

function vertex_sep_bissect_fast!(n, colptr::Array{Cint,1}, rowval::Array{Cint,1}, part::Array{Cint,1})
    nparts = Cint(2)
    one = Cint(1)
    two = Cint(2)
    n2 = Cint(n)
    objval = zeros(Cint,1)
    err = ccall((:METIS_PartGraphRecursive,libmetis), Int32,
                (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                 Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                 Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                Ref(n2), Ref(one), colptr, rowval, C_NULL, C_NULL, C_NULL, Ref(two),
                C_NULL, C_NULL, metis_options, objval, part)
    err == METIS_OK || error("METIS_PartGraphRecursive returned error code $err")
    for i1 = 1:n # update part: those with (0), connected to at least one (1), gets a (2)
        if part[i1] == 1
            continue
        end
        for k0 = colptr[i1]:colptr[i1+1]-1 # 0-based
            k1 = k0 + 1
            j0 = rowval[k1]
            j1 = j0 + 1
            if part[j1] == 1
                part[i1] = 2 # Becomes a sep
                break
            end
        end
    end
end

