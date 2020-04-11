import Base.fetch; fetch(t::Vector) = [fetch(tt) for tt in t]
import Base.length; length(r1::RemoteRef)=length(fetch(r1))
addprocs(8-nprocs()) #Run on 8 processors

chunkrange(i, nc) = floor((i-1)*nc)+1:floor(i*nc)
function chunk(z, n)
    nc = length(z)/n
    y = Vector{typeof(z[1])}[]
    for i=1:n #Number of chunks
        push!(y, (typeof(z[1]))[])
        for j in chunkrange(i, nc)
            push!(y[i], z[j])
        end
    end
    y   
end

function unchunk(y)
    z = (typeof(y[1][1]))[]
    for chunk in y, datum in chunk
        push!(z, datum)
    end
    z
end

@everywhere function prefix!(zz,*)
    z=fetch(zz)
    for i in 2:length(z)
        z[i]=z[i-1]*z[i]
    end
    z
end #serial prefix

function cumsum_p(chunks, op::Function)
    np = nprocs()
    @everywhere @eval p(a,b) = [($op)(a[end], bb) for bb in b] #prefix operator
    #r = pmap(x->prefix!(fetch(x), op), chunks) #??? Why is this slower???
    r = [@spawnat i prefix!((chunks[i]), op) for i=1:np]  #pmap
    @inbounds for i=2:np #Serial reduction
        r[i] = @spawnat i p(fetch(r[i-1]), fetch(r[i]))
    end
    r
end

#Do some timings
n=1024
t=1024

np = nprocs()
nc = n/np #Chunk size

#Generate chunked data
chunks = [@spawnat i [randn(n,n) for j in chunkrange(i, nc)] for i=1:np]

#Prepare serialized data for serial calculation
data = unchunk(fetch(chunks))

tic(); prefix!(data, +); t_ser = toc()
tic(); @sync cumsum_p(chunks, +); t_par = toc() #Caution: race condition bug #4330
@printf("Serial: %.3f sec  Parallel: %.3f sec  speedup: %.3fx (theory=%.1fx)", t_ser, t_par, t_ser/t_par, np/2)