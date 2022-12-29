using KomaMRI, CUDA

#Test parameters
ns = 50_000 #20630
nt = 1_000
#Phantom
T1s = 1000.
T2s = 50.
M0s = 1.
x = zeros(ns)
T1 = T1s * ones(ns)
T2 = T2s * ones(ns)
M0 = M0s * ones(ns)
obj = Phantom(x=x,ρ=M0,T1=T1,T2=T2)
#DiscreteSequnce
T = 50.
Δt = T * ones(nt) / nt
t = [0; cumsum(Δt)[1:end-1]]
O = zeros(nt)
A = zeros(Bool, nt)
B1 = 1im*ones(nt)
seqd = KomaMRI.DiscreteSequence(O,O,O,B1,O,A,t,Δt)
#Magnetization
Mxy = 1im*ones(ns)
Mz = zeros(ns)
M = Mag(Mxy, Mz)

function BlochExcitationKernel(obj, M, seq)
    φ = 0f0
    nz = 1f0
    nxy = 0f0im
    for s ∈ seq
        #Spinor Rotation
        NVTX.@range "Rotation" begin
        #Effective field
        ΔB0 = obj.Δw./(2π*γ) .- s.Δf./γ
        B = (s.Gx.*obj.x .+ s.Gy.*obj.y .+ s.Gz.*obj.z) .+ ΔB0
        φ = -2π*γ*B.*s.Δt
        #Spinor
        α = cos.(φ./2) .- 1im.*nz.*sin.(φ./2)   
        β = -1im.*nxy.*sin.(φ./2)
        #Spinor rot
        M.xy .= 2 .*conj.(α).*β.*M.z .+ conj.(α).^2 .*M.xy .- β.^2 .*conj.(M.xy)
        M.z  .= (abs.(α).^2 .- abs.(β).^2).*M.z .- 2 .*real.(α.*β.*conj.(M.xy))
        end
        #Relaxation
        NVTX.@range "Relaxation" begin
        @. M.xy = M.xy * exp(-s.Δt / obj.T2)
        @. M.z  = M.z  * exp(-s.Δt / obj.T1) + obj.ρ * (1 - exp(-s.Δt / obj.T1))
        end
    end
    return M
end

function print_mxy()
    Mxy_mean = sum(M.xy) / ns; Mz_mean = sum(M.z) / ns
    println(exp(-sum(Δt)/T2s), " ", M0s*(1 - exp(-sum(Δt)/T1s)))
    println(Mxy_mean, " ", Mz_mean)
end
println("##################################")
println("------------ CPU Time ------------")
M.xy = ones(ns); M.z = zeros(ns)
M = BlochExcitationKernel(obj, M, seqd); print_mxy()
M.xy = ones(ns); M.z = zeros(ns)
stat_cpu = Base.@timed BlochExcitationKernel(obj, M, seqd);
println("CPU time = $(stat_cpu.time) s")
println("------------ GPU Time ------------")
M.xy = ones(ns); M.z = zeros(ns)
M = BlochExcitationKernel(KomaMRI.gpu(obj), KomaMRI.gpu(M), KomaMRI.gpu(seqd)); print_mxy()
M.xy = ones(ns); M.z = zeros(ns)
stat_gpu = CUDA.@timed BlochExcitationKernel(KomaMRI.gpu(obj), KomaMRI.gpu(M), KomaMRI.gpu(seqd));
println("GPU time = $(stat_gpu.time) s")
println("SPEEDUP = $(round(stat_cpu.time/stat_gpu.time))")
# Serious profiling
CUDA.reclaim()
CUDA.@profile begin 
    NVTX.@range "Objects to GPU" begin
        obj_gpu  = obj  |> KomaMRI.gpu
        M_gpu    = M    |> KomaMRI.gpu
        seqd_gpu = seqd |> KomaMRI.gpu
    end
    NVTX.@range "BlochExcitationKernel" BlochExcitationKernel(obj_gpu, M_gpu, seqd_gpu)
end
CUDA.reclaim()