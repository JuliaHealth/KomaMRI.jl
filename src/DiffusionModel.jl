## UNDER CONSTRUCTION!
"""
Slab oriented along the x axis.

Bar, L., & Sochen, N. (2015). A spectral framework for NMR signal with restricted diffusion. Concepts in Magnetic Resonance Part A, 44(1), 16–53. doi:10.1002/cmr.a.21326 
Grebenkov, D.S. (2008), Laplacian eigenfunctions in NMR. I. A numerical tool. Concepts Magn. Reson., 32A: 277-301. doi:10.1002/cmr.a.20117
"""
function Planes(L,D=2e-9,M=30)
    #Normalization
    ϵ(m) = m==0 ? 1 : sqrt(2)
    #Eigen values
    λ(m) = π^2*m^2
    Λ =  D/L^2 * diagm([λ(m) for m=0:M])
    #Matrix A
    Ax = L * [m!=n ? ((-1)^(m+n)-1)*ϵ(m)*ϵ(n)*(λ(m)+λ(n))/(λ(m)-λ(n))^2 : 1/2 for m=0:M, n=0:M]
    A = (Ax,0,0) 
    μ = (Λ, A)
end
"""
Infinite cylinder oriented along the z axis.

Bar, L., & Sochen, N. (2015). A spectral framework for NMR signal with restricted diffusion. Concepts in Magnetic Resonance Part A, 44(1), 16–53. doi:10.1002/cmr.a.21326 
Grebenkov, D.S. (2008), Laplacian eigenfunctions in NMR. I. A numerical tool. Concepts Magn. Reson., 32A: 277-301. doi:10.1002/cmr.a.20117
"""
function Cylinder(R,D=2e-9,M=20)
    #J'n(αnk) = 0
    α = [1.841184, 3.054237, 3.831706, 4.201189, 5.317553,
         5.331443, 6.415616, 6.706133, 7.015587, 7.501266,
         8.015237, 8.536316, 8.577836, 9.282396, 9.647422, 
         9.969468, 10.17347, 10.51986, 10.71143, 11.34592]
    n = [1,2,0,3,4, 1,5,2,0,6, 3,1,7,4,8, 2,0,5,9,3]
    #Eigen values
    λ(nk) = nk==0 ? 0 : α[nk]^2
    Λ =  D/R^2 * diagm([λ(nk) for nk=0:M])
    #βnk
    β(nk) = nk==0 ? 1 : sqrt(λ(nk)/(λ(nk)-n[nk]^2))    
    #Some definitions
    cond1(mk,nk) = abs(n[mk] - n[nk]) == 1
    ϵ1(mk,nk) = sqrt( 1 + (n[mk] == 0) + (n[nk] == 0) )
    cond2(mk,nk) = (n[mk] == n[nk]-1) - (n[mk] == n[nk]+1)
    #Matrix A
    Ax = R * [cond1(i,j)*ϵ1(i,j)*β(i)*β(j)*(λ(i)+λ(j)-2*n[i]*n[j])/(λ(i)-λ(j))^2 for i=0:M, j=0:M]
    Ay = 1im*R * [cond2(i,j)*β(i)*β(j)*(λ(i)+λ(j)-2*n[i]*n[j])/(λ(i)-λ(j))^2 for i=0:M, j=0:M]
    A = (Ax,Ay,0)
    μ = (Λ, A)
end

"""
Sphere of radius R.

Bar, L., & Sochen, N. (2015). A spectral framework for NMR signal with restricted diffusion. Concepts in Magnetic Resonance Part A, 44(1), 16–53. doi:10.1002/cmr.a.21326 
Grebenkov, D.S. (2008), Laplacian eigenfunctions in NMR. I. A numerical tool. Concepts Magn. Reson., 32A: 277-301. doi:10.1002/cmr.a.20117
"""
function Sphere(R,D=2e-9,M=20)
    #j'n(αlk) = 0
    α = [2.081576, 3.342094, 4.493409, 4.514100, 5.646704, 
         5.940370, 6.756456, 7.289932, 7.725252, 7.851078, 
         8.583755, 8.934839, 9.205840, 9.840446, 10.01037, 
         10.61386, 10.90412, 11.07021, 11.07942, 11.97273]  
    l = [1,2,0,3,4, 1,5,2,0,6, 3,7,1,4,8, 2,0,5,9,3]
    m = [0,0,1,0,0, 1,0,1,2,0, 1,0,2,1,0, 2,3,1,0,2]
    #Eigen values
    λ(lm) = lm==0 ? 0 : α[lm]^2
    Λ =  D/R^2 * diagm([λ(lm) for lm=0:M])
    #βnk
    β(lm) = lm==0 ? sqrt(3/2) : sqrt((2*l[lm]+1)*λ(lm)/(λ(lm)-l[lm]*(l[lm]+1)))    
    #Some definitions
    δlδm(i,j) = (abs(l[i]-l[j])==1)*(abs(m[i]-m[j])==1)
    ϵ1(mk,nk) = (1+n[mk]+n[nk])/((2*n[mk]+1)*(2*n[nk]+1))
    ϵ2(mk,nk) = n[mk]*(n[nk]+1) + n[nk]*(n[mk]+1) + 1
    #Matrix A
    Ax = R * [cond1(mk,nk)*ϵ1(mk,nk)*β(mk)*β(nk)*(λ(mk)+λ(nk)-ϵ2(mk,nk))/(λ(mk)-λ(nk))^2 for mk=0:M, nk=0:M]
    A = (Ax,Ay,Az)
    μ = (Λ, A)
end

function SignalE(μ, seq)
    𝒊 = 1im;
    M, N = size(seq.GR)
    G = getproperty.(seq.GR,:A)
    δ = getproperty.(seq.GR[1,:],:T)
    # E = [ Π exp( -(Λ + iγ Gn⋅A) ⋅ δn ) ]_{0,0}
    E = *([exp(-(μ[1] .+ 𝒊*2π*γ*.+([μ[2][m]'*G[m,n] for m = 1:M]...))*δ[n]) for n = 1:N]...)[1,1]
end
